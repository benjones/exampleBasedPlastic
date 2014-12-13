
#include "plasticBody.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>

#include <BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h>

#include "eigenTypedefs.h"
#include "readDMAT.h"
#include "readMESH.h"

#include "utils.h"
#include "kernels.h"

#include "range.hpp"
#include "enumerate.hpp"
using benlib::range;
using benlib::enumerate;

#include "json/json.h"

void PlasticBody::loadFromJson(const Json::Value& poi,
	btDiscreteDynamicsWorld& bulletWorld,
	int objectIndex /* for user index var */
  ){
  
  const std::string directory = poi["directory"].asString();
  std::cout << "reading from directory: " << directory << std::endl;

  computeBoneIndices(directory + "/bone_roots.bf");
  numBoneTips = boneIndices.size();

  exampleGraph.load(directory + "/exampleGraph.txt");
  numNodes = exampleGraph.nodes.size();

  RMMatrix3d vertices;
  RMMatrix4i tets;
  RMMatrix3i triangles;

  if(!igl::readMESH(directory + "/mesh.mesh",
	  vertices, tets, triangles)){
	std::cout << "couldn't read tetmesh: " 
			  << directory << "/mesh.mesh" << std::endl;
	exit(1);
  }
  numPhysicsVertices = vertices.rows();
  std::cout << "vertices read size: " << vertices.rows() << ' ' << vertices.cols() << std::endl;
  //assert(vertices.allFinite());
  if(!igl::readDMAT(directory + "/coarseWeights.dmat", boneWeights)){
	std::cout << "couldn't read skinning weights: "
			  << directory << "/coarseWeights.dmat" << std::endl;
	exit(1);
  }
  numRealBones = boneWeights.cols();

  //read in whole object parameters
  density = poi.get("density", 1000).asDouble();
  plasticityImpulseYield = poi.get("plasticityImpulseYield", 1e6).asDouble();
  plasticityImpulseScale = poi.get("plasticityImpulseScale", 0).asDouble();
  plasticityKernelScale = poi.get("plasticityKernelScale", 1.0).asDouble();
  plasticityRate = poi.get("plasticityRate", 1.0).asDouble();
  jacobianAlpha = poi.get("jacobianAlpha", 1.0).asDouble();

  perExampleScale.resize(numNodes);
  auto& pesIn = poi["perExampleScale"];
  if(!pesIn.isNull()){
	for(auto j : range(pesIn.size())){
	  perExampleScale[j] = pesIn[j].asDouble();
	}
  } else {
	perExampleScale.setOnes();
  }
  
  scaleFactor = poi.get("scaleFactor", 1).asDouble();
  restitution = poi.get("restitution", 0.8).asDouble();

  //apply COM offset stuff
  btVector3 offset(0.0,0.0,0.0); //default to no offset
  auto offsetIn = poi["offset"];
  if(!offsetIn.isNull() && offsetIn.isArray()){
	assert(offsetIn.size() == 3);
	offset = btVector3{offsetIn[0].asDouble(),
					   offsetIn[1].asDouble(),
					   offsetIn[2].asDouble()};
  }
	
  auto rotation = btQuaternion::getIdentity();
  auto rotationIn = poi["rotation"];
  if(!rotationIn.isNull()){
	btVector3 axis{rotationIn["axis"][0].asDouble(),
		rotationIn["axis"][1].asDouble(),
		rotationIn["axis"][2].asDouble()};
	rotation = btQuaternion{axis,
							rotationIn["angle"].asDouble()};
	
	std::cout << "rotation: " << rotation << std::endl;
  }
  
  auto constantVelocityIn = poi["constantVelocity"];
  if(!constantVelocityIn.isNull() && 
	constantVelocityIn.isArray() && 
	constantVelocityIn.size() == 3){
	
	hasConstantVelocity = true;
	constantVelocity = btVector3{constantVelocityIn[0].asDouble(),
								 constantVelocityIn[1].asDouble(),
								 constantVelocityIn[2].asDouble()};
  }

  //read in tets in each piece
  //tetmesh.#.txt
  std::string basename = "/tetmesh.";
  for(int i = 0; ; ++i){

	std::string filename = 
	  directory + basename + std::to_string(i) + ".txt";
	std::ifstream ins(filename);
	if(!ins.good()){ break; }
	
	plasticPieces.emplace_back();
	auto& piece = plasticPieces.back();
	piece.scaleFactor = scaleFactor;
	piece.tetmeshVertices = vertices;
	piece.currentBulletVertexPositions = 
	  piece.scaleFactor*piece.tetmeshVertices;
	
	//assert(piece.currentBulletVertexPositions.allFinite());

	//read tets in this piece
	int nTets;
	ins >> nTets;
	piece.tetmeshTets.resize(nTets, 4);
	std::copy(std::istream_iterator<int>(ins),
	  std::istream_iterator<int>(),
	  piece.tetmeshTets.data());

	piece.numPhysicsVertices = vertices.rows();
	piece.computeTriangleFaces(); //and triangles
	piece.computeVertexNeighbors(); //connectivity in this piece
	
	piece.computeActiveVertices();

	//initialize other stuff
	piece.barycentricCoordinates.resize(piece.numPhysicsVertices,
	  numNodes);
	piece.barycentricCoordinates.setZero();
	piece.barycentricCoordinates.col(0).setOnes();

	piece.deltaBarycentricCoordinates.resize(piece.numPhysicsVertices,
	  numNodes);
	piece.deltaBarycentricCoordinates.setZero();
	
	piece.skinMeshVaryingBarycentricCoords(boneWeights,
	  boneIndices, exampleGraph);


	//setup rigid bodies
	piece.btTriMesh = std::unique_ptr<btTriangleIndexVertexArray>{
	  new btTriangleIndexVertexArray{
		static_cast<int>(piece.tetmeshTriangles.rows()),
		piece.tetmeshTriangles.data(),
		3*sizeof(decltype(piece.tetmeshTriangles)::Scalar),
		static_cast<int>(piece.currentBulletVertexPositions.rows()),
		piece.currentBulletVertexPositions.data(),
		3*sizeof(decltype(piece.currentBulletVertexPositions)::Scalar)
	  }};
	piece.bulletShape = std::unique_ptr<btGImpactMeshShape>{
	  new btGImpactMeshShape{piece.btTriMesh.get()}};
	
	piece.motionState = 
	  std::unique_ptr<btDefaultMotionState>{
	  new btDefaultMotionState{}};
	
	piece.mass = 1; //worked before...?
	piece.bulletBody = std::unique_ptr<btRigidBody>{
	  new btRigidBody{piece.mass,
					  piece.motionState.get(),
					  piece.bulletShape.get()}};

	piece.computeMassesAndVolume(density);
  }
  //compute whole object mass and COM:
  mass = std::accumulate(
	  plasticPieces.begin(), plasticPieces.end(),
	  0.0, 
	  [](double acc, const PlasticPiece& piece){
		return acc + piece.mass;
	  });
		
  Eigen::Vector3d com =
	std::accumulate(plasticPieces.begin(),
		plasticPieces.end(),
		Eigen::Vector3d::Zero().eval(),
		[](const Eigen::Vector3d& acc,
			const PlasticPiece& piece){
		  return 
		  (acc + 
			  (piece.tetmeshVertexMasses.asDiagonal()*
				  piece.currentBulletVertexPositions)
			  .colwise().sum().transpose()
			  /piece.mass).eval();
		});
  
  offset -= quatRotate(rotation, eigenToBullet(com));
  
  //set inertia stuff, collision stuff, and add to world
  for(auto& piece : plasticPieces){
	piece.inertiaAligningTransform.setIdentity();
	piece.worldTransform = btTransform{rotation, offset};
	piece.bulletBody->setCenterOfMassTransform(piece.worldTransform);
	piece.saveBulletSnapshot();
	piece.updateBulletProperties();
	piece.saveBulletSnapshot();
	piece.bulletShape->updateBound();
	piece.bulletBody->setRestitution(restitution);
	piece.bulletBody->setUserIndex(objectIndex);
	bulletWorld.addRigidBody(piece.bulletBody.get());
  }
  //setup constraint stuff
  computeConstraints();
  std::cout << "num constraints: " << constraints.size() << std::endl;
  updateConstraints(); //compute the correct positions
  for(auto& c : constraints){
	bulletWorld.addConstraint(std::get<3>(c).get(), true);
  }
  
}

void PlasticBody::computeBoneIndices(const std::string& boneFile){
  std::ifstream ins(boneFile);
  boneIndices.clear();
  
  int wi;
  double trash;

  while(ins){
	ins >> wi;
	for(auto i : range(26)){ ins >> trash;}
	if(ins.good()){
	  boneIndices.push_back(wi);
	}
  }
  
}


void PlasticBody::saveBulletSnapshots(){
  for(auto& piece : plasticPieces){
	piece.saveBulletSnapshot();
  }
}
void PlasticBody::restoreBulletSnapshots(){
  for(auto& piece: plasticPieces){
	piece.restoreBulletSnapshot();
  }
}


Eigen::Vector3d PlasticBody::getDeformationVectorFromImpulse (
	const PlasticPiece& piece,
	const btManifoldPoint& manPoint,
	double dt,
	bool isObject0) const{
  
  
  const btVector3 impulseGlobal = 
	manPoint.m_normalWorldOnB*manPoint.m_appliedImpulse;
  btVector3 displacementGlobal = impulseGlobal*dt/mass;
  if(!isObject0){displacementGlobal *= -1;}

  const auto currentRotation = piece.bulletBody->getCenterOfMassTransform().getRotation()*
	piece.inertiaAligningTransform.getRotation().inverse();
  const btVector3 displacementLocal = 
	quatRotate(currentRotation.inverse(), displacementGlobal);
  
  if(displacementLocal.norm() < plasticityImpulseYield){
	return Eigen::Vector3d::Zero();
  }
  
  const auto normalizedDisplacement = displacementLocal.normalized();
  const auto displacementToApply = plasticityImpulseScale*
	(displacementLocal - plasticityImpulseYield*normalizedDisplacement);
  
  return bulletToEigen(displacementToApply);

}

void PlasticBody::projectImpulsesOntoExampleManifoldLocally(double dt){
  if(plasticityImpulseScale == 0){
	return;
  }

  Eigen::MatrixXd jacobian;
  for(auto& tup : manifoldPoints){
	const auto* body = std::get<0>(tup);
	auto& manPoint = std::get<1>(tup);
	bool isObject0 = std::get<2>(tup);
	
	//find out which piece it is
	auto it = std::find_if(plasticPieces.begin(), plasticPieces.end(),
		[body](const PlasticPiece& piece){
		  return piece.bulletBody.get() == body;
		});
	assert(it != plasticPieces.end());
	auto& piece = *it;
	const auto pieceIndex = std::distance(plasticPieces.begin(), it);

	auto& localPoint = 
	  isObject0 ? manPoint.m_localPointA : manPoint.m_localPointB;
	auto vInd = piece.getNearestVertex(bulletToEigen(localPoint));
	
	auto impulseAtContact = 
	  getDeformationVectorFromImpulse(piece, manPoint, dt, isObject0);
	if(impulseAtContact.squaredNorm() <= 0){continue;} //ignore 0 impulses

	//compute the derivative
	jacobian.resize(3, numNodes);
	jacobian.setZero();
	
	for(auto nodeIndex : range(numNodes)){
	  auto& node = exampleGraph.nodes[nodeIndex];
	  for(auto boneIndex : range(numBoneTips)){
		auto weightIndex = boneIndices[boneIndex];//bones[boneIndex]->get_wi();
		Eigen::AngleAxis<double> nodeRotation{node.transformations[weightIndex].rotation};
		Eigen::AngleAxis<double> currentRotation{piece.perVertexRotations[nodeIndex*numRealBones +
			  weightIndex]};
		
		jacobian.col(nodeIndex) += 
		  boneWeights(vInd, weightIndex)*
		  (node.transformations[boneIndex].translation +
			  nodeRotation.angle()*
			  currentRotation.axis().cross(currentRotation*
				  piece.tetmeshVertices.row(vInd).transpose()));
	  }
	  
	}
	
	
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian, 
										  Eigen::ComputeThinU | 
										  Eigen::ComputeThinV);

	Eigen::VectorXd singularVectorContribution = 
	  svd.matrixV()*impulseAtContact.norm()*(svd.singularValues().asDiagonal()*
		  svd.matrixU().transpose()*
		  impulseAtContact).normalized();
	Eigen::VectorXd jacobianTransposeContribution = 
	  jacobian.transpose()*impulseAtContact;

	Eigen::VectorXd deltaS = jacobianAlpha*singularVectorContribution +
	  (1.0 - jacobianAlpha)*jacobianTransposeContribution;

	//scale it
	deltaS /= std::max(std::fabs(deltaS.maxCoeff()), 1.0);
	
	if(deltaS.squaredNorm() > 0.000001){

	  computeGeodesicDistances(pieceIndex, vInd, 1.0/plasticityKernelScale);
	  
	  //distribute this to evereyone else
	  for(auto& pp : plasticPieces){
		for(auto i : pp.activeVertices){
		  auto scale = 
			Kernels::simpleCubic(plasticityKernelScale*pp.geodesicDistances[i]);
		  pp.deltaBarycentricCoordinates.row(i) += scale*deltaS.transpose();
		}
	  }
	}
  }

  for(auto& pp : plasticPieces){
	pp.barycentricCoordinates += 
	  plasticityRate*pp.deltaBarycentricCoordinates*perExampleScale.asDiagonal();
	pp.deltaBarycentricCoordinates *= 1.0 - plasticityRate;

  
	//clamp and rescale:
	for(auto i : pp.activeVertices){
	  for(auto j : range(numNodes)){
		pp.barycentricCoordinates(i, j) = 
		  std::max(pp.barycentricCoordinates(i,j), 0.0);
	  }
	  const double sum = pp.barycentricCoordinates.row(i).sum();
	  if(fabs(sum) < 0.0001){
		pp.barycentricCoordinates(i, 0) = 1;
	  } else {
		pp.barycentricCoordinates.row(i) /= sum;
	  }
	}
  }

}

void PlasticBody::computeGeodesicDistances(size_t pieceIndex, size_t vInd, double radius){

  auto& mainPiece = plasticPieces[pieceIndex];
  mainPiece.geodesicDistances.assign(numPhysicsVertices,
	  std::numeric_limits<double>::infinity());
  mainPiece.geodesicDistances[vInd] = 0;
  mainPiece.geodesicDistancesPropogate({vInd}, radius);
  
  //propgate to other pieces
  std::vector<size_t> seeds;
  for(auto&& pr : enumerate(plasticPieces)){
	auto otherIndex = pr.first;
	auto& otherPiece = pr.second;
	if(otherIndex != pieceIndex){
	  //everything is infinity
	  otherPiece.geodesicDistances.assign(numPhysicsVertices,
		  std::numeric_limits<double>::infinity());
	  //compute seeds
	  seeds.clear();
	  transform_if(constraints.begin(), constraints.end(),
		  std::back_inserter(seeds),
		  [](const Constraint& c){
			return std::get<2>(c); //return vInd
		  },
		  [pieceIndex, otherIndex](const Constraint& c){
			//if the constraint is between the two pieces we're looking at
			return 
			  (std::get<0>(c) == pieceIndex ||
				  std::get<1>(c) == pieceIndex) &&
			  (std::get<0>(c) == otherIndex ||
				  std::get<1>(c) == otherIndex);
		  });
	  //assign distances
	  for(auto seed : seeds){
		otherPiece.geodesicDistances[seed] = 
		  mainPiece.geodesicDistances[seed];
	  }
	  //propogate
	  otherPiece.geodesicDistancesPropogate(seeds, radius);

	}
  }



}

void PlasticBody::computeConstraints(){
  constraints.clear(); //just in case
  std::vector<size_t> sharedVertices;
  for(auto i : range(plasticPieces.size())){
	for(auto j : range(i+1, plasticPieces.size())){
	  sharedVertices.clear();
	  std::set_intersection(
		  plasticPieces[i].activeVertices.begin(), 
		  plasticPieces[i].activeVertices.end(),
		  plasticPieces[j].activeVertices.begin(), 
		  plasticPieces[j].activeVertices.end(),
		  std::back_inserter(sharedVertices));

	  std::transform(sharedVertices.begin(), sharedVertices.end(),
		  std::back_inserter(constraints),
		  [i,j,this](size_t vInd){ 
			return std::make_tuple(i, j, vInd, 
				std::unique_ptr<btPoint2PointConstraint>{
				  new btPoint2PointConstraint(*plasticPieces[i].bulletBody,
					  *plasticPieces[j].bulletBody,
					  btVector3{}, btVector3{})}// set up the anchors 
			  );});
	}
  }
  
}

void PlasticBody::updateConstraints(){
  for(auto& c : constraints){
	size_t p1, p2, vInd;
	btPoint2PointConstraint* bcon = std::get<3>(c).get();
	bcon->setBreakingImpulseThreshold(50);
	std::tie(p1, p2, vInd, std::ignore) = c;
	if(bcon->isEnabled()){
	  
	  //compute the position on each and update them
	  Eigen::Vector3d pA = plasticPieces[p1].currentBulletVertexPositions.row(vInd).transpose();
	  Eigen::Vector3d pB = plasticPieces[p2].currentBulletVertexPositions.row(vInd).transpose();
	  bcon->setPivotA(eigenToBullet(pA));
	  bcon->setPivotB(eigenToBullet(pB));

	}
  }
}


int PlasticBody::dump(int currentFrame, int objectStart) const{

  static const std::string base = "frames/foo-";
  std::stringstream currentFrameStream;
  currentFrameStream << std::setfill('0') << std::setw(4) << currentFrame;
  std::string currentFrameString = currentFrameStream.str();

  for(auto& piece : plasticPieces){
	std::stringstream objectIndexStream;
	objectIndexStream << std::setfill('0') << std::setw(3) << objectStart;
	++objectStart;

	std::string start = base + objectIndexStream.str();
	std::cout << "writing " << start + '.' + currentFrameString + ".ply" << std::endl;
	piece.dumpPly(start + '.' + currentFrameString + ".ply");
	piece.dumpBcc(start + '.' + currentFrameString + ".bcs");
  }
  return objectStart;
}

void PlasticBody::skinAndUpdate(){
  for(auto& piece : plasticPieces){
	piece.skinMeshVaryingBarycentricCoords(boneWeights,
		boneIndices,
		exampleGraph);
	piece.updateBulletProperties();
  }
}
