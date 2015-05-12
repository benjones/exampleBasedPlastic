
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

#include "world.h"

void PlasticBody::loadFromJson(const Json::Value& poi,
	World& world,
	int objectIndex /* for user index var */
  ){
  auto& bulletWorld = world.bulletWorld;
  
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

  std::cout << "body bounding box:\n" << vertices.colwise().minCoeff() 
			<< '\n' << vertices.colwise().maxCoeff() << std::endl;

  if(!igl::readDMAT(directory + "/coarseWeights.dmat", boneWeights)){
	std::cout << "couldn't read skinning weights: "
			  << directory << "/coarseWeights.dmat" << std::endl;
	exit(1);
  }
  numRealBones = boneWeights.cols();

  //read in whole object parameters
  density = poi.get("density", 1000).asDouble();
  breakingThreshold = poi.get("breakingThreshold", 1e10).asDouble();
  plasticityImpulseYield = poi.get("plasticityImpulseYield", 1e6).asDouble();
  plasticityImpulseScale = poi.get("plasticityImpulseScale", 0).asDouble();
  plasticityKernelScale = poi.get("plasticityKernelScale", 1.0).asDouble();
  plasticityRate = poi.get("plasticityRate", 1.0).asDouble();
  jacobianAlpha = poi.get("jacobianAlpha", 1.0).asDouble();

  useVolumetricCollisions = poi.get("useVolumetricCollisions", true).asBool();


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
  restitutionAddition = poi.get("restitutionAddition", 1.0).asDouble();
  restitutionExponent = poi.get("restitutionExponent", 300).asDouble();

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
	constantVelocityFrames = poi.get("constantVelocityFrames", 60).asInt();
  }



  //setup CL stuff
  hostTranslations.resize(3*exampleGraph.nodes.size()*numBoneTips);
  hostRotations.resize(4*exampleGraph.nodes.size()*numBoneTips);
  for(auto j : range(numBoneTips)){
	for(auto k : range(exampleGraph.nodes.size())){
	  hostTranslations[3*(j*numNodes + k)    ] = 
		  exampleGraph.nodes[k].transformations[j].translation(0);
	  hostTranslations[3*(j*numNodes + k) + 1] = 
		  exampleGraph.nodes[k].transformations[j].translation(1);
	  hostTranslations[3*(j*numNodes + k) + 2] = 
		  exampleGraph.nodes[k].transformations[j].translation(2);

	  hostRotations[4*(j*numNodes + k)    ] =
		exampleGraph.nodes[k].transformations[j].rotation.x();
	  hostRotations[4*(j*numNodes + k) + 1] =
		exampleGraph.nodes[k].transformations[j].rotation.y();
	  hostRotations[4*(j*numNodes + k) + 2] =
		exampleGraph.nodes[k].transformations[j].rotation.z();
	  hostRotations[4*(j*numNodes + k) + 3] =
		exampleGraph.nodes[k].transformations[j].rotation.w();

	}
  }
  deviceTranslations = cl::Buffer(world.queue, 
	  hostTranslations.begin(), hostTranslations.end(),
	  true);
	  
  deviceRotations = cl::Buffer(world.queue,
	  hostRotations.begin(), hostRotations.end(), true);
  

  hostBoneWeights.resize(numRealBones*numPhysicsVertices);
  std::copy(boneWeights.data(), boneWeights.data() + hostBoneWeights.size(), 
	  hostBoneWeights.begin());
  deviceBoneWeights = cl::Buffer(world.queue,
	  hostBoneWeights.begin(), hostBoneWeights.end(), true);

  deviceBoneIndices = cl::Buffer(world.queue,
	  boneIndices.begin(), boneIndices.end(), true);

  hostUnskinnedPositions.assign(vertices.data(),
	  vertices.data() + numPhysicsVertices*3);
  deviceUnskinnedPositions = cl::Buffer(world.queue,
	  hostUnskinnedPositions.begin(), hostUnskinnedPositions.end(), true);
  

  //read in tets in each piece
  //tetmesh.#.txt
  std::string basename = "/tetmesh.";
  for(int i = 0; ; ++i){




	std::string filename = 
	  directory + basename + std::to_string(i) + ".txt";
	std::ifstream ins(filename);
	if(!ins.good() && i > 0){ //support unsliced meshes
	  break; 
	}
	
	plasticPieces.emplace_back();
	auto& piece = plasticPieces.back();

	if(!ins.good()){
	  piece.tetmeshTets = tets;
	}

	piece.initialize(directory, *this,
		vertices, world, i);


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
	piece.updateAabbs();//bulletShape->updateBound();
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
		auto kRange = range(numNodes);
		Eigen::Vector4d currentRotationCoeffs  =
		  std::accumulate(kRange.begin(), kRange.end(),
			  Eigen::Vector4d::Zero().eval(),
			  [this, vInd,boneIndex,&piece](const Eigen::Vector4d& acc, size_t k){
				return acc + piece.barycentricCoordinates(vInd, k)*
				exampleGraph.nodes[k].transformations[boneIndex].rotation.coeffs();
			  });
		Quat currentRotationQuat(currentRotationCoeffs);
		currentRotationQuat.normalize();
		Eigen::AngleAxis<double> currentRotation{currentRotationQuat};
		//Eigen::AngleAxis<double> currentRotation{piece.perVertexRotations[nodeIndex*numRealBones +
		//	  weightIndex]};
		
		jacobian.col(nodeIndex) += 
		  boneWeights(vInd, weightIndex)*
		  (node.transformations[boneIndex].translation +
			  nodeRotation.angle()*
			  nodeRotation.axis().cross(
				  //currentRotation.axis().cross(
				  currentRotation*
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
		pp.framesToSkin = 10;
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
				  std::get<1>(c) == otherIndex) && //and is active
			  std::get<3>(c)->isEnabled();
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
	std::tie(p1, p2, vInd, std::ignore) = c;
	if(bcon->isEnabled()){
	  
	  //compute the position on each and update them
	  Eigen::Vector3d pA = plasticPieces[p1].currentBulletVertexPositions.row(vInd).transpose();
	  Eigen::Vector3d pB = plasticPieces[p2].currentBulletVertexPositions.row(vInd).transpose();
	  bcon->setPivotA(eigenToBullet(pA));
	  bcon->setPivotB(eigenToBullet(pB));// + btVector3{0,0, -0.1});

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
	if(piece.framesToSkin > 0){
	  piece.skinMeshVaryingBarycentricCoords(boneWeights,
		  boneIndices,
		  exampleGraph);
	  piece.updateBulletProperties();
	}
	piece.framesToSkin--;
  }
}

void PlasticBody::skinAndUpdateCL(World& world, cl::Kernel& clKernel){

  for(auto& piece : plasticPieces){
	if(piece.framesToSkin > 0){
	  piece.skinMeshOpenCL(world, *this, clKernel);
	  piece.updateBulletProperties();
	}
	piece.framesToSkin--;
  }
  

}

