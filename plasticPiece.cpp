
#include "plasticPiece.h"
#include "plasticBody.h"
#include "world.h"

#include "range.hpp"
#include "enumerate.hpp"
using benlib::range;
using benlib::enumerate;

#include "utils.h"
#include "plyio.hpp"
#include "exampleGraph.h"

#include <tbb/tbb.h>

#include <array>
#include <numeric>
#include <set>

//#include <BulletCollision/BroadphaseCollision/btDbvt.h>

void PlasticPiece::computeMassesAndVolume(double density){

  volume = 0;
  mass = 0;
  
  tetmeshVertexMasses.resize(numPhysicsVertices);
  tetmeshVertexMasses.setZero();

  Eigen::Matrix3d basis(3,3);
  for(auto tetInd : range(tetmeshTets.rows())){
	basis.row(0) = 
	  scaleFactor*(tetmeshVertices.row(tetmeshTets(tetInd,1)) - 
		tetmeshVertices.row(tetmeshTets(tetInd, 0)));
	basis.row(1) = 
	  scaleFactor*(tetmeshVertices.row(tetmeshTets(tetInd,2)) - 
		tetmeshVertices.row(tetmeshTets(tetInd, 0)));
	basis.row(2) = 
	  scaleFactor*(tetmeshVertices.row(tetmeshTets(tetInd,3)) - 
		tetmeshVertices.row(tetmeshTets(tetInd, 0)));
	
	//shouldn't be necessary, but whatever
	double tetVolume = -basis.determinant()/6.0;
	if(tetVolume < 0){
	  std::cout << "inverted tet!!!!!!" << std::endl;
	  std::cout << " at tet: " << tetInd 
				<< " vol " << tetVolume << std::endl;
	}
	
	if(std::binary_search(cutTets.begin(), cutTets.end(), tetInd)){
	  tetVolume /= 2; //this tet is cut, distribute it's mass to both pieces
	}

	volume += tetVolume;
	double nodeMass = 0.25*tetVolume*density;
	mass += 4*nodeMass;
	for(auto j : {0,1,2,3}){
	  tetmeshVertexMasses(tetmeshTets(tetInd,j)) += nodeMass;
	}
  }
  //assert(tetmeshVertexMasses.allFinite());
}

void PlasticPiece::updateBulletProperties(){
  //the world part of the transform before we do this
  worldTransform = 
	bulletSnapshot.comTransform;

  Eigen::Vector3d centerOfMass = 
	(tetmeshVertexMasses.asDiagonal()*
	 currentBulletVertexPositions).colwise().sum().transpose()/
	mass;

  inertiaTensor.setZero();
  
  for(auto i : activeVertices){
	Eigen::Vector3d r = 
	  currentBulletVertexPositions.row(i).transpose() - 
	  centerOfMass;
	inertiaTensor += tetmeshVertexMasses(i)*
	  (r.squaredNorm()*Eigen::Matrix3d::Identity() -
		  r*r.transpose());
  }
    
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(inertiaTensor);
  Eigen::Matrix3d evecs = solver.eigenvectors();
  Eigen::Vector3d evals = solver.eigenvalues();
  if(evecs.determinant() < 0){
	//make it a rotation, not a refelction
	evecs.col(2) *= -1;
  }


  //translate and then rotate each vertex so that
  //the object is centered at the origin and aligned with the 
  //princial moments of inertia

  for(auto i : activeVertices){
	
	currentBulletVertexPositions.row(i) = 
	  (evecs.transpose()*
	   (currentBulletVertexPositions.row(i).transpose() - 
		centerOfMass)
	   ).transpose();
	
  }
  //update the spliting plane vertices too
  for(auto i : range(splittingPlaneVertices.rows())){
	currentBulletSplittingPlaneVertices.row(i) =
	  (evecs.transpose()*
		  (currentBulletSplittingPlaneVertices.row(i).transpose() -
			  centerOfMass)
		).transpose();
	if(!currentBulletSplittingPlaneVertices.row(i).allFinite()){
	  std::cout << "row i sucks: " 
				<< currentBulletSplittingPlaneVertices.row(i) << std::endl;
	  exit(1);
	}
  }
  



  //update the tet vertices
  if(useVolumetricCollisions){
	for(auto&& pr : enumerate(bulletTets)){
	  for(auto i : range(4)){
		Eigen::Vector3d ePos = 
		  currentBulletVertexPositions.row(
			  tetmeshTets(collisionTetIndices[pr.first], i));
		pr.second.m_vertices[i] = eigenToBullet(ePos);
	  }
	}
  } else {
	bulkMesh->postUpdate();
	bulkMesh->updateBound();
  }


  if(fractureMesh){
	fractureMesh->postUpdate();
	fractureMesh->updateBound();
  }
  updateAabbs();
  

  bulletBody->setActivationState(DISABLE_DEACTIVATION);
  
  inertiaAligningTransform = btTransform{eigenToBullet(evecs),
  										 eigenToBullet(centerOfMass)};

  bulletBody->setCenterOfMassTransform(worldTransform*
									   inertiaAligningTransform);

  bulletBody->setMassProps(mass, eigenToBullet(evals));
  bulletBody->updateInertiaTensor(); //rotate it
  bulletBody->setDamping(0, 0.01);
}

int PlasticPiece::getNearestVertex(const Eigen::Vector3d& localPoint) const {

  int best = 0;
  double sqrdist = std::numeric_limits<double>::max();
  for(auto i : activeVertices){
	double thisDist = (localPoint - currentBulletVertexPositions.row(i).transpose()).squaredNorm();
	if(thisDist < sqrdist){ best = i; sqrdist = thisDist;}
  }
  return best;

  return *std::min_element(activeVertices.begin(), activeVertices.end(),
	[this, &localPoint](int a, int b){
	  return 
		(localPoint - 
		  currentBulletVertexPositions.row(a).transpose()).squaredNorm()
		< (localPoint - 
		  currentBulletVertexPositions.row(b).transpose()).squaredNorm();
	});
}


void PlasticPiece::saveBulletSnapshot(){
  bulletSnapshot.linearVelocity = bulletBody->getLinearVelocity();
  bulletSnapshot.angularMomentum = bulletBody->getInvInertiaTensorWorld().inverse()*
	bulletBody->getAngularVelocity();
  //keep the inertia part out of it
  bulletSnapshot.comTransform = bulletBody->getCenterOfMassTransform()*
	inertiaAligningTransform.inverse();
}
 
void PlasticPiece::restoreBulletSnapshot(){
  bulletBody->setLinearVelocity(bulletSnapshot.linearVelocity);
  bulletBody->setAngularVelocity(bulletBody->getInvInertiaTensorWorld()*bulletSnapshot.angularMomentum);
  bulletBody->setCenterOfMassTransform(bulletSnapshot.comTransform*inertiaAligningTransform);

}

void PlasticPiece::skinMeshVaryingBarycentricCoords(
  const Eigen::MatrixXd& boneWeights,
  const std::vector<int>& boneIndices,
  const ExampleGraph& exampleGraph){
  
  //const auto numRealBones = boneWeights.cols();
  const auto numBoneTips = boneIndices.size();
  const auto numNodes = exampleGraph.nodes.size();

  //perVertexTranslations.resize(numPhysicsVertices*numRealBones);
  //perVertexRotations.resize(numPhysicsVertices*numRealBones);

  //should be a no-op, but could cause issues, i guess?
  //currentBulletVertexPositions.resize(numPhysicsVertices, 3);
  currentBulletVertexPositions.setZero();


  //todo just use activeVertices here...
  tbb::parallel_for(
	tbb::blocked_range<size_t>(0, numPhysicsVertices),
	[&](const tbb::blocked_range<size_t>& r){
	  for(auto i = r.begin(); i != r.end(); ++i){
		for(auto j : range(numBoneTips)){
		  if(boneIndices[j] < 0){continue;} //bones[j]->get_wi() < 0) //not a real bone
		  
		  //const auto kRange = range(numNodes);
		  Eigen::Vector3d translation = Eigen::Vector3d::Zero();
		  Eigen::Vector4d rotationCoeffs = Eigen::Vector4d::Zero();
		  for(auto k : range(numNodes)){
			translation += barycentricCoordinates(i,k)*
			  exampleGraph.nodes[k].transformations[j].translation;
			rotationCoeffs += barycentricCoordinates(i,k)*
			  exampleGraph.nodes[k].transformations[k].rotation.coeffs();
		  }

/*
		  //compute bone transform
		  const Eigen::Vector3d translation = 
			std::accumulate(kRange.begin(), kRange.end(),
			  Eigen::Vector3d::Zero().eval(),
			  [this, i,j,&exampleGraph](const Eigen::Vector3d& acc, size_t k){
				return (acc + barycentricCoordinates(i, k)*
				  exampleGraph.nodes[k].transformations[j].translation).eval();
			  });
		  const Eigen::Vector4d rotationCoeffs =
			std::accumulate(kRange.begin(), kRange.end(),
			  Eigen::Vector4d::Zero().eval(),
			  [this, i,j,&exampleGraph](const Eigen::Vector4d& acc, size_t k){
				return acc + barycentricCoordinates(i, k)*
				exampleGraph.nodes[k].transformations[j].rotation.coeffs();
			  });
*/
		  Quat rotation(rotationCoeffs);
		  rotation.normalize();
		  
		  //rotation = Quat::Identity();

		  const auto weightIndex = boneIndices[j];//bones[j]->get_wi();
		  currentBulletVertexPositions.row(i) +=
			scaleFactor*boneWeights(i,weightIndex)*
			(rotation*(tetmeshVertices.row(i).transpose()) + 
			  translation).transpose();
		  
		  //perVertexTranslations[i*numRealBones + weightIndex] = translation;
		  //perVertexRotations[i*numRealBones + weightIndex] = rotation;
		  
		}
	  }
	});
  //assert(currentBulletVertexPositions.allFinite());

  //skin the clipping plane vertices
  for(auto i : range(splittingPlaneVertices.rows())){
	auto& inds = std::get<0>(clippingPlaneTetInfo[i]);
	auto& bcs = std::get<1>(clippingPlaneTetInfo[i]);
	currentBulletSplittingPlaneVertices.row(i) =
	  bcs(0)*currentBulletVertexPositions.row(inds(0)) +
	  bcs(1)*currentBulletVertexPositions.row(inds(1)) +
	  bcs(2)*currentBulletVertexPositions.row(inds(2)) +
	  bcs(3)*currentBulletVertexPositions.row(inds(3));
  }


}

RMMatrix3i PlasticPiece::computeTriangleFaces(
	const std::vector<size_t>& tetIndices) const{
  
  //use sorted vector instead of set if speed is an issue?
  std::set<std::tuple<int,int,int>> seenFaces;
  int triIndices[4][3] = { {0,1,2},{0,1,3},{0,2,3},{1,2,3}};
  auto makeTri = [](int a, int b, int c) -> std::tuple<int,int,int>{
	int max = std::max({a, b, c});
	int min = std::min({a, b, c});
	int mid = ((a != max) && ( a != min)) ? a :
	( ((b != max) && ( b != min)) ? b : c);
	return std::tuple<int,int,int>(min,mid,max);
  };
  //loop only over collision tets so we don't output overlapping tets
  for(auto i : tetIndices){ //range(tetmeshTets.rows())){
	for(auto j : range(4)){
	  auto tri = makeTri(tetmeshTets(i, triIndices[j][0]),
		tetmeshTets(i, triIndices[j][1]),
		tetmeshTets(i, triIndices[j][2]));
	  auto it = seenFaces.find(tri);
	  if(it == seenFaces.end()){
		seenFaces.insert(tri);
	  } else {
		seenFaces.erase(it);
	  }
	}
  }
  
  //std::cout << "mesh file tris: " << tetmeshTriangles.rows() << std::endl;
  //std::cout << "computed " << seenFaces.size() << std::endl;
  
  //orientation is ignored by bullet, so don't reorient stuff
  RMMatrix3i ret(seenFaces.size(), 3);

  for(auto&& pr : enumerate(seenFaces)){
	ret(pr.first, 0) = std::get<0>(pr.second);
	ret(pr.first, 1) = std::get<1>(pr.second);
	ret(pr.first, 2) = std::get<2>(pr.second);
  }
  return ret;
}

RMMatrix3i PlasticPiece::computeAllTriangles(
	const std::vector<size_t>& tetIndices) const{
  
  //use sorted vector instead of set if speed is an issue?
  std::set<std::tuple<int,int,int>> seenFaces;
  int triIndices[4][3] = { {0,1,2},{0,1,3},{0,2,3},{1,2,3}};
  auto makeTri = [](int a, int b, int c) -> std::tuple<int,int,int>{
	int max = std::max({a, b, c});
	int min = std::min({a, b, c});
	int mid = ((a != max) && ( a != min)) ? a :
	( ((b != max) && ( b != min)) ? b : c);
	return std::tuple<int,int,int>(min,mid,max);
  };
  //loop only over collision tets so we don't output overlapping tets
  for(auto i : tetIndices){ //range(tetmeshTets.rows())){
	for(auto j : range(4)){
	  auto tri = makeTri(tetmeshTets(i, triIndices[j][0]),
		  tetmeshTets(i, triIndices[j][1]),
		  tetmeshTets(i, triIndices[j][2]));
	  seenFaces.insert(tri);  //always insert
	}
  }
  
  //std::cout << "mesh file tris: " << tetmeshTriangles.rows() << std::endl;
  //std::cout << "computed " << seenFaces.size() << std::endl;
  
  //orientation is ignored by bullet, so don't reorient stuff
  RMMatrix3i ret(seenFaces.size(), 3);

  for(auto&& pr : enumerate(seenFaces)){
	ret(pr.first, 0) = std::get<0>(pr.second);
	ret(pr.first, 1) = std::get<1>(pr.second);
	ret(pr.first, 2) = std::get<2>(pr.second);
  }
  return ret;
}



void PlasticPiece::computeVertexNeighbors(){

  vertexNeighbors.resize(numPhysicsVertices);
  neighborDistances.resize(numPhysicsVertices);
  for(auto row : range(tetmeshTets.rows())){
	for(auto i : range(4)){
	  vertexNeighbors[tetmeshTets(row, i)].push_back(tetmeshTets(row, (i + 1)%4));
	  vertexNeighbors[tetmeshTets(row, i)].push_back(tetmeshTets(row, (i + 2)%4));
	  vertexNeighbors[tetmeshTets(row, i)].push_back(tetmeshTets(row, (i + 3)%4));
	}
  }

  //trash duplicates
  for(auto vInd : range(numPhysicsVertices)){
	auto& vlist = vertexNeighbors[vInd];
	std::sort(vlist.begin(), vlist.end());
	vlist.erase(std::unique(vlist.begin(), vlist.end()),
				vlist.end());
	
	neighborDistances[vInd].resize(vlist.size());
	std::transform(vlist.begin(), vlist.end(),
				   neighborDistances[vInd].begin(),
				   [this, vInd](size_t nInd){
					 return (tetmeshVertices.row(vInd) - tetmeshVertices.row(nInd)).norm();
				   });
  }
}


void PlasticPiece::geodesicDistancesPropogate(const std::vector<size_t>& seeds, double radius){

  const size_t passes = 4;
  std::deque<size_t> queue;
  for(auto pass : range(passes)){
	(void)pass; //don't use it

	queue.clear();
	for(auto seed : seeds){
	  queue.insert(queue.end(), vertexNeighbors[seed].begin(), vertexNeighbors[seed].end());
	}
	while(!queue.empty()){
	  const size_t currentVertex = queue.front();
	  queue.pop_front();
	  const double oldDist = geodesicDistances[currentVertex];
	  auto neighborRange = benlib::range(vertexNeighbors[currentVertex].size());
	  auto minIt = 
		std::min_element(neighborRange.begin(),
		  neighborRange.end(),
		  [this, currentVertex](size_t ind1, size_t ind2){
			return (geodesicDistances[vertexNeighbors[currentVertex][ind1]] + 
			  neighborDistances[currentVertex][ind1]) <
			(geodesicDistances[vertexNeighbors[currentVertex][ind2]] +
			  neighborDistances[currentVertex][ind2]);
		  });
	  double bestDist = geodesicDistances[vertexNeighbors[currentVertex][*minIt]] + 
		neighborDistances[currentVertex][*minIt];
	  if(bestDist < oldDist &&
		bestDist < radius){
		geodesicDistances[currentVertex] = bestDist;
		queue.insert(queue.end(), 
		  vertexNeighbors[currentVertex].begin(), 
		  vertexNeighbors[currentVertex].end());
	  }
	}
  }
}

void PlasticPiece::computeActiveVertices(){
  activeVertices.resize(4*tetmeshTets.rows());
  std::copy(tetmeshTets.data(), tetmeshTets.data() + 4*tetmeshTets.rows(),
	activeVertices.begin());
  
  std::sort(activeVertices.begin(), activeVertices.end());
  activeVertices.erase(std::unique(activeVertices.begin(),
	  activeVertices.end()),
	activeVertices.end());

  activeVertices.shrink_to_fit();
}

void PlasticPiece::dumpPly(const std::string& filename) const{
  std::ofstream outs(filename);
  
  auto bulletTrans = bulletBody->getCenterOfMassTransform();
  auto eigenRot = bulletToEigen(bulletTrans.getBasis());
  Eigen::Vector3d eigenTrans{bulletToEigen(bulletTrans.getOrigin())};
  RMMatrix3f transformedVertices(currentBulletVertexPositions.rows() + 
	  currentBulletSplittingPlaneVertices.rows(), 3);
  for(auto i : range(currentBulletVertexPositions.rows())){
	Eigen::Vector3d transformedVertex = 
	  eigenRot*currentBulletVertexPositions.row(i).transpose() +
	  eigenTrans;
	Eigen::Vector3f transFloat = transformedVertex.template cast<float>();
	transformedVertices.row(i) = transFloat.transpose();
  }

  assert(currentBulletVertexPositions.allFinite());
  for(auto i : range(currentBulletSplittingPlaneVertices.rows())){
	Eigen::Vector3d transformedVertex = 
	  eigenRot*currentBulletSplittingPlaneVertices.row(i).transpose() +
	  eigenTrans;
	Eigen::Vector3f transFloat = transformedVertex.template cast<float>();
	transformedVertices.row(i + currentBulletVertexPositions.rows()) = transFloat.transpose();
  }
  
  RMMatrix3i combinedTriangles(tetmeshTriangles.rows() +
	  splittingPlaneTriangles.rows(), 3);

  combinedTriangles.block(0,0,tetmeshTriangles.rows(),3) = tetmeshTriangles;

  combinedTriangles.block(tetmeshTriangles.rows(), 0, 
	  splittingPlaneTriangles.rows(), 3) = splittingPlaneTriangles;
  
  combinedTriangles.block(tetmeshTriangles.rows(), 0,
	  splittingPlaneTriangles.rows(), 3).array() += currentBulletVertexPositions.rows();
  
  assert(transformedVertices.allFinite());

  RMMatrix3f colors(tetmeshVertices.rows() + splittingPlaneVertices.rows(), 3);
  colors.block(0, 0, tetmeshVertices.rows(), 3) = 
	tetmeshVertices.template cast<float>();
  colors.block(tetmeshVertices.rows(), 0,
	  splittingPlaneVertices.rows(), 3) =
	splittingPlaneVertices.template cast<float>();

  //normalize so it's in a nice RGB space
  colors.array() -= colors.minCoeff();
  colors.array() /= colors.maxCoeff();

  assert(fabs(colors.minCoeff()) < 1e-5);
  assert(fabs(1 - colors.maxCoeff()) < 1e-5);
  assert(colors.allFinite());

  std::ofstream plyStream(filename);
  writePLYWithColor(plyStream, transformedVertices, 
	  colors,
	  combinedTriangles);
  //writePLY(plyStream, splittingPlaneVertices, splittingPlaneTriangles);
  plyStream.close();
  
}


void PlasticPiece::dumpBcc(const std::string& filename) const{
  Eigen::MatrixXd paddedBarycentricCoordinates =
	Eigen::MatrixXd::Zero(barycentricCoordinates.rows() +
		splittingPlaneVertices.rows(),
		barycentricCoordinates.cols());

  paddedBarycentricCoordinates.block(0,0,
	  barycentricCoordinates.rows(), barycentricCoordinates.cols()) = 
	barycentricCoordinates;
  paddedBarycentricCoordinates.block(barycentricCoordinates.rows(), 0, 
	  splittingPlaneVertices.rows(), 1).array() = 1;

  writeMatrixBinary(filename, paddedBarycentricCoordinates);
  
}

void PlasticPiece::updateAabbs(){

  btTransform id;
  id.setIdentity();
  if(useVolumetricCollisions){
	for(auto&& pr : enumerate(bulletTets)){
	  auto& tet = pr.second;
	  tet.m_numVertices = 4;
	  tet.recalcLocalAabb();
	  bulletShape->updateChildTransform(pr.first, id, false);
	}
  } else {
	//update bulkMesh's info
	bulletShape->updateChildTransform( 0, id, false);
  }
  if(fractureMesh){
	bulletShape->updateChildTransform(
		bulletShape->getNumChildShapes() -1, id, false);
  }
  //auto* dbvt = bulletShape->getDynamicAabbTree();
  //dbvt->optimizeIncremental(4); //??? how many passes?
  bulletShape->recalculateLocalAabb();
  
}



void PlasticPiece::initializeNoDynamics(const std::string& directory, 
	const PlasticBody& parent,
	const RMMatrix3d& verticesIn){
  
  framesToSkin = 10000;
  scaleFactor = 1;
  useVolumetricCollisions = false;
  tetmeshVertices = verticesIn;
  currentBulletVertexPositions = tetmeshVertices;
  numPhysicsVertices = verticesIn.rows();
  computeVertexNeighbors();
  computeActiveVertices();
  barycentricCoordinates.resize(numPhysicsVertices, parent.numNodes);
  barycentricCoordinates.col(0).setOnes();
  
  deltaBarycentricCoordinates.resize(numPhysicsVertices, parent.numNodes);
  deltaBarycentricCoordinates.setZero();

  computeCollisionTetIndices();
  tetmeshTriangles = computeTriangleFaces(collisionTetIndices);

  
  skinMeshVaryingBarycentricCoords(
	  parent.boneWeights,
	  parent.boneIndices,
	  parent.exampleGraph);
  

}

void PlasticPiece::initialize(const std::string& directory,
	const PlasticBody& parent,
	const RMMatrix3d& verticesIn, 
	World& world,
	int pieceNumber){

  framesToSkin = 5; //start skinning at the beginning of the sim

  scaleFactor = parent.scaleFactor;
  tetmeshVertices = verticesIn;
  currentBulletVertexPositions = 
	scaleFactor*tetmeshVertices;

  useVolumetricCollisions = parent.useVolumetricCollisions;

  //assert(piece.currentBulletVertexPositions.allFinite());

  //read tets in this piece
  {
	
	std::ifstream ins(
		directory + "/tetmesh." + std::to_string(pieceNumber) + ".txt");
	if(ins.good()){
	  int nTets;
	  ins >> nTets;
	  tetmeshTets.resize(nTets, 4);
	  std::copy(std::istream_iterator<int>(ins),
		  std::istream_iterator<int>(),
		  tetmeshTets.data());
	} //else the parent already set this to all the tets
	
  }
  {
	std::ifstream ins(
		directory + "/cutTets." + std::to_string(pieceNumber) + ".txt");
	if(!ins.good()){
	  std::cout << "couldn't open " 
				<< directory + "/cutTets." + std::to_string(pieceNumber) + ".txt" << std::endl;
	} else{
	  int nCutTets;
	  ins >> nCutTets;
	  std::cout << "n cut tets: " << nCutTets << std::endl;
	  cutTets.resize(nCutTets);
	  std::copy(std::istream_iterator<size_t>(ins),
		  std::istream_iterator<size_t>(),
		  cutTets.begin());
	  std::cout << "back cut tet: " << cutTets.back() << std::endl;
	}
  }


  {
	std::ifstream ins(
		directory + "/planes.txt");
	if(ins.good()){

	  for(std::string line; std::getline(ins, line); ){
		Eigen::Vector3d normal;
		double offset;
		std::istringstream lineStream(line);
		lineStream >> normal.x() >> normal.y() >> normal.z() >> offset;
		//if(linestream.good()){
		//		  std::cout << "good: " << std::endl;
		  planes.emplace_back(normal, offset);
		  //}
	  }
	}
  }
  std::cout << "read " << planes.size() << "planes: \n";
  for(auto& plane : planes){
	std::cout << plane.normal << '\t' <<plane.offset << std::endl;
  }
  numPhysicsVertices = verticesIn.rows();
  computeCollisionTetIndices();
  tetmeshTriangles = computeTriangleFaces(collisionTetIndices); //and triangles
  computeCutTetGeometry();

  computeVertexNeighbors(); //connectivity in this piece
	
  computeActiveVertices();


  //setup rigid bodies
  bulletShape = std::unique_ptr<btCompoundShape>{
	new btCompoundShape{true}};
  auto idTransform = btTransform{};  
  if(useVolumetricCollisions){
  
	bulletTets.resize(collisionTetIndices.size());
	for(auto& tet : bulletTets){ tet.setMargin(0.001);}
	//actually fill this stuff in later...

	idTransform.setIdentity();
	for(auto i : range(bulletTets.size())){
	  bulletShape->addChildShape(idTransform, &(bulletTets[i]));
	}
	bulletShape->setMargin(0.001);
  } else {
	
	bulkTriMesh = std::unique_ptr<btTriangleIndexVertexArray>{
	  new btTriangleIndexVertexArray{
		static_cast<int>(tetmeshTriangles.rows()),
		tetmeshTriangles.data(),
		3*sizeof(decltype(tetmeshTriangles)::Scalar),
		static_cast<int>(currentBulletVertexPositions.rows()),
		currentBulletVertexPositions.data(),
		3*sizeof(decltype(currentBulletVertexPositions)::Scalar)
	  }};

	bulkMesh = std::unique_ptr<btGImpactMeshShape>{
	  new btGImpactMeshShape{bulkTriMesh.get()}};

	bulletShape->addChildShape(idTransform, bulkMesh.get());

  } 
  

  //load the plane ply
  {
	std::ifstream ins(
		directory + "/clippingTriangles." + std::to_string(pieceNumber) + ".ply");
	if(ins.good()){
	  RMMatrix3f splitVerticesIn;
	  RMMatrix3i splittingPlaneTrianglesIn;
	  readPLY(ins, splitVerticesIn, splittingPlaneTrianglesIn);
	  
	  //include the cut tet geometry also
	  splittingPlaneVertices.resize(
		  splitVerticesIn.rows() + cutGeomVertices.rows(), 3);
	  splittingPlaneTriangles.resize(
		  splittingPlaneTrianglesIn.rows() + cutGeomIndices.rows(), 3);
	  

	  splittingPlaneVertices.block(0,0, splitVerticesIn.rows(), 3) = 
		  splitVerticesIn.template cast<double>();
	  splittingPlaneVertices.block(splitVerticesIn.rows(), 0, 
		  cutGeomVertices.rows(), 3) = cutGeomVertices;

	  cutGeomIndices.array() += splitVerticesIn.rows(); //get the vertex offsets right
	  splittingPlaneTriangles.block(0,0,splittingPlaneTrianglesIn.rows(), 3) =
		splittingPlaneTrianglesIn;
	  splittingPlaneTriangles.block(splittingPlaneTrianglesIn.rows(), 0,
		  cutGeomIndices.rows(), 3) = cutGeomIndices;
	  

	  currentBulletSplittingPlaneVertices = 
		scaleFactor*splittingPlaneVertices;
	  assert(splittingPlaneVertices.allFinite());
	  assert(currentBulletSplittingPlaneVertices.allFinite());
	} else {
	  std::cout << "object has no clipping triangles" << std::endl;
	}
  }

  if(splittingPlaneVertices.rows() > 0){
	btTriMesh = std::unique_ptr<btTriangleIndexVertexArray>{
	  new btTriangleIndexVertexArray{
		static_cast<int>(splittingPlaneTriangles.rows()),
		splittingPlaneTriangles.data(),
		3*sizeof(decltype(splittingPlaneTriangles)::Scalar),
		static_cast<int>(currentBulletSplittingPlaneVertices.rows()),
		currentBulletSplittingPlaneVertices.data(),
		3*sizeof(decltype(currentBulletSplittingPlaneVertices)::Scalar)
	  }};
	
	fractureMesh = std::unique_ptr<btGImpactMeshShape>{
	  new btGImpactMeshShape{btTriMesh.get()}};
	
	bulletShape->addChildShape(idTransform, fractureMesh.get());
	
	computeClippingPlaneTetInfo();
	std::cout << "added fracture mesh" << std::endl;
  } 


  //initialize other stuff
  barycentricCoordinates.resize(numPhysicsVertices,
	  parent.numNodes);
  barycentricCoordinates.setZero();
  barycentricCoordinates.col(0).setOnes();
  
  deltaBarycentricCoordinates.resize(numPhysicsVertices,
	  parent.numNodes);
  deltaBarycentricCoordinates.setZero();
  
  skinMeshVaryingBarycentricCoords(
	  parent.boneWeights,
	  parent.boneIndices, 
	  parent.exampleGraph);
  



  motionState = 
	std::unique_ptr<btDefaultMotionState>{
	new btDefaultMotionState{}};
	
  mass = 1; //worked before...?
  bulletBody = std::unique_ptr<btRigidBody>{
	new btRigidBody{mass,
					motionState.get(),
					bulletShape.get()}};

  computeMassesAndVolume(parent.density);


  //put stuff on the GPU
  hostBarycentricCoordinates.resize(
	  barycentricCoordinates.rows()*barycentricCoordinates.cols());
  deviceBarycentricCoordinates = cl::Buffer(world.context,
	  CL_MEM_READ_ONLY, 
	  hostBarycentricCoordinates.size()*sizeof(float));
  
  hostSkinnedPositions.resize(3*numPhysicsVertices);
  deviceSkinnedPositions = cl::Buffer(world.context,
	  CL_MEM_READ_WRITE,
	  hostSkinnedPositions.size()*sizeof(float));


  //const auto numRealBones = parent.boneWeights.cols();
  //hostPerVertexTranslations.resize(3*numRealBones*numPhysicsVertices);
  //devicePerVertexTranslations = cl::Buffer(world.context,
  //  CL_MEM_WRITE_ONLY,
  //  hostPerVertexTranslations.size()*sizeof(float));
  //hostPerVertexRotations.resize(4*numRealBones*numPhysicsVertices);
  //devicePerVertexRotations = cl::Buffer(world.context,
  //  CL_MEM_WRITE_ONLY,
  //  hostPerVertexRotations.size()*sizeof(float));




}

void PlasticPiece::computeCollisionTetIndices(){
  collisionTetIndices.resize(tetmeshTets.rows() - cutTets.size());
  auto rangeNum = range(tetmeshTets.rows());
  //save the ones that aren't cut
  std::copy_if(rangeNum.begin(), rangeNum.end(),
	  collisionTetIndices.begin(),
	  [this](size_t a){
		return 
		  !std::binary_search(cutTets.begin(), cutTets.end(), a);
	  });
}


void PlasticPiece::computeClippingPlaneTetInfo(){
  clippingPlaneTetInfo.resize(splittingPlaneVertices.rows());
  Eigen::Matrix<double,3,4> mat;
  for(auto i : range(splittingPlaneVertices.rows())){

	clippingPlaneTetInfo[i] = std::make_tuple(
		Eigen::Vector4i(-1,-1,-1,-1), Eigen::Vector4d(0,0,0,0));
	//std::cout << "checking vertex: " << i << std::endl;
	Eigen::Vector3d thisVertex = splittingPlaneVertices.row(i).transpose();
	//std::cout << "coords: " << thisVertex << std::endl;

	int bestTet;
	double bestError = std::numeric_limits<double>::infinity();
	Eigen::Vector3d bestBcs;

	for(auto j : range(tetmeshTets.rows())){//cutTets){
	  for(auto k : range(4)){
		mat.col(k) = tetmeshVertices.row(tetmeshTets(j, k)).transpose();
	  }
	  Eigen::Vector3d slop{0.1,0.1, 0.1};
	  Eigen::Vector3d lower = mat.rowwise().minCoeff() - slop;
	  Eigen::Vector3d upper = mat.rowwise().maxCoeff() + slop;



	  if((lower.array() <= thisVertex.array()).all() &&
		  (upper.array() >= thisVertex.array()).all()){
		//std::cout << "passed bb check" << std::endl;
		//std::cout << "tet bb: " << lower << upper << std::endl;
		Eigen::Matrix3d bMatrix;
		bMatrix.col(0)= mat.col(0) - mat.col(3);
		bMatrix.col(1)= mat.col(1) - mat.col(3);
		bMatrix.col(2)= mat.col(2) - mat.col(3);

		Eigen::Vector3d barycentricCoords = bMatrix.inverse()*
		  (thisVertex - mat.col(3));
		//std::cout << "barycentric: " << barycentricCoords << std::endl;
		if((barycentricCoords.array() > -0.05).all() &&
			barycentricCoords.sum() < 1.15){
		  //std::cout << "found containing tet " << i << std::endl;
		  //in this tet
		  clippingPlaneTetInfo[i] = std::make_tuple(
			  Eigen::Vector4i(tetmeshTets(j,0), tetmeshTets(j,1), 
				  tetmeshTets(j,2), tetmeshTets(j,3)),
			  Eigen::Vector4d(
				  barycentricCoords(0),
				  barycentricCoords(1),
				  barycentricCoords(2),
				  1.0 - barycentricCoords.sum()));
			  
		  break; //stop checking tets
		} else {
		  double thisError = 
			std::max( barycentricCoords.minCoeff() < 0 ? fabs(barycentricCoords.minCoeff()) : 0,
				barycentricCoords.sum() > 1 ? barycentricCoords.sum() -1 : 0);
		  if(thisError < bestError){
			bestError = thisError;
			bestTet = j;
			bestBcs = barycentricCoords;
		  }
		}
	  }
	}
	if(std::get<0>(clippingPlaneTetInfo[i]).minCoeff() < 0){
	  std::cout << "vertex " << i << " not actually in any tets, using closest with bcs: " << bestBcs 
				<< ", " << 1.0 - bestBcs.sum() << std::endl;
	  clippingPlaneTetInfo[i] = std::make_tuple(
		  Eigen::Vector4i(tetmeshTets(bestTet, 0),
			  tetmeshTets(bestTet, 1),
			  tetmeshTets(bestTet, 2),
			  tetmeshTets(bestTet, 3)),
		  Eigen::Vector4d(
			  bestBcs(0),
			  bestBcs(1),
			  bestBcs(2),
			  1.0 - bestBcs.sum()));
	}
  }
}

Eigen::Vector3d PlasticPiece::intersectPlane(size_t ind1, size_t ind2, 
	const PlasticPiece::Plane& plane) const{
  
  //std::cout << "v1: " << tetmeshVertices.row(ind1) << std::endl;
  //std::cout << "v2: " << tetmeshVertices.row(ind2) << std::endl;

  const double p1dn = tetmeshVertices.row(ind1).dot(plane.normal);
  const double p2dn = tetmeshVertices.row(ind2).dot(plane.normal);
  const auto t = (plane.offset - p2dn)/(p1dn - p2dn);

  //std::cout << t << std::endl;

  return t*tetmeshVertices.row(ind1) + 
	(1 - t)*tetmeshVertices.row(ind2);

}

void PlasticPiece::computeCutTetGeometry(){
  if(cutTets.empty()){
	return;
  }
  auto cutTetFaces = computeTriangleFaces(cutTets);
  /*std::ofstream plyOut("cutPieceGeom.ply");
  RMMatrix3f verticesFloat = tetmeshVertices.cast<float>();
  RMMatrix3i allTris = computeAllTriangles(cutTets);
  writePLY(plyOut, verticesFloat, allTris);
  */
  /*{
	RMMatrix3f tmp = tetmeshVertices.template cast<float>();
	std::ofstream outs("cutTets.ply");
	writePLY(outs, tmp, cutTetFaces);
  }
  */

  std::cout << "cut tet faces: " << cutTetFaces.rows() << std::endl;
  std::vector<double> cutGeomVerticesVec;
  std::vector<int> cutGeomIndicesVec;
  
  for(auto i : range(cutTetFaces.rows())){
	//find which plane cuts

	for(const auto& plane : planes){
	  
	  //vertices of any collisionTetIndices tets should are on the "right"
	  //side of the plane
	  bool correctSide = 
		tetmeshVertices.row(tetmeshTets(collisionTetIndices[0],0)).dot(
			plane.normal) > plane.offset;

	  std::array<bool,3> sides; //which side of the plane are things on
	  for(auto j : range(3)){
		sides[j] = tetmeshVertices.row(cutTetFaces(i, j)).dot(plane.normal) > 
		  plane.offset;
	  }
	  
	  if(sides[0] == sides[1] &&
		  sides[0] == sides[2]){
		
		//not cut, don't do anything
		
	  } else {
		//triangle is cut
		//assume it's OK to duplicate vertices

		/*std::cout << "cutting triangle with vertices:\n" 
				  << tetmeshVertices.row(cutTetFaces(i,0)) << '\n'
				  << tetmeshVertices.row(cutTetFaces(i,1)) << '\n'
				  << tetmeshVertices.row(cutTetFaces(i,2)) << std::endl;;
		*/
		auto numCorrect = std::count(sides.begin(), sides.end(), correctSide);
		if(numCorrect == 1){
		  //std::cout << "one new triangle: " << std::endl;

		  //tessalate into a single triangle
		  auto whichCorrect = (sides[0] == correctSide) ? 0 : 
			((sides[1] == correctSide) ? 1 : 2);
		  

		  Eigen::Vector3d i1 = 
			intersectPlane(cutTetFaces(i, whichCorrect),
				cutTetFaces(i, (whichCorrect +1)%3),
				plane);
		  Eigen::Vector3d i2 =
			intersectPlane(cutTetFaces(i, whichCorrect),
				cutTetFaces(i, (whichCorrect +2)%3),
				plane);
		  
		  //std::cout << "new vertices: " << i1 << '\n' << i2 << std::endl;


		  int startIndex = cutGeomVerticesVec.size()/3;
		  cutGeomVerticesVec.insert(cutGeomVerticesVec.end(),
			  tetmeshVertices.row(cutTetFaces(i, whichCorrect)).data(),
			  tetmeshVertices.row(cutTetFaces(i, whichCorrect)).data() + 3);
		  cutGeomVerticesVec.insert(cutGeomVerticesVec.end(),
			  i1.data(), i1.data() +3);
		  cutGeomVerticesVec.insert(cutGeomVerticesVec.end(),
			  i2.data(), i2.data() +3);
		  
		  cutGeomIndicesVec.insert(cutGeomIndicesVec.end(),
			  {startIndex, startIndex +1, startIndex +2});



		} else if(numCorrect == 2){
		  //std::cout << "two triangles: " << std::endl;
		  //tesselate into two triangles
		  auto whichIncorrect = (sides[0] != correctSide) ? 0 : 
			( (sides[1] != correctSide) ? 1 : 2);
		  

		  Eigen::Vector3d i1 = 
			intersectPlane(
				cutTetFaces(i, whichIncorrect),
				cutTetFaces(i, (whichIncorrect +1)%3),
				plane);
		  Eigen::Vector3d i2 = 
			intersectPlane(
				cutTetFaces(i, whichIncorrect),
				cutTetFaces(i, (whichIncorrect +2)%3),
				plane);

		  //std::cout << "new vertices: " << i1 << '\n' << i2 << std::endl;

		  //Quadrangulating is trickier than I thought...
		  //Approach here: 
		  //t1 is original 2 vertices plus whichever
		  //i1 or i2 projects closer to original1
		  //t2 is new points and original2
		  Eigen::Vector3d o1 = 
			tetmeshVertices.row(cutTetFaces(i, (whichIncorrect +1)%3));
		  Eigen::Vector3d o2 = 
			tetmeshVertices.row(cutTetFaces(i, (whichIncorrect +2)%3));

		  Eigen::Vector3d existingEdge = 
			o2 - o1;
		  double i1dot = existingEdge.dot(i1 - o1);
		  double i2dot = existingEdge.dot(i2 - o1);
		  
		  int startIndex = cutGeomVerticesVec.size()/3;
		  cutGeomVerticesVec.insert(cutGeomVerticesVec.end(),
			  o1.data(), o1.data() +3);
		  cutGeomVerticesVec.insert(cutGeomVerticesVec.end(),
			  o2.data(), o2.data() +3);
		  cutGeomVerticesVec.insert(cutGeomVerticesVec.end(),
			  i1.data(), i1.data() +3);
		  cutGeomVerticesVec.insert(cutGeomVerticesVec.end(),
			  i2.data(), i2.data() +3);

		  if(i1dot < i2dot){
			//t1 is o1, o2, i1
			cutGeomIndicesVec.insert(cutGeomIndicesVec.end(),
				{startIndex, startIndex + 1, startIndex + 2});
			
		  } else {
			//t1 is o1, o2, i2
			cutGeomIndicesVec.insert(cutGeomIndicesVec.end(),
				{startIndex, startIndex + 1, startIndex + 3});
		  }
		  //t2 is always i1, i2, o2
		  cutGeomIndicesVec.insert(cutGeomIndicesVec.end(),
			  {startIndex + 1, startIndex + 2, startIndex + 3});
		} else {
		  //none or 3?  Something's wrong
		  assert(false);
		}
	  }
	}
  }
  cutGeomVertices.resize(cutGeomVerticesVec.size()/3, 3);
  cutGeomIndices.resize(cutGeomIndicesVec.size()/3, 3);
  std::copy(cutGeomVerticesVec.begin(), cutGeomVerticesVec.end(),
	  cutGeomVertices.data());
  std::copy(cutGeomIndicesVec.begin(), cutGeomIndicesVec.end(),
	  cutGeomIndices.data());
  std::cout << "new border triangles: " << cutGeomIndices.rows() << std::endl;
  /*  {
	std::ofstream outs("brandNewTriangles.ply");
	RMMatrix3f tmp = cutGeomVertices.template cast<float>();
	writePLY(outs, tmp, cutGeomIndices);
  }
  */
  //exit(1);
}


void PlasticPiece::skinMeshOpenCL(World& world, PlasticBody& parent, cl::Kernel& clKernel){
  
  const auto numRealBones = parent.boneWeights.cols();
  const auto numBoneTips = parent.boneIndices.size();
  const auto numNodes = parent.exampleGraph.nodes.size();

  /*  std::cout << "skinned before" << std::endl;
  for(auto i : range(10)){
	std::cout << currentBulletVertexPositions.data()[i] << ' ';
  }
  std::cout << std::endl;
*/

  std::copy(barycentricCoordinates.data(), 
	  barycentricCoordinates.data() + hostBarycentricCoordinates.size(),
	  hostBarycentricCoordinates.begin());
  cl::copy(world.queue, hostBarycentricCoordinates.begin(), hostBarycentricCoordinates.end(),
	  deviceBarycentricCoordinates);


  auto kernelInvoc = cl::make_kernel<int, int, int, int, float,
									 cl::Buffer&,cl::Buffer&,cl::Buffer&,
									 cl::Buffer&,cl::Buffer&,cl::Buffer&,
									 //cl::Buffer&,cl::Buffer&,
									 cl::Buffer&>(clKernel);
  
  
  cl::NDRange clRange(numPhysicsVertices);
  auto event = kernelInvoc(cl::EnqueueArgs(world.queue, clRange),
	  numRealBones, numBoneTips, numNodes,
	  numPhysicsVertices, 
	  parent.scaleFactor,
	  parent.deviceUnskinnedPositions,
	  parent.deviceTranslations,
	  parent.deviceRotations,
	  deviceBarycentricCoordinates,
	  parent.deviceBoneIndices,
	  parent.deviceBoneWeights,
	  deviceSkinnedPositions//,
	  //devicePerVertexTranslations,
	  //devicePerVertexRotations
	);

  //  std::cout << "waiting on CL" << std::endl;
  event.wait(); //skinning done now
  

  cl::copy(world.queue, deviceSkinnedPositions, hostSkinnedPositions.begin(), hostSkinnedPositions.end());
  std::copy(hostSkinnedPositions.begin(), hostSkinnedPositions.end(),
	  currentBulletVertexPositions.data());

  //skin the cut vertices
  for(auto i : range(splittingPlaneVertices.rows())){
	auto& inds = std::get<0>(clippingPlaneTetInfo[i]);
	auto& bcs = std::get<1>(clippingPlaneTetInfo[i]);
	currentBulletSplittingPlaneVertices.row(i) =
	  bcs(0)*currentBulletVertexPositions.row(inds(0)) +
	  bcs(1)*currentBulletVertexPositions.row(inds(1)) +
	  bcs(2)*currentBulletVertexPositions.row(inds(2)) +
	  bcs(3)*currentBulletVertexPositions.row(inds(3));
  }


  /*
  cl::copy(world.queue, devicePerVertexTranslations, 
	  hostPerVertexTranslations.begin(), hostPerVertexTranslations.end());
  std::copy(hostPerVertexTranslations.begin(), hostPerVertexTranslations.end(),
	  perVertexTranslations.front().data());

  cl::copy(world.queue, devicePerVertexRotations,
	  hostPerVertexRotations.begin(), hostPerVertexRotations.end());
  std::copy(hostPerVertexRotations.begin(), hostPerVertexRotations.end(),
	  perVertexRotations.front().coeffs().data());
  */

  /*std::cout << "done with CL" << std::endl;

  std::cout << "skinned after" << std::endl;
  for(auto i : range(10)){
	std::cout << currentBulletVertexPositions.data()[i] << ' ';
  }
  std::cout << std::endl;
  */

  //should be a no-op, but could cause issues, i guess?
  //currentBulletVertexPositions.resize(numPhysicsVertices, 3);
  //currentBulletVertexPositions.setZero();
  
  /*
  //todo just use activeVertices here...
  for(auto i : range(numPhysicsVertices)){ 
	for (auto j : range(numBoneTips)) {
	  if (boneIndices[j] < 0) { continue; } //bones[j]->get_wi() < 0) //not a real bone
	    
	  const auto kRange = range(numNodes);
	  //compute bone transform
	  Eigen::Vector3d translation = Eigen::Vector3d::Zero().eval(); //TODO: Ron - accumulate
	  for (auto k : kRange) {
		translation += barycentricCoordinates(i, k) * exampleGraph.nodes[k].transformations[j].translation;
	  }
	    
	  Eigen::Vector4d rotationCoeffs = Eigen::Vector4d::Zero().eval(); //TODO: Ron - accumulate
	  for (auto k : kRange){
		rotationCoeffs += barycentricCoordinates(i, k) * 
		  exampleGraph.nodes[k].transformations[j].rotation.coeffs();
	  }
	  Quat rotation(rotationCoeffs);
	  rotation.normalize();
	  
	  const auto weightIndex = boneIndices[j];
	    
	  currentBulletVertexPositions.row(i) +=
		scaleFactor*boneWeights(i, weightIndex)*
		(rotation*(tetmeshVertices.row(i).transpose()) +
			translation).transpose();
	    
	  perVertexTranslations[i*numRealBones + weightIndex] = translation;
	  perVertexRotations[i*numRealBones + weightIndex] = rotation;
	    
	}
  }
  //assert(currentBulletVertexPositions.allFinite());
  
  //skin the clipping plane vertices
  for (auto i : range(splittingPlaneVertices.rows())) {
	auto& inds = std::get<0>(clippingPlaneTetInfo[i]);
	auto& bcs = std::get<1>(clippingPlaneTetInfo[i]);
	currentBulletSplittingPlaneVertices.row(i) =
	  bcs(0)*currentBulletVertexPositions.row(inds(0)) +
	  bcs(1)*currentBulletVertexPositions.row(inds(1)) +
	  bcs(2)*currentBulletVertexPositions.row(inds(2)) +
	  bcs(3)*currentBulletVertexPositions.row(inds(3));
  }
    */
}

