
#include "plasticPiece.h"
#include "range.hpp"
#include "enumerate.hpp"
using benlib::range;
using benlib::enumerate;

#include "utils.h"
#include "plyio.hpp"
#include "exampleGraph.h"

#include <tbb/tbb.h>

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

  //update the tet vertices
  for(auto&& pr : enumerate(bulletTets)){
	for(auto i : range(4)){
	  Eigen::Vector3d ePos = currentBulletVertexPositions.row(tetmeshTets(pr.first, i));
	  pr.second.m_vertices[i] = eigenToBullet(ePos);
	}
  }
  updateAabbs();
  
  //  bulletShape->postUpdate();
  //  bulletShape->updateBound();
  bulletBody->setActivationState(DISABLE_DEACTIVATION);
  
  inertiaAligningTransform = btTransform{eigenToBullet(evecs),
  										 eigenToBullet(centerOfMass)};

  bulletBody->setCenterOfMassTransform(worldTransform*
									   inertiaAligningTransform);

  bulletBody->setMassProps(mass, eigenToBullet(evals));
  bulletBody->updateInertiaTensor(); //rotate it
}

int PlasticPiece::getNearestVertex(const Eigen::Vector3d& localPoint) const {

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
  
  const auto numRealBones = boneWeights.cols();
  const auto numBoneTips = boneIndices.size();
  const auto numNodes = exampleGraph.nodes.size();

  perVertexTranslations.resize(numPhysicsVertices*numRealBones);
  perVertexRotations.resize(numPhysicsVertices*numRealBones);

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
		  
		  const auto kRange = range(numNodes);
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
		  Quat rotation(rotationCoeffs);
		  rotation.normalize();
		  
		  const auto weightIndex = boneIndices[j];//bones[j]->get_wi();
		  
		  currentBulletVertexPositions.row(i) +=
			scaleFactor*boneWeights(i,weightIndex)*
			(rotation*(tetmeshVertices.row(i).transpose()) + 
			  translation).transpose();
		  
		  perVertexTranslations[i*numRealBones + weightIndex] = translation;
		  perVertexRotations[i*numRealBones + weightIndex] = rotation;
		  
		}
	  }
	});
  //assert(currentBulletVertexPositions.allFinite());
}

void PlasticPiece::computeTriangleFaces(){
  
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
  
  for(auto i : range(tetmeshTets.rows())){
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
  
  std::cout << "mesh file tris: " << tetmeshTriangles.rows() << std::endl;
  std::cout << "computed: " << seenFaces.size() << std::endl;
  
  //orientation is ignored by bullet, so don't reorient stuff
  tetmeshTriangles.resize(seenFaces.size(), 3);

  for(auto pr : enumerate(seenFaces)){
	tetmeshTriangles(pr.first, 0) = std::get<0>(pr.second);
	tetmeshTriangles(pr.first, 1) = std::get<1>(pr.second);
	tetmeshTriangles(pr.first, 2) = std::get<2>(pr.second);
  }
  
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
  RMMatrix3f transformedVertices(currentBulletVertexPositions.rows(), 3);
  for(auto i : range(transformedVertices.rows())){
	Eigen::Vector3d transformedVertex = 
	  eigenRot*currentBulletVertexPositions.row(i).transpose() +
	  eigenTrans;
	Eigen::Vector3f transFloat = transformedVertex.template cast<float>();
	transformedVertices.row(i) = transFloat.transpose();
	
  }
  
  std::ofstream plyStream(filename);
  writePLY(plyStream, transformedVertices, tetmeshTriangles);
  plyStream.close();
  
}


void PlasticPiece::dumpBcc(const std::string& filename) const{
  writeMatrixBinary(filename, barycentricCoordinates);
}

void PlasticPiece::updateAabbs(){

  btTransform id;
  id.setIdentity();
  for(auto&& pr : enumerate(bulletTets)){
	auto& tet = pr.second;
	tet.m_numVertices = 4;
	tet.recalcLocalAabb();
	bulletShape->updateChildTransform(pr.first, id, false);
  }
  //auto* dbvt = bulletShape->getDynamicAabbTree();
  //dbvt->optimizeIncremental(4); //??? how many passes?
  bulletShape->recalculateLocalAabb();
  
}
