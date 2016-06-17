
#include "plasticPieceSpheres.h"
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


void PlasticPieceSpheres::updateBulletProperties(){

  //do transform bullshit
  btTransform principal;
  btVector3 momentOfInertia;
  bulletShape->calculatePrincipalAxisTransform(
	  sphereMasses.data(), principal, momentOfInertia);

  for(auto i : range(bulletShape->getNumChildShapes())){
	btTransform newTransform =
	  principal.inverse()*bulletShape->getChildTransform(i);
	bulletShape->updateChildTransform(i, newTransform);
  }

  bulletBody->setCenterOfMassTransform(bulletBody->getCenterOfMassTransform()*principal);
  
  bulletBody->setMassProps(mass, momentOfInertia);

  bulletBody->updateInertiaTensor(); //rotate it
  bulletBody->setDamping(0, 0.01);  
  bulletBody->setActivationState(DISABLE_DEACTIVATION);
  
  


}


void PlasticPieceSpheres::computeVertexNeighbors(){

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


void PlasticPieceSpheres::geodesicDistancesPropogate(const std::vector<size_t>& seeds, double radius){

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

void PlasticPieceSpheres::dumpBcc(const std::string& filename) const{
  writeMatrixBinary(filename, barycentricCoordinates);
  
}

void PlasticPieceSpheres::dumpSpheres(const std::string& filename) const {
  std::ofstream outs(filename);
  auto numSpheres = bulletShape->getNumChildShapes();
  auto worldTransform = bulletBody->getCenterOfMassTransform();
  outs << numSpheres << std::endl;
  for(auto i : range(numSpheres)){
	auto position = (worldTransform*
		bulletShape->getChildTransform(i)).getOrigin();
	outs << position.x() << ' ' << position.y() << ' ' << position.z() << ' '
		 << bulletSpheres[i].getRadius() << std::endl;
	
  }
}





void PlasticPieceSpheres::initialize(const std::string& directory,
	const PlasticBody& parent,
	const RMMatrix3d& verticesIn, 
	World& world,
	int pieceNumber){

  framesToSkin = 5; //start skinning at the beginning of the sim

  scaleFactor = parent.scaleFactor;
  tetmeshVertices = verticesIn;

  numPhysicsVertices = verticesIn.rows();

  computeVertexNeighbors(); //connectivity in this piece
	
  //setup rigid bodies
  bulletShape = std::unique_ptr<btCompoundShape>{
	new btCompoundShape{true}};
  auto transform = btTransform{};

  //load spheres
  {
	std::ifstream spheresIn(directory + "/spheres.sph");
	size_t nSpheres;
	spheresIn >> nSpheres;
	bulletSpheres.reserve(nSpheres);
	sphereMasses.reserve(nSpheres);

	double x, y, z, r;
	double pi4_3 = M_PI*4/3;
	while(bulletSpheres.size() < nSpheres){
	  spheresIn >> x >> y >> z >> r;

	  bulletSpheres.emplace_back(scaleFactor*r);
	  sphereMasses.push_back(pi4_3*pow(scaleFactor*r, 3));
	  
	  transform.setOrigin(
		  btVector3{scaleFactor*x, scaleFactor*y, scaleFactor*z});
	  bulletShape->addChildShape(transform, &(bulletSpheres.back()));
	}

  }
  

  //initialize other stuff
  barycentricCoordinates.resize(numPhysicsVertices,
	  parent.numNodes);
  barycentricCoordinates.setZero();
  barycentricCoordinates.col(0).setOnes();
  
  deltaBarycentricCoordinates.resize(numPhysicsVertices,
	  parent.numNodes);
  deltaBarycentricCoordinates.setZero();
  


  motionState = 
	std::unique_ptr<btDefaultMotionState>{
	new btDefaultMotionState{}};
	
  mass = 1; //worked before...?
  bulletBody = std::unique_ptr<btRigidBody>{
	new btRigidBody{mass,
					motionState.get(),
					bulletShape.get()}};

  computeMassesAndVolume(parent.density);
  //set moment of inertia, etc
  updateBulletProperties();



}

