#pragma once

#include "eigenTypedefs.h"
#include <vector>
#include "Quat.h"

#include <LinearMath/btTransform.h>
#include <LinearMath/btDefaultMotionState.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>


class ExampleGraph;

class PlasticPiece{
public:

  //move only type because of the unique_ptrs
  //also, we should never actually need to copy them

  
  PlasticPiece(){}
  ~PlasticPiece() = default;
  PlasticPiece(PlasticPiece&&) = default;
  PlasticPiece& operator=(PlasticPiece&&) = default;
  PlasticPiece(const PlasticPiece&) = delete;
  PlasticPiece& operator=(const PlasticPiece&) = delete;
  
  void computeMassesAndVolume(double density);

  void updateBulletProperties();

  int getNearestVertex(const Eigen::Vector3d& localPoint) const;

  void saveBulletSnapshot();
  void restoreBulletSnapshot();

  void skinMeshVaryingBarycentricCoords(const Eigen::MatrixXd& boneWeights,
										const std::vector<int>& boneIndices,
										const ExampleGraph& exampleGraph);


  void computeTriangleFaces();

  //neighbors within this piece.  
  //We can just union with vertices
  //in other pieces starting at vertices with the same index
  void computeVertexNeighbors();

  //propogate out from the seed points
  //there can be more than one if the real seed is in a different piece
  void geodesicDistancesPropogate(const std::vector<size_t>& seeds,
	double radius);
  

  struct BulletSnapshot{
	btVector3 linearVelocity, angularMomentum;
	btTransform comTransform;
  };
  BulletSnapshot bulletSnapshot;
  
  
  Eigen::MatrixXd barycentricCoordinates; //per vertex
  Eigen::MatrixXd deltaBarycentricCoordinates;


  std::vector<Eigen::Vector3d> perVertexTranslations; 
  std::vector<Quat> perVertexRotations;
  //should call world.removeRigidBody....
  std::unique_ptr<btRigidBody> bulletBody;
  std::unique_ptr<btGImpactMeshShape> bulletShape;
  std::unique_ptr<btDefaultMotionState> motionState;
  std::unique_ptr<btTriangleIndexVertexArray> btTriMesh; //only used for trimesh shapes

  RMMatrix3d tetmeshVertices;
  RMMatrix4i tetmeshTets;
  RMMatrix3i tetmeshTriangles; 
  Eigen::VectorXd tetmeshVertexMasses;  //for computing inertia tensor, etc

  int numPhysicsVertices; //include the ones that aren't actually part of any
  //tets in this piece

  RMMatrix3d currentBulletVertexPositions;

  double mass, volume, scaleFactor;

  Eigen::Matrix3d inertiaTensor;
  btTransform inertiaAligningTransform, worldTransform;

  std::vector<std::vector<size_t>> vertexNeighbors;
  std::vector<std::vector<double>> neighborDistances;
  std::vector<double> geodesicDistances; //to a given point

  void computeActiveVertices();
  std::vector<size_t> activeVertices; //vertices that are actually part of tets
  
  void dumpPly(const std::string& filename) const;
  //void dumpImpulses(const std::string& filename) const; //skip for now
  void dumpBcc(const std::string& filename) const;

};
