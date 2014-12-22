#pragma once

#include "eigenTypedefs.h"
#include <vector>
#include "Quat.h"

#include <LinearMath/btTransform.h>
#include <LinearMath/btDefaultMotionState.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <BulletCollision/CollisionShapes/btCompoundShape.h>
#include <BulletCollision/CollisionShapes/btTetrahedronShape.h>

using TetrahedronShape = btBU_Simplex1to4;


class ExampleGraph;
class PlasticBody;

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
  
  void initialize(const std::string& directory,
	  const PlasticBody& parent,
	  const RMMatrix3d& verticesIn,
	  int pieceNumber);


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

  std::unique_ptr<btCompoundShape> bulletShape;
  std::vector<TetrahedronShape> bulletTets;
  std::unique_ptr<btDefaultMotionState> motionState;

  //for the planar fracture mesh
  std::unique_ptr<btTriangleIndexVertexArray> btTriMesh; //only used for trimesh shapes
  std::unique_ptr<btGImpactMeshShape> fractureMesh;
  


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
  void computeCollisionTetIndices();
  std::vector<size_t> cutTets; //read from slicer output file
  std::vector<size_t> collisionTetIndices; 

  //indices of the vertices in the containing tet,
  //and the corresponding barycentric coordinates
  void computeClippingPlaneTetInfo();
  std::vector<std::tuple<Eigen::Vector4i,Eigen::Vector4d>> clippingPlaneTetInfo;
  
  

  RMMatrix3d splittingPlaneVertices;
  RMMatrix3d currentBulletSplittingPlaneVertices;
  RMMatrix3i splittingPlaneTriangles;


  
  void dumpPly(const std::string& filename) const;
  //void dumpImpulses(const std::string& filename) const; //skip for now
  void dumpBcc(const std::string& filename) const;

  void updateAabbs();

};
