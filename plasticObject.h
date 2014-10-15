#pragma once

#include <memory>
#include <exampleGraph.h>
#include "egTraverser.h"
#include <Skeleton.cpp> //templates...
#include <Bone.h>

#include <unsupported/Eigen/SparseExtra>

#include <LinearMath/btDefaultMotionState.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <BulletCollision/CollisionShapes/btCompoundShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h>



using RMMatrix3d = Eigen::Matrix<double,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3f = Eigen::Matrix<float ,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix4i = Eigen::Matrix<int,   Eigen::Dynamic, 4, Eigen::RowMajor>;

class PlasticObject{
public:


  struct BulletSnapshot{
	btVector3 linearVelocity, angularMomentum;
	btTransform comTransform;
  };
  BulletSnapshot bulletSnapshot;

  PlasticObject(){}
  ~PlasticObject() = default;
  PlasticObject(PlasticObject&&) = default;
  PlasticObject(const PlasticObject&) = delete;

  
  void loadFromFiles(std::string directory);
  
  void computeMassesAndVolume();
  
  void updateBulletProperties(const RMMatrix3d& vertices,
							  const RMMatrix4i& tets);
  

  void deformBasedOnImpulses(btPersistentManifold* man, bool isObject0);
  
  bool deformBasedOnImpulseLocal(btPersistentManifold* man, bool isObject0);
  
  bool projectImpulsesOntoExampleManifold();
  
  bool projectImpulsesOntoExampleManifoldLocally();


  btVector3 getDeformationVectorFromImpulse(const btManifoldPoint& manPoint,
											bool isObject0);
  
  btVector3 getLocalDeformationVectorFromImpulse(const btManifoldPoint& manPoint,
												 bool isObject0);
  


  Eigen::Vector3d getBarycentricCoordinates(btVector3 localPoint, int triangleIndex);
  int getNearestVertex(Eigen::Vector3d localPoint);
  

  void computeTransformationDerivatives(const EGPosition& egPosition, 
										const std::vector<int>& indices,
										Eigen::MatrixXd& deriv);
  

  void skinMesh(btDiscreteDynamicsWorld& bulletWorld);

  void dump(std::string filename);

  void dumpVerticesBinary(std::string filename) const;
  void dumpBarycentricCoords(std::string filename) const;


  void saveBulletSnapshot();
  void restoreBulletSnapshot();

  void skinMeshVaryingBarycentricCoords();


  ExampleGraph exampleGraph;
  EGTraverser egTraverser; //ignore current position stuff, but useful for shortest paths

  Eigen::MatrixXd barycentricCoordinates; //per vertex
  Eigen::MatrixXd deltaBarycentricCoordinates;
  RMMatrix3d desiredDeformations; //per vertex

  //nVertex x nBones, row-major matrix of the translation/rotation of each bone
  //store this so that we can use it when computing derivatives for the projection
  std::vector<Vec3> perVertexTranslations; 
  std::vector<Quat> perVertexRotations;

  std::unique_ptr<btRigidBody> bulletBody;
  std::unique_ptr<btGImpactMeshShape> bulletShape;
  
  std::unique_ptr<btDefaultMotionState> motionState;
  
  
  std::unique_ptr<btTriangleIndexVertexArray> btTriMesh; //only used for trimesh shapes
  //std::vector<double> trimeshVertices;
  //std::vector<int> trimeshIndices;
  
  RMMatrix3d trimeshVertices, tetmeshVertices;
  size_t numPhysicsVertices, numBoneTips, numRealBones, numNodes;
  RMMatrix3i trimeshTriangles, tetmeshTriangles;
  RMMatrix4i tetmeshTets;
  
  RMMatrix3d currentBulletVertexPositions;
  RMMatrix3d localImpulseBasedOffsets;

  double scaleFactor;//scale the object's size by this much
  bool deformedThisFrame = false;

  double dt;
  double density;
  double mass;
  double volume;
  Eigen::Matrix3d inertiaTensor;
  
  Eigen::VectorXd tetmeshVertexMasses;
  
  std::unique_ptr<Skeleton<Bone>> skeleton;
  std::vector<Bone*> bones;
  
  Eigen::MatrixXd boneWeights;
  Eigen::MatrixXd currentTransformMatrix;
  
  //an extra term to add to the translation of each handle based on
  //collision impacts
  Eigen::MatrixXd collisionHandleAdjustments;
  //only deform once things cross this point
  double plasticityImpulseYield;
  //scale*(mag(impulse) - yield) gets distributed to handles
  double plasticityImpulseScale;
  double plasticityKernelScale;
  
  double localPlasticityImpulseYield;
  double localPlasticityImpulseScale;


  bool hasConstantVelocity = false;
  btVector3 constantVelocity;

  //second tells us if po is object0 in the contact
  std::vector<std::pair<btManifoldPoint, bool>> manifoldPoints;

  Eigen::MatrixXd LBSMatrix;

  btTransform inertiaAligningTransform; //align the skinned position so that it's 
  //centered at its COM and aligned with the principal axes of inertia
  
  btTransform worldTransform; //the transform bullet gets is the product of these 2

  Eigen::SparseMatrix<double> cotangentLaplacian; //for smoothing barycentric weights

  void computeGeodesicDistances(size_t vInd, double radius);
  std::vector<double> geodesicDistances; //geodesic distances to a given point

  void computeVertexNeighbors();
  std::vector<std::vector<size_t>> vertexNeighbors;
  std::vector<std::vector<double>> neighborDistances;

};
