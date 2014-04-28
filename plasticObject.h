#pragma once

#include <memory>
#include <exampleGraph.h>
#include <Skeleton.cpp> //templates...
#include <Bone.h>

#include <LinearMath/btDefaultMotionState.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <BulletCollision/CollisionShapes/btCompoundShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>


using RMMatrix3d = Eigen::Matrix<double,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix4i = Eigen::Matrix<int,   Eigen::Dynamic, 4, Eigen::RowMajor>;

class PlasticObject{
public:
  PlasticObject(){}
  
  void writeOutputFiles(std::string filename);
  
  void loadFromFiles(std::string directory);
  
  void computeMassesAndVolume();
  
  void updateBulletProperties(const RMMatrix3d& vertices,
							  const RMMatrix4i& tets);
  

  void dump(std::string filename);

  ExampleGraph exampleGraph;
  std::unique_ptr<btRigidBody> bulletBody;
  std::unique_ptr<btGImpactMeshShape> bulletShape;
  std::unique_ptr<btDefaultMotionState> motionState;
  
  
  std::unique_ptr<btTriangleIndexVertexArray> btTriMesh; //only used for trimesh shapes
  //std::vector<double> trimeshVertices;
  //std::vector<int> trimeshIndices;
  
  RMMatrix3d trimeshVertices, tetmeshVertices;
  RMMatrix3i trimeshTriangles, tetmeshTriangles;
  RMMatrix4i tetmeshTets;
  
  RMMatrix3d currentBulletVertexPositions;


  double density;
  double mass;
  double volume;
  Eigen::Matrix3d inertiaTensor;
  
  std::vector<double> tetmeshVertexMasses;
  
  Skeleton<Bone> skeleton;
  std::vector<Bone*> boneRoots;
  
  Eigen::MatrixXd boneWeights;
  
  
  
};
