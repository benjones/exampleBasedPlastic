#pragma once


#include "plasticPiece.h"
#include "exampleGraph.h"

#include <BulletDynamics/ConstraintSolver/btPoint2PointConstraint.h>


/*

  holds many fracturable plastic pieces, responsible for spreading
  deformation and stuff to them.

 */
namespace Json{ class Value;}

class World;
class PlasticBody{

public:
  
  PlasticBody(){}
  ~PlasticBody() = default;
  PlasticBody(PlasticBody&&) = default;
  PlasticBody& operator=(PlasticBody&& other) = default;
  PlasticBody(const PlasticBody& other) = delete;
  PlasticBody& operator=(const PlasticBody& other) = delete;


  void loadFromJson(const Json::Value& poi,
	  World& world,
	int objectIndex);
  
  void loadFromJsonNoDynamics(const Json::Value& poi);

  void projectImpulsesOntoExampleManifoldLocally(double dt);
  Eigen::VectorXd projectSingleImpulse( 
	  const PlasticPiece& piece, const Eigen::Vector3d& impulseAtContact, int vInd) const;
  Eigen::Vector3d getDeformationVectorFromImpulse(
	  const PlasticPiece& piece,
	  const btManifoldPoint& manPoint,
	  double dt,
	  bool isObject0) const;
  
  std::vector<PlasticPiece> plasticPieces;
  
  //return objectStart + numPiecesWritten
  int dump(int currentFrame, int objectStart) const;
  
  void saveBulletSnapshots();
  void restoreBulletSnapshots();
  
  void computeBoneIndices(const std::string& boneFile);
  std::vector<int> boneIndices;
  
  //geodesic distance stuff
  void computeGeodesicDistances(size_t pieceIndex, size_t vInd, double radius);

  ExampleGraph exampleGraph;
  Eigen::MatrixXd boneWeights;
  Eigen::VectorXd perExampleScale;
  int numNodes, numBoneTips, numRealBones, numPhysicsVertices;

  double density, plasticityImpulseYield, plasticityImpulseScale,
	plasticityKernelScale, plasticityRate, jacobianAlpha, scaleFactor,
	restitution, mass, breakingThreshold, restitutionAddition, restitutionExponent;

  bool useVolumetricCollisions;

  bool hasConstantVelocity = false;
  int constantVelocityFrames;
  btVector3 constantVelocity;
  
  std::vector<std::tuple<const btRigidBody*, btManifoldPoint, bool>> manifoldPoints;


  //                            p1      p2      vInd   constraint
  using Constraint = 
	std::tuple<size_t, size_t, size_t, 
			   std::unique_ptr<btPoint2PointConstraint>>;
  //should call world.removeConstraint in deleter...

  void computeConstraints();
  void updateConstraints();
  std::vector<Constraint> constraints;


  void skinAndUpdate();

  void skinAndUpdateCL(World& world, cl::Kernel& clKernel);

  //CL stuff
  std::vector<float> hostTranslations;
  std::vector<float> hostRotations;
  std::vector<float> hostBoneWeights;
  std::vector<float> hostUnskinnedPositions;

  cl::Buffer deviceTranslations, deviceRotations, 
	deviceBoneWeights, deviceBoneIndices, deviceUnskinnedPositions;

};
