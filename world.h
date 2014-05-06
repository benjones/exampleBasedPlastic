#pragma once
#include <vector>
#include <list>
#include <BulletCollision/BroadphaseCollision/btBroadphaseInterface.h>
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
#include <BulletCollision/BroadphaseCollision/btDispatcher.h>
#include <BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolver.h>
#include <BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h>

#include "rigidBody.h"
#include "fem.H"
#include "plasticObject.h"

namespace Json{ class Value;}

class World{
  public:
  
  World(std::string filename);
  World(const World& other) = delete;  //no moving or copying
  World(World&& other) = delete;

  World& operator=(World&& other) = delete;
  World& operator=(const World& other) = delete;


  void loadPlasticObjects(const Json::Value& root);

  void computeConstraints();
  void countConstraints();

  void massScale(double *x, double *y);
  void inverseMassScale(double *x, double *y);
  void solve();
  void solveMinres();
  void project(double *in);
  void checkVector(double *in);
  void project(double *in, RigidBody &rb, int rbIndex, CouplingConstraint &c);
  void timeStep();
  void timeStepRigidCollisions();

  void timeStepDynamicSprites();
  void deformBasedOnImpulses();

  void dumpFrame();

  void computeFemVelocities();
  void updateFemPositions();

  void computeCrossProductMatrices();

  void mulSpringMatrix(double *in, double *out);
  void applySpringPreconditioner(double *in, double *out);

  void constraintCullingCallback(btBroadphasePair& collisionPair, 
								 btCollisionDispatcher& _dispatcher, 
								 const btDispatcherInfo& dispatcherInfo);

  //multiply the constaint matrix by in and put the result in out.
  //this is the lower left block of the system.
  //make sure that in and out is correctly allocated already
  void mulJV(double* in, double* out);

  //out += JT*In.  This is the upper right block of the system
  void mulJTLambda(double* in, double* out);

  //set out to be the rigid mass/inertia matrix times in
  //out = M_Rigid*In
  void mulRigidMassMatrix(double* in, double* out);

  //out = M^-1 I^-1 * in
  void applyRigidPreconditioner(double* in, double* out);

  //out += regularizerAlpha*in;
  void applyRegularizer(double* in, double* out);
  //out += 1/regularizerAlpha*in
  void applyRegularizerPreconditioner(double* in, double* out);

  std::unique_ptr<btBroadphaseInterface> broadphaseInterface;  
  std::unique_ptr<btCollisionConfiguration> collisionConfiguration;
  std::unique_ptr<btCollisionDispatcher> dispatcher;

  //maybe using the MLCP solver makes sense, since it's supposed to give higher fidelity results
  std::unique_ptr<btSequentialImpulseConstraintSolver> bulletSolver;

  btDiscreteDynamicsWorld bulletWorld;

  bool ground;
  double groundHeight;
  double dt;
  btVector3 gravity;
  double friction;
  
  double regularizerAlpha;
  size_t totalNumConstraints;

  std::vector<RigidBody> rigidBodies;
  std::vector<FemObject> femObjects;
  std::vector<PlasticObject> plasticObjects;

  unsigned int nrbdof, nfemdof;

  int currentFrame;
 
	bool RIGIDS;
	bool FEMS;
};
