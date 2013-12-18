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



class World{
  public:
  
  World(std::string filename);
  World(const World& other) = delete;  //move only, no copying
  World(World&& other) = default;

  void computeConstraints();
  void countConstraints();

  void timeStep();
  void dumpFrame();

  void computeFemVelocities();
  void updateFemPositions();

  void computeCrossProductMatrices();

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
  std::unique_ptr<btDispatcher> dispatcher;

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


  int currentFrame;

};
