
#pragma once

#include <eigen3/Eigen/Eigen>

class World;
class RigidBody;

//solve for the impules necessary to keep FEM nodes attached to the rigid body
class CouplingConstraintSolver{
  public:
  
  //compute impulses to make the velocities at the FEM nodes and the rigid body surface match
  void solveVelocityConstraints(World& world, RigidBody& rigidBody);
  
  Eigen::MatrixXd constraintMatrix;
  Eigen::VectorXd rhs, impulseSolutionVector;

};
