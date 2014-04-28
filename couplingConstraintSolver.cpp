#include "couplingConstraintSolver.h"
#include "world.h"
#include "rigidBody.h"

#include "cppitertools/enumerate.hpp"

using iter::enumerate;




void CouplingConstraintSolver::solveVelocityConstraints(World& world, RigidBody& rigidBody){

  if(rigidBody.constraints.empty()){
    return;
  }
  constraintMatrix.resize(3*rigidBody.constraints.size(),
			  3*rigidBody.constraints.size());

  rhs.resize(3*rigidBody.constraints.size());

  
  Eigen::Matrix3d eigenInvInertiaTensor;
  auto& bit = rigidBody.bulletBody->getInvInertiaTensorWorld();
  eigenInvInertiaTensor << 
    bit[0][0], bit[0][1], bit[0][2],
    bit[1][0], bit[1][1], bit[1][2],
    bit[2][0], bit[2][1], bit[2][2];
    


  for(auto c1 : enumerate(rigidBody.constraints)){
    auto& fem1 = world.femObjects[c1.element.femIndex];
    //diagonal part for change in FEM node velocity due to impulse
    constraintMatrix(3*c1.index    , 3*c1.index    ) = fem1.massInv[c1.element.nodeIndex];
    constraintMatrix(3*c1.index + 1, 3*c1.index + 1) = fem1.massInv[c1.element.nodeIndex];
    constraintMatrix(3*c1.index + 2, 3*c1.index + 2) = fem1.massInv[c1.element.nodeIndex];

    //rhs = difference in velocities at the constraint points
    auto rigidVelocity1 = rigidBody.bulletBody->getVelocityInLocalPoint(c1.element.localPosition);
    rhs(3*c1.index    ) = rigidVelocity1.x() - fem1.vel[c1.index].x();
    rhs(3*c1.index + 1) = rigidVelocity1.y() - fem1.vel[c1.index].y();
    rhs(3*c1.index + 2) = rigidVelocity1.z() - fem1.vel[c1.index].z();

    //this part = -1/m_rI_33 + rel_i*I^-1 r_j* (* = crossprodut matrix)
    for(auto c2 : enumerate(rigidBody.constraints)){
      
      constraintMatrix.block(3*c1.index, 3*c2.index, 3, 3) =
	-rigidBody.bulletBody->getInvMass()*Eigen::Matrix3d::Identity() +
	(crossProductMatrix(c1.element.localPosition)*
	 eigenInvInertiaTensor*
	 crossProductMatrix(c2.element.localPosition));      
    }
  }

  
  impulseSolutionVector = constraintMatrix.colPivHouseholderQr().solve(rhs);
  
  //apply impulses.  Positive on femNodes, negative on the RB
  for(auto c : enumerate(rigidBody.constraints)){
    world.femObjects[c.element.femIndex].vel[c.element.nodeIndex] += 
      world.femObjects[c.element.femIndex].massInv[c.element.nodeIndex]*
      SlVector3(impulseSolutionVector(3*c.index),
		impulseSolutionVector(3*c.index +1),
		impulseSolutionVector(3*c.index +2));

    rigidBody.bulletBody->applyImpulse(btVector3(impulseSolutionVector(3*c.index),
						 impulseSolutionVector(3*c.index + 1),
						 impulseSolutionVector(3*c.index + 2)),
				       c.element.localPosition);

  }
  //std::cout << "rhs: " << rhs << std::endl;
  //std::cout << "impulse: " << impulseSolutionVector << std::endl;

}
