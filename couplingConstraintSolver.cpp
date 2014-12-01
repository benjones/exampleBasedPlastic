#include "couplingConstraintSolver.h"
#include "world.h"
#include "rigidBody.h"

//#include "cppitertools/enumerate.hpp"
//using iter::enumerate;
#include "enumerate.hpp"
using benlib::enumerate;



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
    


  for(auto&& c1 : enumerate(rigidBody.constraints)){
    auto& fem1 = world.femObjects[c1.second.femIndex];
    //diagonal part for change in FEM node velocity due to impulse
    constraintMatrix(3*c1.first    , 3*c1.first    ) = fem1.massInv[c1.second.nodeIndex];
    constraintMatrix(3*c1.first + 1, 3*c1.first + 1) = fem1.massInv[c1.second.nodeIndex];
    constraintMatrix(3*c1.first + 2, 3*c1.first + 2) = fem1.massInv[c1.second.nodeIndex];

    //rhs = difference in velocities at the constraint points
    auto rigidVelocity1 = rigidBody.bulletBody->getVelocityInLocalPoint(c1.second.localPosition);
    rhs(3*c1.first    ) = rigidVelocity1.x() - fem1.vel[c1.first].x();
    rhs(3*c1.first + 1) = rigidVelocity1.y() - fem1.vel[c1.first].y();
    rhs(3*c1.first + 2) = rigidVelocity1.z() - fem1.vel[c1.first].z();

    //this part = -1/m_rI_33 + rel_i*I^-1 r_j* (* = crossprodut matrix)
    for(auto&& c2 : enumerate(rigidBody.constraints)){
      
      constraintMatrix.block(3*c1.first, 3*c2.first, 3, 3) =
	-rigidBody.bulletBody->getInvMass()*Eigen::Matrix3d::Identity() +
	(crossProductMatrix(c1.second.localPosition)*
	 eigenInvInertiaTensor*
	 crossProductMatrix(c2.second.localPosition));      
    }
  }

  
  impulseSolutionVector = constraintMatrix.colPivHouseholderQr().solve(rhs);
  
  //apply impulses.  Positive on femNodes, negative on the RB
  for(auto&& c : enumerate(rigidBody.constraints)){
    world.femObjects[c.second.femIndex].vel[c.second.nodeIndex] += 
      world.femObjects[c.second.femIndex].massInv[c.second.nodeIndex]*
      SlVector3(impulseSolutionVector(3*c.first),
		impulseSolutionVector(3*c.first +1),
		impulseSolutionVector(3*c.first +2));

    rigidBody.bulletBody->applyImpulse(btVector3(impulseSolutionVector(3*c.first),
						 impulseSolutionVector(3*c.first + 1),
						 impulseSolutionVector(3*c.first + 2)),
				       c.second.localPosition);

  }
  //std::cout << "rhs: " << rhs << std::endl;
  //std::cout << "impulse: " << impulseSolutionVector << std::endl;

}
