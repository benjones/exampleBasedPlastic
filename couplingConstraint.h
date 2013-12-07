#pragma once

#include <LinearMath/btVector3.h>

struct CouplingConstraint{

  size_t rigidBodyIndex; //rigid body the constraint belongs to
  btVector3 localPosition; //local position in the rigid body
  size_t femIndex; //FEM body the constraint belongs to
  size_t nodeIndex; //node index in that FEM
};
