#pragma once
#include <eigen3/Eigen/Eigen>
#include <LinearMath/btVector3.h>

struct CouplingConstraint{

  btVector3 localPosition; //local position in the rigid body
  size_t femIndex; //FEM body the constraint belongs to
  size_t nodeIndex; //node index in that FEM
  size_t rbIndex;
  Eigen::Matrix3d crossProductMatrix; //updated each time we solve.
	double weight;
  //cross product matrix of currentRigidRotate(localPositition)

};

inline Eigen::Matrix3d crossProductMatrix(const btVector3& v){
  
  Eigen::Matrix3d ret;
  ret << 
    0,      -v.z(), v.y(),
    v.z(),  0,      -v.x(),
    -v.y(), v.x(),  0;
    
  return ret;
}

inline void inPlaceMult(double* in, double* out, const Eigen::Matrix3d& mat){
  
  out[0] += in[0]*mat(0,0) + in[1]*mat(0,1) + in[2]*mat(0,2);
  out[1] += in[0]*mat(1,0) + in[1]*mat(1,1) + in[2]*mat(1,2);
  out[2] += in[0]*mat(2,0) + in[1]*mat(2,1) + in[2]*mat(2,2);

}

inline void inPlaceMultTranspose(double* in, double* out, const Eigen::Matrix3d& mat){
  out[0] += in[0]*mat(0,0) + in[1]*mat(1,0) + in[2]*mat(2,0);
  out[1] += in[0]*mat(0,1) + in[1]*mat(1,1) + in[2]*mat(2,1);
  out[2] += in[0]*mat(0,2) + in[1]*mat(1,2) + in[2]*mat(2,2);

}

inline void bulletOverwriteMultiply(double* in, double* out, const btMatrix3x3& mat){
  
  out[0] = in[0]*mat[0][0] + in[1]*mat[0][1] + in[2]*mat[0][2];
  out[1] = in[0]*mat[1][0] + in[1]*mat[1][1] + in[2]*mat[1][2];
  out[2] = in[0]*mat[2][0] + in[1]*mat[2][1] + in[2]*mat[2][2];

}
