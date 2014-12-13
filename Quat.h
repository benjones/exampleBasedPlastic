#pragma once
#include "Eigen/Geometry"


// typedef for quaternions so it's easier to switch between float and double if
// need be

typedef Eigen::Quaterniond Quat;


inline Quat qlerp(const std::vector<Quat>& quats, std::vector<double> weights){
  assert(!quats.empty());

  Eigen::Matrix<double,4,1> coeffs = weights[0]*quats[0].coeffs();
  for(size_t i = 1; i < quats.size(); ++i){
	if(quats[0].dot(quats[i]) >= 0){
	  coeffs += weights[i]*quats[i].coeffs();
	} else {
	  coeffs -= weights[i]*quats[i].coeffs();
	}
  }
  Quat ret(coeffs);
  ret.normalize();
  return ret;
}


