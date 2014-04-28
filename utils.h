#pragma once

#include <LinearMath/btVector3.h>
#include <LinearMath/btQuaternion.h>
#include <LinearMath/btTransform.h>
#include <Eigen/Dense>

template <typename T>
inline T sqr(const T& t){
  return t*t;
}

inline std::ostream& operator<<(std::ostream& outs, const btVector3& v){
  
  outs << '[' << v.x() << ' ' << v.y() << ' ' << v.z() << ']';
  return outs;
}

inline std::ostream& operator<<(std::ostream& outs, const btQuaternion& q){
  outs << '[' << q.getAxis() << ", " << q.getAngle() << ']';
  return outs;
}

inline std::ostream& operator<<(std::ostream& outs, const btMatrix3x3& t){
  outs << '[' << t[0][0] << ' ' << t[0][1] << ' ' << t[0][2] << "]\n"
	   << '[' << t[1][0] << ' ' << t[1][1] << ' ' << t[1][2] << "]\n"
	   << '[' << t[2][0] << ' ' << t[2][1] << ' ' << t[2][2] << ']';

  return outs;
}


inline std::ostream& operator<<(std::ostream& outs, const btTransform& t){
  outs << '[' << t.getOrigin() << ", \n" << t.getBasis() << ']';
  return outs;
}


inline Eigen::Matrix3d bulletToEigen(const btMatrix3x3& mat){
  Eigen::Matrix3d ret;
  ret << 
	mat[0][0], mat[0][1], mat[0][2],
	mat[1][0], mat[1][1], mat[1][2],
	mat[2][0], mat[2][1], mat[2][2];
  return ret;
}


inline Eigen::Vector3d bulletToEigen(const btVector3& vec){
  return Eigen::Vector3d{vec[0], vec[1], vec[2]};
}
