#pragma once

#include <LinearMath/btVector3.h>
#include <LinearMath/btQuaternion.h>
#include <LinearMath/btTransform.h>
#include <Eigen/Dense>
#include <vector>
#include <fstream>

#include "range.hpp"
using benlib::range;

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

template<typename T>
inline std::ostream& operator<<(std::ostream& outs, const std::vector<T>& t){
  for(auto&& u : t){
	outs << u << ' ';
  }
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

inline btMatrix3x3 eigenToBullet(const Eigen::Matrix3d& mat){
  return btMatrix3x3{mat(0,0), mat(0,1), mat(0,2),
	  mat(1,0), mat(1,1), mat(1, 2),
	  mat(2,0), mat(2,1), mat(2,2)};
}

inline btVector3 eigenToBullet(const Eigen::Vector3d& v){
  return btVector3{v(0), v(1), v(2)};
}


template <typename T>
inline std::pair<T,T> makeSortedPair(T a, T b){
  return (a < b) ?
	std::pair<T,T>(a, b) :
	std::pair<T,T>(b,a);
}

inline std::vector<double> eigenToStd(const Eigen::VectorXd& v){
  return std::vector<double>{v.data(), v.data() + v.rows()};
}

inline Eigen::VectorXd stdToEigen(const std::vector<double>& v){
  Eigen::VectorXd ret(v.size());
  std::copy(v.begin(), v.end(), ret.data());
  return ret;
}

template <typename TType>
void writeEle(std::string filename, const TType& triangles){
  std::ofstream outs(filename);
  assert(triangles.cols() == 3);
  assert(outs);
  outs << triangles.rows() << ' ' << 3 << ' ' << "0\n";
  for(auto i : range(triangles.rows())){
	outs << i << ' ' << triangles(i,0) << ' ' << triangles(i, 1) << ' ' << triangles(i, 2) << '\n';
  }
  
}

template<typename T>
void checkNans(const T * const data, size_t size){
  assert(!std::any_of(data, data + size, [](const T& t){return std::isnan(t);}));
}

template <typename T>
inline void checkNans(const std::vector<T>& vec){
  checkNans(vec.data(), vec.size());
}
template <typename EigenType>
inline void checkNans(const EigenType& mat){
  checkNans(mat.data(), mat.rows()*mat.cols());
}


template <typename EigenType>
void writeMatrixBinary(const std::string& filename, const EigenType& mat){
  
  std::ofstream outs(filename);
  size_t rows = mat.rows();
  size_t cols = mat.cols();
  outs.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
  outs.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
  outs.write(reinterpret_cast<const char*>(mat.data()), rows*cols*sizeof(double));

}

inline Eigen::MatrixXd readMatrixBinary(const std::string& filename){
  std::ifstream ins(filename);
  size_t rows, cols;
  ins.read(reinterpret_cast<char*>(&rows), sizeof(rows));
  ins.read(reinterpret_cast<char*>(&cols), sizeof(cols));
  if(!ins.good()){
	return Eigen::MatrixXd(0,0);
  }
  Eigen::MatrixXd ret(rows, cols);
  ins.read(reinterpret_cast<char*>(ret.data()), rows*cols*sizeof(double));
  return ret;
}
