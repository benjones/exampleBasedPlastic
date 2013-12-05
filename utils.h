#pragma once

#include <LinearMath/btVector3.h>

template <typename T>
inline T sqr(const T& t){
  return t*t;
}

inline std::ostream& operator<<(std::ostream& outs, const btVector3& v){

  outs << '[' << v.x() << ' ' << v.y() << ' ' << v.z() << ']';
  return outs;
}
