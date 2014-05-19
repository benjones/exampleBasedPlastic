#pragma once


namespace Kernels{

  /* f(0) = 1, f'(0) = 0, f(1) = 0, f'(1) = 0 f(>1) = 0
   */
  double simpleCubic(double t){
	
	return t > 1 ? 0 : 2*t*t*t - 3*t*t + 1;

  }

}
