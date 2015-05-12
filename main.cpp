
#include "world.h"

#include <iostream>
//#include <fenv.h>
//#include <time.h>
//#include <float.h>

#include "utils.h"

//#define __CL_ENABLE_EXCEPTIONS
//in makefile
#include "cl.hpp"
#include <tbb/tbb.h>


int main(int argc, char *argv[]) {

  
  


  //FemSimulator simulator;
  if(argc < 2){
    std::cout << "usage: ./fracture inputFiles/inputFile.json" << std::endl;
    exit(1);
  }

  //tbb::task_scheduler_init init(4);

  //  try{
	World world(argv[1]);
	
	double framerate = 1/60.0;//60.0;
	double time = 0;
	
	world.dumpFrame();
	double nextFrameTime = 0;
	while (time < world.duration){
	  time += world.dt;
	  if(time > nextFrameTime){
		world.dumpFrame();
		nextFrameTime += framerate;
	  }
	  world.timeStepDynamicSprites();
	}
	world.profiler.dump<std::chrono::seconds>(std::cout);
	//  } catch (cl::Error& err){
	//	std::cout << "cl error: " << err.what() << std::endl << err.err() << std::endl;;
	//  }
}
