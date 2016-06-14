
#include "world.h"

#include <iostream>
#include <fstream>
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

  tbb::task_scheduler_init init(6);

  //  try{
	double framerate = (argc == 2) ? 1/60.0 : atof(argv[2]);//60.0;
	std::cout << "recip framerate: " << 1.0/framerate << std::endl;

	World world(argv[1]);
	
	
	double time = 0;
	
	std::ofstream timingOut("timingInfo.txt");

	world.dumpFrame();
	//world.currentFrame++;
	double nextFrameTime = 0;
	int frameCount = 0;
	while (time < world.duration){
	  time += world.dt;
	  if(time > nextFrameTime){
		world.dumpFrame();
		//++world.currentFrame;

		//world.profiler.dump<std::chrono::milliseconds>(timingOut);
		//		world.profiler.clear();

		nextFrameTime += framerate;
		std::cout << "done with frame: " << ++frameCount << std::endl;
	  }
	  world.timeStepDynamicSprites();
	  
	}
	world.profiler.dump<std::chrono::milliseconds>(std::cout);
	std::cout << world.profiler.totalTime<std::chrono::milliseconds>() << std::endl;
	//  } catch (cl::Error& err){
	//	std::cout << "cl error: " << err.what() << std::endl << err.err() << std::endl;;
	//  }
}
