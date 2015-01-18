
#include "world.h"

#include <iostream>
//#include <fenv.h>
//#include <time.h>
//#include <float.h>

#include "utils.h"


int main(int argc, char *argv[]) {
  //FemSimulator simulator;
  if(argc < 2){
    std::cout << "usage: ./fracture inputFiles/inputFile.json" << std::endl;
    exit(1);
  }
  World world(argv[1]);

  //	double dt = 1/900.0 + 1e-12; 
  //	SlVector3 gravity(0.0,0.0,-9.8);
  double framerate = 1/60.0;//60.0;
  //double framerate = 0.00005960464478;
  //double framerate = 0.00011920928956;
	//	bool ground = true;
	//	unsigned int frame = 0;
	//	char framestring[80];
	//	sprintf(framestring, "frames/foo-%%02i.%%04i.obj");
	//	double friction = 0.5;
  double time = 0;

	//	simulator.nobjects = 1;
	//	simulator.objects = new FemObject[simulator.nobjects];
	//	FemObject &o = simulator.objects[0];
	//	o.load(argv[1]);

	//	for (unsigned int i=0; i<o.nv; i++) {
	//		o.pos[i][2] += 1;
	//	}
	
	//simulator.initCollisions();

  //while (time < 0.2){//5) {
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
}
