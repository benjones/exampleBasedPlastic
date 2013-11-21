#include "fem.H"
#include <iostream>
#include <fenv.h>
#include <time.h>
#include <float.h>

int main(int argc, char *argv[]) {
	FemSimulator simulator;

	double dt = 1/900.0 + 1e-12; 
	SlVector3 gravity(0.0,0.0,-9.8);
	double framerate = 1/30.0;
	bool ground = true;
	unsigned int frame = 0;
	double timeSinceLastFrame = 100;
	char framestring[80];
	sprintf(framestring, "frames/foo-%%02i.%%04i.obj");
	double friction = 0.5;
	double time = 0;

	simulator.nobjects = 1;
	simulator.objects = new FemObject[simulator.nobjects];
	FemObject &o = simulator.objects[0];
	o.load(argv[1]);

	for (unsigned int i=0; i<o.nv; i++) {
		o.pos[i][2] += 1;
	}
	
	simulator.initCollisions();

	while (time < 5) {
		time += dt;
		timeSinceLastFrame += dt;

		if (timeSinceLastFrame > framerate) {
			timeSinceLastFrame = 0;
			char fname[80];
#ifndef SUPPRESS_OBJ_OUTPUT
			for (unsigned int i=0; i<simulator.nobjects; i++) {
				sprintf(fname, framestring, i, frame);
				std::cout<<"writing frame "<<fname<<std::endl;
				simulator.objects[i].dumpObj(fname);
			}
#endif
			
			frame++;
#ifdef TIMING
			if (frame > 0)
				{
					std::cout<<totalTime/frame<<" "<<computeForcesTime/frame<<" "<<" "<<computeMatrixTime/frame<<" "<<solveTime/frame<<" "<<collisionTime/frame<<std::endl;
				}
#endif
		}

	
		simulator.computeVelocities(dt, gravity);
		
		if (ground) {
			for (unsigned int i=0; i<simulator.nobjects; i++) {
				FemObject &o = simulator.objects[i];
				for (unsigned int j=0; j<o.nv; j++) {
					double x = o.pos[j][2]+dt*o.vel[j][2];
					if (x < 0.0) {
						double impulse = -x / dt;
						double nm=0.0,m = sqrt(sqr(o.vel[j][0]) + sqr(o.vel[j][1]));
						if (m > 0) nm = std::max<double>(m - impulse*friction, 0.0) / m;
						o.vel[j].set(nm*o.vel[j][0], nm*o.vel[j][1], o.vel[j][2]+impulse);
					}
				}
			}
		}
		
		simulator.handleSelfCollisions(dt);
		simulator.updatePositions(dt);
	}
}
