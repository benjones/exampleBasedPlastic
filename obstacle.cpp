#include "obstacle.H"
#include "fem.H"

void SphereObstacle::collide(FemObject &o, double dt) {
	double r2 = radius*radius;
	for (unsigned int i=0; i<o.nv; i++) {
		SlVector3 x = o.pos[i] + dt*o.vel[i];
		SlVector3 d = x - center;
		double m2 = sqrMag(d);
		if (m2 < r2 && m2 > 1e-10) {	
			double m = sqrt(m2);
			d/=m;
			SlVector3 impulse = ((radius-m) / dt)*d;
			o.vel[i] += impulse;
			std::cout<<"collision"<<std::endl;
		}
	}
}

void SphereObstacle::step(double dt) {
	vel[2] -= 9.8*dt;
	center += dt*vel;
	if (center[2]-radius < 0.0) {
		vel[2] = -0.9*vel[2];
	}
}
