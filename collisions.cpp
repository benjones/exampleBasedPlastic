#include <iostream>
#include "fem.H"
#include <sys/time.h>
#include <dynamicsurface.h>

void FemSimulator::initCollisions() {
	unsigned int nv=0,totalcv=0,totalv=0,totalt=0;
	std::vector<Vec3d> vertex_positions;
	std::vector<Vec3ui> triangles;
	std::vector<double> masses;
	for (unsigned int i=0; i<nobjects; i++) {
		totalcv += objects[i].nv;
		//std::cout<<objects[i].nv<<std::endl;
	}
	bool *used = new bool[totalcv];
	for (unsigned int i=0; i<totalcv; i++) used[i] = false;

	//std::cout<<totalcv<<std::endl;
	
	for (unsigned int i=0; i<nobjects; i++) {
		for (unsigned int j=0; j<objects[i].ntris; j++) {
			used[objects[i].tris[j][0]+nv] = true;
			used[objects[i].tris[j][1]+nv] = true;
			used[objects[i].tris[j][2]+nv] = true;
		}
		nv += objects[i].nv;
		totalt += objects[i].ntris;
	}

	//std::cout<<totalcv<<std::endl;
	int *newIndex = new int[totalcv];
	nv = 0;
	totalv = 0;
	for (unsigned int i=0; i<nobjects; i++) {
		for (unsigned int j=0; j<objects[i].nv; j++) {
			if (used[j+nv]) newIndex[j+nv] = totalv++;
			else newIndex[j+nv] = -1;
		}
		nv += objects[i].nv;
	}
	
	//std::cout<<totalv<<std::endl;
	collisionSurfaceToVertices = new SlInt2[totalv];
	vertex_positions.resize(totalv);
	triangles.resize(totalt);
	masses.resize(totalv);

	totalv = 0;
	nv = 0;
	for (unsigned int i=0; i<nobjects; i++) {
		for (unsigned int j=0; j<objects[i].nv; j++) {
			if (used[j+nv]) collisionSurfaceToVertices[totalv++].set(i,j);
		}
		nv += objects[i].nv;
	}

	//std::cout<<"a "<<totalv<<std::endl;
	totalt = 0;
	nv = 0;
	for (unsigned int i=0; i<nobjects; i++) {
		for (unsigned int j=0; j<objects[i].ntris; j++) {
			triangles[totalt][0] = newIndex[objects[i].tris[j][0]+nv];
			triangles[totalt][1] = newIndex[objects[i].tris[j][1]+nv];
			triangles[totalt][2] = newIndex[objects[i].tris[j][2]+nv];
			//std::cout<<newIndex[objects[i].tris[j][0]+nv]<<" "<<newIndex[objects[i].tris[j][1]+nv]<<" "<<newIndex[objects[i].tris[j][2]+nv]<<" "<<std::endl;
			//std::cout<<triangles[totalt][0]<<" "<<triangles[totalt][1]<<" "<<triangles[totalt][2]<<std::endl;
			totalt++;
		}
		nv += objects[i].nv;
	}
	
	//std::cout<<totalv<<std::endl;
	for (unsigned int i=0; i<totalv; i++) {
		SlVector3 &x = objects[collisionSurfaceToVertices[i][0]].pos[collisionSurfaceToVertices[i][1]];
		vertex_positions[i] = Vec3d(x[0], x[1], x[2]);
		masses[i] = objects[collisionSurfaceToVertices[i][0]].mass[collisionSurfaceToVertices[i][1]];
	}

	//std::cout<<totalv<<std::endl;
	collisionSurface = new DynamicSurface(vertex_positions, triangles, masses, 1e-4, true, false);
	collisionSurface->m_mesh.update_connectivity(totalv);
	collisionSurface->m_velocities.resize(totalv);

	nCollisionSurfaceVertices = totalv;

	delete [] newIndex;
	delete [] used;
}

void FemSimulator::handleSelfCollisions(double dt) {
	timeval startTime, endTime;
	gettimeofday(&startTime, NULL);

	if (nCollisionSurfaceVertices == 0) return;
	for (unsigned int i=0; i<nCollisionSurfaceVertices; i++) {
		SlVector3 &x = objects[collisionSurfaceToVertices[i][0]].pos[collisionSurfaceToVertices[i][1]];
		collisionSurface->m_positions[i] = Vec3d(x[0], x[1], x[2]);
		SlVector3 &v = objects[collisionSurfaceToVertices[i][0]].vel[collisionSurfaceToVertices[i][1]];
		collisionSurface->m_velocities[i] = Vec3d(v[0], v[1], v[2]);
	}
	collisionSurface->integrate(dt);

	for (unsigned int i=0; i<nCollisionSurfaceVertices; i++) {
		Vec3d &x = collisionSurface->m_positions[i];
		SlVector3 y(x[0], x[1], x[2]);
		//SlVector3 & v = objects[collisionSurfaceToVertices[i][0]].vel[collisionSurfaceToVertices[i][1]];
		objects[collisionSurfaceToVertices[i][0]].vel[collisionSurfaceToVertices[i][1]] = 
			(y - objects[collisionSurfaceToVertices[i][0]].pos[collisionSurfaceToVertices[i][1]]) / dt;
	}

	gettimeofday(&endTime, NULL);
	collisionTime += (endTime.tv_sec - startTime.tv_sec) +
		(endTime.tv_usec-startTime.tv_usec)*1.0e-6;
}
