
#include <iostream>
#include <fstream>

#ifdef __APPLE__
#include<GLUT/glut.h>
#include<OpenGL/gl.h>
#include<OpenGL/glu.h>
#else
#include<GL/glut.h>
#include<GL/gl.h>
#include<GL/glu.h>
#endif

#include "Eigen/Dense"

#include "json/json.h"
#include "plasticBody.h"
#include "kernels.h"

#include "range.hpp"
using benlib::range;

using RMMatrix3d = Eigen::Matrix<double,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3f = Eigen::Matrix<float ,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor>;



struct Camera{
  Eigen::Vector3d center, eye, up;
  
  int width, height;
  
};

Camera camera;
Eigen::Vector2i mouse;
PlasticBody plasticBody;

struct ImpulseInfo{
  Eigen::Vector3d initialImpulse, tangent;
  int vertexIndex;
};

ImpulseInfo impulseInfo;

int currentFrame = 0;


void keyInput(unsigned char key, int x, int y);
void specialInput(int key, int x, int y);
void reshape(int width, int height);
void mouseClicks(int button, int state, int x, int y);
void mouseMove(int x, int y);
void displayFrame();
void idleFunc();


Eigen::Vector3d getVectorForFrame(const Eigen::Vector3d v, double maxScale, int timestep);


void drawTriangles(bool changeColor);
Eigen::Vector3d get3dMousePos(int x, int y);

Eigen::Vector3d get3dMousePos(int x, int y)
{
  GLint viewport[4];
  GLdouble modelview[16];
  GLdouble projection[16];
  GLdouble winX, winY;
  GLfloat winZ;
  GLdouble posX, posY, posZ;
 
  glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
  glGetDoublev( GL_PROJECTION_MATRIX, projection );
  glGetIntegerv( GL_VIEWPORT, viewport );
 
  winX = static_cast<double>(x);
  winY = static_cast<double>(viewport[3]) - static_cast<double>(y);
  glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
 
  gluUnProject( winX, winY, winZ, 
	  modelview, projection, viewport, 
	  &posX, &posY, &posZ);
 
  return {posX, posY, posZ};
}


int main(int argc, char** argv){

  const std::string usage = "./impulseTester <inputJsonFile>";
  if(argc != 2){
	std::cout << usage << std::endl;
	return 1;
  }

  std::ifstream ins(argv[1]);
  Json::Value root;
  Json::Reader reader;
  if(!reader.parse(ins, root)){
	std::cout << "couldn't read json file " << argv[1] << '\n'
			  << reader.getFormattedErrorMessages() << std::endl;
	return 1;
  }

  auto& pbsIn = root["plasticBodies"];
  if(!pbsIn.size()){
	std::cout << "no plastic bodies" << std::endl;
	return 1;
  }

  auto& pbIn = pbsIn[0];
  if(pbsIn.size() > 1){
	std::cout << "multiple bodies defined, using the first" << std::endl;
  }



  plasticBody.loadFromJsonNoDynamics(pbIn);

  PlasticPiece& plasticPiece = plasticBody.plasticPieces.front();
  
  //impulse stuff
  auto& impulseIn = root["testImpulse"];
  impulseInfo.vertexIndex = impulseIn["vertexIndex"].asInt();
  //compute normal at that point
  
  Eigen::Vector3d normalSoFar = Eigen::Vector3d::Zero();
  for(auto i : range(plasticPiece.tetmeshTriangles.rows())){
	if((plasticPiece.tetmeshTriangles.row(i).array() == 
			impulseInfo.vertexIndex).any()){
	  Eigen::Vector3i tri = plasticPiece.tetmeshTriangles.row(i);
	  Eigen::Vector3d v1 = plasticPiece.currentBulletVertexPositions.row(tri(0));
	  Eigen::Vector3d v2 = plasticPiece.currentBulletVertexPositions.row(tri(1));
	  Eigen::Vector3d v3 = plasticPiece.currentBulletVertexPositions.row(tri(2));
	  
	  Eigen::Vector3d thisNormal = (v1 - v2).cross(v1 - v3);
	  if(thisNormal.dot(normalSoFar) < 0){
		thisNormal *= -1;
	  }
	  normalSoFar += thisNormal;
	}
  }
  normalSoFar.normalize();
  impulseInfo.initialImpulse = plasticBody.plasticityImpulseYield*normalSoFar;
  if(impulseIn.get("invertNormal", false).asBool()){
	impulseInfo.initialImpulse *= -1;
  }

  
  //get a tangent vector	auto 
  Eigen::Vector3d start = impulseInfo.initialImpulse.normalized();
  Eigen::Vector3d c1 = start.cross(Eigen::Vector3d(1,0,0));
  while(c1.norm() < 1e-2){
	c1 = start.cross(Eigen::Vector3d::Random());
  }
  Eigen::Vector3d tangent = c1.cross(start);
  tangent.normalize();
  impulseInfo.tangent = tangent;

  Eigen::Vector3d bbMin = plasticPiece.currentBulletVertexPositions.colwise().minCoeff();
  Eigen::Vector3d bbMax = plasticPiece.currentBulletVertexPositions.colwise().maxCoeff();

  double width = (bbMax - bbMin).norm();
  
  camera.center = (bbMax + bbMin)/2;
  camera.up = Eigen::Vector3d(0,1,0);
  camera.eye = camera.center + width*Eigen::Vector3d(1,0,0);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
  glutInitWindowSize(800,800);
  
  glutCreateWindow("Mesh Viewer");
  camera.width = glutGet(GLUT_SCREEN_WIDTH);
  camera.height = glutGet(GLUT_SCREEN_HEIGHT);
  glutPositionWindow(camera.width/4, camera.height/5);

  glutDisplayFunc(displayFrame);
  glutKeyboardFunc(keyInput);
  glutSpecialFunc(specialInput);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouseClicks);
  glutMotionFunc(mouseMove);
  glutIdleFunc(idleFunc);

  glClearColor(0,0,0,0);
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glViewport(0,0,camera.width, camera.height);
  glutMainLoop();

  return 0;
}

void reshape(int width, int height){
  glViewport(0,0,width, height);
  camera.width = width;
  camera.height = height;
  glutPostRedisplay();
}

void drawTriangles(bool changeColor){
  RMMatrix3d drawableColors(5, 3);
  drawableColors << 1, 0, 0,
	0, 1, 0,
	0, 0, 1,
	0, 0, 0,
	1, 1, 1;
  
  
  glBegin(GL_TRIANGLES);
  //std::cout << "using beginning triangles: " << (oj.tris.size() ? "no" : "yes") << std::endl;
  const auto& piece = plasticBody.plasticPieces.front();
  const auto& tris = piece.tetmeshTriangles;
  const auto& vertices = piece.currentBulletVertexPositions;
  for(size_t i = 0; i < tris.rows(); ++i){
	Eigen::Vector3d v1 = vertices.row(tris(i,0)).transpose();
	Eigen::Vector3d v2 = vertices.row(tris(i,1)).transpose();
	Eigen::Vector3d v3 = vertices.row(tris(i,2)).transpose();
	
	//hack.  Fix to use more than 3 barycentric coordinates
	if(changeColor && piece.barycentricCoordinates.rows() > 0){
	  Eigen::Vector3d color = Eigen::Vector3d::Zero();
	  for(auto j : range(piece.barycentricCoordinates.cols())){
		color += piece.barycentricCoordinates(tris(i, 0),j)*drawableColors.row(j);
	  }
	  glColor3dv(color.data());
	  glVertex3dv(v1.data());
	  
	  color = Eigen::Vector3d::Zero();
	  for(auto j : range(piece.barycentricCoordinates.cols())){
		color += piece.barycentricCoordinates(tris(i, 1),j)*drawableColors.row(j);
	  }
	  glColor3dv(color.data());
	  
	  glVertex3dv(v2.data());
	  
	  color = Eigen::Vector3d::Zero();
	  for(auto j : range(piece.barycentricCoordinates.cols())){
		color += piece.barycentricCoordinates(tris(i, 2),j)*drawableColors.row(j);
	  }
	  glColor3dv(color.data());
	  glVertex3dv(v3.data());
	  
	} else {
	  glVertex3dv(v1.data());
	  glVertex3dv(v2.data());
	  glVertex3dv(v3.data());
	}
	//	glEnd();
    
  }      
  
  glEnd();
}


void displayFrame(){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


  glDisable(GL_CULL_FACE); //triangles aren't wound properly...
  glCullFace(GL_BACK);


  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(camera.eye(0), camera.eye(1), camera.eye(2), 
			camera.center(0), camera.center(1), camera.center(2),
			camera.up(0), camera.up(1), camera.up(2));
  

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40, (double)(camera.width)/camera.height, .5, 500);


  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glMatrixMode(GL_MODELVIEW);
  glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
  glPolygonOffset(1.0, 1.0);
  drawTriangles(true);
  glDisable(GL_POLYGON_OFFSET_FILL);
  
  glColor3d(1, 1, 1); // Color for your polygon border
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  drawTriangles(false);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  

  Eigen::Vector3d startPoint = 
	plasticBody.plasticPieces.front().currentBulletVertexPositions.row(
		impulseInfo.vertexIndex);
  Eigen::Vector3d endPoint = startPoint -
	5*getVectorForFrame(impulseInfo.initialImpulse,
		plasticBody.plasticityImpulseYield +
		100*plasticBody.plasticityImpulseScale,
		currentFrame);
  
  glDisable(GL_DEPTH_TEST);
  glColor3d(0,1,0);
  glLineWidth(3);
  glBegin(GL_LINES);
  glVertex3dv(startPoint.data());
  glVertex3dv(endPoint.data());
  glEnd();
  glEnable(GL_DEPTH_TEST);
  glLineWidth(1);
  glFlush();
  glutSwapBuffers();

  

}


void keyInput(unsigned char key, int x, int y) {
  switch (key) {
  case 27: // ESC
	exit(0);
	break;
  case 'r':
	currentFrame = 0;
	std::cout << "current frame reset" << std::endl;
	break;
  default:
	; //do nothing
  }
}

void specialInput(int key, int x, int y){
  if(key == GLUT_KEY_UP){
	camera.eye += 0.1*(camera.center - camera.eye);
    glutPostRedisplay();
  } else if(key == GLUT_KEY_DOWN){
	camera.eye -= 0.1*(camera.center - camera.eye);
    glutPostRedisplay();
  }
}


void mouseClicks(int button, int state, int x, int y){
  mouse.x() = x;
  mouse.y() = y;

  auto mouse3D = get3dMousePos(x, y);
  auto nearest = plasticBody.plasticPieces.front().getNearestVertex(mouse3D);

  std::cout << "nearest vertex: " << nearest << std::endl;
}

void mouseMove(int x, int y){
  Eigen::Vector3d eyeToCenter = camera.center - camera.eye;
  double dist = eyeToCenter.norm();
  Eigen::Vector3d xDirection = eyeToCenter.cross(camera.up);
  xDirection.normalize();
  Eigen::Vector3d yDirection = xDirection.cross(eyeToCenter);
  yDirection.normalize();

  double dx = x - mouse.x();
  double dy = y - mouse.y();
  mouse.x() = x;
  mouse.y() = y;

  eyeToCenter -= 
	(-50*dx)/(.5*camera.width ) * xDirection + 
	( 50*dy)/(.5*camera.height) * yDirection;


  eyeToCenter.normalize();
  eyeToCenter = eyeToCenter*dist;


  camera.eye = camera.center - eyeToCenter;
  glutPostRedisplay();
  
}

//scale, then spiral out
Eigen::Vector3d getVectorForFrame(const Eigen::Vector3d v, double maxScale, int timestep){

  const double nScaleTimesteps = 120;
  const double nSpiralTimesteps = 120;

  timestep = timestep % static_cast<int>(nScaleTimesteps + nSpiralTimesteps);

  if(timestep < nScaleTimesteps){
	return maxScale*std::min(timestep/nScaleTimesteps, 1.0)*v;
  } else {
	Eigen::Vector3d start = maxScale*v;

	auto fakeTimestep = timestep - nScaleTimesteps;
	std::cout << "faketimestep: " << fakeTimestep << std::endl;
	Eigen::AngleAxis<double> tangentRot(M_PI_2*fakeTimestep/nSpiralTimesteps, impulseInfo.tangent);
	Eigen::AngleAxis<double> normalRot(2*M_PI*fakeTimestep/nSpiralTimesteps, v.normalized());
	
	return normalRot*(tangentRot*start);
	
  }
}

void idleFunc(){
  
  auto impulse = getVectorForFrame(
	  impulseInfo.initialImpulse,
	  plasticBody.plasticityImpulseYield +
	  100*plasticBody.plasticityImpulseScale,
	  currentFrame);
  auto deltaS = 
	plasticBody.projectSingleImpulse(plasticBody.plasticPieces.front(),
		impulse,
		impulseInfo.vertexIndex);
		
  std::cout << "impulse: " << impulse << '\n'
			<< "deltaS: " << deltaS << std::endl;
  auto& pp = plasticBody.plasticPieces.front();


  pp.geodesicDistances.assign(pp.currentBulletVertexPositions.rows(),
	  std::numeric_limits<double>::infinity());  
  pp.geodesicDistances[impulseInfo.vertexIndex] = 0;
  pp.geodesicDistancesPropogate(
	  {static_cast<size_t>(impulseInfo.vertexIndex)},
	  1.0/plasticBody.plasticityKernelScale);


  pp.deltaBarycentricCoordinates.resize(
	  pp.barycentricCoordinates.rows(),
	  pp.barycentricCoordinates.cols());
  pp.deltaBarycentricCoordinates.setZero();
  for(auto i : pp.activeVertices){
	auto scale = 
	  Kernels::simpleCubic(
		  plasticBody.plasticityKernelScale*pp.geodesicDistances[i]);
	pp.deltaBarycentricCoordinates.row(i) += scale*deltaS.transpose();
  }  
  pp.barycentricCoordinates.setZero();
  pp.barycentricCoordinates.col(0).setOnes();
  
  pp.barycentricCoordinates += pp.deltaBarycentricCoordinates;

  
  //clamp and rescale:
  for(auto i : pp.activeVertices){
	for(auto j : range(plasticBody.exampleGraph.nodes.size())){
	  pp.barycentricCoordinates(i, j) = 
		std::max(pp.barycentricCoordinates(i,j), 0.0);
	}
	const double sum = pp.barycentricCoordinates.row(i).sum();
	if(fabs(sum) < 0.0001){
	  pp.barycentricCoordinates(i, 0) = 1;
	} else {
	  pp.barycentricCoordinates.row(i) /= sum;
	}
  }
  
  pp.skinMeshVaryingBarycentricCoords(
	  plasticBody.boneWeights,
	  plasticBody.boneIndices,
	  plasticBody.exampleGraph);

  ++currentFrame;
  glutPostRedisplay();
}
