#include<iostream>
#include<vector>
#include<fstream>
#include <cstdlib>
#include <cassert>
#include "Eigen/Dense"
//#include "Common/slVector.H"
//#include "Common/slUtil.H"

#ifdef __APPLE__
#include<GLUT/glut.h>
#include<OpenGL/gl.h>
#include<OpenGL/glu.h>
#else
#include<GL/glut.h>
#include<GL/gl.h>
#include<GL/glu.h>
#endif


//#include "cppitertools/enumerate.hpp"
//using iter::enumerate;
#include "range.hpp"
#include "enumerate.hpp"
using benlib::range;
using benlib::enumerate;

#include "plyIO.hpp"
#include "utils.h"

using RMMatrix3d = Eigen::Matrix<double,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3f = Eigen::Matrix<float ,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor>;
using Vector3 = Eigen::Vector3f;

struct objStruct {
  RMMatrix3f pts;
  RMMatrix3i tris;
  Eigen::MatrixXd barycentricCoordinates;
  //  std::vector<SlVector3> pts;
  //  std::vector<SlTri> tris;
  
};

// struct particleInfo
// {
//     std::vector<double> material;
//     std::vector<double> positions;
//     std::vector<double> color;
// };




double randDouble(){
  return (double)(rand())/((double)(RAND_MAX));
}

std::string makeFilename(const std::string &format, unsigned frame);
void readFrame(const std::string &plyFile, size_t frameNumber);
void displayFrame();
void keyInput(unsigned char key, int x, int y);
void specialInput(int key, int x, int y);
void reshape(int width, int height);
void mouseClicks(int button, int state, int x, int y);
void mouseMove(int x, int y);

void writeMitsuba();
void drawImpulsesForFrame();

//std::vector<particleInfo> globalPInfo; //Gotta because of glut
std::vector<std::vector<objStruct>> globalObjs;
Vector3 center;//for gluLookat
Vector3 eye;
Vector3 upVector;
unsigned currentFrame = 0;
std::string objFormat;
std::string barycentricFormat;
std::string impulseFormat;

int windWidth, windHeight;
double zoomFactor = 1.0;

float lightPos1[4] = {0.0f, 20.0f, 90.0f, 0.0f};
float lightPos2[4] = {3.0f, 20.0f, -30.0f, 0.0f};
float lightColor[4] = {1.0f, 1.0f, 1.0f, 1.0f};
float lightAmbient[4] = {.2f, .2f, .2f, 1.0f};

bool useBarycentricColors;
bool drawImpulses;

int main(int argc, char** argv){
  if(argc < 2){
    std::cout << "Usage: ./openglViewer plyFormatString [barycentricCoordinatesString] [impulse files string]\n"
              << "(eg frames/foo-%03d.%04d.ply)\n";
    exit(1);
  }
  objFormat = argv[1];
  useBarycentricColors = argc > 2;
  if(useBarycentricColors){
	barycentricFormat = argv[2];
  }
  drawImpulses = argc > 3;
  if(drawImpulses){
	impulseFormat = argv[3];
  }
  readFrame(objFormat, 0);
  



  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
  glutInitWindowSize(800,800);
  windWidth = windHeight = 800;

  glutCreateWindow("Mesh Viewer");
  const auto screenWidth = glutGet(GLUT_SCREEN_WIDTH);
  const auto screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
  glutPositionWindow(screenWidth/4, screenHeight/5);
  glutDisplayFunc(displayFrame);
  glutKeyboardFunc(keyInput);
  glutSpecialFunc(specialInput);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouseClicks);
  glutMotionFunc(mouseMove);

  glClearColor(0,0,0,0);
  glEnable(GL_DEPTH_TEST);

  center = Vector3{0, 1, 0};
  eye = Vector3{0, 1, 4};
  upVector = Vector3{0, 1, 0};

  glutMainLoop();
}

//number of particles subject to change

void readFrame(const std::string &plyFile, size_t frameNumber){

  globalObjs.emplace_back();
  for( auto i = 0; ; i++){
    char filename[1024];
    sprintf(filename, plyFile.c_str(), i, frameNumber);
    std::ifstream test(filename);
	std::cout << "reading " << filename << std::endl;
    if(test.good()){
	  
	  std::cout << "read " << filename << std::endl;
	  
	  globalObjs.back().emplace_back();
	  auto& oj = globalObjs.back().back();
	  readPLY(test, oj.pts, oj.tris);
	  
	  //SlVector3 uc, lc;
	  //objStruct oj;
	  //readObjFile(objName.c_str(), oj.pts, oj.tris, uc, lc);
	  std::cout << "num tris: " << oj.tris.rows() << std::endl;


	  if(useBarycentricColors){
		sprintf(filename, barycentricFormat.c_str(), i, frameNumber);
		oj.barycentricCoordinates = readMatrixBinary(filename);
	  }

	}
	/*else{
	  const auto& binName = std::string(filename) + std::string(".bin");
	  std::ifstream ins(binName, std::ios_base::in | std::ios_base::binary);
	  if(ins.good()){
		std::cout << "reading " << binName << std::endl;
		size_t numVerts;
		ins.read(reinterpret_cast<char*>(&numVerts),
				 sizeof(numVerts));
		assert(numVerts = globalObjs.front()[i].pts.size());
		globalObjs.back().emplace_back();
		globalObjs.back().back().pts.resize(numVerts);
		ins.read(reinterpret_cast<char*>(globalObjs.back().back().pts.data()),
				 3*numVerts*sizeof(double));

				 }*/
	else{
	  break;
	}
	//	std::cout << "pt 1: " << globalObjs.back().back().pts[0] << std::endl;
	//	std::cout << "pt 2: " << globalObjs.back().back().pts[1] << std::endl;
  }
  std::cout << globalObjs.back().size() << " objs this frame " << std::endl;
}

				  


void drawTriangles(bool changeColor){
  RMMatrix3d drawableColors(5, 3);
  drawableColors << 1, 0, 0,
	0, 1, 0,
	0, 0, 1,
	0, 0, 0,
	1, 1, 1;


  glBegin(GL_TRIANGLES);
  for(auto&& e : enumerate(globalObjs[currentFrame])){
	auto& oj = e.second;
	if(changeColor){
	  auto greenChannel = static_cast<double>(e.first)/globalObjs[currentFrame].size();
		
		glColor4d(.5,greenChannel,.5,1);
	}
	//std::cout << "using beginning triangles: " << (oj.tris.size() ? "no" : "yes") << std::endl;
	const auto& tris = oj.tris;
    for(size_t i = 0; i < tris.rows(); ++i){
      Vector3 v1 = oj.pts.row(tris(i,0)).transpose();
	  Vector3 v2 = oj.pts.row(tris(i,1)).transpose();
	  Vector3 v3 = oj.pts.row(tris(i,2)).transpose();
	  //		v2 = oj.pts[tris[i].indices[1]],
	  //		v3 = oj.pts[tris[i].indices[2]];
      //glBegin(GL_LINE_LOOP);
      Vector3 normal =(v2 - v1).cross(v3 -v1);
      normal.normalize();
      glNormal3fv(normal.data());
	  //hack.  Fix to use more than 3 barycentric coordinates
	  if(changeColor && useBarycentricColors && oj.barycentricCoordinates.rows() > 0){
		Eigen::Vector3d color = Eigen::Vector3d::Zero();
		for(auto j : range(oj.barycentricCoordinates.cols())){
		  color += oj.barycentricCoordinates(tris(i, 0),j)*drawableColors.row(j);
		}
		glColor3dv(color.data());
		glVertex3fv(v1.data());
		color = Eigen::Vector3d::Zero();
		for(auto j : range(oj.barycentricCoordinates.cols())){
		  color += oj.barycentricCoordinates(tris(i, 1),j)*drawableColors.row(j);
		}
		glColor3dv(color.data());

		glVertex3fv(v2.data());

		color = Eigen::Vector3d::Zero();
		for(auto j : range(oj.barycentricCoordinates.cols())){
		  color += oj.barycentricCoordinates(tris(i, 2),j)*drawableColors.row(j);
		}
		glColor3dv(color.data());
		glVertex3fv(v3.data());

	  } else {
		glVertex3fv(v1.data());
		glVertex3fv(v2.data());
		glVertex3fv(v3.data());
	  }
      //	glEnd();
      
    }      
  }
  glEnd();
}

void drawPoints(){


  glPointSize(5);
  glDisable(GL_DEPTH_TEST);
  glBegin(GL_POINTS);
  for(auto&& e : enumerate(globalObjs[currentFrame])){
	auto redChannel = static_cast<double>(e.first)/globalObjs[currentFrame].size();
	glColor3d(redChannel, 1, 1);
	auto& oj = e.second;
	for(auto i : range(oj.pts.rows())){
	  glVertex3fv(oj.pts.row(i).data());
	}
  }
  glEnd();
  

  glPointSize(8);

  glBegin(GL_POINTS);

  glColor3d(1, 0, 0);
  auto& oj = globalObjs[currentFrame].back();
  for(auto i : range(oj.pts.rows())){
	glVertex3fv(oj.pts.row(i).data());
  }
  
  glEnd();
  glEnable(GL_DEPTH_TEST);

}


void displayFrame(){

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glShadeModel(GL_SMOOTH);
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos1);
  glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.0f);
  glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0f);  
  glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.0f);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
  glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);

  glLightfv(GL_LIGHT1, GL_POSITION, lightPos2);
  glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 1.0f);
  glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.0f);  
  glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.0f);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor);
  glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient);


  glDisable(GL_CULL_FACE); //triangles aren't wound properly...
  glCullFace(GL_BACK);

  //glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(eye(0), eye(1), eye(2), 
			center(0), center(1), center(2),
			upVector(0), upVector(1), upVector(2));
  

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40*zoomFactor, (double)(windWidth)/windHeight, .5, 1000);

    //glLineWidth (2.0);
    //glBegin(GL_LINES);
    //glColor3f(1,0,0); // X axis is red.
    //glVertex3f(0, 0, 0);
    //glVertex3f(2, 0, 0); 
    //glColor3f(0,1,0); // Y axis is green.
    //glVertex3f(0, 0, 0);
    //glVertex3f(0, 2, 0);
    //glColor3f(0,0,1); // z axis is blue.
    //glVertex3f(0, 0, 0);
    //glVertex3f(0, 0, 2); 
    //glEnd();

  /*  glColor4d(.6,.6,.6, 1.0f);
    glBegin(GL_QUADS);
    glVertex3d(-10,-1.5,10);
    glVertex3d(-10,-1.5, -10);
    glVertex3d(10, -1.5, -10);
    glVertex3d(10, -1.5, 10);
    glEnd();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glTranslated(-1.3, -1.5, -0.4);
    glutSolidSphere(0.6, 20, 20);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.5, -1.5, -0.4);
    glutSolidSphere(0.8, 20, 20);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-0.5, -2, 1.3);
    glutSolidSphere(1.2, 20, 20);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.3, -1.5, 0.7);
    glutSolidSphere(0.5, 20, 20);
    glPopMatrix();
  */
    //for (unsigned i = 0; i < globalPInfo.size(); ++i) {

    //draw obj

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

	//drawPoints();

	if(drawImpulses){
	  drawImpulsesForFrame();
	}

    glFlush();
    glutSwapBuffers();
	//writeMitsuba();
	glutSetWindowTitle((std::string("Mesh viewer: Frame ") + std::to_string(currentFrame)).c_str());
}

void keyInput(unsigned char key, int x, int y) {
    switch (key) {
    case 27: // ESC
        exit(0);

    case '0':
        currentFrame = 1;
        specialInput(GLUT_KEY_LEFT, 0, 0);
        break;
    }
}

std::string makeFilename(const std::string &format, unsigned frame)
{
    char result[512];
    sprintf(result, format.c_str(), frame);
    return result;
}

void specialInput(int key, int x, int y){
  if(key == GLUT_KEY_RIGHT){
    std::cout << "right pressed" << std::endl;
	int framesToRead = (glutGetModifiers() & GLUT_ACTIVE_SHIFT) ? 10 : 1;
	for(auto i : range(framesToRead)){
	  (void)i;//unused
	  currentFrame++;
	  std::cout << "currentFrame: " <<currentFrame << std::endl;
	  //if(globalPInfo[0].positions.size() <= (currentFrame)*3){
	  if(globalObjs.size() <= currentFrame){
		readFrame(objFormat, currentFrame);     
	  }
	}
    glutPostRedisplay();
  } else if (key == GLUT_KEY_LEFT){
    std::cout << "left pressed" << std::endl;
	int framesToBacktrack = (glutGetModifiers() & GLUT_ACTIVE_SHIFT) ? 10 : 1;
	for(auto i : range(framesToBacktrack)){
	  (void)i; //unused
		if(currentFrame != 0)
		  --currentFrame;
		std::cout << "currentFrame: " <<currentFrame << std::endl;
	  }
    glutPostRedisplay();
  } else if(key == GLUT_KEY_UP){
	//    zoomFactor /= 1.1;
	eye = eye + 0.1*(center - eye);
    glutPostRedisplay();
  } else if(key == GLUT_KEY_DOWN){
    //zoomFactor *= 1.1;
	eye = eye - 0.1*(center - eye);
    glutPostRedisplay();
  }
}

void reshape(int width, int height){
  glViewport(0,0,width, height);
  windWidth = width;
  windHeight = height;
  glutPostRedisplay();
}

int lastMouseX, lastMouseY;
void mouseClicks(int button, int state, int x, int y){
  lastMouseX = x; lastMouseY = y;
}

void mouseMove(int x, int y){
  Vector3 eyeToCenter = center - eye;
  double dist = eyeToCenter.norm();
  Vector3 xDirection = eyeToCenter.cross(upVector);
  xDirection.normalize();
  Vector3 yDirection = xDirection.cross(eyeToCenter);
  yDirection.normalize();

  double dx = x - lastMouseX;
  double dy = y - lastMouseY;
  lastMouseX = x;
  lastMouseY = y;
  
  
  double sideScale = eyeToCenter.norm();

  eyeToCenter -= (-sideScale*dx)/(.5*windWidth) * xDirection + (sideScale*dy)/(.5*windHeight) * yDirection;


  eyeToCenter.normalize();
  eyeToCenter = eyeToCenter*dist;


  eye = center - eyeToCenter;
  glutPostRedisplay();
  
}

/*void writeMitsuba(){

  char fname[1024];
  sprintf(fname, "frames/mitsubaFrame_%04d.xml", currentFrame);

  std::ofstream outs(fname);

  std::cout << "writing: " << fname << std::endl;

  outs << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
	"<scene version=\"0.5.0\">\n";

  for(size_t i = 0; ; ++i){
	sprintf(fname, objFormat.c_str(), i, currentFrame);
	std::ifstream test(fname);
	if(!test.good()){break;}

	outs << "<shape type=\"obj\">\n"
	  "<string name=\"filename\" value=\"" 
		 <<fname 
		 << "\" />\n"
	  "</shape>\n";

	
  }
  outs << "</scene>" << std::endl;



  }*/


void drawImpulsesForFrame(){
  glDisable(GL_DEPTH_TEST);
  glColor3d(1, 1, 0);
  glBegin(GL_LINES);
  for(auto i : range(globalObjs[currentFrame].size())){
	char filename[1024];
	sprintf(filename, impulseFormat.c_str(), i, currentFrame);
	std::ifstream ins(filename);
	
	if(ins.good()){
	  Vector3 p1, p2;
	  ins >> p1.x() >> p1.y() >> p1.z()
		  >> p2.x() >> p2.y() >> p2.z();
	  while(ins.good()){
		glVertex3fv(p1.data());
		glVertex3fv(p2.data());
		ins >> p1.x() >> p1.y() >> p1.z()
			>> p2.x() >> p2.y() >> p2.z();
	  }
	}
  }
  glEnd();
  glEnable(GL_DEPTH_TEST);
}
