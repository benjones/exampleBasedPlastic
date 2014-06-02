#include<iostream>
#include<vector>
#include<fstream>
#include <cstdlib>
#include <cassert>
#include "Common/slVector.H"
#include "Common/slUtil.H"

#ifdef __APPLE__
#include<GLUT/glut.h>
#include<OpenGL/gl.h>
#include<OpenGL/glu.h>
#else
#include<GL/glut.h>
#include<GL/gl.h>
#include<GL/glu.h>
#endif

#include "cppitertools/enumerate.hpp"
using iter::enumerate;

struct objStruct {
  std::vector<SlVector3> pts;
  std::vector<SlTri> tris;
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
void readFrame(const std::string &objFile, size_t frameNumber);
void displayFrame();
void keyInput(unsigned char key, int x, int y);
void specialInput(int key, int x, int y);
void reshape(int width, int height);
void mouseClicks(int button, int state, int x, int y);
void mouseMove(int x, int y);

void writeMitsuba();


//std::vector<particleInfo> globalPInfo; //Gotta because of glut
std::vector<std::vector<objStruct>> globalObjs;
SlVector3 center;//for gluLookat
SlVector3 eye;
SlVector3 upVector;
unsigned currentFrame = 0;
std::string objFormat;

int windWidth, windHeight;
double zoomFactor = 1.0;

float lightPos1[4] = {0.0f, 20.0f, 90.0f, 0.0f};
float lightPos2[4] = {3.0f, 20.0f, -30.0f, 0.0f};
float lightColor[4] = {1.0f, 1.0f, 1.0f, 1.0f};
float lightAmbient[4] = {.2f, .2f, .2f, 1.0f};

int main(int argc, char** argv){
  if(argc < 2){
    std::cout << "Usage: ./openglViewer objFormatString\n"
              << "(eg frames/foo-%03d.%04d)\n";
    exit(1);
  }
  objFormat = argv[1];
  readFrame(objFormat, 0);
  


  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
  glutInitWindowSize(800,800);
  windWidth = windHeight = 800;

  glutCreateWindow("Mesh Viewer");
  glutDisplayFunc(displayFrame);
  glutKeyboardFunc(keyInput);
  glutSpecialFunc(specialInput);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouseClicks);
  glutMotionFunc(mouseMove);

  glClearColor(0,0,0,0);
  glEnable(GL_DEPTH_TEST);

  center[0] = center[2] = 0;
  center[1] = 1; 

  eye[0] = 0;
  eye[1] = 1;
  eye[2] = 4;

  upVector[0] = upVector[2] = 0;
  upVector[1] = 1;  

  glutMainLoop();
}

//number of particles subject to change

void readFrame(const std::string &objFile, size_t frameNumber){

  globalObjs.emplace_back();
  for( auto i = 0; ; i++){
    char filename[1024];
    sprintf(filename, objFile.c_str(), i, frameNumber);
    
    std::ifstream test(filename);
    if(!test.good()){
      break;
    }
    std::cout << "reading " << filename << std::endl;
    
    SlVector3 uc, lc;
    objStruct oj;
    readObjFile(filename, oj.pts, oj.tris, uc, lc);
    std::cout << "num tris: " << oj.tris.size() << std::endl;
    globalObjs.back().push_back(std::move(oj));


  }
  std::cout << globalObjs.back().size() << " objs this frame " << std::endl;
}

void drawTriangles(bool changeColor){
  glBegin(GL_TRIANGLES);
  for(auto e : enumerate(globalObjs[currentFrame])){
	auto& oj = e.element;
	if(changeColor){
	  auto greenChannel = static_cast<double>(e.index)/globalObjs[currentFrame].size();
		
		glColor4d(.5,greenChannel,.5,1);
	}
    for(size_t i = 0; i < oj.tris.size(); ++i){
      SlVector3 v1 = oj.pts[oj.tris[i].indices[0]],
	v2 = oj.pts[oj.tris[i].indices[1]],
	v3 = oj.pts[oj.tris[i].indices[2]];
      //glBegin(GL_LINE_LOOP);
      SlVector3 normal = cross(v2 - v1, v3 -v1);
      normalize(normal);
      glNormal3d(normal[0], normal[1], normal[2]);
      glVertex3d(v1[0], v1[1], v1[2]);
      glVertex3d(v2[0], v2[1], v2[2]);
      glVertex3d(v3[0], v3[1], v3[2]);
      //	glEnd();
      
    }      
  }
  glEnd();
}

void drawPoints(){


  glPointSize(5);
  glDisable(GL_DEPTH_TEST);
  glBegin(GL_POINTS);
  for(auto e : enumerate(globalObjs[currentFrame])){
	auto redChannel = static_cast<double>(e.index)/globalObjs[currentFrame].size();
	glColor3d(redChannel, 1, 1);
	auto& oj = e.element;
	for(auto& v : oj.pts){
	  glVertex3d(v[0], v[1], v[2]);
	}
  }
  glEnd();
  

  glPointSize(8);

  glBegin(GL_POINTS);

  glColor3d(1, 0, 0);
  auto& oj = globalObjs[currentFrame].back();
  for(auto& v : oj.pts){
	glVertex3d(v[0], v[1], v[2]);
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


  glEnable(GL_CULL_FACE);
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
  gluLookAt(eye[0], eye[1], eye[2], center[0], center[1], center[2],
	    upVector[0], upVector[1], upVector[2]);


  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40*zoomFactor, (double)(windWidth)/windHeight, .5, 200);

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

    glFlush();
    glutSwapBuffers();
	//writeMitsuba();
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
    currentFrame++;
    std::cout << "currentFrame: " <<currentFrame << std::endl;
    //if(globalPInfo[0].positions.size() <= (currentFrame)*3){
    if(globalObjs.size() <= currentFrame){
      readFrame(objFormat, currentFrame);     
    }
    glutPostRedisplay();
  } else if (key == GLUT_KEY_LEFT){
    std::cout << "left pressed" << std::endl;
    if(currentFrame != 0)
      --currentFrame;
    std::cout << "currentFrame: " <<currentFrame << std::endl;
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
  SlVector3 eyeToCenter = center - eye;
  double dist = mag(eyeToCenter);
  SlVector3 xDirection = cross(eyeToCenter, upVector);
  normalize(xDirection);
  SlVector3 yDirection = cross(xDirection, eyeToCenter);
  normalize(yDirection);

  double dx = x - lastMouseX;
  double dy = y - lastMouseY;
  lastMouseX = x;
  lastMouseY = y;

  eyeToCenter -= (-4*dx)/(.5*windWidth) * xDirection + (4*dy)/(.5*windHeight) * yDirection;


  normalize(eyeToCenter);
  eyeToCenter = eyeToCenter*dist;


  eye = center - eyeToCenter;
  glutPostRedisplay();
  
}

void writeMitsuba(){

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



}
