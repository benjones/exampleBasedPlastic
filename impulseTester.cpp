
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#ifdef __APPLE__
#include<GLUT/glut.h>
#include<OpenGL/gl.h>
#include<OpenGL/glu.h>
#else
#include<GL/glut.h>
#include<GL/gl.h>
#include<GL/glu.h>
#endif


#include <png.h>
#include <cstdio>

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
  double maxScale;
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


Eigen::Vector3d getVectorForFrame(const Eigen::Vector3d v,  int timestep);


void drawTriangles(bool changeColor);

std::pair<Eigen::Vector3d, Eigen::Vector3d> inline getTangents(Eigen::Vector3d vec){
  
  vec.normalize();
  Eigen::Vector3d v1 = vec.cross(Eigen::Vector3d(1,0,0));
  while(v1.norm() < 1e-2){
	v1 = vec.cross(Eigen::Vector3d::Random());
  }
  Eigen::Vector3d t1 = vec.cross(v1).normalized();
  Eigen::Vector3d t2 = t1.cross(vec).normalized();
  return {t1, t2};

}

void drawArrow(const Eigen::Vector3d& direction, const Eigen::Vector3d& pointPosition);

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


bool dumpMitsuba;
std::string mitsubaBaseName;

void dumpMitsubaFrame();
void dumpPng();
int main(int argc, char** argv){

  const std::string usage = "./impulseTester <inputJsonFile> [impulseFrames/]";
  if(argc < 2){
	std::cout << usage << std::endl;
	return 1;
  }
  dumpMitsuba = argc > 2;
  if(dumpMitsuba){
	mitsubaBaseName = std::string(argv[2]);
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
  impulseInfo.initialImpulse = //plasticBody.plasticityImpulseYield*
	normalSoFar;
  if(impulseIn.get("invertNormal", false).asBool()){
	impulseInfo.initialImpulse *= -1;
  }

  impulseInfo.maxScale = 
	impulseIn.get("maxScale", 1.0).asDouble();
  
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
  
  if(root["eye"].size() == 3){
	camera.eye = Eigen::Vector3d(
		root["eye"][0].asDouble(),
		root["eye"][1].asDouble(),
		root["eye"][2].asDouble());
  } else {
	camera.eye = camera.center + width*Eigen::Vector3d(1,0,0);
  }

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
  glutInitWindowSize(960,540);
  
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
  drawableColors << 
	1, 0, 0,
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
  

  Eigen::Vector3d frameVec = getVectorForFrame(impulseInfo.initialImpulse, currentFrame);
  
  Eigen::Vector3d startPoint = 
	plasticBody.plasticPieces.front().currentBulletVertexPositions.row(
		impulseInfo.vertexIndex);
  //Eigen::Vector3d endPoint = startPoint - 5*frameVec;
  
  glDisable(GL_DEPTH_TEST);
  /*glColor3d(0,1,0);
  glLineWidth(3);
  glBegin(GL_LINES);
  glVertex3dv(startPoint.data());
  glVertex3dv(endPoint.data());
  glEnd();

  glLineWidth(1);
  */
  glColor3f(1,1,1);
  drawArrow(frameVec, startPoint);
  glEnable(GL_DEPTH_TEST);  

  glFlush();
  glutSwapBuffers();

  

}

void drawArrow(const Eigen::Vector3d& direction, const Eigen::Vector3d& pointPosition){
  
  auto tangents = getTangents(direction);
  double length = (direction.norm()/impulseInfo.maxScale);
  double tipHeight = length/10.0;
  double tipWidth = tipHeight;
  double cylinderWidth = tipWidth*0.5;

  Eigen::Vector3d dNorm = direction.normalized();
  
  glBegin(GL_TRIANGLE_FAN);
  glVertex3dv(pointPosition.data());
  int nDivs = 20;
  for(auto i : range(nDivs + 1)){
	Eigen::Vector3d pos = pointPosition - tipHeight*dNorm +
	  cos(i*2.0*M_PI/nDivs)*tipWidth*tangents.first +
	  sin(i*2.0*M_PI/nDivs)*tipWidth*tangents.second;
	glVertex3dv(pos.data());
  }
  glEnd();

  //draw cylinder
  glBegin(GL_QUAD_STRIP);
  for(auto i : range(nDivs + 1)){
	Eigen::Vector3d p1 = pointPosition - tipHeight*dNorm +
	  cos(i*2.0*M_PI/nDivs)*cylinderWidth*tangents.first +
	  sin(i*2.0*M_PI/nDivs)*cylinderWidth*tangents.second;
	Eigen::Vector3d p2 = p1 - length*dNorm;
	glVertex3dv(p1.data());
	glVertex3dv(p2.data());
  }
  glEnd();
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
  std::cout << camera.eye << std::endl;
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


const double nScaleTimesteps = 120;
const double nSpiralTimesteps = 120;

Eigen::Vector3d getVectorForFrame(const Eigen::Vector3d v, int timestep){

  //timestep = timestep % static_cast<int>(nScaleTimesteps + nSpiralTimesteps);

  if(timestep < nScaleTimesteps){
	return impulseInfo.maxScale*std::min(timestep/nScaleTimesteps, 1.0)*v;
  } else {
	Eigen::Vector3d start = impulseInfo.maxScale*v;

	auto fakeTimestep = timestep - nScaleTimesteps;
	//std::cout << "faketimestep: " << fakeTimestep << std::endl;
	Eigen::AngleAxis<double> tangentRot(
		0.5*M_PI_2*std::min(fakeTimestep, nSpiralTimesteps)/nSpiralTimesteps, impulseInfo.tangent);
	Eigen::AngleAxis<double> normalRot(2*M_PI*fakeTimestep/nSpiralTimesteps, v.normalized());
	
	return normalRot*(tangentRot*start);
	
  }
}

void idleFunc(){
  
  auto impulse = getVectorForFrame(
	  impulseInfo.initialImpulse,
	  currentFrame);
  auto deltaS = 
	plasticBody.projectSingleImpulse(plasticBody.plasticPieces.front(),
		impulse,
		impulseInfo.vertexIndex);
		
  //std::cout << "impulse: " << impulse << '\n'
  //<< "deltaS: " << deltaS << std::endl;
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
  if(dumpMitsuba && currentFrame < (nScaleTimesteps + 2*nSpiralTimesteps) ){
	//dumpPng();
	dumpMitsubaFrame();
  }
  glutPostRedisplay();
}

void dumpMitsubaFrame(){
  std::ostringstream paddedFrameStream;
  paddedFrameStream << std::setw(4)  << std::setfill('0') << currentFrame;
  std::string frameString = paddedFrameStream.str();
  //std::cout << frameString  << std::endl;

  const std::string meshFilename = std::string("mesh_") + frameString + ".ply";
  plasticBody.plasticPieces.front().dumpColoredShapeOnly(mitsubaBaseName + "/" + meshFilename);

  std::ofstream outs(mitsubaBaseName + "/mitsubaFrame_" + paddedFrameStream.str() + ".xml");
  
  outs << "<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n"
	   "<scene version=\"0.5.0\">\n"
	"<integrator type=\"direct\" ><boolean name=\"hideEmitters\" value=\"true\" /></integrator>\n"
	"<sensor type=\"perspective\">\n"
	"<float name=\"farClip\" value=\"1241.42\"/>\n"
	"<float name=\"focusDistance\" value=\"10.5293\"/>\n"
	"<float name=\"fov\" value=\"40\"/>\n"
	"<string name=\"fovAxis\" value=\"y\"/>\n"
	"<float name=\"nearClip\" value=\"2.41421\"/>\n"
	"<transform name=\"toWorld\">"
	
	"<lookat target=\"" << camera.center.x() << ',' << camera.center.y() << ',' << camera.center.z() << 
	"\" origin=\"" << camera.eye.x() << ',' << camera.eye.y() << ',' << camera.eye.z() << 
	"\" up=\"0, 1, 0\"/>\n"
	"</transform>\n <sampler type=\"sobol\"> <integer name=\"sampleCount\" value=\"128\"/> </sampler>\n"
	"<film type=\"hdrfilm\">\n"
	"<integer name=\"height\" value=\"540\"/> <integer name=\"width\" value=\"960\"/>\n"
	"<boolean name=\"banner\" value=\"false\"/> <rfilter type=\"gaussian\"/> </film></sensor>\n";

  outs << "<shape type=\"ply\" ><string name=\"filename\" value=\""
	   << meshFilename << "\" />\n"
	"<boolean name=\"faceNormals\" value=\"true\" /><boolean name=\"srgb\" value=\"true\" />\n"
	"<bsdf type=\"mixturebsdf\"><string name=\"weights\" value=\"0.5, 0.5\"/>\n"
	"<bsdf type=\"twosided\"><bsdf type=\"diffuse\" >\n"
	"<texture name=\"reflectance\" type=\"vertexColors\" /></bsdf></bsdf>\n"
	"<bsdf type=\"twosided\"><bsdf type=\"diffuse\">\n"
	"<texture name=\"reflectance\" type=\"wireframe\" >\n"
	"<spectrum name=\"interiorColor\" value=\"1.0\" /><spectrum name=\"color1\" value=\"0.0\" />\n"
	"</texture></bsdf></bsdf></bsdf>\n"
	"</shape>\n";
	   
  //arrow
  Eigen::Vector3d vertexPosition = 
	plasticBody.plasticPieces.front().currentBulletVertexPositions.row(
		impulseInfo.vertexIndex);
  
  auto arrow = 
	getVectorForFrame(impulseInfo.initialImpulse, 
		currentFrame < nScaleTimesteps ? nScaleTimesteps -1 : currentFrame);
  auto axis = Eigen::Vector3d(0,1,0).cross(arrow);

  //  if(axis.norm() < 1e-2){
  //	axis = Eigen::Vector3d::Random().cross(arrow)
  //  }

  //std::cout << "axis: " << axis << std::endl;
  axis.normalize();
  //std::cout << "axis normalized: " << axis << std::endl;

  double angle = acos(std::max(-1.0, std::min(1.0, Eigen::Vector3d(0,1,0).dot(arrow.normalized()))));
  //std::cout << "angle: " << angle << std::endl;

  double scaleValue = 10*getVectorForFrame(impulseInfo.initialImpulse, currentFrame).norm();
  
  outs << "<shape type=\"obj\" ><string name=\"filename\" value=\"arrow.obj\" /><bsdf type=\"diffuse\" >\n"
	"<rgb name=\"reflectance\" value=\"#ffffff\" /></bsdf>\n"
	"<transform name=\"toWorld\" >\n"
	"<scale y=\"" << scaleValue << "\" />\n"
	"<rotate x=\"" << axis.x() << "\" y=\"" << axis.y() << "\" z=\"" << axis.z() << "\" angle=\""
	   << angle*180.0/M_PI << "\" />\n"
	"<translate x=\"" << vertexPosition.x() << "\" y=\"" << vertexPosition.y() << "\" z=\""
	   << vertexPosition.z() << "\" />\n"
	"</transform></shape>\n";


  outs << "<emitter type=\"constant\" ><spectrum name=\"radiance\" value=\"1\" /></emitter>\n";
  //outs << "<emitter type=\"point\"><spectrum name=\"radiance\" value=\"1\" />\n"
  //	"<point name=\"position\" x=\"" << camera.eye.x() << "\" y=\"" 
  //	   << camera.eye.y() << "\" z=\"" << camera.eye.z() << "\" />\n"
  //	"</emitter>";
  outs << "</scene>" << std::endl;

  
}
void dumpPng(){
  std::ostringstream paddedFrameStream;
  paddedFrameStream << std::setw(4)  << std::setfill('0') << currentFrame;
  std::string frameString = paddedFrameStream.str();
  
  std::string filename = mitsubaBaseName + "/frame_" + frameString + ".png";

  FILE* fp = fopen(filename.c_str(), "wb");
  if(!fp){
    std::cout << "couldn't open png file to write" << std::endl;
  }

  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 
						NULL, NULL, NULL);
  if (!png_ptr)
    return;
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr){
    png_destroy_write_struct(&png_ptr,
			     (png_infopp)NULL);
    return;
  }
  
  png_init_io(png_ptr, fp);

  
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  
  //read opengl stuff
  GLubyte* glPixelData = new GLubyte[viewport[2]*viewport[3]*4];
  glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3],
	       GL_RGBA, GL_UNSIGNED_BYTE, glPixelData);

  png_set_IHDR(png_ptr, info_ptr, viewport[2], viewport[3],
	       8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_bytep* row_pointers = new png_bytep[viewport[3]];
  for(int i = 0; i < viewport[3]; ++i){
    row_pointers[i] = glPixelData + 4*(viewport[3] - i -1)*viewport[2];
  }
  png_write_info(png_ptr, info_ptr);
  
  png_write_image(png_ptr, row_pointers);

  png_write_end(png_ptr, NULL);


  delete [] glPixelData;
  delete [] row_pointers;
  fclose(fp);


}
