#-----------------------------------------
#Basic Stuff -----------------------------

CC          = clang++
cc          = clang++


FLAGS=-Wall -DTIMING -Wno-c++11-extensions -std=c++11 -stdlib=libc++\
-DBT_USE_DOUBLE_PRECISION -Wno-deprecated-declarations \
-DBULLET_TRIANGLE_COLLISION

#-----------------------------------------
#Optimization ----------------------------
#OPT   = -O3 -g #-flto for more speed
 #ignore glut deprecation warnings
#-DPCUBE -DINVENTOR_DEFINED
OPT = -O3 -g 

#OPT   = -O3 -g -Wall -DTIMING -Wno-c++11-extensions -std=c++11 -DBT_USE_DOUBLE_PRECISION  -I/usr/include/c++/4.8 -I/usr/include/x86_64-linux-gnu/c++/4.8 #-DPCUBE -DINVENTOR_DEFINED
#OPT   = -g -Wall



#-----------------------------------------
#-----------------------------------------

TARGETS = fracture

OBJECTS =   main.o world.o jsoncpp.o \
rigidBody.o  exampleGraph.o \
plasticPiece.o plasticBody.o
#egTraverser.o plasticObject.o
#fem.o grip.o globalMatrix.o obstacle.o  couplingConstraintSolver.o kdTree.o
HEADERS = *.h *.H *.hpp
#collisions.o, punt on this

#-----------------------------------------

LIBS = -lm  -L./bullet-2.82-r2704/build/src/LinearMath \
-L./bullet-2.82-r2704/build/src/BulletCollision \
-L./bullet-2.82-r2704/build/src/BulletDynamics \
-llapack -lblas -lBulletDynamics \
-lBulletCollision -lLinearMath \
-framework OpenGL -L/usr/local/lib -ltbb #./eltopo/eltopo3d/libeltopo_release.a
#-L./Common
#-L./bullet-2.82-r2704/build/Extras/HACD 
#-lslcommon
#-lHACD

INCS = -I./bullet-2.82-r2704/src/ -I./bullet-2.82-r2704/Extras -I/usr/local/include 
#-I./eltopo/eltopo3d -I./eltopo/common
#-I./Common

EIGEN_INCLUDE= #use our local copy
#FAST_INCLUDE = -I../fastSkinning
#FAST_LIB = -L../fastSkinning -lskinning

IGL_INCLUDE=-I../libigl/include -DIGL_HEADER_ONLY
#IGL_LIB=-L../libigl/lib -ligl -liglmosek -ligltetgen

#MOSEK=/Users/ben/mosek
#MOSEKPLATFORM=osx64x86
#MOSEK_INC=-I$(MOSEK)/7/tools/platform/$(MOSEKPLATFORM)/h
#MOSEK_LIB=-lmosek64 -L$(MOSEK)/7/tools/platform/$(MOSEKPLATFORM)/bin #-liglmosek -lmosek64 -liomp5 -lpthread

#ANT_LIB=-L../libigl/external/AntTweakBar/lib -lAntTweakBar -framework AppKit

#CETRA_LIB=-L../cetra/ -lcetra
#ELTOPO_LIB=-L../cetra/eltopo/eltopo3d -leltopo_release


CCOPTS = $(OPT) $(FLAGS) $(INCS) $(EIGEN_INCLUDE) 
#$(IGL_INCLUDE) 
#$(FAST_INCLUDE) -pg
LDOPTS = $(OPT) 

#$(INCS) $(IGL_LIB) $(FAST_LIB) $(MOSEK_LIB) $(ANT_LIB) $(CETRA_LIB) $(ELTOPO_LIB)#-pg

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS) openglViewer 
#decomposeLaplacian


clean: 
	/bin/rm -fv *.o $(TARGETS) openglViewer

#-----------------------------------------
#-----------------------------------------

fracture: $(OBJECTS) $(HEADERS)
	$(CC) $(LDOPTS) -o fracture $(OBJECTS) $(LIBS) 


openglViewer: openglViewer.cpp
	$(CC) $(CCOPTS) -o openglViewer openglViewer.cpp  \
-framework OpenGL -framework GLUT
#-I./Common -L./Common -lslCommon 
decomposeLaplacian: decomposeLaplacian.cpp
	$(cc) $(CCOPTS) -o decomposeLaplacian decomposeLaplacian.cpp $(EIGEN_INCLUDE) $(IGL_INCLUDE)

remakeObjs: remakeObjs.cpp
	$(CC) $(CCOPTS) -o remakeObjs remakeObjs.cpp $(EIGEN_INCLUDE) $(IGL_INCLUDE)


#-----------------------------------------

.cpp.o: 
	$(CC) $(CCOPTS) -c $< -o $@ 

.C.o: 
	$(CC) $(CCOPTS) -c $< -o $@

.c.o: 
	$(cc) $(CCOPTS) -c $< -o $@

.o: $(OBJECTS)
	$(CC) $(LDOPTS) $(OBJS) $(OBJECTS) $< $(LIBS) -o $@

.C: $(OBJECTS)
	$(CC) $(LDOPTS) $(OBJS) $(OBJECTS) $< $(LIBS) $(FOR_LIB) -o $@

#-----------------------------------------
#-----------------------------------------













