#-----------------------------------------
#Basic Stuff -----------------------------
CC          = clang++ #`which c++` #g++
cc          = clang #`which c++` #gcc

#-----------------------------------------
#Optimization ----------------------------
OPT   = -O3 -g -Wall -DTIMING -Wno-c++11-extensions -std=c++11 -DBT_USE_DOUBLE_PRECISION  -I/usr/include/c++/4.8 -I/usr/include/x86_64-linux-gnu/c++/4.8 #-DPCUBE -DINVENTOR_DEFINED
#OPT = -O0 -g -Wall -Wno-c++11-extensions -std=c++11 -stdlib=libc++ -DBT_USE_DOUBLE_PRECISION
#OPT   = -g -Wall

#-----------------------------------------
#-----------------------------------------

TARGETS = fracture
OBJECTS =  fem.o grip.o main.o globalMatrix.o obstacle.o world.o jsoncpp.o rigidBody.o couplingConstraintSolver.o kdTree.o #collisions.o, punt on this

#-----------------------------------------

LIBS = -lm -L./Common -L./bullet-2.82-r2704/build/src/LinearMath -L./bullet-2.82-r2704/build/src/BulletCollision -L./bullet-2.82-r2704/build/src/BulletDynamics -lslcommon -llapack -lblas -lBulletDynamics -lBulletCollision -lLinearMath  #./eltopo/eltopo3d/libeltopo_release.a
INCS = -I./Common -I./bullet-2.82-r2704/src/ -I/usr/local/include #-I./eltopo/eltopo3d -I./eltopo/common

CCOPTS = $(OPT) $(INCS) #-pg
LDOPTS = $(OPT) $(INCS) #-pg

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS) openglViewer


clean: 
	/bin/rm -fv *.o $(TARGETS) openglViewer

#-----------------------------------------
#-----------------------------------------

fracture: $(OBJECTS)
	$(CC) $(LDOPTS) -o fracture $(OBJECTS) $(LIBS) 


openglViewer: openglViewer.cpp
	$(CC) $(OPT) -o openglViewer openglViewer.cpp -I./Common -L./Common -lslCommon -framework OpenGL -framework GLUT
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













