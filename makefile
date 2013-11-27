#-----------------------------------------
#Basic Stuff -----------------------------
CC          = `which c++` #g++
cc          = `which c++` #gcc

#-----------------------------------------
#Optimization ----------------------------
OPT   = -O3 -g -Wall -DTIMING -Wno-c++11-extensions -std=c++11 -stdlib=libc++ -DBT_USE_DOUBLE_PRECISION #-DPCUBE -DINVENTOR_DEFINED
#OPT   = -g -Wall

#-----------------------------------------
#-----------------------------------------

TARGETS = fracture
OBJECTS =  fem.o grip.o main.o globalMatrix.o obstacle.o world.o jsoncpp.o rigidBody.o #collisions.o, punt on this

#-----------------------------------------

LIBS = -lm -L./Common -L./bullet-2.82-r2704/build/src/LinearMath -L./bullet-2.82-r2704/build/src/BulletCollision -L./bullet-2.82-r2704/build/src/BulletDynamics -lslcommon -llapack -lblas -lBulletDynamics -lBulletCollision -lLinearMath  #./eltopo/eltopo3d/libeltopo_release.a
INCS = -I/usr/local/include -I./Common -I./bullet-2.82-r2704/src/ #-I./eltopo/eltopo3d -I./eltopo/common

CCOPTS = $(OPT) $(DEBUG) $(INCS) #-pg
LDOPTS = $(OPT) $(DEBUG) $(INCS) #-pg

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)


clean: 
	/bin/rm -f *.o $(TARGETS)

#-----------------------------------------
#-----------------------------------------

fracture: $(OBJECTS)
	$(CC) $(LDOPTS) -o fracture $(OBJECTS) $(LIBS) 

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













