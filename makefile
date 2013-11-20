#-----------------------------------------
#Basic Stuff -----------------------------
CC          = g++
cc          = gcc

#-----------------------------------------
#Optimization ----------------------------
#OPT   = -O3 -fast -g -DPCUBE -DINVENTOR_DEFINED -Wall
OPT   = -g -Wall

#-----------------------------------------
#-----------------------------------------

TARGETS = fracture
OBJECTS = collisions.o fem.o grip.o main.o globalMatrix.o obstacle.o

#-----------------------------------------

LIBS = -lm -L/usr/local/lib -L../../Library/Common -lslcommon 	../../Misc/eltopo-old/eltopo3d/libeltopo_release.a -llapack -lblas
INCS = -I/usr/local/include -I../../Library/Common -I../../Misc/eltopo-old/eltopo3d -I../../Misc/eltopo-old/common

CCOPTS = $(OPT) $(DEBUG) $(INCS) -pg
LDOPTS = $(OPT) $(DEBUG) $(INCS) -pg

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)


clean: 
	/bin/rm -f *.o $(TARGETS)

#-----------------------------------------
#-----------------------------------------

fracture: $(OBJECTS)
	g++ $(LDOPTS) -o fracture $(OBJECTS) $(LIBS) 

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













