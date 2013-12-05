#-----------------------------------------
#Basic Stuff -----------------------------
CC          = `which c++` #g++
cc          = `which c++` #gcc

OS = $(shell uname -s)

#-----------------------------------------
#Optimization ----------------------------
#Dans OPT flags
# -march=native was not available in gcc 4.2 which is the latest version on the mac
ifeq ($(OS),Darwin)
OPT	= -O3 -fomit-frame-pointer -mtune=native -finline-functions -funroll-loops -finline -fno-strict-aliasing -g -stdlib=libc++ -std=c++11
else
OPT	= -O3 -fomit-frame-pointer -mtune=native -march=native -finline-functions -funroll-loops -finline -fno-strict-aliasing -g
endif


#-----------------------------------------
# Debug mode -----------------------------
DEBUG = -g -pg
#DEBUG = -pg
#DEBUG = 

#-----------------------------------------
#-----------------------------------------

DEBUG_TARGET = libslcommonD.a
DEBUG_OBJECTS = slVectorD.o slMatrixD.o slRotationsD.o slArrayD.o slUtilD.o slBoundsD.o

RELEASE_TARGET = libslcommon.a
RELEASE_OBJECTS = slVector.o slMatrix.o slRotations.o slArray.o slUtil.o slBounds.o
#-----------------------------------------

default: $(RELEASE_TARGET) $(DEBUG_TARGET)

all: $(RELEASE_TARGET) $(DEBUG_TARGET)

debug: $(DEBUG_TARGET)

release: $(RELEASE_TARGET)

clean: 
	/bin/rm -f *.o $(DEBUG_TARGET) $(RELEASE_TARGET)

#-----------------------------------------
#-----------------------------------------

libslcommon.a: $(RELEASE_OBJECTS) 
	ar cr libslcommon.a $(RELEASE_OBJECTS) 
	ranlib libslcommon.a

libslcommonD.a: $(DEBUG_OBJECTS)
	ar cr libslcommonD.a $(DEBUG_OBJECTS)
	ranlib libslcommonD.a

#-----------------------------------------


%D.o: %.cpp 
	$(CC) -c -DDEBUG $(DEBUG) ${INCLUDEFLAGS} -o $@ $< 

%D.o: %.c 
	$(cc) -c -DDEBUG $(DEBUG) ${INCLUDEFLAGS} -o $@ $< 

%.o: %.cpp
	$(CC) -c -DNDEBUG $(OPT) ${INCLUDEFLAGS} -o $@ $< 

%.o: %.c
	$(cc) -c -DNDEBUG $(OPT) ${INCLUDEFLAGS} -o $@ $<




#-----------------------------------------
#-----------------------------------------













