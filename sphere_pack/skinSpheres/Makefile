
CC=clang++
OPTS= -O3 -Wall -std=c++1y -stdlib=libc++ -g -I/usr/local/include/eigen3

all: skinSpheres

skinSpheres: skinSpheres.o
	$(CC) $(OPTS) -o skinSpheres skinSpheres.o

%.o: %.cpp
	$(CC) $(OPTS) $(INCLUDE_DIRS) -c  -o $@ $<

clean:
	rm -v *.o skinSpheres
