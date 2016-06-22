
CC=clang++
OPTS= -O3 -Wall -std=c++1y -stdlib=libc++ -g -I/usr/local/include/eigen3

all: acg

acg: acg.o slUtil.o slVector.o
	$(CC) $(OPTS) -o acg acg.o slUtil.o slVector.o

%.o: %.cpp
	$(CC) $(OPTS) $(INCLUDE_DIRS) -c  -o $@ $<

clean:
	rm -v *.o acg
