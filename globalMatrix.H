#ifndef __GLOBALMATRIX_H__
#define __GLOBALMATRIX_H__

#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <vector>
#include <algorithm>
#include <slMatrix.H>

class GlobalMatrix {
public:
	unsigned int dim;
	unsigned int *colIndices; 
	unsigned int *cols; 
	unsigned int *nblocks; //number of blocks currently in each col (on or below the diagonal)
	SlMatrix3x3 *blocks; //the blocks
	unsigned int *neighborIndices;
	unsigned int totalblocks;
	SlMatrix3x3 zero;
	unsigned int nonzeros;
	
	GlobalMatrix():dim(0),colIndices(NULL),nblocks(NULL),blocks(NULL),zero(0.0){}

	void output();

	GlobalMatrix(unsigned int nvertices, unsigned int ntets, const Tet *tets);

	SlMatrix3x3 &find(unsigned int r, unsigned int c) {
		unsigned int i=colIndices[r];
		unsigned int j=colIndices[r+1]-1;
		unsigned int k = (i+j)/2;
		while (i<j) {
			//std::cout<<i<<" "<<j<<" "<<k<<" "<<r<<" "<<c<<" "<<std::endl;
			if (cols[k] == c) break;
			if (cols[k] < c) i = k+1;
			else j = k;
			k = (i+j)/2;
		}
		//if (cols[k] != c) std::cout<<r<<" "<<c<<" "<<k<<" "<<cols[k]<<" "<<colIndices[r]<<std::endl<<blocks[k]<<std::endl;
		//if (cols[k] != c) {std::cout<<"not found"<<std::endl; return zero;}
		return blocks[k];
	}

	//SlMatrix3x3 &find(unsigned int r, unsigned int c);

	void mvmul(const SlVector3 *x, SlVector3 *y);

	void setPrecon(SlVector3 *precon);

	void clear();
};


#endif
