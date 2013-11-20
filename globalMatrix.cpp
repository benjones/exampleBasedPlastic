#include "fem.H"
#include "globalMatrix.H"

void GlobalMatrix::mvmul(const SlVector3 *x, SlVector3 *y) {
	SlVector3 *tmp = y;
	for(unsigned int i = 0; i < dim; i++, tmp++) {
		(*tmp) = 0.0;
	}

	for (unsigned int r=0; r<dim; r++) {
		unsigned int *c = cols + colIndices[r];
		SlMatrix3x3 *m = blocks + colIndices[r];
		for (unsigned int j=0; j<nblocks[r]; j++, c++, m++) {
			if ((*c) == r) {
				y[r] += (*m) * (x[r]);
			} else {
				y[r] += (*m) * (x[*c]);
				y[*c] += transmult(*m,x[r]);
			}
		}
	}
}

void GlobalMatrix::output() {
		for (unsigned int i=0; i<totalblocks; i++) {
			std::cout<<blocks[i]<<std::endl;
		}
	}

GlobalMatrix::GlobalMatrix(unsigned int nvertices, unsigned int ntets, const Tet *tets) {
		dim = nvertices;
		colIndices = new unsigned int[dim+1];
		nblocks = new unsigned int[dim];
		nonzeros = 0;
		
		std::vector<std::vector<unsigned int> > vneighbors;
		vneighbors.resize(dim);

		Tet const *tet = tets;
		for (unsigned int t=0; t<ntets; t++, tet++) {
			for (unsigned int i=0; i<4; i++) {
				unsigned int vi = (*tet)(i);
				for (unsigned int j=i; j<4; j++) {
					unsigned int vj = (*tet)(j);
					bool found = false;
					for (unsigned int k=0; k<vneighbors[vi].size() && !found; k++)
						if (vneighbors[vi][k] == vj) found = true;
					if (!found) {
						vneighbors[vi].push_back(vj);
						if (vi != vj) vneighbors[vj].push_back(vi);
					}
				}
			}
		}

		totalblocks = 0;
		for (unsigned int i=0; i<vneighbors.size(); i++) {
			std::sort(vneighbors[i].begin(), vneighbors[i].end());
			unsigned int rowblocks = 0;
			for (unsigned int j=0; j<vneighbors[i].size(); j++) {
				if (i >= vneighbors[i][j]) {
					totalblocks++;
					rowblocks++;
				}
			}
			nblocks[i] = rowblocks;
		}

		blocks = new SlMatrix3x3[totalblocks]();
		cols = new unsigned int[totalblocks];
		colIndices[0] = 0;


		for (unsigned int i=0; i<vneighbors.size(); i++) {
			if (vneighbors[i].size() == 0) {std::cout<<"panic"<<std::endl;}
			for (unsigned int j=0; j<vneighbors[i].size(); j++) {
				if (vneighbors[i][j] <= i) {
					if (i < nvertices) {
						cols[colIndices[i]+j] = vneighbors[i][j];
						nonzeros += 9;
					}
					if (j+1 >= vneighbors[i].size() || vneighbors[i][j+1] > i) {
						colIndices[i+1] = colIndices[i] + j + 1;
					}
				}
			}
		}
		
		vneighbors.resize(0);
		zero = 0.0;
}

void GlobalMatrix::setPrecon(SlVector3 *precon) {
	for (unsigned int r=0; r<dim; r++) {
		unsigned int *c = cols + colIndices[r];
		SlMatrix3x3 *m = blocks + colIndices[r];
		for (unsigned int j=0; j<nblocks[r]; j++, c++, m++) {		
			if ((*c) == r) {
				precon[r].set(1.0/((*m)(0,0)), 1.0/((*m)(1,1)), 1.0/((*m)(2,2)));
			}
		}
	}
}

void GlobalMatrix::clear() {
	for (unsigned int i=0; i<totalblocks; i++) blocks[i] = 0.0;
}
