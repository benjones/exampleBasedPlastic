//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// Simple Spring Mass System
// -- Basic Mass points
//
// Primary Author: James F. O'Brien (obrienj@cc.gatech.edu)
//
// (C) Copyright James F. O'Brien, 1995, 1999
// (C) Copyright Georgia Institute of Technology, 1995, 1999
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// RCS Revision History
//
// $Log: smBounds.C,v $
// Revision 1.4  2006/12/23 02:29:22  adamwb
// integrated chris's changes
//
// Revision 1.3  2006/12/22 03:29:05  cwojtan
// Chris added compatibility for WINDOWS operating system,
// because SOMEBODY has been neglecting poor little windows for too long.
//
// Revision 1.2  2006/11/29 01:29:04  adamwb
// checkin
//
// Revision 1.1.1.1  2003/03/17 10:03:24  adamb
// Initial Revision
//
// Revision 3.2  2001/03/12 02:53:58  job
// Changed SmRayHitBound, old version was buggy.
//
// Revision 3.1  2001/03/12 02:52:45  job
// Added Self Traversal
//
// Revision 3.0  1999/03/11  22:29:14  obrienj
// Move from O32 ABI to n32
//
// Revision 2.10  1998/12/20  02:00:04  obrienj
// Added traversal for spheres
// New traversal methods that cache BBox intersections
//
// Revision 2.9  1998/01/19  20:48:08  obrienj
// Fixed bug for building trees with degenerate centers.
//
// Revision 2.8  1997/10/25  22:33:35  obrienj
// Added new versions of SmUpdateTree.
//
// Revision 2.7  1997/09/23  19:32:10  obrienj
// Traversal routines return number of bound overlaps rather than 1/0.
//
// Revision 2.6  1997/05/30  20:09:35  obrienj
// Bound tree member for plane traversal.
// Minor change to check in ray traversal.
//
// Revision 2.5  1997/02/12  20:38:25  obrienj
// Added test of ray against boundBox.
//
// Revision 2.4  1997/02/04  19:06:03  obrienj
// Added tree based ray intersect w/ face set
// Ray traversal added.
//
// Revision 2.3  1996/10/04  04:17:13  obrienj
// Modified to compile under IRIX 6.2 ( mips2 ).
//
// Revision 2.2  1996/09/24  06:34:36  obrienj
// Minor optimizations.
//
// Revision 2.1  1996/09/20  08:29:56  obrienj
// Updatetree calles for SmVertexLists and SmMassPointLists.
//
// Revision 2.0  1996/09/13  03:22:34  obrienj
// Update bound call for rigid bodies with transform pos and vel.
//
// Revision 1.5  1996/09/11  18:08:48  obrienj
// Changed SmOverlapEngine to SmBoundOverlapServer
// Bound trees nolonger stored as a part of the overlapServer.
//
// Revision 1.1  1996/09/10  01:41:39  obrienj
// Initial revision
//
//
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#include "slIO.H"
#include "slBounds.H"
using namespace std;

//-------------------------------------------------------------------
#if 0
inline static
istream &eatChar(char c,istream &buf) {
  char r;
  buf >> r;
  if (strncasecmp(&c,&r,1)) {
    buf.clear(buf.rdstate() | ios::failbit);
  }
  return buf;
}


inline static
istream &eatStr(const char *s,istream &buf) {
  while (*s != '\0') {
    eatChar(*s,buf);
    s++;
  }
  return buf;
}
#endif

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
//class SlBoundBox


istream &operator>>(istream &buf, SlBoundBox &bb) {
  ios::fmtflags orgFlags = buf.setf(ios::skipws);
  eatChar('{',buf);
  buf >> bb.min;
  eatChar(',',buf);
  buf >> bb.max;
  eatChar('}',buf);
  buf.flags(orgFlags);
  return buf;
}
  
ostream &operator<<(ostream &buf, const SlBoundBox &bb) {
  buf << '{';
  buf << bb.min;
  buf << ',';
  buf << bb.max;
  buf << '}';
  return buf;
}
 

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// SlBoundTree


SlBoundTree::SlBoundTree() {
  _numInnerLevels  = 0;
  _innerLevelSizes = NULL;
  _leafLevelSize   = 0;
  _innerLevels = NULL;
  _leafLevel   = NULL;
}


SlBoundTree::SlBoundTree(const SlBoundTree &that){
  _numInnerLevels  = 0;
  _innerLevelSizes = NULL;
  _leafLevelSize   = 0;
  _innerLevels = NULL;
  _leafLevel   = NULL;

  (*this) = that;
}

SlBoundTree::~SlBoundTree(){
  register unsigned short i;

  for(i=0;i<_numInnerLevels;i++) 
    delete [] _innerLevels[i];

  delete [] _innerLevelSizes;
  delete [] _innerLevels;
  delete [] _leafLevel;
}
  


SlBoundTree &SlBoundTree::operator=(const SlBoundTree &that) {
  setNumInnerLevels(that._numInnerLevels);
  setLeafLevelSize (that._leafLevelSize );
  
  register unsigned int   j;
  register unsigned short i;
  for(i=0;i<_numInnerLevels;i++) {
    setInnerLevelSize(i,that._innerLevelSizes[i]);
    for(j=0;j<_innerLevelSizes[i];j++) {
      _innerLevels[i][j] = that._innerLevels[i][j];
    }
  }

  for(j=0;j<_leafLevelSize;j++) {
    _leafLevel[j] = that._leafLevel[j];
  }

  return (*this);
}


void SlBoundTree::updateLevel(unsigned short level) {
  register unsigned int p,pl=_innerLevelSizes[level];

  for(p=0;p<pl;p++) {

    register innerNode &parent = _innerLevels[level][p];

    if (parent.childCount == 0) {
      parent.bound.min =  (DBL_MAX);
      parent.bound.max = -(DBL_MAX);
    }else{
      register unsigned short cLevel = parent.childLevel;
      register unsigned int   cStart = parent.childStart;
      register unsigned int   cStop  = parent.childCount + cStart;
      register unsigned int   cMark  ;
			
      if (cLevel==_numInnerLevels) {
				// leaf level
				parent.bound = _leafLevel[cStart].bound;
				for(cMark=cStart+1;cMark<cStop;cMark++) {
					parent.bound += _leafLevel[cMark].bound;
				}
      }else{
				// inner level
				parent.bound = _innerLevels[cLevel][cStart].bound;
				for(cMark=cStart+1;cMark<cStop;cMark++) {
					parent.bound += _innerLevels[cLevel][cMark].bound;
				}
      }
    }
  }
}
	  

void SlBoundTree::setNumInnerLevels(unsigned short numLevels) {
  register unsigned int i;
  for(i=0;i<_numInnerLevels;i++) {
    delete [] _innerLevels[i];
  }
  delete [] _innerLevels;
  delete [] _innerLevelSizes;

  _innerLevels     = new SlBoundTree::innerNode*[numLevels];
  _innerLevelSizes = new unsigned int[numLevels];
  _numInnerLevels  = numLevels;

  for(i=0;i<_numInnerLevels;i++) {
    _innerLevels    [i] = NULL;
    _innerLevelSizes[i] = 0   ;
  }
}

  
void SlBoundTree::setInnerLevelSize(unsigned short level    , 
				    unsigned int   levelSize) {
  if (level >= _numInnerLevels) {
    cerr << "SlBoundTree::setInnerLevelSize range error.\n" << flush;
    abort();
  }
  delete [] _innerLevels[level];
  
   _innerLevels    [level] = new SlBoundTree::innerNode[levelSize];
   _innerLevelSizes[level] = levelSize;
}
  

void SlBoundTree::setLeafLevelSize(unsigned int levelSize) {
  delete [] _leafLevel;
  _leafLevel     = new SlBoundTree::leafNode [levelSize];
  _leafLevelSize = levelSize;
}
				   

//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------


istream &operator>>(istream &buf, SlBoundTree::innerNode &n) {
  ios::fmtflags orgFlags = buf.setf(ios::skipws);
  eatChar('{',buf);
  buf >> n.bound;
  eatChar(',',buf);
  buf >> n.childLevel;
  eatChar(',',buf);
  buf >> n.childStart;
  eatChar(',',buf);
  buf >> n.childCount;
  eatChar('}',buf);
  buf.flags(orgFlags);
  return buf;
}  


ostream &operator<<(ostream &buf,const SlBoundTree::innerNode &n) {
  buf << '{';
  buf << n.bound;
  buf << ',';
  buf << n.childLevel;
  buf << ',';
  buf << n.childStart;
  buf << ',';
  buf << n.childCount;
  buf << '}';
  return buf;
}  


//----------------------------------------------------------


istream &operator>>(istream &buf, SlBoundTree::leafNode &n) {
  ios::fmtflags orgFlags = buf.setf(ios::skipws);
  eatChar('{',buf);
  buf >> n.bound;
  eatChar(',',buf);
  buf >> n.mapIndex;
  eatChar('}',buf);
  buf.flags(orgFlags);
  return buf;
}  

ostream &operator<<(ostream &buf,const SlBoundTree::leafNode &n) {
  buf << '{';
  buf << n.bound;
  buf << ',';
  buf << n.mapIndex;
  buf << '}';
  return buf;
}


//----------------------------------------------------------


istream &operator>>(istream &buf, SlBoundTree &t) {
  ios::fmtflags orgFlags = buf.setf(ios::skipws);

  unsigned short newNIL = 0,i;
  unsigned int   newLLS = 0,j;
  unsigned int   newILS = 0  ;

  eatStr("{BoundTree:",buf);

  buf >> newNIL;
  t.setNumInnerLevels(newNIL);

  eatChar('[',buf);
  for(i=0;i<newNIL;i++) {
    buf >> newILS;
    t.setInnerLevelSize(i,newILS);
  }
  eatChar(']',buf);
  
  buf >> newLLS;
  t.setLeafLevelSize(newLLS);
  eatChar(';',buf);
  
  for(i=0;i<newNIL;i++) {
    for(j=0;j<t.innerLevelSize(i);j++) {
      buf >> t.constructInnerLevel(i)[j];
    }
  }

  for(j=0;j<t.leafLevelSize();j++) {
    buf >> t.constructLeafLevel()[j];
  }
 
  eatChar('}',buf);

  buf.flags(orgFlags);
  return buf;
}


ostream &operator<<(ostream &buf,const SlBoundTree &t) {

  unsigned short i;
  unsigned int   j;

  buf << "{ BoundTree: ";
  
  buf << t.numInnerLevels();

  buf << "\n[ ";
  
  for(i=0;i<t.numInnerLevels();i++) {
    buf << t.innerLevelSize(i);
    buf << " ";
  }
  buf << " ] ";

  buf << t.leafLevelSize() << ";\n";
  
  for(i=0;i<t.numInnerLevels();i++) {
    for(j=0;j<t.innerLevelSize(i);j++) {
      buf << t(i,j) << ' ';
    }
    buf << '\n';
  }

  for(j=0;j<t.leafLevelSize();j++) {
    buf << t[j] << ' ';
  }
 
  buf << "\n}\n";

  return buf;
}


//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------

SlVector3 *SlCreateCenterList(const SlInt3 *tris, int nt, const SlVector3 *pts) {

  register unsigned int  t;
  register SlVector3 *retList = new SlVector3[nt];

  for(t=0; t<nt; t++) {
    retList[t] = (pts[tris[t][0]] +
		   pts[tris[t][1]] +
		   pts[tris[t][2]]);
  }
  return retList;
}

    
//----------------------------------------------------------
//----------------------------------------------------------

void SlUpdateTree(const SlInt3 *tris, const SlVector3 *pos, const SlVector3 *vel,
									SlBoundTree &btree, double dt) {
  SlBoundTree::leafNode *leaves = btree.beginLeafUpdate();

  register unsigned int l,numL = btree.leafLevelSize();

  for(l=0;l<numL;l++) {
    
    register unsigned int mapIndex = leaves->mapIndex;

    register const SlVector3 &pos0 = pos[tris[mapIndex][0]];
    register const SlVector3 &pos1 = pos[tris[mapIndex][1]];
    register const SlVector3 &pos2 = pos[tris[mapIndex][2]];

    leaves->bound.max = pos0;
    leaves->bound.min = pos0;

    leaves->bound.max.maxSet(pos1);
    leaves->bound.min.minSet(pos1);

    leaves->bound.max.maxSet(pos2);
    leaves->bound.min.minSet(pos2);
    
    if (dt !=0.0) {
      SlVector3 del0 = pos0 + dt * vel[tris[mapIndex][0]];
      SlVector3 del1 = pos1 + dt * vel[tris[mapIndex][1]];
      SlVector3 del2 = pos2 + dt * vel[tris[mapIndex][2]];
     
      leaves->bound.max.maxSet(del0);
      leaves->bound.min.minSet(del0);

      leaves->bound.max.maxSet(del1);
      leaves->bound.min.minSet(del1);

      leaves->bound.max.maxSet(del2);
      leaves->bound.min.minSet(del2);
    }
    leaves++;
  }

  btree.endLeafUpdate();
}

//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------



static 
unsigned int findEnd(short *kList) {
  register unsigned int ret;
  for(ret=0;kList[ret]<0;ret++);
  return ret;
}

//----------------------------------------------------------

static
void findCV(const SlVector3 *cList,
						unsigned int    *iList,
						unsigned int     segB ,
						unsigned int     segE ,
						SlVector3       &cnt  ,
						SlVector3       &var) {
  
  register unsigned int i;

  cnt = 0.0;
  var = 0.0;

  if ((segE-segB) <= 0) return;

  for(i=segB;i<segE;i++) {
    cnt += cList[iList[i]];
  }

  cnt /= (segE-segB);

  for(i=segB;i<segE;i++) {
    register SlVector3 tmp = cnt - cList[iList[i]];
    var[0] += tmp[0] * tmp[0];
    var[1] += tmp[1] * tmp[1];
    var[2] += tmp[2] * tmp[2];
  }
}  

//----------------------------------------------------------

static 
int splitLevel(const SlVector3 *cList,
							 unsigned int    *iList,
							 short           *kList,
							 unsigned int     leafGroupSize,
							 unsigned int     level) {
  
  unsigned int chunkEnd   = 0;
  unsigned int chunkBegin = 0;

  unsigned int b,e;

  SlVector3 cnt;
  SlVector3 var;

  int split;

  int numSplits = 0;

  while(kList[chunkEnd] != 0) {
    
    chunkBegin = chunkEnd;
    chunkEnd   = chunkBegin+1+findEnd(&(kList[chunkBegin+1]));

    if ((chunkEnd - chunkBegin) > leafGroupSize) {

      findCV(cList,iList,chunkBegin,chunkEnd,cnt,var);

      if      ((var[0] > var[1]) && (var[0] >= var[2])) split = 0;
      else if ((var[1] > var[0]) && (var[1] >= var[2])) split = 1;
      else if ((var[2] > var[0]) && (var[2] >= var[1])) split = 2;
      else split = 2;

      if (split >= 0) {
      
	b = chunkBegin;
	e = chunkEnd-1;

	while(b<=e) {
	  register int testB = cnt[split] >= cList[iList[b]][split];
	  register int testE = cnt[split] <= cList[iList[e]][split];
	  if ((!testB) && (!testE)) {
	    register unsigned int tmp = iList[b];
	    iList[b] = iList[e];
	    iList[e] = tmp     ;
	    b++;
	    e--;
	  }else{
	    if (testB) b ++;
	    if (testE) e --;
	  }
	}	  
	
	// If group can not split, chop group in half.
	if ((b == chunkBegin) || (e == chunkEnd-1)) {
	  b = (chunkBegin + chunkEnd-1)/2;
	  e = b -1;
	}

	if ((b>0) && (b<chunkEnd)) {
	  kList[b] = level;
	  numSplits ++;
	}
      }
    }
  }
  return numSplits;
}

//----------------------------------------------------------

static
void buildLeafLevel(unsigned int *iList  ,
										SlBoundTree  &bTree  ,
										unsigned int numItems) {
  bTree.setLeafLevelSize(numItems);
  SlBoundTree::leafNode *leaves = bTree.constructLeafLevel();
  unsigned int i;
  for(i=0;i<numItems;i++) {
    leaves[i].mapIndex = iList[i];
  }
}

//----------------------------------------------------------

static 
void buildInnerLevels(short        *kList,
											SlBoundTree  &bTree,
											unsigned int  numItems,
											unsigned int  numLevel,
											unsigned int  totalSplits) {

  int i,l;

  kList ++;

  register unsigned int childStart;
  register unsigned int nodePlace ;

  unsigned int numBreaks = totalSplits;
  unsigned int listLen   = numItems   ;

  SlBoundTree::innerNode *nodes;

  bTree.setNumInnerLevels(numLevel);
  for(l=numLevel-1;l>=0;l--) {
    bTree.setInnerLevelSize(l,numBreaks);
    nodes = bTree.constructInnerLevel(l);

    nodePlace  = 0;
    childStart = 0;

    for(i=0;i<listLen;i++) {
      if (kList[i]>=0) {
				nodes[nodePlace].childLevel = l+1;
				nodes[nodePlace].childStart = childStart;
				nodes[nodePlace].childCount = i+1-childStart;
				
				if (kList[i]==l) {
					kList[nodePlace] = -1;
					numBreaks--;
				}else{
					kList[nodePlace] = kList[i];
				}
				nodePlace++;
				childStart = i+1;	  
      }
    }
    listLen = nodePlace;
  }
}
  
//----------------------------------------------------------

void SlBuildTree(const SlVector3   *cList, unsigned int numItems    ,
								 SlBoundTree &bTree, unsigned int leafGroupSize) {

  // cList is the list of centers
  // iList is a list of indexs into cList
  // kList break list

  unsigned int *iList = new unsigned int[numItems  ];
  short        *kList = new short       [numItems+1];

  //--------------------------------------------------------

  register unsigned int i;

  //--------------------------------------------------------
    
  for(i=0;i<numItems;i++) {
    iList[i] = i ;
    kList[i] = -1;
  }
  kList[numItems] = 0;

  //--------------------------------------------------------

  unsigned int level       = 1;
  unsigned int totalSplits = 1;
  unsigned int splits         ;

  while(splits = splitLevel(cList,iList,kList,leafGroupSize,level)) {
    totalSplits += splits;
    level++;
  }

  buildLeafLevel(iList,bTree,numItems);

  buildInnerLevels(kList,bTree,numItems,level,totalSplits);

  //--------------------------------------------------------

  delete []kList;
  delete []iList;
}

//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------

// SlBoundOverlapServer

SlBoundOverlapServer::SlBoundOverlapServer() :
  _chunkList   (NULL),
  _freeListHead(NULL),
  _pendListHead(NULL),
  _cacheListHead(NULL),
  _isInUse     (   0),
  _cacheUse    (SlBoundOverlapServer::CNotReady) {
  addChunk();
}

//----------------------------------------------------------

SlBoundOverlapServer::~SlBoundOverlapServer() {
  blowChunks();
}

//----------------------------------------------------------

int SlBoundOverlapServer::isReady() const {
  return (!isInUse());
}

int SlBoundOverlapServer::isInUse() const {
  return _isInUse;
}

//----------------------------------------------------------

void SlBoundOverlapServer::addChunk() {
  register SlBoundOverlapServer::overlapChunk *newChunk 
    = new SlBoundOverlapServer::overlapChunk;

  newChunk->chunk[newChunk->size-1].next = _freeListHead;
  _freeListHead = &(newChunk->chunk[0]);
  
  newChunk->next = _chunkList;
  _chunkList = newChunk;
}
 
//----------------------------------------------------------

void SlBoundOverlapServer::blowChunks(){
  while(_chunkList!=NULL) {
    register SlBoundOverlapServer::overlapChunk *tmp = _chunkList;
    _chunkList = _chunkList->next;
    delete tmp;
  }

  _pendListHead  = NULL;
  _freeListHead  = NULL;
  _cacheListHead = NULL;
}

//----------------------------------------------------------

void SlBoundOverlapServer::dumpPending() {
  if (_pendListHead != NULL) {
    SlBoundOverlapServer::overlap *tail = _pendListHead;
    while (tail->next != NULL) {
      tail = tail->next;
    }
    tail -> next = _freeListHead;
    _freeListHead = _pendListHead;
    _pendListHead = NULL;
  }
}
    
//----------------------------------------------------------

void SlBoundOverlapServer::clearCache() {
  if (_cacheListHead != NULL) {
    SlBoundOverlapServer::overlap *tail = _cacheListHead;
    while (tail->next != NULL) {
      tail = tail->next;
    }
    tail -> next = _freeListHead;
    _freeListHead = _cacheListHead;
    _cacheListHead = NULL;
  }
}

//----------------------------------------------------------
//----------------------------------------------------------


// Dummy functions to send to doTraversal to that it caches leaf hits


static 
int LeafFlag(const SlBoundOverlapServer  &,
						 const SlBoundTree::leafNode &, void *,
						 const SlBoundTree::leafNode &, void *) { 
  abort(); return 0; 
}

static 
int SelfFlag(const SlBoundOverlapServer  &,
						 const SlBoundTree::leafNode &, 
						 const SlBoundTree::leafNode &, void *) { 
  abort(); return 0; 
}


//----------------------------------------------------------

int SlBoundOverlapServer::cacheTraversal(const SlBoundTree &treeOne, 
																				 const SlBoundTree &treeTwo) {
  _cacheUse = CReady;
  int ret = doTraversal(LeafFlag,treeOne,NULL,treeTwo,NULL);
  _isInUse = 1;
  return ret;
}
		   
int SlBoundOverlapServer::cacheSelfTraversal(const SlBoundTree &tree) {
  _cacheUse = CReady;
  int ret = doSelfTraversal(SelfFlag,tree,NULL);
  _isInUse = 1;
  return ret;
}
		   
//----------------------------------------------------------

void SlBoundOverlapServer::doneWithCache() {
  if ((_cacheUse == CNotReady) || (!_isInUse)) {
    cerr << "SlBoundOverlapServer::doneWithCache :: Error: Cache was not ready...\n";
    abort();
  }

  _cacheUse = CNotReady;
  _isInUse  = 0;
  
  clearCache();
}

//----------------------------------------------------------

int SlBoundOverlapServer::traverseCache(SlLeafOverlapCB cb, 
					const SlBoundTree &treeOne, void *uprmOne, 
					const SlBoundTree &treeTwo, void *uprmTwo) {

  if ((_cacheUse != CReady) || (!_isInUse)) {
    cerr << "SlBoundOverlapServer::traverseCache :: Error: Cache was not ready...\n";
    abort();
  }

  overlap ** pos = &_cacheListHead;

  int result = 0;

  while ((pos != NULL) && ((*pos) != NULL)) {
    register unsigned short l1 = (*pos)->level1;
    register unsigned short l2 = (*pos)->level2;
    register unsigned int   n1 = (*pos)-> node1;
    register unsigned int   n2 = (*pos)-> node2;
    
    if ((treeOne.numInnerLevels() != l1) ||
				(treeTwo.numInnerLevels() != l2)) {
      cerr << "SlBoundOverlapServer::traverseCache :: Error bad trees.\n";
      abort();
    }
		
    if (SlCheckInterfere(treeOne[n1].bound,treeTwo[n2].bound)) {
      result ++;
      if (cb==NULL) {
				pos = NULL;
      }else{
				int cbrv = (*cb)(*this,treeOne[n1],uprmOne,treeTwo[n2],uprmTwo);
				if (cbrv <  0) removeFromCache(pos);
				if (cbrv == 0) pos = NULL;
      }
    }
    if (pos != NULL) pos = &((*pos)->next);
  }
  return result;
}
     
//----------------------------------------------------------

int SlBoundOverlapServer::traverseSelfCache(SlSelfOverlapCB cb, 
					    const SlBoundTree &tree, void *uprm) {

  if ((_cacheUse != CReady) || (!_isInUse)) {
    cerr << "SlBoundOverlapServer::traverseSelfCache :: Error: Cache was not ready...\n";
    abort();
  }

  overlap ** pos = &_cacheListHead;

  int result = 0;

  while ((pos != NULL) && ((*pos) != NULL)) {
    register unsigned short l1 = (*pos)->level1;
    register unsigned short l2 = (*pos)->level2;
    register unsigned int   n1 = (*pos)-> node1;
    register unsigned int   n2 = (*pos)-> node2;
    
    if ((tree.numInnerLevels() != l1) ||
				(tree.numInnerLevels() != l2)) {
      cerr << "SlBoundOverlapServer::traverseSelfCache :: Error bad trees.\n";
      abort();
    }

    if (SlCheckInterfere(tree[n1].bound,tree[n2].bound)) {
      result ++;
      if (cb==NULL) {
	pos = NULL;
      }else{
	int cbrv = (*cb)(*this,tree[n1],tree[n2],uprm);
	if (cbrv <  0) removeFromCache(pos);
	if (cbrv == 0) pos = NULL;
      }
    }
    if (pos != NULL) pos = &((*pos)->next);
  }
  return result;
}
     

//-------------------------------------------------------------------
//-------------------------------------------------------------------


// This macro is used to add items to pending list.
// it assumes that a,aMax,lA,b,limB,lB are set.

#define ADD_TOO_PEND_LIST {			\
  for(a=aMin;a<aMax;a++) {			\
    for(b=bMin;b<bMax;b++) {			\
     pushPend(1);				\
      _pendListHead->level1 = lA;		\
      _pendListHead->level2 = lB;		\
      _pendListHead-> node1 =  a;		\
      _pendListHead-> node2 =  b;		\
    }						\
  }						\
}

//----------------------------------------------------------
//----------------------------------------------------------

int SlBoundOverlapServer::doTraversal(SlLeafOverlapCB cb, 
				      const SlBoundTree &treeOne, void *uprmOne, 
				      const SlBoundTree &treeTwo, void *uprmTwo) {

  int result = 0;
  register unsigned int a,aMax,aMin,lA;
  register unsigned int b,bMax,bMin,lB;
  //--------------------------------------------------------
	
  if (!isReady()) {
    cerr << "SlBoundOverlapServer::doTraversal :: Error: Traversal attempted on not ready server.\n";
    abort();
  }
	
  _isInUse = 1;

  //--------------------------------------------------------
  // Fisrt initalize pending list
	
  if (treeOne.numInnerLevels() == 0) {
    aMax = treeOne.leafLevelSize();
  }else{
    aMax = treeOne.innerLevelSize(0);
  }
  if (treeTwo.numInnerLevels() == 0) {
    bMax = treeTwo.leafLevelSize();
  }else{
    bMax = treeTwo.innerLevelSize(0);
  }
  aMin = 0; lA = 0;
  bMin = 0; lB = 0;
	
  ADD_TOO_PEND_LIST;
  
  //--------------------------------------------------------
  // Loop until the pending list is empty.
  
  while(popPend()) {
    
    register unsigned short l1 = _freeListHead->level1;
    register unsigned short l2 = _freeListHead->level2;
    register unsigned int   n1 = _freeListHead-> node1;
    register unsigned int   n2 = _freeListHead-> node2;
    if (treeOne.numInnerLevels() == l1) {
      if (treeTwo.numInnerLevels() == l2) {  
				//--------------------------------------------------
				// Leaf,Leaf
				if (SlCheckInterfere(treeOne[n1].bound,treeTwo[n2].bound)) {
					result ++;
					if (cb == LeafFlag) {
						addToCache();
					}else{
						if ((cb==NULL)||
								(!((*cb)(*this,treeOne[n1],uprmOne,treeTwo[n2],uprmTwo)))) {
							dumpPending();
						}
					}
				}
      }else{
				//--------------------------------------------------
				// Leaf,Inner
				if (SlCheckInterfere(treeOne[n1].bound,treeTwo(l2,n2).bound)) {
					lA   = l1; 
					aMin = n1; 
					aMax = n1+1;
					lB   = treeTwo(l2,n2).childLevel     ;
					bMin = treeTwo(l2,n2).childStart     ;
					bMax = treeTwo(l2,n2).childCount+bMin;
					ADD_TOO_PEND_LIST;
				}
      }
    }else{
      if (treeTwo.numInnerLevels() == l2) {
				//--------------------------------------------------
				// Inner,Leaf
				if (SlCheckInterfere(treeOne(l1,n1).bound,treeTwo[n2].bound)) {
					lA   = treeOne(l1,n1).childLevel     ;
					aMin = treeOne(l1,n1).childStart     ;
					aMax = treeOne(l1,n1).childCount+aMin;
					lB   = l2; 
					bMin = n2; 
					bMax = n2+1;
					ADD_TOO_PEND_LIST;
				}
      }else{
				//--------------------------------------------------
				// Inner,Inner
				if (SlCheckInterfere(treeOne(l1,n1).bound,treeTwo(l2,n2).bound)) {
					lA   = treeOne(l1,n1).childLevel     ;
					aMin = treeOne(l1,n1).childStart     ;
					aMax = treeOne(l1,n1).childCount+aMin;
					lB   = treeTwo(l2,n2).childLevel     ;
					bMin = treeTwo(l2,n2).childStart     ;
					bMax = treeTwo(l2,n2).childCount+bMin;
					ADD_TOO_PEND_LIST;
				}
      }
    }
  }
  _isInUse = 0;
  return result;
}
  

//----------------------------------------------------------
//----------------------------------------------------------

int SlBoundOverlapServer::doSelfTraversal(SlSelfOverlapCB cb, 
																					const SlBoundTree &tree, void *uprm) {
  int result = 0;
  register unsigned int a,aMax,aMin,lA;
  register unsigned int b,bMax,bMin,lB;
  //--------------------------------------------------------

  if (!isReady()) {
    cerr << "SlBoundOverlapServer::doSelfTraversal :: Error: Traversal attempted on not ready server.\n";
    abort();
  }

  _isInUse = 1;

  //--------------------------------------------------------
  // Fisrt initalize pending list

  if (tree.numInnerLevels() == 0) {
    aMax = tree.leafLevelSize();
    bMax = tree.leafLevelSize();
  }else{
    aMax = tree.innerLevelSize(0);
    bMax = tree.innerLevelSize(0);
  }
  aMin = 0; lA = 0;
  bMin = 0; lB = 0;

  ADD_TOO_PEND_LIST;
  
  //--------------------------------------------------------
  // Loop until the pending list is empty.
  
  while(popPend()) {
    
    register unsigned short l1 = _freeListHead->level1;
    register unsigned short l2 = _freeListHead->level2;
    register unsigned int   n1 = _freeListHead-> node1;
    register unsigned int   n2 = _freeListHead-> node2;

    if ((l1 != l2) || (n2>=n1)) {
      if (tree.numInnerLevels() == l1) {
				if (tree.numInnerLevels() == l2) {  
					//--------------------------------------------------
					// Leaf,Leaf
					if (n1 != n2) {
						if (SlCheckInterfere(tree[n1].bound,tree[n2].bound)) {
							result ++;
							if (cb == SelfFlag) {
								addToCache();
							}else{
								if ((cb==NULL)||
										(!((*cb)(*this,tree[n1],tree[n2],uprm)))) {
									dumpPending();
								}
							}
						}
					}
				}else{
					//--------------------------------------------------
					// Leaf,Inner
					if (SlCheckInterfere(tree[n1].bound,tree(l2,n2).bound)) {
						lA   = l1; 
						aMin = n1; 
						aMax = n1+1;
						lB   = tree(l2,n2).childLevel     ;
						bMin = tree(l2,n2).childStart     ;
						bMax = tree(l2,n2).childCount+bMin;
						ADD_TOO_PEND_LIST;
					}
				}
      }else{
				if (tree.numInnerLevels() == l2) {
					//--------------------------------------------------
					// Inner,Leaf
					if (SlCheckInterfere(tree(l1,n1).bound,tree[n2].bound)) {
						lA   = tree(l1,n1).childLevel     ;
						aMin = tree(l1,n1).childStart     ;
						aMax = tree(l1,n1).childCount+aMin;
						lB   = l2; 
						bMin = n2; 
						bMax = n2+1;
						ADD_TOO_PEND_LIST;
					}
				}else{
					//--------------------------------------------------
					// Inner,Inner
					if (SlCheckInterfere(tree(l1,n1).bound,tree(l2,n2).bound)) {
						lA   = tree(l1,n1).childLevel     ;
						aMin = tree(l1,n1).childStart     ;
						aMax = tree(l1,n1).childCount+aMin;
						lB   = tree(l2,n2).childLevel     ;
						bMin = tree(l2,n2).childStart     ;
						bMax = tree(l2,n2).childCount+bMin;
						ADD_TOO_PEND_LIST;
					}
				}
      }
    }
  }  
  _isInUse = 0;
  return result;
}


