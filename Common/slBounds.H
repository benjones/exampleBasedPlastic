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
//
// $Log: smBounds.H,v $
// Revision 1.3  2006/09/22 21:24:29  adamwb
// Revisions for the goop project
//
// Revision 1.2  2004/03/10 07:41:38  adamb
// fixing warnings in shell simulator
//
// Revision 1.1.1.1  2003/03/17 10:03:25  adamb
// Initial Revision
//
// Revision 3.1  2001/03/12 02:52:45  job
// Added Self Traversal
//
// Revision 3.0  1999/03/11  22:29:14  obrienj
// Move from O32 ABI to n32
//
// Revision 2.8  1998/12/20  02:00:04  obrienj
// Added traversal for spheres
// New traversal methods that cache BBox intersections
//
// Revision 2.7  1998/02/04  22:08:47  obrienj
// Changed calls to minComp and maxComp to min and max.
//
// Revision 2.6  1997/10/25  22:32:56  obrienj
// Added new versions of SmUpdateTree(...)
//
// Revision 2.5  1997/09/23  19:32:47  obrienj
// Traversal routines return number of bound overlaps rather than 1/0.
//
// Revision 2.4  1997/05/30  20:07:39  obrienj
// Add plane vs bound and half space vs bound tests
// Add plane traversal of bound tree.
// Member function for access of the root bound.
//
// Revision 2.3  1997/02/12  20:38:51  obrienj
// Added test of ray against boundBox.
//
// Revision 2.2  1997/02/04  19:06:03  obrienj
// Added tree based ray intersect w/ face set
// Ray traversal added.
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
// Revision 1.1  1996/09/10  01:43:09  obrienj
// Initial revision
//
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------


#ifndef SlBoundsDEFINED
#define SlBoundsDEFINED

//-------------------------------------------------------------------
//-------------------------------------------------------------------

#include <cstdlib>
#include <limits.h>
#include <iostream>
#include <cfloat>
#include "slVector.H"
#include "slUtil.H"

//-------------------------------------------------------------------
// Classes from this header:

class SlBoundBox;
class SlBoundTree;
class SlBoundCashe;

class SlBoundOverlapServer;

//-------------------------------------------------------------------

#include "slVector.H"

//-------------------------------------------------------------------

#define SlDefaultLeafGroupSize ( 8    )
#define SlDefaultOverChunkSize ( 1000 )

//-------------------------------------------------------------------
//-------------------------------------------------------------------

inline int SlCheckInterfere(const SlBoundBox &a,const SlBoundBox &b); 
// returns true if 'a' overlaps the area of 'b', else false.

//-------------------------------------------------------------------
 
SlVector3 *SlCreateCenterList(const SlInt3 *tris, int nt, const SlVector3 *pts);
// These routines create a list of vectors froim their first
// argument. This list of vectors holds the center of each 
// object in the list (first arg). This list can then be used
// to build a SmBoundTree.
// The user should delete[] the return value when finished
// using it (ie after building the tree).


void SlBuildTree(const SlVector3   *cList, unsigned int numItems,
								 SlBoundTree &bTree, unsigned int leafGroupSize = SlDefaultLeafGroupSize);
// Builds a boundtree which does a spatial sort of
// its leaves based on the vector locations in cList.
// Does not modify cList.
// Bounds are uninitalized, see SmUpdateTree.

//-------------------------------------------------------------------

void SlUpdateTree(const SlInt3 *tris, const SlVector3 *pos,
									const SlVector3 *vel, SlBoundTree &bTree, double dt);
// Updates the bounds in a tree after the points in pList|vList have moved.
// Also inializes the bounds for a newly constructed tree.
// The bound boxes have extra "slop" space added to them
// based on deltaT and the velocity of a point.  Depending
// on the order of the calls with respect to the integration,
// a positive or negative deltaT should be used.

//-------------------------------------------------------------------

class SlBoundBox {
  // An axis-aligned bounding box.

public:

  //----------------------------------------------------------
  // Constructors...
  
  inline SlBoundBox(                      );                     
  inline SlBoundBox(const SlBoundBox &that);                     
  inline SlBoundBox(const SlVector3 &mins,const SlVector3 &maxs);
    
  //----------------------------------------------------------
  // Overloaded operators
  
  inline SlBoundBox &operator= (const SlBoundBox &that)     ;
  inline SlBoundBox &operator+=(const SlBoundBox &that)     ;
  inline int         operator==(const SlBoundBox &that)const;

  //----------------------------------------------------------
  // Accessable data members.

  SlVector3 min;
  SlVector3 max;
};


inline SlBoundBox boundUnion(const SlBoundBox &a, const SlBoundBox &b);
inline SlBoundBox operator+ (const SlBoundBox &a, const SlBoundBox &b);
// Union operator

std::istream &operator>>(std::istream &buf,      SlBoundBox &bb);
std::ostream &operator<<(std::ostream &buf,const SlBoundBox &bb);


//-------------------------------------------------------------------
//-------------------------------------------------------------------

class SlBoundTree {
  // Stores a bounding tree built from a list of
  // objects.  This class does not actually store
  // the triangles (or whatever) so it only has
  // meaning when considered with the object from
  // which it was build.
  // 
  // A tree is of depth "n" if it has n+1 levels.
  // There are n inner levels and 1 leaf level.
  // The inner levels are levels 0..(n-1) and
  // the leaf level is level n. The root level 
  // is inner level 0.
  //
  // Inner nodes hold a child level field, a
  // start child field, and a numChildren field.
  //
  // Leaf nodes holds and index which holds the
  // position index in a list of the things from
  // which the BoundTree was built.

public:

  //-------------------------------------------------------------------
  // internal classes.

  class innerNode {
  public:
    SlBoundBox     bound;
    unsigned short childLevel;
    unsigned int   childStart;
    unsigned short childCount;

    // Operators:
    inline innerNode           (                     ); 
    inline innerNode           (const innerNode &that); 
    inline innerNode &operator=(const innerNode &that); 
  };

  class leafNode {
  public:
    SlBoundBox     bound;
    unsigned int   mapIndex;

    // Operators:
    inline leafNode           (                    ); 
    inline leafNode           (const leafNode &that); 
    inline leafNode &operator=(const leafNode &that); 
  };

  //----------------------------------------------------------
  // Constructors & destructor...

  SlBoundTree(                       );
  SlBoundTree(const SlBoundTree &that);
 ~SlBoundTree(                       );
  
  //----------------------------------------------------------
  // Size querries
  
  inline unsigned short numInnerLevels()                   const; 

  inline unsigned int innerLevelSize(unsigned short level) const; 
  inline unsigned int  leafLevelSize(                    ) const; 

  //----------------------------------------------------------
  // Access operators
  
  inline const  leafNode &operator[](unsigned int   index) const; 
  inline const innerNode &operator()(unsigned short level, 
																		 unsigned int   index) const; 

  //----------------------------------------------------------

  SlBoundTree &operator=(const SlBoundTree &that);

  //----------------------------------------------------------
  // Root bound access

  inline const SlBoundBox rootBound() const;
  // If there is no root node for some reason, then 
  // an empty bound ((0,0,0),(0,0,0)) will be returned.


  //----------------------------------------------------------
  // Level Access

  inline const  leafNode * leafLevel(                    ) const; 
  inline const innerNode *innerLevel(unsigned short level) const; 

  //----------------------------------------------------------
  // Level update

  inline void updateTree ();              
  inline void updateAbove(unsigned short level); 
  void        updateLevel(unsigned short level);  
  // The update tree will update all inner nodes.
  // updateAbove will update all innerLevels starting
  // at level 'level' to level 0.
  // updateLevel will only update the innerlevel 'level'.

  inline leafNode *beginLeafUpdate()                 ; 
  inline void        endLeafUpdate(){ updateTree(); }; 


//protected:

  //----------------------------------------------------------
  // Size operations
  // Using these calls invalidate the tree. They
  // are for use durring construction or the tree.

  void setNumInnerLevels(unsigned short numLevels);
  
  void setInnerLevelSize(unsigned short level, unsigned int levelSize);
  void  setLeafLevelSize(                      unsigned int levelSize);
  
  inline leafNode  *constructLeafLevel (                    );
  inline innerNode *constructInnerLevel(unsigned short level);
  
protected:

  //----------------------------------------------------------
  // Data members

  unsigned short  _numInnerLevels;
  unsigned int   *_innerLevelSizes;
  unsigned int    _leafLevelSize;

  innerNode **_innerLevels;
  leafNode   *_leafLevel  ;
};

//-------------------------------------------------------------------

std::istream &operator>>(std::istream &buf,      SlBoundTree::innerNode &n);
std::ostream &operator<<(std::ostream &buf,const SlBoundTree::innerNode &n);

std::istream &operator>>(std::istream &buf,      SlBoundTree::leafNode &n);
std::ostream &operator<<(std::ostream &buf,const SlBoundTree::leafNode &n);

std::istream &operator>>(std::istream &buf,      SlBoundTree       &t);
std::ostream &operator<<(std::ostream &buf,const SlBoundTree       &t);



//-------------------------------------------------------------------
//-------------------------------------------------------------------

typedef int (* SlLeafOverlapCB)(const SlBoundOverlapServer  &oe,
				const SlBoundTree::leafNode &indexOne, void *uprmOne,
				const SlBoundTree::leafNode &indexTwo, void *uprmTwo);


typedef int (* SlSelfOverlapCB)(const SlBoundOverlapServer  &oe,
				const SlBoundTree::leafNode &indexOne, 
				const SlBoundTree::leafNode &indexTwo, void *uprm);

class SlBoundOverlapServer {
  // This class provides a set of methods for preforming travrsals
  // of bound tree pairs searching for overlaping bound boxes.
  // Note that while both trees may be the SAME tree (ie to find self
  // overlap) this is not an efficent way to do things.
  // This object uses an internal structure to track
  // its recursion.  For effency reasons, the structures
  // are alowcated in blocks so that the server can avoid the
  // overhead of memory allocation durring traversal.  If 
  // many traversals will be done (ie once a time step) it
  // it more effecient to reuse the server rather than constantly
  // creating and disposing of them.
  
public:

  //----------------------------------------------------------
  // Constructors...
  
  SlBoundOverlapServer();
 ~SlBoundOverlapServer();
 
  //----------------------------------------------------------
  // Traversal related calls

  int isReady() const;
  // return true if the overlapServer is in a state valid to
  // start a traversal.

  int isInUse() const;
  // returns true if the overlapServer is in the midst
  // of doing a traversal.


  //----------------------------------------------------------
  // Actuall traversal calls

  int doTraversal(SlLeafOverlapCB cb, 
		  const SlBoundTree &treeOne, void *uprmOne, 
		  const SlBoundTree &treeTwo, void *uprmTwo);
  // Does a traversal looking for overlaping bound boxes at the leaf
  // level.  Each time overlapping leaves are found, cb will be called.  
  // If cb returns  true, the traversal will continue.  The return
  // from this call will be the number of overlaping bounds  found
  // at the leaf levels. If cb is NULL then the traversal will be
  // stoped at the first overlaping leaf pair.

  //----------------------------------------------------------

  int doSelfTraversal(SlSelfOverlapCB cb, 
		      const SlBoundTree &tree, void *uprmOne); 
  // Does a self traversal looking for overlaping bound boxes at the
  // leaf level.  Each time overlapping leaves are found, cb will be
  // called.  If cb returns true, the traversal will continue.  The
  // return from this call will be the number of overlaping bounds
  // found at the leaf levels. If cb is NULL then the traversal will
  // be stoped at the first overlaping leaf pair.

  //----------------------------------------------------------
  //  Traversal cache 

  int cacheTraversal      (const SlBoundTree &treeOne,
													 const SlBoundTree &treeTwo);
  int cacheSelfTraversal  (const SlBoundTree &tree);
			
  // These calls are like the doTraversal version except they bult an
  // internal list of the overlaping bounds.  Once a traversal has
  // been done, the appropriate call to traverseCache can be used many
  // times.  The return value and other behavior is the same as the
  // equivalent doTraversal call.

  int traverseCache(SlLeafOverlapCB cb, 
										const SlBoundTree &treeOne, void *uprmOne, 
										const SlBoundTree &treeTwo, void *uprmTwo);
  int traverseSelfCache(SlSelfOverlapCB cb, 
												const SlBoundTree &tree, void *uprm);

  // Like an equivalent call to doTraversal, but on the leaf results
  // of a previous call to cacheTraversal.  The return value of the
  // callback is interpreted slightly differently.  As with the
  // doTraversal method, a positive value means continue and a zero
  // value will stop.  Additionaly, a negative value can be returned
  // to indicate that the traversaal should continue, but that the
  // current pair could be removed from the cache.

  void doneWithCache();
  // Clears out the caches traversal and prepares the server for a
  // call to cacheTraversal or doTraversal.


  //----------------------------------------------------------

protected:

  class overlap {
  public:
    unsigned short level1,level2;
    unsigned int    node1, node2;
    SlBoundOverlapServer::overlap *next;
  };

  class overlapChunk {
  public:
    inline  overlapChunk(unsigned int size = SlDefaultOverChunkSize);
    inline ~overlapChunk(                                          );
    SlBoundOverlapServer::overlapChunk *next ;
    SlBoundOverlapServer::overlap      *chunk;
    unsigned int                   size ;
  };
    
private:

  SlBoundOverlapServer::overlapChunk * _chunkList;
  SlBoundOverlapServer::overlap      * _freeListHead;
  SlBoundOverlapServer::overlap      * _pendListHead;
  SlBoundOverlapServer::overlap      * _cacheListHead;
  
  int  _isInUse;

  enum CacheUse { CNotReady , CReady , CReadyPlane , CReadyRay , CReadySphere } _cacheUse;

protected:
  void addChunk();
  void blowChunks();  // hehe, aka freeChunks

  void dumpPending();

  inline int pushPend(int allowGrowth = 1);
  inline int  popPend(                   );
  
  inline int addToCache(int allowGrowth = 1);
  inline void removeFromCache(overlap **node);

  void clearCache();
};

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------

// Inline member functions

inline
int SlCheckInterfere(const SlBoundBox &a,const SlBoundBox &b){
  return !(((a.min[0] > b.max[0]) || (b.min[0] > a.max[0])) ||
					 ((a.min[1] > b.max[1]) || (b.min[1] > a.max[1])) ||
					 ((a.min[2] > b.max[2]) || (b.min[2] > a.max[2])));
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

inline 
SlBoundBox::SlBoundBox() :
  min( (DBL_MAX)),
  max(-(DBL_MAX)) {
}

inline 
SlBoundBox::SlBoundBox(const SlBoundBox &that) :
  min(that.min),
  max(that.max) {
}

inline 
SlBoundBox::SlBoundBox(const SlVector3 &mins,const SlVector3 &maxs) :
  min(mins),
  max(maxs) {
}

inline 
SlBoundBox &SlBoundBox::operator=(const SlBoundBox &that) {
  min = that.min;
  max = that.max;
  return (*this);
}

inline
SlBoundBox &SlBoundBox::operator+=(const SlBoundBox &that) {
  min.minSet(that.min);
  max.maxSet(that.max);
  return (*this);
}


inline
int SlBoundBox::operator==(const SlBoundBox &that) const {
  return (min == that.min) && (max == that.max);
}
    
//-------------------------------------------------------------------


inline
SlBoundBox boundUnion(const SlBoundBox &a, const SlBoundBox &b){
  return SlBoundBox(min(a.min,b.min),
		    max(a.max,b.max));
}

inline
SlBoundBox operator+(const SlBoundBox &a, const SlBoundBox &b) {
  return boundUnion(a,b);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------

inline
SlBoundTree::innerNode::innerNode() :
  childCount(0) {
}


inline
SlBoundTree::innerNode::innerNode(const SlBoundTree::innerNode &that) :
  bound(that.bound),
  childLevel(that.childLevel),
  childStart(that.childStart),
  childCount(that.childCount) {
}
 

inline
SlBoundTree::innerNode &
SlBoundTree::innerNode::operator=(const SlBoundTree::innerNode &that) {
  bound = that.bound;
  childLevel = that.childLevel;
  childCount = that.childCount;
  childStart = that.childStart;
  return (*this);
}


inline
SlBoundTree::leafNode::leafNode() {
}


inline   
SlBoundTree::leafNode::leafNode(const SlBoundTree::leafNode &that) :
  bound(that.bound),
  mapIndex(that.mapIndex) {
}


inline   
SlBoundTree::leafNode &
SlBoundTree::leafNode::operator=(const SlBoundTree::leafNode &that) {
  bound = that.bound;
  mapIndex = that.mapIndex;
  return (*this);
}


inline 
unsigned short SlBoundTree::numInnerLevels() const {
  return _numInnerLevels;
}
  

inline 
unsigned int SlBoundTree::innerLevelSize(unsigned short level) const {
  return _innerLevelSizes[level];
}
  

inline
unsigned int SlBoundTree::leafLevelSize() const {
  return _leafLevelSize;
}

inline 
const SlBoundTree::leafNode &SlBoundTree::operator[](unsigned int index) const {
  return _leafLevel[index];
}


inline
const SlBoundTree::innerNode &SlBoundTree::operator()(unsigned short level, 
						      unsigned int   index) const {
  return _innerLevels[level][index];
}

inline 
const SlBoundBox SlBoundTree::rootBound() const {
  if (numInnerLevels() == 0) {
    if (leafLevelSize() == 1) {
      return (*this)[0].bound;
    }
  }else{
    if (innerLevelSize(0) == 1) {
      return (*this)(0,0).bound;
    }
  }
  return SlBoundBox(SlVector3(0,0,0),SlVector3(0,0,0));
}


inline
const SlBoundTree::leafNode *SlBoundTree::leafLevel() const {
  return _leafLevel;
}

inline       
const SlBoundTree::innerNode *SlBoundTree::innerLevel(unsigned short level) const {
  return _innerLevels[level];
}


inline
void SlBoundTree::updateTree() {
  updateAbove(_numInnerLevels-1);
}

inline    
void SlBoundTree::updateAbove(unsigned short level) {
  register unsigned short i;
  for(i=level;;i--) {
    updateLevel(i);
    if (i==0) break;
  }
}


inline
SlBoundTree::leafNode *SlBoundTree::beginLeafUpdate() {
  return _leafLevel;
}

inline
SlBoundTree::leafNode *SlBoundTree::constructLeafLevel() {
  return _leafLevel;
}

inline
SlBoundTree::innerNode *SlBoundTree::constructInnerLevel(unsigned short level) {
  return _innerLevels[level];
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

inline SlBoundOverlapServer::overlapChunk::overlapChunk(unsigned int size) {
  this->next  = NULL;
  this->size  = size;
  chunk = new SlBoundOverlapServer::overlap[size];

  register unsigned int i;

  for(i=1;i<size;i++) {
    chunk[i-1].next = &(chunk[i]);
  }
  chunk[size-1].next = NULL;
}

inline SlBoundOverlapServer::overlapChunk::~overlapChunk() {
  delete [] chunk;
}


int SlBoundOverlapServer::pushPend(int allowGrowth) {
  register int ret;

  if (( _freeListHead == NULL ) && (allowGrowth)) {
		std::cout<<"adding chunk"<<std::endl;
    addChunk();
  }

  if ( _freeListHead == NULL ) {
    ret = 0;
  }else{
    SlBoundOverlapServer::overlap *hold = _freeListHead;
    _freeListHead = _freeListHead->next;

    hold->next = _pendListHead;
    _pendListHead = hold;

    ret = 1;
  }

  return ret;
}
  
int SlBoundOverlapServer::popPend() {
  register int ret;
  if ( _pendListHead == NULL ) {
    ret = 0;
  }else{
    SlBoundOverlapServer::overlap *hold = _pendListHead;
    _pendListHead = _pendListHead->next;

    hold->next = _freeListHead;
    _freeListHead = hold;

    ret = 1;
  }
  return ret;
}

int SlBoundOverlapServer::addToCache(int allowGrowth) {
  register int ret;

  if (( _freeListHead == NULL ) && (allowGrowth)) {
    addChunk();
  }

  if ( _freeListHead == NULL ) {
    ret = 0;
  }else{
    SlBoundOverlapServer::overlap *hold = _freeListHead;
    _freeListHead = _freeListHead->next;

    hold->next = _cacheListHead;
    _cacheListHead = hold;

    ret = 1;
  }
  return ret;
}

void SlBoundOverlapServer::removeFromCache(overlap **node) {
  
  if ((node == NULL) && ((*node) != NULL)) {
    SlBoundOverlapServer::overlap * hold = *node;
    (*node) = hold->next;
    hold->next = _freeListHead;
    _freeListHead = hold;
  }
}
 
//----------------------------------------------------------
//----------------------------------------------------------

#endif
