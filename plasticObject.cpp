#include "plasticObject.h"

#include <iostream>
#include <cstdlib>

#include <skeleton.h>
#include <read_BF.h>
#include <lbs_matrix.h>
#include <gather_transformations.h>
#include <columnize.h>

#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readMESH.h>

#include "utils.h"

#include "cppitertools/range.hpp"
using iter::range;

void PlasticObject::loadFromFiles(std::string directory){
  
  if(!read_BF((directory + "/bone_roots.bf").c_str(), 
			  &skeleton, skeleton.roots)){
	std::cout << "couldn't read bone roots: " << std::endl;
	exit(1);
  }

  if(!igl::readDMAT(directory + "/weights.dmat", boneWeights)){
	std::cout << "couldn't read dmat file: " << std::endl;
	exit(1);
  }
  
  exampleGraph.load(directory + "/exampleGraph.txt");
  exampleGraph.parent = &skeleton;
  egTraverser = EGTraverser{&exampleGraph};
  
  if(!igl::readOBJ(directory + "/mesh.obj", 
				   trimeshVertices,
				   trimeshTriangles)){
	std::cout << "couldn't read obj mesh" << std::endl;
	exit(1);
  }
  if(!igl::readMESH(directory + "/mesh.mesh",
					tetmeshVertices,
					tetmeshTets,
					tetmeshTriangles)){
	std::cout << "couldn't read tetmesh" << std::endl;
	exit(1);
  }
  std::cout << "tetmesh has: " << tetmeshVertices.rows() << " verts" << std::endl
			<< "and " << tetmeshTets.rows() << " tets" << std::endl;

  if(!igl::readDMAT(directory + "/coarseWeights.dmat", boneWeights)){
	std::cout << "couldn't read skinning weights" << std::endl;
	exit(1);
  }
  
  //just the translational bit for each handle
  collisionHandleAdjustments = Eigen::MatrixXd::Zero(3, boneWeights.cols());

  lbs_matrix(tetmeshVertices, boneWeights, LBSMatrix);


}


void PlasticObject::computeMassesAndVolume(){
  
  volume = 0;
  mass = 0;
  //compute object centerOfMass
  tetmeshVertexMasses.resize(tetmeshVertices.rows());
  tetmeshVertexMasses.setZero();

  Eigen::Matrix3d basis(3,3);
  for(auto tetInd : range(tetmeshTets.rows())){
	basis.row(0) = tetmeshVertices.row(tetmeshTets(tetInd,1)) - 
	  tetmeshVertices.row(tetmeshTets(tetInd, 0));
	basis.row(1) = tetmeshVertices.row(tetmeshTets(tetInd,2)) - 
	  tetmeshVertices.row(tetmeshTets(tetInd, 0));
	basis.row(2) = tetmeshVertices.row(tetmeshTets(tetInd,3)) - 
	  tetmeshVertices.row(tetmeshTets(tetInd, 0));

	//shouldn't be necessary, but whatever
	double tetVolume = fabs(basis.determinant())/6.0;
	volume += tetVolume;
	double nodeMass = 0.25*tetVolume*density;
	mass += 4*nodeMass;
	for(auto j : {0,1,2,3}){
	  tetmeshVertexMasses(tetmeshTets(tetInd,j)) += nodeMass;
	}
  }
  
}


void PlasticObject::updateBulletProperties(const RMMatrix3d& vertices,
										   const RMMatrix4i& tets){
  
  //the world part of the transform before we do this
  worldTransform = bulletBody->getCenterOfMassTransform()*
	inertiaAligningTransform.inverse();


  Eigen::Vector3d centerOfMass = 
	(tetmeshVertexMasses.asDiagonal()*vertices).colwise().sum().transpose()/
	mass;
  //  std::cout << " com: " << centerOfMass << std::endl;

  inertiaTensor.setZero();
  
  for(auto i : range(vertices.rows())){
	Eigen::Vector3d r = 
	  vertices.row(i).transpose() - 
	  centerOfMass;
	inertiaTensor += tetmeshVertexMasses(i)*
	  (r.squaredNorm()*Eigen::Matrix3d::Identity() -
	   r*r.transpose());
  }
  
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(inertiaTensor);
  Eigen::Matrix3d evecs = solver.eigenvectors();
  Eigen::Vector3d evals = solver.eigenvalues();
  if(evecs.determinant() < 0){
	//make it a rotation, not a refelction
	evecs.col(2) *= -1;
  }

  //  std::cout << "rotation: " << evecs << std::endl;
  //    std::cout << "evals: " << evals << std::endl;

  //translate and then rotate each vertex so that
  //the object is centered at the origin and aligned with the 
  //princial moments of inertia

  for(auto i : range(vertices.rows())){
	
	currentBulletVertexPositions.row(i) = (evecs.transpose()*
											   (vertices.row(i).transpose() - 
												centerOfMass)
											   ).transpose();
	
  }
  
  inertiaAligningTransform = btTransform{eigenToBullet(evecs),
  										 eigenToBullet(centerOfMass)};

  bulletBody->setCenterOfMassTransform(worldTransform*
									   inertiaAligningTransform);

  /*bulletBody->setCenterOfMassTransform(btTransform{
	  btMatrix3x3{
		evecs(0,0),evecs(0,1),evecs(0,2),
		  evecs(1,0),evecs(1,1),evecs(1,2),
		  evecs(2,0),evecs(2,1),evecs(2,2)},
		btVector3{centerOfMass(0), centerOfMass(1), centerOfMass(2)}
	});
  */
  bulletBody->setMassProps(mass, eigenToBullet(evals));
  
}

void PlasticObject::dump(std::string filename){
  //dump coarse mesh from bullet I guess
  auto bulletTrans = bulletBody->getCenterOfMassTransform();
  auto eigenRot = bulletToEigen(bulletTrans.getBasis());
  Eigen::Vector3d eigenTrans{bulletToEigen(bulletTrans.getOrigin())};
  RMMatrix3d transformedVertices(currentBulletVertexPositions);
  for(auto i : range(transformedVertices.rows())){
	transformedVertices.row(i) 
	  = (eigenRot*transformedVertices.row(i).transpose() +
		 eigenTrans).transpose();
	
  }
  igl::writeOBJ(filename,
				transformedVertices,
				tetmeshTriangles);
  
  
}


void PlasticObject::skinMesh(){

  exampleGraph.setInterpolatedTransform(egTraverser.currentPosition.simplex,
										egTraverser.currentPosition.coords);

  
  currentTransformMatrix.resize(3, 4*boneWeights.cols());
  if(!gather_transformations(skeleton.roots, false, 
							 currentTransformMatrix)){
	std::cout << "couldn't gather transforms" << std::endl;
  }

  //add the collisionHandleAdjustments to this
  for(auto i : range(collisionHandleAdjustments.cols())){
	currentTransformMatrix.col(4*i + 3) += collisionHandleAdjustments.col(i);
  }


  //do as the FAST skinners do
  
  Eigen::MatrixXd TCol;
  columnize<double, 
			Eigen::Dynamic>(currentTransformMatrix.block(0,0,3, 
														 currentTransformMatrix.cols()).eval(),
							currentTransformMatrix.cols()/4, 2, TCol);


  Eigen::MatrixXd U3 = LBSMatrix*TCol.cast<double>().eval();
  
  currentBulletVertexPositions.col(0) = U3.block(0,0, 
												 currentBulletVertexPositions.rows(),1);
  currentBulletVertexPositions.col(1) = U3.block(currentBulletVertexPositions.rows(), 0,
												 currentBulletVertexPositions.rows(), 1);
  currentBulletVertexPositions.col(2) = U3.block(2*currentBulletVertexPositions.rows(), 0,
												 currentBulletVertexPositions.rows(), 1);
  
}


void PlasticObject::deformBasedOnImpulses(btPersistentManifold* man, bool isObject0){
  for(auto j : range(man->getNumContacts())){
	auto& manPoint = man->getContactPoint(j);
	
	btVector3 impulseGlobal = manPoint.m_normalWorldOnB*manPoint.m_appliedImpulse;
	if(!isObject0){impulseGlobal *= -1;}
	
	//rotate this into the unskinned space
	
	auto currentRotation = bulletBody->getCenterOfMassTransform().getRotation();
	btVector3 impulseLocal = quatRotate(currentRotation.inverse(),impulseGlobal);

	if(impulseLocal.norm() < plasticityImpulseYield){
	  continue; //not enough impulse here to deform anything
	}

	auto triangleIndex = isObject0 ? manPoint.m_index0 : manPoint.m_index1;
	auto localPoint = isObject0 ? manPoint.m_localPointA : manPoint.m_localPointB;
	

	Eigen::VectorXd weightsAtContact;
	if(triangleIndex < 0){
	  auto vIndex = getNearestVertex(bulletToEigen(localPoint));
	  weightsAtContact = boneWeights.row(vIndex).transpose();
	} else {
	  auto bcCoords = getBarycentricCoordinates(localPoint, triangleIndex);
	  auto tri = tetmeshTriangles.row(triangleIndex);
	  weightsAtContact = (bcCoords(0)*boneWeights.row(tri(0)) +
						  bcCoords(1)*boneWeights.row(tri(1)) +
						  bcCoords(2)*boneWeights.row(tri(2))).transpose();
	}

	//	std::cout << "weigths at contact: " << std::endl
	//			  << weightsAtContact << std::endl;
	auto normalizedImpulse = impulseLocal.normalized();
	//	std::cout << "unnormalized impulse: " << impulseLocal << std::endl;
	//	std::cout << "normalized impulse: " << normalizedImpulse << std::endl;
	auto impulseToApply = bulletToEigen(plasticityImpulseScale*
										(impulseLocal - 
										 plasticityImpulseYield*normalizedImpulse));
	
	//add the outer product (scale impulse by the weights at contact)
	collisionHandleAdjustments += impulseToApply*weightsAtContact.transpose();
	

  }
}


Eigen::Vector3d PlasticObject::getBarycentricCoordinates(btVector3 localPoint, int triangleIndex){
  //adapted from check_point_triangle_proximity from eltopo's common lib
  
  
  auto tri = tetmeshTriangles.row(triangleIndex);
  Eigen::Vector3d x13 = (currentBulletVertexPositions.row(tri(0)) -
						 currentBulletVertexPositions.row(tri(2))).transpose();
  double r00 = x13.norm() + 1e-30; //use eltopo tolerances
  x13 /= r00;
  
  Eigen::Vector3d x23 = (currentBulletVertexPositions.row(tri(1)) -
			  currentBulletVertexPositions.row(tri(2))).transpose();
  double r01 = x23.dot(x13);
  x23-= r01*x13;

  double r11 = x23.norm()+1e-30;
  x23 /= r11;

  Eigen::Vector3d x03 = bulletToEigen(localPoint) -
	currentBulletVertexPositions.row(tri(2)).transpose();
  
  Eigen::Vector3d ret;
  ret(1) = (x23.dot(x03))/r11;
  ret(0) = ((x13.dot(x03)) - r01*ret(1))/r00;
  ret(2) = 1 - ret(0) - ret(1);

  if(ret.minCoeff() < 0){
	std::cout << "warning, point is outside of the triangle.  Might be bad news" << std::endl;
  }

  return ret;
}


int PlasticObject::getNearestVertex(Eigen::Vector3d localPoint){

  //eigen doesn't have nice iterators so I can't use stdlib algorithms :(
  
  int best = 0;
  double bestDist = (localPoint - currentBulletVertexPositions.row(0).transpose()).squaredNorm();
  for(auto i : range(1l, currentBulletVertexPositions.rows())){
	auto currDist = (localPoint  - currentBulletVertexPositions.row(i).transpose()).squaredNorm();
	if(currDist < bestDist){
	  bestDist = currDist;
	  best = i;
	}
  }
  return best;
}
