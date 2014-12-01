#include "plasticObject.h"

#include "kernels.h"

#include <iostream>
#include <cstdlib>
#include <numeric> //accumulate
#include <deque>
#include <set>

#include <tbb/tbb.h>

#include <skeleton.h>
#include <read_BF.h>
//#include <lbs_matrix.h>
//#include <gather_transformations.h>
//#include <columnize.h>

#include <igl/readDMAT.h>
//#include <igl/readOBJ.h>
//#include <igl/writeOBJ.h>
#include <igl/readMESH.h>
//#include <igl/cotmatrix.h>
//#include <igl/adjacency_matrix.h>
//#include <igl/sum.h>
//#include <igl/diag.h>

#include "plyIO.hpp"

//#include <Eigen/QR>

#include "utils.h"

#include "range.hpp"
using benlib::range;

#include "enumerate.hpp"
using benlib::enumerate;

void PlasticObject::loadFromFiles(std::string directory){
  
  skeleton = std::unique_ptr<Skeleton<Bone>>{new Skeleton<Bone>{}};
  if(!read_BF((directory + "/bone_roots.bf").c_str(), 
			  skeleton.get(), skeleton->roots)){
	std::cout << "couldn't read bone roots: " << std::endl;
	exit(1);
  }

  /*  if(!igl::readDMAT(directory + "/weights.dmat", boneWeights)){
	std::cout << "couldn't read dmat file: " << std::endl;
	exit(1);
  }
  */
  exampleGraph.parent = skeleton.get();
  
  bones = gather_bones(skeleton->roots);
  numBoneTips = bones.size();
  std::cout << "numBoneTips" << numBoneTips << std::endl;

  exampleGraph.load(directory + "/exampleGraph.txt");
  numNodes = exampleGraph.nodes.size();

  //egTraverser = EGTraverser{&exampleGraph};
  

  std::cout << "eg parent: " << exampleGraph.parent << " roots " 
			<< exampleGraph.parent->roots << std::endl;

  //we're not using the high res mesh for anything
  /*if(!igl::readOBJ(directory + "/mesh.obj", 
				   trimeshVertices,
				   trimeshTriangles)){
	std::cout << "couldn't read obj mesh" << std::endl;
	exit(1);
	}*/
  if(!igl::readMESH(directory + "/mesh.mesh",
					tetmeshVertices,
					tetmeshTets,
					tetmeshTriangles)){
	std::cout << "couldn't read tetmesh" << std::endl;
	exit(1);
  }
  numPhysicsVertices = tetmeshVertices.rows();
  std::cout << "tetmesh has: " << numPhysicsVertices << " verts" << std::endl
			<< "and " << tetmeshTets.rows() << " tets" << std::endl;

  
  computeTriangleFaces();

  computeVertexNeighbors();

  if(!igl::readDMAT(directory + "/coarseWeights.dmat", boneWeights)){
	std::cout << "couldn't read skinning weights" << std::endl;
	exit(1);
  }
  
  numRealBones = boneWeights.cols();

  //sanity check
  assert (numRealBones == 
		  std::count_if(bones.begin(), bones.end(), 
						[](const Bone* b){
						  return b->get_wi() >= 0;
						}));

  //just the translational bit for each handle
  collisionHandleAdjustments = Eigen::MatrixXd::Zero(3, boneWeights.cols());

  //lbs_matrix(tetmeshVertices, boneWeights, LBSMatrix);

  //igl::cotmatrix(tetmeshVertices, tetmeshTets, cotangentLaplacian);
  /*Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(tetmeshTets, A);
  Eigen::SparseVector<double> Asum;
  igl::sum(A, 1, Asum);
  Eigen::SparseMatrix<double> Adiag;
  igl::diag(Asum, Adiag);
  cotangentLaplacian = 2*Adiag + A;//A - Adiag;

  for(auto i : range(cotangentLaplacian.cols())){
	cotangentLaplacian.col(i) /= cotangentLaplacian.col(i).sum();
  }
  //  std::cout << cotangentLaplacian.col(0) << std::endl;
  //  exit(1);
  */
  barycentricCoordinates = Eigen::MatrixXd::Zero(numPhysicsVertices, numNodes);
  barycentricCoordinates.col(0).setOnes();// = 1; //set everything to be in rest pose

  deltaBarycentricCoordinates.resize(numPhysicsVertices, numNodes);
  deltaBarycentricCoordinates.setZero();

  perExampleScale = Eigen::VectorXd::Ones(numNodes);

  skinMeshVaryingBarycentricCoords();//precompute rotations and stuff

}


void PlasticObject::computeMassesAndVolume(){
  
  volume = 0;
  mass = 0;
  //compute object centerOfMass
  tetmeshVertexMasses.resize(numPhysicsVertices);
  tetmeshVertexMasses.setZero();

  Eigen::Matrix3d basis(3,3);
  for(auto tetInd : range(tetmeshTets.rows())){
	basis.row(0) = scaleFactor*(tetmeshVertices.row(tetmeshTets(tetInd,1)) - 
								tetmeshVertices.row(tetmeshTets(tetInd, 0)));
	basis.row(1) = scaleFactor*(tetmeshVertices.row(tetmeshTets(tetInd,2)) - 
								tetmeshVertices.row(tetmeshTets(tetInd, 0)));
	basis.row(2) = scaleFactor*(tetmeshVertices.row(tetmeshTets(tetInd,3)) - 
								tetmeshVertices.row(tetmeshTets(tetInd, 0)));

	//shouldn't be necessary, but whatever
	double tetVolume = -basis.determinant()/6.0;
	if(tetVolume < 0){
	  std::cout << "inverted tet!!!!!!" << std::endl;
	  std::cout << " at tet: " << tetInd << " vol " << tetVolume << std::endl;
	}
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
  worldTransform = 
	//bulletBody->getCenterOfMassTransform()*
	bulletSnapshot.comTransform;//*
  //inertiaAligningTransform.inverse();


  Eigen::Vector3d centerOfMass = 
	(tetmeshVertexMasses.asDiagonal()*vertices).colwise().sum().transpose()/
	mass;

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


  //translate and then rotate each vertex so that
  //the object is centered at the origin and aligned with the 
  //princial moments of inertia

  for(auto i : range(vertices.rows())){
	
	currentBulletVertexPositions.row(i) = (evecs.transpose()*
											   (vertices.row(i).transpose() - 
												centerOfMass)
											   ).transpose();
	
  }
  
  bulletShape->postUpdate();
  bulletShape->updateBound();
  bulletBody->setActivationState(DISABLE_DEACTIVATION);
  

  inertiaAligningTransform = btTransform{eigenToBullet(evecs),
  										 eigenToBullet(centerOfMass)};

  bulletBody->setCenterOfMassTransform(worldTransform*
									   inertiaAligningTransform);

  bulletBody->setMassProps(mass, eigenToBullet(evals));
  bulletBody->updateInertiaTensor(); //roate it
  
}

void PlasticObject::dump(std::string filename){
  //dump coarse mesh from bullet I guess
  auto bulletTrans = bulletBody->getCenterOfMassTransform();
  auto eigenRot = bulletToEigen(bulletTrans.getBasis());
  Eigen::Vector3d eigenTrans{bulletToEigen(bulletTrans.getOrigin())};
  RMMatrix3f transformedVertices(currentBulletVertexPositions.rows(), 3);
  for(auto i : range(transformedVertices.rows())){
	Eigen::Vector3d transformedVertex = 
	  eigenRot*currentBulletVertexPositions.row(i).transpose() +
	  eigenTrans;
	Eigen::Vector3f transFloat = transformedVertex.template cast<float>();
	transformedVertices.row(i) = transFloat.transpose();
	
  }
  
  std::ofstream plyStream(filename);
  writePLY(plyStream, transformedVertices, tetmeshTriangles);
  plyStream.close();

  /*
  RMMatrix3d sanityVertices;
  RMMatrix3i sanityTriangles;
  std::ifstream plyIns(filename + std::string(".ply"));
  readPLY(plyIns, sanityVertices, sanityTriangles);
  assert(sanityVertices.rows() == transformedVertices.rows());
  assert(sanityTriangles.rows() == tetmeshTriangles.rows());
  */
  
}

void PlasticObject::dumpBarycentricCoords(std::string filename) const{
  writeMatrixBinary(filename, barycentricCoordinates);
}

void PlasticObject::dumpVerticesBinary(std::string filename) const {
  auto bulletTrans = bulletBody->getCenterOfMassTransform();
  auto eigenRot = bulletToEigen(bulletTrans.getBasis());
  Eigen::Vector3d eigenTrans{bulletToEigen(bulletTrans.getOrigin())};
  RMMatrix3d transformedVertices(currentBulletVertexPositions);
  for(auto i : range(transformedVertices.rows())){
	transformedVertices.row(i) 
	  = (eigenRot*transformedVertices.row(i).transpose() +
		 eigenTrans).transpose();
	
  }


  std::ofstream outs(filename, std::ios_base::out | std::ios_base::binary);
  
  size_t nVerts = transformedVertices.rows();
  outs.write(reinterpret_cast<const char*>(&nVerts), sizeof(nVerts));
  outs.write(reinterpret_cast<const char*>(transformedVertices.data()),
			 3*nVerts*sizeof(double));

}




void PlasticObject::skinMesh(){//btDiscreteDynamicsWorld& bulletWorld){
  
  
  //btTransform trans;
  //trans.setIdentity();
  //btVector3 aabbMin, aabbMax;
  //bulletShape->getAabb(trans, aabbMin,aabbMax);
  /*std::cout << "before: " //<< trans << std::endl
			<< "aabbmin" <<aabbMin << std::endl
			<< "aabbmax" <<aabbMax << std::endl;
  */

  /*
  exampleGraph.setInterpolatedTransform(egTraverser.currentPosition.simplex,
										egTraverser.currentPosition.coords);

  
  currentTransformMatrix.resize(3, 4*boneWeights.cols());
  if(!gather_transformations(skeleton->roots, false, 
							 currentTransformMatrix)){
	std::cout << "couldn't gather transforms" << std::endl;
	return;
  }

  //add the collisionHandleAdjustments to this
  for(auto i : range(collisionHandleAdjustments.cols())){
	currentTransformMatrix.col(4*i + 3) += collisionHandleAdjustments.col(i);
  }

  //std::cout << "current trans matrix: " << currentTransformMatrix << std::endl;

  //do as the FAST skinners do
  
  Eigen::MatrixXd TCol;
  columnize<double, 
			Eigen::Dynamic>(currentTransformMatrix.block(0,0,3, 
														 currentTransformMatrix.cols()).eval(),
							currentTransformMatrix.cols()/4, 2, TCol);


  Eigen::MatrixXd U3 = scaleFactor*(LBSMatrix*TCol.cast<double>().eval());
  
  
  auto* dataBefore = currentBulletVertexPositions.data();
												  

  currentBulletVertexPositions.col(0) = U3.block(0,0, 
												 currentBulletVertexPositions.rows(),1).eval();
  currentBulletVertexPositions.col(1) = U3.block(currentBulletVertexPositions.rows(), 0,
												 currentBulletVertexPositions.rows(), 1);
  currentBulletVertexPositions.col(2) = U3.block(2*currentBulletVertexPositions.rows(), 0,
												 currentBulletVertexPositions.rows(), 1);


  if(isnan(currentBulletVertexPositions(0,0))){
	std::cout << "post skin: " << currentBulletVertexPositions.row(0) << std::endl;
	exit(1);
  }
  //std::cout << "data after: " << currentBulletVertexPositions.data() << std::endl;
  if(dataBefore != currentBulletVertexPositions.data()){
	std::cout << "data changed" << std::endl;
	exit(1);
	}*/

  /*  
	  set everything to the same barycentric coordinate
Eigen::VectorXd basicRow = Eigen::VectorXd::Zero(numNodes);
  for(auto i : range(egTraverser.currentPosition.coords.size())){
	basicRow(exampleGraph.simplices[egTraverser.currentPosition.simplex][i]) =
	  egTraverser.currentPosition.coords[i];
  }
  for(auto i : range(barycentricCoordinates.rows())){
	barycentricCoordinates.row(i) = basicRow.transpose();
  }
  */
  skinMeshVaryingBarycentricCoords();

  //add the extra offsets
  currentBulletVertexPositions += localImpulseBasedOffsets;




}


/*void PlasticObject::deformBasedOnImpulses(btPersistentManifold* man, bool isObject0){
  for(auto j : range(man->getNumContacts())){
	auto& manPoint = man->getContactPoint(j);
	
	auto impulseToApply = bulletToEigen(getDeformationVectorFromImpulse(manPoint, isObject0));


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
*/

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
  //  std::cout << "nearest vertex: " << best << " dist: " << bestDist << std::endl;

  return best;
}


bool PlasticObject::deformBasedOnImpulseLocal(btPersistentManifold* man, 
											  bool isObject0){

  bool deformed = false;
  
  for(auto j : range(man->getNumContacts())){
	auto& manPoint = man->getContactPoint(j);
	
	auto triangleIndex = isObject0 ? 
	  manPoint.m_index0 : manPoint.m_index1;
	auto localPoint = isObject0 ? 
	  manPoint.m_localPointA : manPoint.m_localPointB;
	//	  std::cout << "local point: " << localPoint << std::endl;
	auto impulseToApply = 
	  bulletToEigen(getLocalDeformationVectorFromImpulse(manPoint, 
														 isObject0));
	if(impulseToApply.squaredNorm() <= 0){continue;} //ignore 0 impulses
	deformed = true;
	Eigen::VectorXd weightsAtContact;

	auto vIndex = getNearestVertex(bulletToEigen(localPoint));
	weightsAtContact = boneWeights.row(vIndex).transpose();
	
	for(auto i : range(numPhysicsVertices)){
	  auto weightSpaceDistance = (weightsAtContact - 
								  boneWeights.row(i).transpose()).norm();
	  
	  auto scale = Kernels::simpleCubic(weightSpaceDistance/0.5);
	  localImpulseBasedOffsets.row(i) += 
		scale*impulseToApply.transpose();
	  //		if(scale >0){
	  //		  std::cout << "adding to vertex " << i << " " << scale*impulseToApply << std::endl;
	  //		}
	  
	}
	
	
  }
  return deformed;

}

btVector3 PlasticObject::getDeformationVectorFromImpulse(const btManifoldPoint& manPoint,
														 bool isObject0){

  btVector3 impulseGlobal = manPoint.m_normalWorldOnB*manPoint.m_appliedImpulse;

  btVector3 displacementGlobal = impulseGlobal*dt/mass;
  if(!isObject0){displacementGlobal *= -1;}
  //if(displacementGlobal.norm() > 0){std::cout << "global impulse: " << displacementGlobal << std::endl;}
  //rotate this into the unskinned space
  
  auto currentRotation = bulletBody->getCenterOfMassTransform().getRotation()*
	inertiaAligningTransform.getRotation().inverse();
  btVector3 displacementLocal = quatRotate(currentRotation.inverse(),displacementGlobal);
  
  if(displacementLocal.norm() < plasticityImpulseYield){
	return btVector3(0,0,0);
  }

  //	std::cout << "weigths at contact: " << std::endl
  //			  << weightsAtContact << std::endl;
  auto normalizedDisplacement = displacementLocal.normalized();
  //	std::cout << "unnormalized impulse: " << impulseLocal << std::endl;
  //	std::cout << "normalized impulse: " << normalizedImpulse << std::endl;
  auto displacementToApply = plasticityImpulseScale*
	(displacementLocal - 
	 plasticityImpulseYield*normalizedDisplacement);
  
  return displacementToApply;
}

btVector3 PlasticObject::getLocalDeformationVectorFromImpulse(const btManifoldPoint& manPoint,
														 bool isObject0){

  btVector3 impulseGlobal = manPoint.m_normalWorldOnB*manPoint.m_appliedImpulse;

  btVector3 displacementGlobal = impulseGlobal*dt/mass;
  if(!isObject0){displacementGlobal *= -1;}
  //if(displacementGlobal.norm() > 0){std::cout << "global impulse: " << displacementGlobal << std::endl;}
  //rotate this into the unskinned space
  
  auto currentRotation = bulletBody->getCenterOfMassTransform().getRotation()*
	inertiaAligningTransform.getRotation().inverse();
  btVector3 displacementLocal = quatRotate(currentRotation.inverse(),displacementGlobal);
  

  if(displacementLocal.norm() < localPlasticityImpulseYield){
	return btVector3(0,0,0);
  }

  //std::cout << "displacement local norm: " << displacementLocal.norm() << std::endl;

  //	std::cout << "weigths at contact: " << std::endl
  //			  << weightsAtContact << std::endl;
  auto normalizedDisplacement = displacementLocal.normalized();
  //	std::cout << "unnormalized impulse: " << impulseLocal << std::endl;
  //	std::cout << "normalized impulse: " << normalizedImpulse << std::endl;
  auto displacementToApply = localPlasticityImpulseScale*
	(displacementLocal - 
	 localPlasticityImpulseYield*normalizedDisplacement);
  
  return displacementToApply;
}

bool PlasticObject::projectImpulsesOntoExampleManifoldLocally(){
  if(plasticityImpulseScale == 0){
	return false;
  }
  bool deformed = false;
  //desiredDeformations.resize(numPhysicsVertices, 3);
  //desiredDeformations.setZero();
  //  deltaBarycentricCoordinates.resize(numPhysicsVertices, numNodes);
  //  deltaBarycentricCoordinates.setZero();


  Eigen::VectorXd weightsAtContact; //save on allocations
  Eigen::MatrixXd jacobian;
  /*if(manifoldPoints.size() > 0){
	std::cout << "processing " << manifoldPoints.size() << " contacts for object\n";
	}*/
  for(auto& pr : manifoldPoints){
	auto& manPoint = pr.first;
	bool isObject0 = pr.second;
	
	auto& localPoint = isObject0 ? manPoint.m_localPointA : manPoint.m_localPointB;
	auto vInd = getNearestVertex(bulletToEigen(localPoint));
	
	auto impulseAtContact = 
	  bulletToEigen(getDeformationVectorFromImpulse(manPoint, isObject0));
	if(impulseAtContact.squaredNorm() <= 0){continue;} //ignore 0 impulses
	//std::cout << "contact impulse SN: " << impulseAtContact.squaredNorm() << std::endl;
	deformed = true;

	weightsAtContact = boneWeights.row(vInd).transpose();

	//compute the derivative
	jacobian.resize(3, numNodes);
	jacobian.setZero();
	
	for(auto nodeIndex : range(numNodes)){
	  auto& node = exampleGraph.nodes[nodeIndex];
	  for(auto boneIndex : range(numBoneTips)){
		auto weightIndex = bones[boneIndex]->get_wi();
		Eigen::AngleAxis<double> nodeRotation{node.transformations[weightIndex].rotation};
		Eigen::AngleAxis<double> currentRotation{perVertexRotations[nodeIndex*numRealBones +
																	weightIndex]};
		
		//		std::cout << "node rotaiton: " << nodeRotation.angle() << ", " << nodeRotation.axis() << std::endl;
		//		std::cout << "current rotation: " << currentRotation.angle() << ", " << currentRotation.axis()  << std::endl;
		jacobian.col(nodeIndex) += 
		  boneWeights(vInd, weightIndex)*
		  (node.transformations[boneIndex].translation +
		   nodeRotation.angle()*
		   currentRotation.axis().cross(currentRotation*
										tetmeshVertices.row(vInd).transpose()));
		   }
	  
	}
	
	//	std::cout << "jacobian: " << jacobian << std::endl;
	//scale the columns
	//std::cout << "col norms: " << jacobian.colwise().norm() << std::endl;

	//jacobian = (jacobian*perExampleScale.asDiagonal().inverse()).eval();
	
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian, 
										  Eigen::ComputeThinU | 
										  Eigen::ComputeThinV);
	//	std::cout << "singular values: " << svd.singularValues() << std::endl;
	//	Eigen::VectorXd modifiedSingularValues = sqrt(svd.singularValues().array());

	Eigen::VectorXd singularVectorContribution = 
	  //(svd.matrixV()*svd.singularValues())*impulseAtContact.norm();
	  svd.matrixV()*impulseAtContact.norm()*(svd.singularValues().asDiagonal()*
											 svd.matrixU().transpose()*
											 impulseAtContact).normalized();
	Eigen::VectorXd jacobianTransposeContribution = 
	  jacobian.transpose()*impulseAtContact;

	Eigen::VectorXd deltaS = jacobianAlpha*singularVectorContribution +
	  (1.0 - jacobianAlpha)*jacobianTransposeContribution;
	//	  modifiedSingularValues.asDiagonal()*
	//	  svd.matrixU().transpose()*
	//	  impulseAtContact;

	  //(1.0*jacobian.colwise().norm().maxCoeff())*jacobian.transpose()*impulseAtContact;
	
	//std::cout << "deltaS " << deltaS << std::endl;
	/*	if((impulseAtContact - jacobian*deltaS).norm() > impulseAtContact.norm()){
	  std::cout << "norm grew from " << impulseAtContact.norm() 
				<< " to " << (impulseAtContact - jacobian*deltaS).norm() << std::endl;
	  std::cout << "i - jj^t norm: " << (Eigen::MatrixXd::Identity(jacobian.rows(), jacobian.rows()) -
										 jacobian*jacobian.transpose()).norm() << std::endl;

	}
	*/
	//	deltaS = jacobian.jacobiSvd(Eigen::ComputeThinU | 
	//								Eigen::ComputeThinV).solve(impulseAtContact -
	//													   jacobian*deltaS);

  //	std::cout << "deltaS both: " << deltaS << std::endl;
	
	//scale it
	deltaS /= std::max(std::fabs(deltaS.maxCoeff()), 1.0);
	
	if(deltaS.squaredNorm() > 0.000001){

	  computeGeodesicDistances(vInd, 1.0/plasticityKernelScale);
	  
	  //distribute this to evereyone else
	  for(auto i : range(numPhysicsVertices)){
		//		auto weightSpaceDistance =
		//		  (weightsAtContact - boneWeights.row(i).transpose()).norm();
		//		auto scale = Kernels::simpleCubic(weightSpaceDistance*plasticityKernelScale);
		auto scale = Kernels::simpleCubic(plasticityKernelScale*geodesicDistances[i]);
		deltaBarycentricCoordinates.row(i) += scale*deltaS.transpose();
	  }
	}
  }

  //hit the delta with the laplacian to smooth it
  /*for(auto i : range(8)){
	deltaBarycentricCoordinates = (cotangentLaplacian.transpose()*deltaBarycentricCoordinates).eval();
	//	deltaBarycentricCoordinates -= 0.005*(cotangentLaplacian*deltaBarycentricCoordinates).eval();
	//	if(deltaBarycentricCoordinates.row(i).squaredNorm() > 0.001){
	//	  std::cout << "smoothing iter: " << i << ' ' << deltaBarycentricCoordinates.row(0) << std::endl;
	//	}
  }
  */

  barycentricCoordinates += plasticityRate*deltaBarycentricCoordinates*perExampleScale.asDiagonal();
  deltaBarycentricCoordinates *= 1.0 - plasticityRate;

  //clamp and rescale:
  for(auto i : range(numPhysicsVertices)){
	for(auto j : range(numNodes)){
	  barycentricCoordinates(i, j) = std::max(barycentricCoordinates(i,j), 0.0);
	}
	double sum = barycentricCoordinates.row(i).sum();
	if(fabs(sum) < 0.0001){
	  barycentricCoordinates(i, 0) = 1;
	} else {
	  barycentricCoordinates.row(i) /= sum;
	}
  }
  //  std::cout << "row 0: " << std::endl;
  //  std::cout << barycentricCoordinates.row(0) << std::endl;
  //  if(barycentricCoordinates(0,2) > 0.1) {exit(1);}
  checkNans(barycentricCoordinates);

  return deformed;
  
}

/*bool PlasticObject::projectImpulsesOntoExampleManifold(){
  
  if(plasticityImpulseScale == 0){
	return false;
  }
  
  std::vector<int> vertexIndices;
  std::vector<btVector3> impulses;
  
  bool deformed = false;
  for(auto& pr : manifoldPoints){
	auto& manPoint = pr.first;
	bool isObject0 = pr.second;

	//auto triangleIndex = isObject0 ? manPoint.m_index0 : manPoint.m_index1;
	auto localPoint = isObject0 ? manPoint.m_localPointA : manPoint.m_localPointB;
	
	//always use closest point, even it's actually on a triangle
	//if(triangleIndex < 0){
	auto vInd = getNearestVertex(bulletToEigen(localPoint));
	auto it = std::find(vertexIndices.begin(), vertexIndices.end(),
						vInd);
	
	
	if(it == vertexIndices.end()){
	  vertexIndices.push_back(vInd);
	  impulses.push_back(getDeformationVectorFromImpulse(manPoint, isObject0));
	} else {
	  impulses[std::distance(vertexIndices.begin(), it)] +=
		getDeformationVectorFromImpulse(manPoint, isObject0);
	}
	//} else {
	
	  
	//}
  }
  
  //copy impulses into an Eigen::vector to do a least squares solve
  Eigen::VectorXd rhs(3*impulses.size());
  for(auto i : range(impulses.size())){
	rhs(3*i    ) = impulses[i].x();
	rhs(3*i + 1) = impulses[i].y();
	rhs(3*i + 2) = impulses[i].z();
  }
  
  //compute derivates wrt to skinning weights

  //these are partial (entry of Transformation matrix) wrt skinning weights
  //this is n*3x simplex.size() matrix
  Eigen::MatrixXd jacobian;
  computeTransformationDerivatives(egTraverser.currentPosition,
								   vertexIndices,
								   jacobian);

  //do a least squares solve for the change in S
  Eigen::VectorXd deltaS = jacobian.colPivHouseholderQr().solve(rhs);
  if(deltaS.array().abs().sum() > 0.0001){
	deformed = true;
	std::cout << "delta S" << deltaS << std::endl;
	
	//EGTraverser::zeroWeightsSum(deltaS);
	//std::cout << "zero summed: " << deltaS << std::endl;
	
	egTraverser.currentPosition = EGTraverser::addBcVector(egTraverser.currentPosition,
														   deltaS);
	std::cout <<egTraverser.currentPosition.coords << std::endl;
  }
  return deformed;
  }*/

/*void PlasticObject::computeTransformationDerivatives(const EGPosition& egPosition, 
													 const std::vector<int>& indices,
													 Eigen::MatrixXd& deriv){
  
  std::vector<double> currentWeights(numNodes, 0.0);
  auto& simplex = exampleGraph.simplices[egPosition.simplex];
  for(auto i : range(simplex.size())){
	currentWeights[simplex[i]] = egPosition.coords[i];
  }
  


  deriv = Eigen::MatrixXd::Zero(3*indices.size(), simplex.size());
  
  //loop over handles
  for(auto hInd : range(boneWeights.cols())){
	
	auto currentInterpedTransform = 
	  exampleGraph.interpolateTransform(hInd, currentWeights);
	
	Eigen::AngleAxis<double> currentRotation{currentInterpedTransform.rotation};

	
	//loop over vertices that we care about
	for(auto vInd : range(indices.size())){
	  
	  //loop over simplex nodes
	  for(auto sInd : range(simplex.size())){
		
		auto& node = exampleGraph.nodes[simplex[sInd]];
		Eigen::AngleAxis<double> nodeRotation{node.transformations[hInd].rotation};		
		//add in the translational part
		deriv.block(vInd*3, sInd, 3, 1) += 
		  boneWeights(vInd, hInd)*node.transformations[hInd].translation;
		
		//and the rotational part...
		//I think this is well approximated by
		//w_i*angle(Q(sInd))*QTotal.axis x QTotal.rotate(v)
		
		deriv.block(vInd*3, sInd, 3, 1) +=
		  boneWeights(vInd, hInd)*
		  nodeRotation.angle()*
		  currentRotation.axis().cross(currentRotation*
									   tetmeshVertices.row(vInd).transpose());
	  }
	}
  }

  }*/

void PlasticObject::saveBulletSnapshot(){
  bulletSnapshot.linearVelocity = bulletBody->getLinearVelocity();
  bulletSnapshot.angularMomentum = bulletBody->getInvInertiaTensorWorld().inverse()*
	bulletBody->getAngularVelocity();
  //keep the inertia part out of it
  bulletSnapshot.comTransform = bulletBody->getCenterOfMassTransform()*
	inertiaAligningTransform.inverse();
}

void PlasticObject::restoreBulletSnapshot(){
  bulletBody->setLinearVelocity(bulletSnapshot.linearVelocity);
  bulletBody->setAngularVelocity(bulletBody->getInvInertiaTensorWorld()*bulletSnapshot.angularMomentum);
  bulletBody->setCenterOfMassTransform(bulletSnapshot.comTransform*inertiaAligningTransform);
}


void PlasticObject::skinMeshVaryingBarycentricCoords(){
  
  perVertexTranslations.resize(numPhysicsVertices*numRealBones);
  perVertexRotations.resize(numPhysicsVertices*numRealBones);

  currentBulletVertexPositions.resize(numPhysicsVertices, 3);
  currentBulletVertexPositions.setZero();

//  for(auto i : range(numPhysicsVertices)){
  tbb::parallel_for(tbb::blocked_range<size_t>(0, numPhysicsVertices),
					[this](const tbb::blocked_range<size_t>& r){
					  for(auto i = r.begin(); i != r.end(); ++i){
						for(auto j : range(numBoneTips)){
						  if(bones[j]->get_wi() < 0){ continue;} //not a real bone
						  
						  const auto kRange = range(numNodes);
						  //compute bone transform
						  const Vec3 translation = 
							std::accumulate(kRange.begin(), kRange.end(),
											Vec3::Zero().eval(),
											[this, i,j](const Vec3& acc, size_t k){
											  return (acc + barycentricCoordinates(i, k)*
													  exampleGraph.nodes[k].transformations[j].translation).eval();
											});
						  const Eigen::Vector4d rotationCoeffs =
							std::accumulate(kRange.begin(), kRange.end(),
											Eigen::Vector4d::Zero().eval(),
											[this, i,j](const Eigen::Vector4d& acc, size_t k){
											  return acc + barycentricCoordinates(i, k)*
											  exampleGraph.nodes[k].transformations[j].rotation.coeffs();
											});
						  Quat rotation(rotationCoeffs);
						  rotation.normalize();
	  
						  const auto weightIndex = bones[j]->get_wi();
						  
						  currentBulletVertexPositions.row(i) +=
							scaleFactor*boneWeights(i,weightIndex)*
							(rotation*(tetmeshVertices.row(i).transpose()) + 
							 translation).transpose();
						  
						  perVertexTranslations[i*numRealBones + weightIndex] = translation;
						  perVertexRotations[i*numRealBones + weightIndex] = rotation;
						  
						}
					  }
					});
}

void PlasticObject::computeGeodesicDistances(size_t vInd, double radius){

  geodesicDistances.assign(numPhysicsVertices, std::numeric_limits<double>::infinity());
  
  geodesicDistances[vInd] = 0;
  
  const size_t passes = 4;
  std::deque<size_t> queue;
  for(auto pass : range(passes)){
	queue.clear();
	queue.insert(queue.end(), vertexNeighbors[vInd].begin(), vertexNeighbors[vInd].end());
	while(!queue.empty()){
	  size_t currentVertex = queue.front();
	  queue.pop_front();
	  double oldDist = geodesicDistances[currentVertex];
	  auto neighborRange = benlib::range(vertexNeighbors[currentVertex].size());
	  auto minIt = 
		std::min_element(neighborRange.begin(),
						 neighborRange.end(),
						 [this, currentVertex](size_t ind1, size_t ind2){
						   return (geodesicDistances[vertexNeighbors[currentVertex][ind1]] + 
								   neighborDistances[currentVertex][ind1]) <
						   (geodesicDistances[vertexNeighbors[currentVertex][ind2]] +
							neighborDistances[currentVertex][ind2]);
						 });
	  double bestDist = geodesicDistances[vertexNeighbors[currentVertex][*minIt]] + 
		neighborDistances[currentVertex][*minIt];
	  if(bestDist < oldDist &&
		 bestDist < radius){
		geodesicDistances[currentVertex] = bestDist;
		queue.insert(queue.end(), 
					 vertexNeighbors[currentVertex].begin(), 
					 vertexNeighbors[currentVertex].end());
	  }
	}
  }
}

void PlasticObject::computeVertexNeighbors(){

  vertexNeighbors.resize(numPhysicsVertices);
  neighborDistances.resize(numPhysicsVertices);
  for(auto row : range(tetmeshTets.rows())){
	for(auto i : range(4)){
	  vertexNeighbors[tetmeshTets(row, i)].push_back(tetmeshTets(row, (i + 1)%4));
	  vertexNeighbors[tetmeshTets(row, i)].push_back(tetmeshTets(row, (i + 2)%4));
	  vertexNeighbors[tetmeshTets(row, i)].push_back(tetmeshTets(row, (i + 3)%4));
	}
  }

  //trash duplicates
  for(auto vInd : range(numPhysicsVertices)){
	auto& vlist = vertexNeighbors[vInd];
	std::sort(vlist.begin(), vlist.end());
	vlist.erase(std::unique(vlist.begin(), vlist.end()),
				vlist.end());
	
	neighborDistances[vInd].resize(vlist.size());
	std::transform(vlist.begin(), vlist.end(),
				   neighborDistances[vInd].begin(),
				   [this, vInd](size_t nInd){
					 return (tetmeshVertices.row(vInd) - tetmeshVertices.row(nInd)).norm();
				   });
  }
}


void PlasticObject::dumpImpulses(std::string filename) const{
  if(manifoldPoints.size()){
	std::ofstream outs(filename);
	for(auto& mp : manifoldPoints){
	  auto& pos = mp.second ? mp.first.getPositionWorldOnA() : mp.first.getPositionWorldOnB();
	  auto imp = 0.2*(mp.second ? 1.0 : -1.0)*mp.first.m_normalWorldOnB;
	  outs << pos.x() << ' ' << pos.y() << ' ' << pos.z() << ' '
		   << pos.x() + imp.x() << ' ' << pos.y() + imp.y() << ' ' << pos.z() + imp.z() << '\n';
	}
  }
}


void PlasticObject::computeTriangleFaces(){
  
  //use sorted vector instead of set if speed is an issue?
  std::set<std::tuple<int,int,int>> seenFaces;
  int triIndices[4][3] = { {0,1,2},{0,1,3},{0,2,3},{1,2,3}};
  auto makeTri = [](int a, int b, int c) -> std::tuple<int,int,int>{
	int max = std::max({a, b, c});
	int min = std::min({a, b, c});
	int mid = ((a != max) && ( a != min)) ? a :
	( ((b != max) && ( b != min)) ? b : c);
	return std::tuple<int,int,int>(min,mid,max);
  };
  
  for(auto i : range(tetmeshTets.rows())){
	for(auto j : range(4)){
	  auto tri = makeTri(tetmeshTets(i, triIndices[j][0]),
						 tetmeshTets(i, triIndices[j][1]),
						 tetmeshTets(i, triIndices[j][2]));
	  auto it = seenFaces.find(tri);
	  if(it == seenFaces.end()){
		seenFaces.insert(tri);
	  } else {
		seenFaces.erase(it);
	  }
	}
  }
  
  std::cout << "mesh file tris: " << tetmeshTriangles.rows() << std::endl;
  std::cout << "computed: " << seenFaces.size() << std::endl;

  //orientation is ignored by bullet, so don't reorient stuff
  tetmeshTriangles.resize(seenFaces.size(), 3);

  for(auto pr : enumerate(seenFaces)){
	tetmeshTriangles(pr.first, 0) = std::get<0>(pr.second);
	tetmeshTriangles(pr.first, 1) = std::get<1>(pr.second);
	tetmeshTriangles(pr.first, 2) = std::get<2>(pr.second);
  }

}
