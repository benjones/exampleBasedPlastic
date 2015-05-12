

__kernel void skinVertexVarying(
	int numRealBones,
	int numBoneTips,
	int numNodes,
	int numVertices,
	float scaleFactor,
	const __global  float* unskinnedPositions,
	const __global float* translations, //all examples for bone j, then all examples for j+1, etc
	const __global float* rotations,
	const __global float* barycentricCoordinates, //column major, floats
	const __global int* boneIndices,
	const __global float* boneWeights, //column major, nVerts x nBones
	__global float* skinnedPositions//,
	//__global float* perVertexTranslations,
	//__global float* perVertexRotations
  ){

  int i = get_global_id(0);
  if(i < numVertices){

	float3 output = (float3)(0,0,0);
	float3 unskinnedPosition = 
	  (float3)(unskinnedPositions[3*i], unskinnedPositions[3*i + 1], unskinnedPositions[3*i + 2]);
  
		  
	for (int j = 0; j < numBoneTips; ++j) {
	  if (boneIndices[j] < 0) { continue; } //bones[j]->get_wi() < 0) //not a real bone
	  
	  float3 interpolatedTranslation = (float3)(0,0,0);
	  float4 interpolatedRotation = (float4)(0,0,0,0);
	  for(int k = 0; k < numNodes; ++k){
		float bc = barycentricCoordinates[i*numNodes + k];
		float3 readTranslation = (float3)(
			translations[3*(j*numNodes + k)    ],
			translations[3*(j*numNodes + k) + 1],
			translations[3*(j*numNodes + k) + 2]);

		float4 readRotation = (float4)(
			rotations[4*(j*numNodes + k)],
			rotations[4*(j*numNodes + k) + 1],
			rotations[4*(j*numNodes + k) + 2],
			rotations[4*(j*numNodes + k) + 3]);
		interpolatedTranslation += bc*readTranslation;
		interpolatedRotation += bc*readRotation;
	  }
	  
	  normalize(interpolatedRotation);
	  //interpolatedRotation = (float4)(0,0,0,1);

	  
	  int weightIndex = boneIndices[j];
	  float boneWeight = boneWeights[numVertices*weightIndex + i];
	  
	  //rotate by the quaternion
	  float3 t = 2 * cross(interpolatedRotation.xyz, unskinnedPosition);
	  float3 rotatedVertex  = unskinnedPosition + 
		interpolatedRotation.w * t + cross(interpolatedRotation.xyz, t);

	  /*	  float3 u = interpolatedRotation.xyz;
	  float s = interpolatedRotation.w;
	  float3 rotatedVertex = 2.0f* dot(u,unskinnedPosition)*u +
		(s*s - dot(u,u))* unskinnedPosition +
		2.0f*s*cross(u, unskinnedPosition);
	  */


	  output += scaleFactor*boneWeight*
		(rotatedVertex + 
			interpolatedTranslation);

	  /*
	  perVertexTranslations[3*(i*numRealBones + weightIndex) ] = interpolatedTranslation.x;
	  perVertexTranslations[3*(i*numRealBones + weightIndex) + 1] = interpolatedTranslation.y;
	  perVertexTranslations[3*(i*numRealBones + weightIndex) + 2] = interpolatedTranslation.z;
	  
	  perVertexRotations[4*(i*numRealBones + weightIndex)] = interpolatedRotation.x;
	  perVertexRotations[4*(i*numRealBones + weightIndex) + 1] = interpolatedRotation.y;
	  perVertexRotations[4*(i*numRealBones + weightIndex) + 2] = interpolatedRotation.z;
	  perVertexRotations[4*(i*numRealBones + weightIndex) + 3] = interpolatedRotation.w;
	  */
	  }
	//float3 finalPosition = scaleFactor*unskinnedPosition;
	skinnedPositions[3*i] = output.x;//finalPosition.x;//output;
	skinnedPositions[3*i + 1] = output.y;//finalPosition.y;//output;
	skinnedPositions[3*i + 2] = output.z;//finalPosition.z;//output;
  }
}
	
