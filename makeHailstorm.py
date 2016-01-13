#!/usr/bin/python

import random

print '''{
  "ground" : true,
  "dt" : 0.004,
  "duration" : 10,
  "bulletObjects" : ['''

numBalls = 10

balls = []
for ballIndex in range(numBalls):
    balls.append( '''
    {
    "mass" : 20000, "shape" : "sphere", "radius" : 0.4,
    "offset" : [ %s, %s, %s], 
    "shadingXml" : "
    <bsdf type=\\"twosided\\">
    <bsdf type=\\"phong\\" >
    
    <rgb name=\\"diffuseReflectance\\" value=\\"#1111ff\\" />
    <spectrum name=\\"specularReflectance\\" value=\\"0.1\\" />
    
    </bsdf>
    </bsdf>
    "

    }
        ''' % (random.uniform(-0.5, 0.5), 6 + 0.5*ballIndex, random.uniform(-2, 2)))


print ",".join(balls) , "],"

print '''	"plasticBodies" : [
		{
			"directory" : "inputFiles/AdamsCrownVic/",
			"offset" : [0,0.75,0],
			"scaleFactor" : 1.0,
			"restitution": 0.5,
			"plasticityImpulseScale" : 2000,
			"plasticityImpulseYield" : 0.9e-2,
			"plasticityKernelScale" : 0.5,
			"plasticityRate" : 0.01,
			"density" : 6000,
			"useVolumetricCollisions" : false,

	"shadingXml" : "

<bsdf type=\\"twosided\\">
<bsdf type=\\"phong\\" >

<rgb name=\\"diffuseReflectance\\" value=\\"#1111ff\\" />
<spectrum name=\\"specularReflectance\\" value=\\"0.1\\" />

</bsdf>
</bsdf>
"


		}	],
	"groundFriction" : 0.6,
	"groundShadingXml" : "
<bsdf type=\\"twosided\\" >

<bsdf type=\\"roughdiffuse\\">
  <float name=\\"alpha\\" value=\\"0.7\\" />
  <rgb name=\\"reflectance\\" value=\\"#444444\\" />
</bsdf>
</bsdf>
"

}
'''
                  
