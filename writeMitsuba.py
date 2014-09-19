#!/usr/bin/python

import sys
import os.path

usage="writeMitsub.py <formatString> [extraFile]"

if len(sys.argv) < 2:
    exit(usage)

formatString = sys.argv[1]

directory, base = os.path.split(formatString)


extraStuff = ""
if len(sys.argv) > 2:
    extraStuff = open(sys.argv[2]).read()


currentFrame = 0
while True:
    
    objectNumber = 0
    if not os.path.exists(os.path.join(directory, base % (objectNumber, currentFrame))):
        break
        
    outfileName = os.path.join(directory, "mitsubaFrame_%04d.xml" % currentFrame)
    print "writing ", outfileName
    outfile = open(outfileName, 'w')

    outfile.write('<?xml version="1.0" encoding="utf-8"?>\n<scene version="0.5.0">\n')
    outfile.write(extraStuff)

    while True:
        if os.path.exists(os.path.join(directory, base % (objectNumber, currentFrame))):
            if objectNumber == 0:
                outfile.write('<shape type="ply">\n<string name="filename" value="%s" />\n<bsdf type="twosided"><bsdf type="diffuse"><srgb name="reflectance" value="#eeeeee"/></bsdf></bsdf></shape>\n' % 
                              (base % (objectNumber, currentFrame)))
            else:
                outfile.write('<shape type="ply">\n<string name="filename" value="%s" />\n<bsdf type="twosided"><bsdf type="diffuse"><srgb name="reflectance" value="#aaaaaa"/></bsdf></bsdf></shape>\n' % 
                #outfile.write('<shape type="ply">\n<string name="filename" value="%s" />\n<bsdf type="twosided"><bsdf type="diffuse"><texture type="wireframe" name="reflectance" ><srgb name="interiorColor" value="#aaaaaa"/></texture></bsdf></bsdf></shape>\n' % 
                              (base % (objectNumber, currentFrame)))
            objectNumber += 1
        else:
            currentFrame += 1
            outfile.write("</scene>\n")
            outfile.close()
            break
    

