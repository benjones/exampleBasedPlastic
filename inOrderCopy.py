#!/usr/bin/python

import sys
import os.path
import os
import re

usage = "inOrderCopy.py <local dir> <remote dir> <file format string with 1 * and 1 %0d's>"


if len(sys.argv) != 4:
    exit(usage)


localDir = sys.argv[1]
remoteDir = sys.argv[2]
fileFormatString = sys.argv[3]

currentFrame = 0
while True:
    
    filename = fileFormatString % (currentFrame)
    print filename
    #if(not os.path.exists(os.path.join(localDir, filename))):
    #    print "no more objects this frame"
    #    break
    
    firstFilename = re.sub('\\*', '000', filename)
    print "first filename", firstFilename
    if(not os.path.exists(os.path.join(localDir, firstFilename))):
       print "no more frames"
       break

    os.system("rsync -tzvv %s %s" % (os.path.join(localDir, filename),
                                     "benjones@brain.cs.utah.edu:" + remoteDir))
                                        #os.path.join(remoteDir, filename)))


    currentFrame += 1




