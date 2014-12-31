#!/bin/zsh
echo $#
if ! [[ $# -eq 2 ]] ; then
	echo "usage: ./upateCuts <directory for fracture> <directory for tetmeshSlicer>"
else
	echo "removing"
	rm -v $1/clippingTriangles.*.ply
	rm -v $1/cutTets.*.txt
	rm -v $1/tetmesh.*.txt
	rm -v $1/planes.txt
	echo "copying"
	cp -v $2/clippingTriangles.*.ply $1
	cp -v $2/cutTets.*.txt $1
	cp -v $2/tetmesh.*.txt $1
	cp -v $2/planes.txt $1

fi
