#
# Makefile for R3 library.
#

CXX=g++
CXXFLAGS=-Wall -I. -g -DUSE_JPEG

ifneq (,$(findstring Darwin,$(shell uname)))
  AR=libtool -static -o
  CXXFLAGS := $(CXXFLAGS) -mmacosx-version-min=10.8
else
  AR=ar -rcs
endif

libR3.a: R3Distance.o R3Matrix.o R3Mesh.o R3Sphere.o R3Cylinder.o R3Cone.o R3Box.o R3Circle.o R3Plane.o R3Segment.o R3Ray.o R3Line.o R3Vector.o R3Point.o
	rm -f $@
	$(AR) $@ $^

R3Distance.o R3Matrix.o R3Mesh.o R3Sphere.o R3Cylinder.o R3Cone.o R3Box.o R3Circle.o R3Plane.o R3Segment.o R3Ray.o R3Line.o R3Vector.o R3Point.o: R3.h

clean:
	rm -f *.o libR3.a



