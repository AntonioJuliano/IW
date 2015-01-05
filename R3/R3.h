// Include files for R3 package
#ifndef R3_INCLUDED
#define R3_INCLUDED



// Basic include files 

#if defined(_MSC_VER)
# pragma warning(disable:4996)
# pragma warning(disable:4244)
#endif
#ifdef _WIN32
# define NOMINMAX
# include <windows.h>
#endif
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cassert>
#include <cstring>
#include <cctype>
#include <cmath>
#include <climits>
#include <algorithm>
#include <vector>
using namespace std;



// OpenGL include files 

#if defined(__APPLE__)
# include <OpenGL/gl.h>
# include <OpenGL/glu.h>
#else
# include <GL/gl.h>
# include <GL/glu.h>
#endif



// R2 include files 

#include "../R2/R2.h"



// Constant declarations

#ifndef R3_X
# define R3_X 0
# define R3_Y 1
# define R3_Z 2
#endif

#ifndef M_PI
# define  M_PI  3.14159265358979323846
#endif



// Class declarations

class R3Point;
class R3Vector;
class R3Line;
class R3Ray;
class R3Segment;
class R3Plane;
class R3Circle;
class R3Box;
class R3Sphere;
class R3Cylinder;
class R3Cone;
class R3Sphere;
class R3Matrix;



// Class include files

#include "R3Point.h"
#include "R3Vector.h"
#include "R3Line.h"
#include "R3Ray.h"
#include "R3Segment.h"
#include "R3Plane.h"
#include "R3Circle.h"
#include "R3Box.h"
#include "R3Cylinder.h"
#include "R3Cone.h"
#include "R3Sphere.h"
#include "R3Mesh.h"
#include "R3Matrix.h"



// Utility include files

#include "R3Distance.h"



#endif



