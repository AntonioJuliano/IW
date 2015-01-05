// Include files 
#include "R3/R3.h"
#include "R3Scene.h"
#include <math.h>
#include "time.h"
#include "limits.h"

// Random number generator
static double 
RandomNumber(void) 
{
#if defined(_WIN32)
  int r1 = rand();
  double r2 = ((double) rand()) / ((double) (RAND_MAX + 1));
  return (r1 + r2) / ((double) (RAND_MAX + 1));
#else
  return drand48();
#endif
}


// Constructor
R3Scene::R3Scene(void)
  : gravity(0,0,0),
    bbox(R3null_box),
    // background(109./255., 192./255., 247./255., 1),
    background(0, 0, 0, 1),
    ambient(0,0,0,1)
{
  // Initialize score
  score = 0;

  // Setup sound
  // sound = Sound();
  // sound.setListeningPosition(R3Point(0, .5, 8), R3Vector(0, 0, -1), R3Vector(0, 1, 0));

  // Setup default camera
  R3Point scene_center = bbox.Centroid();
  camera.eye = R3Point(0, .5, 8);
  camera.towards = R3Vector(0, 0, -1);
  camera.up = R3Vector(0, 1, 0);
  camera.right = R3Vector(1, 0, 0);
  camera.xfov = 0.25;
  camera.yfov = 0.25;
  camera.neardist = 0.01;
  camera.fardist = 100;

  // Create root node
  root = new R3Node();
  root->parent = NULL;
  root->transformation = R3identity_matrix;
  root->material = NULL;
  root->shape = NULL;
  root->bbox = R3null_box;

  // Create array of materials
  vector<R3Material *> materials;

  // Create first directional light
  R3Light *light = new R3Light();
  R3Vector direction(0,0,-1);
  direction.Normalize();
  light->type = R3_DIRECTIONAL_LIGHT;
  light->color = R3Rgb(.8,.8,.8,1);
  light->position = R3Point(0, 0, 0);
  light->direction = direction;
  light->radius = 0;
  light->constant_attenuation = 0;
  light->linear_attenuation = 0;
  light->quadratic_attenuation = 0;
  light->angle_attenuation = 0;
  light->angle_cutoff = M_PI;
  lights.push_back(light);

  enemy.Read("enemy.off");
}
