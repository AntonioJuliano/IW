//////////////////////////////////////////////////
// INCLUDES AND DEFINES
//////////////////////////////////////////////////

#ifndef CH_INCLUDED
#define CH_INCLUDED

#include "R3/R3.h"

#define PROJ_RADIUS 0.01
#define SHOOT_DISTANCE 70



//////////////////////////////////////////////////
// STRUCTS
//////////////////////////////////////////////////

struct R3Scene;
struct R3Node;
struct R3Shape;
struct R3Intersection;

struct Projectile {
  R3Point position;
  R3Vector velocity;
  double lifetime;
  double start;
  double dmg;
  double radius;
};

struct R3Intersection {
  R3Intersection(void){
    hit = false;
  }
  R3Intersection(double _t, R3Point _position, R3Vector _normal){
    t = _t;
    position = _position;
    normal = _normal;
  }
  bool hit;
  R3Node *node;
  R3Point position;
  R3Vector normal;
  double t;
};



//////////////////////////////////////////////////
// CHARACTER DEFINITION
//////////////////////////////////////////////////

class Character
{
public:
  Character(void);
  Character(R3Scene* scene, R3Point pos, double r, double hp, double speed, double gun_rate, double gun_speed, double lifetime, double dmg);
  ~Character(void);
  R3Point GetPosition();
  void MainLoop(double current_time, double delta_time);
  R3Intersection Intersect(R3Ray ray);
  void UpdateHP(double delta);
  double max_hp;
  double hp;
  R3Vector velocity;
  R3Point destination;

  //move self location based on scene's protagonist character/gravity/collision detection
  virtual void UpdateSelf(double current_time, double delta_time);

  //update projectile location 
  virtual void UpdateProjectiles(double current_time, double delta_time);

  //draw self and projectiles
  virtual void Draw();

  R3Scene* scene;
  R3Shape* model;

  double speed;
  double gun_rate;
  double gun_speed;
  double lifetime;
  double dmg;
  vector<Projectile*> projectiles;
};



//////////////////////////////////////////////////
// MAIN CHARACTER DEFINITION
//////////////////////////////////////////////////

class MainCharacter : public Character 
{
public:
  MainCharacter(R3Scene* scene, R3Point pos, double r, double hp, double speed, double gun_rate, double gun_speed, double lifetime, double dmg);
  void Move(R3Vector delta);
  void MoveTo(R3Point newp);
  void ShootProjectile(R3Vector dir);
  virtual void UpdateSelf(double current_time, double delta_time);
  virtual void UpdateProjectiles(double current_time, double delta_time);
  virtual void Draw();
};

#endif
