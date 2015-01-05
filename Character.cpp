//////////////////////////////////////////////////
// INCLUDES AND DEFINES
//////////////////////////////////////////////////

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "Character.h"

#ifdef _WIN32
#  include <windows.h>
#else
#  include <sys/time.h>
#endif

#define MU 0.0001



//////////////////////////////////////////////////
// RANDOM NUMBER GENERATOR
//////////////////////////////////////////////////

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



//////////////////////////////////////////////////
// SYSTEM TIME
//////////////////////////////////////////////////

static double GetTime(void)
{
#ifdef _WIN32
  // Return number of seconds since start of execution
  static int first = 1;
  static LARGE_INTEGER timefreq;
  static LARGE_INTEGER start_timevalue;

  // Check if this is the first time
  if (first) {
    // Initialize first time
    QueryPerformanceFrequency(&timefreq);
    QueryPerformanceCounter(&start_timevalue);
    first = 0;
    return 0;
  }
  else {
    // Return time since start
    LARGE_INTEGER current_timevalue;
    QueryPerformanceCounter(&current_timevalue);
    return ((double) current_timevalue.QuadPart - 
            (double) start_timevalue.QuadPart) / 
            (double) timefreq.QuadPart;
  }
#else
  // Return number of seconds since start of execution
  static int first = 1;
  static struct timeval start_timevalue;

  // Check if this is the first time
  if (first) {
    // Initialize first time
    gettimeofday(&start_timevalue, NULL);
    first = 0;
    return 0;
  }
  else {
    // Return time since start
    struct timeval current_timevalue;
    gettimeofday(&current_timevalue, NULL);
    int secs = current_timevalue.tv_sec - start_timevalue.tv_sec;
    int usecs = current_timevalue.tv_usec - start_timevalue.tv_usec;
    return (double) (secs + 1.0E-6F * usecs);
  }
#endif
}


//////////////////////////////////////////////////
// CHARACTER CONSTRUCTORS
//////////////////////////////////////////////////

// Default constructor
Character::Character(void)
{
}

// Main Character constructor
Character::
Character(R3Scene* scene, R3Point pos, double r, double hp, double speed, double gun_rate, double gun_speed, double lifetime, double dmg)
  : scene(scene),
    max_hp(hp),
    hp(hp),
    speed(speed),
    gun_rate(gun_rate),
    gun_speed(gun_speed),
    lifetime(lifetime),
    dmg(dmg)
{ 
  model = new R3Shape();
  model->sphere = new R3Sphere(pos, r);
  model->mesh = new R3Mesh();

  model->mesh->Read("Duke.off");

  model->type = R3_MESH_SHAPE;
}

// Destructor
Character::~Character(void)
{
  for(int i = projectiles.size()-1; !projectiles.empty(); i--)
  {
    delete projectiles[i];
    projectiles.pop_back();
  }
}



//////////////////////////////////////////////////
// CHARACTER FUNCTIONS
//////////////////////////////////////////////////

// Get position of character
R3Point Character::GetPosition()
{
  R3Point p(0,0,0);
    p = model->sphere->Center();
  return p;
}

// Main loop of character that updates the character and its projectiles
void Character::MainLoop(double current_time, double delta_time)
{
  UpdateSelf(current_time, delta_time);
  UpdateProjectiles(current_time, delta_time);
}

// Gets intersection between a character and a ray
R3Intersection Character::Intersect(R3Ray ray)
{
  R3Intersection intersect;
  R3Sphere sphere = *(model->sphere);
  intersect.hit = false;

  R3Vector L = sphere.Center() - ray.Start();
  double tca = L.Dot(ray.Vector());

  if(tca < MU) 
    return intersect;

  double d2 = L.Dot(L) - tca*tca;
  if(d2 > sphere.Radius()*sphere.Radius()) 
    return intersect;

  double thc = sqrt(sphere.Radius()*sphere.Radius() - d2);
  double t1 = tca - thc;
  double t2 = tca + thc;
  
  intersect.t = min(t1, t2);
  if(intersect.t < 0)
    intersect.t = max(t1, t2);

  if(intersect.t > 0)
  {
    intersect.hit = true;
    intersect.position = ray.Point(intersect.t);
    intersect.normal = intersect.position - sphere.Center();
    intersect.normal.Normalize();
  }

  return intersect;
}

// Updates the HP by delta
void Character::UpdateHP(double delta)
{ 
  if(hp+delta>=0)
    hp += delta;
  else
    hp = 0;
}

// Updates the character
void Character::UpdateSelf(double, double delta_time)
{
  R3Point mcPos = scene->mainChar->GetPosition();
  R3Point c = model->sphere->Center();
  double r = model->sphere->Radius();
  double distToPlayer = R3Distance(c, mcPos);
  R3Vector toPlayer = mcPos - c;
  toPlayer.Normalize();

  /*if(!(scene->mainChar->hp > 0))
  {*/
    velocity = R3Vector(0,0,0);
  //}
  if(distToPlayer <= TERRAIN_DENSITY * 12.0) // too close
  {
    velocity = -toPlayer;
  }
  else if (distToPlayer <= TERRAIN_DENSITY * 30.0) // kind of close
  {
    double rand = RandomNumber();
    if(rand <= 0.01){
      R3Vector randVector = R3Vector(RandomNumber(), RandomNumber(), 0);
      randVector.Normalize();
      randVector *= (12 + 18 * RandomNumber());
      destination = mcPos + randVector;
    }
    velocity = destination - c;
  }
  else // far away
    velocity = toPlayer;

  velocity.SetZ(0);
  velocity.Normalize();
  model->sphere->Translate(velocity * speed * delta_time);
  if (!scene->LocationValid(GetPosition().X(), GetPosition().Y()))
    model->sphere->Translate(-velocity * speed * delta_time);

  c = model->sphere->Center();
  c.SetZ(scene->getRealZValue(c)+r);
  model->sphere->Reposition(c);

  //shoot gun
  r = RandomNumber();
  double dist = (c - mcPos).Length();
  if(dist < SHOOT_DISTANCE && r <= delta_time * gun_rate )
  {
    scene->sound.playEnemyShoot(c);
    Projectile* p = new Projectile();
    p->lifetime = lifetime;
    p->start = GetTime();
    p->position = model->sphere->Center();
    R3Vector v = (mcPos) - p->position;
    v.Normalize();
    v *= .3;
    p->position.Translate(v + 0.55*.4*(R3Vector(0,0,0.1) + (R3Vector(0,0,-1) % v)));
    v = (mcPos) - p->position;
    v.Normalize();
    v *= gun_speed;
    p->velocity = v;
    p->dmg = dmg;
    p->radius = PROJ_RADIUS;
    projectiles.push_back(p);
  }
  
}

//////////////////////////////////////////////////
// PROJECTILE FUNCTIONS
//////////////////////////////////////////////////

// Copy projectile p1 into projectile p2
static void cpy(Projectile* p1, Projectile* p2)
{
  p1->dmg = p2->dmg;
  p1->lifetime = p2->lifetime;
  p1->position = p2->position;
  p1->radius = p2->radius;
  p1->start = p2->start;
  p1->velocity = p2->velocity;
}

// Update the characters projectiles
void Character::UpdateProjectiles(double, double delta_time)
{
  //for each projectile, intersect each other character
  //you hit something, you decrease its health and delete this projectile
  //if you hit terrain, delete this projectile
  //if lifetime runs out delete this projectile

  int n = projectiles.size();
  vector<R3Point*> newp;
  vector<int> toDelete;

  for(int i = 0; i < n; ++i)
  {
    Projectile* p = projectiles[i];

    //run out of life
    if((p->lifetime > 0 && GetTime() - p->start > p->lifetime) || R3Distance(p->position, model->sphere->Center()) > 100)
    {
      toDelete.push_back(i);
      newp.push_back(new R3Point(0,0,0));
      continue;
    }

    //out-of-bounds
    double x_upper = TERRAIN_DENSITY * (TERRAIN_X_SIZE-1);
    double y_upper = TERRAIN_DENSITY * (TERRAIN_Y_SIZE-1);
    R2Point np(p->position.X(), p->position.Y());
    if(np.X() < PROJ_RADIUS || np.Y() < PROJ_RADIUS || np.X() > x_upper - PROJ_RADIUS || np.Y() > y_upper - PROJ_RADIUS)
    {
      toDelete.push_back(i);
      newp.push_back(new R3Point(0,0,0));
      continue;
    }

    //hit other chars and do dmg
    bool del = false;
    R3Ray r(p->position, p->velocity);
    Character* ch = scene->mainChar;
    R3Intersection intersect = ch->Intersect(r);
    if(ch->hp > 0)
    {
      if(intersect.hit && intersect.t <= delta_time * p->velocity.Length() 
        || R3Distance(p->position+delta_time*p->velocity, 
          ch->model->sphere->Center()) < ch->model->sphere->Radius())
      {
        ch->UpdateHP(-p->dmg);
        scene->sound.playPlayerOw();
        if(ch->hp <= 0)
          scene->sound.playGameOver();
        del = true;
      }
    }

    //hit terrain
    if(scene->getRealZValue(p->position) >= p->position.Z())
      del = true;

    //integration or del
    if(del)
    {
      toDelete.push_back(i);
      newp.push_back(new R3Point(0,0,0));
      continue;
    }
    else
    {
      //no drag
      newp.push_back(new R3Point(p->position + delta_time * p->velocity));
    }
  }

  //replace
  for(int i = 0; i < n; ++i)
  {
    projectiles[i]->position = *(newp[i]);
  }

  //delete
  if(toDelete.size() > 0)
  {
    sort(toDelete.begin(), toDelete.end(), std::greater<int>());
    for(unsigned int i = 0; i < toDelete.size(); ++i)
    {
      Projectile* temp = projectiles[toDelete[i]];
      cpy(projectiles[toDelete[i]], projectiles[projectiles.size()-1]);
      cpy(projectiles[projectiles.size()-1], temp);
      projectiles.pop_back();
    }
  }
}

// Draw the character and its projectiles
void Character::Draw()
{
  //dont draw dead chars
  if(hp <= 0)
    return;

  float no_mat[] = {0.0f, 0.0f, 0.0f, 1.0f};
  float full_mat[] = {0.5f, 0.5f, 0.5f, 1.0f};
  float mat_ambient_color[] = {0.8f, 0.0f, 0.0f, 1.0f};

  if(model->type == R3_SPHERE_SHAPE)
  {
    float hp_p = hp/max_hp;
    float mat_ambient[] = {1.1f-hp_p, 1.1f-hp_p, 1.1f-hp_p, 1.0f};
    glPushAttrib(GL_LIGHTING_BIT);
    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_ambient);
    glMaterialfv(GL_FRONT, GL_SPECULAR, full_mat);
    glMaterialf(GL_FRONT, GL_SHININESS, 2);
    glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);
    model->sphere->Draw();
    glPopAttrib();
  }
  if (model->type == R3_MESH_SHAPE)
  {
    glPushMatrix();
    R2Point center(GetPosition().X(), GetPosition().Y());
    R2Point main_center(scene->mainChar->GetPosition().X(), 
      scene->mainChar->GetPosition().Y());
    R2Vector a(0, -1);
    R2Vector b = center - main_center;
    b.Normalize();
    double theta = acos(a.Dot(b));
    theta *= 57.2957795;

    if (b.X() < 0)
      theta = 360 - theta;

    glTranslatef(center.X(), center.Y(), 
      scene->getRealZValue(center.X(), center.Y()));
    glRotatef(theta, 0, 0, 1);
    glRotatef(180, 0, 0, 1);
    glRotatef(90, 1, 0, 0);
    glScalef(.1f, .1f, .1f);

    model->mesh->Draw();
    glPopMatrix();
  }

  //draw projectiles as spheres
  int n = projectiles.size();
  for(int i = 0; i < n; ++i)
  {
    // Draw sphere
    R3Point center = projectiles[i]->position;
    //if(R3Distance(center, eye) > 75 && (eye-center).Dot(projectiles[i]->velocity) < 0) {continue;}
    double radius = projectiles[i]->radius;
    glPushMatrix();
    glTranslated(center[0], center[1], center[2]);
    static GLUquadricObj *glu_sphere = gluNewQuadric();
    gluQuadricTexture(glu_sphere, GL_TRUE);
    gluQuadricNormals(glu_sphere, (GLenum) GLU_SMOOTH);
    gluQuadricDrawStyle(glu_sphere, (GLenum) GLU_FILL);
    glPushAttrib(GL_LIGHTING_BIT);
    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient_color);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_ambient_color);
    glMaterialfv(GL_FRONT, GL_SPECULAR, no_mat);
    glMaterialf(GL_FRONT, GL_SHININESS, 0);
    glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);
    gluSphere(glu_sphere, radius, 32, 32);
    glPopAttrib();
    glPopMatrix();
  }
}

//////////////////////////////////////////////////
// MAIN CHARACTER
//////////////////////////////////////////////////

// Default constructor
MainCharacter::MainCharacter(R3Scene* scene, R3Point pos, double r, double hp, double speed, double gun_rate, double gun_speed, double lifetime, double dmg)
  : Character(scene, pos, r, hp, speed, gun_rate, gun_speed, lifetime, dmg)
{
  model->type = R3_SPHERE_SHAPE;
}

// Move the main character by a vector
void MainCharacter::Move(R3Vector delta)
{
  if(model->type != R3_SPHERE_SHAPE)
    return;
  
  double x_upper = TERRAIN_DENSITY * (TERRAIN_X_SIZE-1);
  double y_upper = TERRAIN_DENSITY * (TERRAIN_Y_SIZE-1);
  R3Sphere* sphere = model->sphere;
  R2Point np((sphere->Center() + delta).X(), (sphere->Center() + delta).Y());
  
  if(np.X() < sphere->Radius() || np.Y() < sphere->Radius() || np.X() > x_upper - sphere->Radius() || np.Y() > y_upper - sphere->Radius())
    return;

  delta.SetZ(0);
  sphere->Translate(delta);
}

// Move the main character to a point
void MainCharacter::MoveTo(R3Point newp)
{
  if(model->type != R3_SPHERE_SHAPE)
    return;

  model->sphere->Reposition(R3Point(newp.X(), newp.Y(), model->sphere->Radius() + scene->getRealZValue(newp)));
}

// Shoot a new projectile
void MainCharacter::ShootProjectile(R3Vector dir)
{
  Projectile* p = new Projectile();
  scene->sound.playPlayerShoot();
  p->lifetime = lifetime;
  p->start = GetTime();
  p->position = model->sphere->Center();
  R3Vector v = dir;
  v.Normalize();
  v *= model->sphere->Radius();
  p->position.Translate(v);
  v.Normalize();
  v *= gun_speed;
  p->velocity = v;
  p->dmg = dmg;
  p->radius = PROJ_RADIUS;
  projectiles.push_back(p);
}

// Update the main character
void MainCharacter::UpdateSelf(double, double)
{
  if(model->type == R3_SPHERE_SHAPE)
  {
    //drop us to terrain height
    R3Point c = model->sphere->Center();
    c.SetZ(scene->getRealZValue(c)+model->sphere->Radius());
    model->sphere->Reposition(c);
  }
}

// Update the projectiles of the main character
void MainCharacter::UpdateProjectiles(double, double delta_time)
{
  int n = projectiles.size();
  vector<R3Point*> newp;
  vector<int> toDelete;

  for(int i = 0; i < n; ++i)
  {
    Projectile* p = projectiles[i];

    //run out of life
    if(p->lifetime > 0 && GetTime() - p->start > p->lifetime)
    {
      toDelete.push_back(i);
      newp.push_back(new R3Point(0,0,0));
      continue;
    }

    //out-of-bounds
    double x_upper = TERRAIN_DENSITY * (TERRAIN_X_SIZE-1);
    double y_upper = TERRAIN_DENSITY * (TERRAIN_Y_SIZE-1);
    R2Point np(p->position.X(), p->position.Y());
    if(np.X() < PROJ_RADIUS || np.Y() < PROJ_RADIUS || np.X() > x_upper - PROJ_RADIUS || np.Y() > y_upper - PROJ_RADIUS)
    {
      toDelete.push_back(i);
      newp.push_back(new R3Point(0,0,0));
      continue;
    }

    //hit other chars and do dmg
    R3Ray r(p->position, p->velocity);
    int nc = scene->characters.size();
    bool del = false;
    for(int j = 0; j < nc; j++)
    {
      Character* ch = scene->characters[j];
      R3Intersection intersect = ch->Intersect(r);
      if(intersect.hit && intersect.t <= delta_time * p->velocity.Length())
      {
        ch->UpdateHP(-p->dmg);
        scene->sound.playEnemyOw(ch->GetPosition());
        del = true;
        break;
      }
      else if(R3Distance(p->position+delta_time*p->velocity, ch->model->sphere->Center()) < ch->model->sphere->Radius())
      {
        ch->UpdateHP(-p->dmg);
        scene->sound.playEnemyOw(ch->GetPosition());
        del = true;
        break;
      }
    }

    //hit terrain
    if(scene->getRealZValue(p->position) >= p->position.Z())
      del = true;

    //integration or del
    if(del)
    {
      toDelete.push_back(i);
      newp.push_back(new R3Point(0,0,0));
      continue;
    }
    else
    {
      newp.push_back(new R3Point(p->position + delta_time * p->velocity));
    }
  }

  //replace
  for(int i = 0; i < n; ++i)
  {
    projectiles[i]->position = *(newp[i]);
  }

  //delete
  if(toDelete.size() > 0)
  {
    sort(toDelete.begin(), toDelete.end(), std::greater<int>());
    for(unsigned int i = 0; i < toDelete.size(); ++i)
    {
      Projectile* temp = projectiles[toDelete[i]];
      cpy(projectiles[toDelete[i]], projectiles[projectiles.size()-1]);
      cpy(projectiles[projectiles.size()-1], temp);
      projectiles.pop_back();
    }
  }
}

// Draw the main character and its projectiles
void MainCharacter::Draw()
{
  float no_mat[] = {0.0f, 0.0f, 0.0f, 1.0f};
  float mat_ambient_color[] = {0.0f, 0.0f, 0.8f, 1.0f};
  //draw projectiles as spheres
  int n = projectiles.size();
  for(int i = 0; i < n; ++i)
  {
    // Draw sphere
    R3Point center = projectiles[i]->position;
    //if(R3Distance(center, eye) > 75 && (eye-center).Dot(projectiles[i]->velocity) < 0) {continue; }
    double radius = projectiles[i]->radius;
    glPushMatrix();
    glTranslated(center[0], center[1], center[2]);
    static GLUquadricObj *glu_sphere = gluNewQuadric();
    gluQuadricTexture(glu_sphere, GL_TRUE);
    gluQuadricNormals(glu_sphere, (GLenum) GLU_SMOOTH);
    gluQuadricDrawStyle(glu_sphere, (GLenum) GLU_FILL);
    glPushAttrib(GL_LIGHTING_BIT);
    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient_color);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_ambient_color);
    glMaterialfv(GL_FRONT, GL_SPECULAR, no_mat);
    glMaterialf(GL_FRONT, GL_SHININESS, 0);
    glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);
    gluSphere(glu_sphere, radius, 32, 32);
    glPopAttrib();
    glPopMatrix();
  }
}