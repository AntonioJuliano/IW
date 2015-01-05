// Source file for the particle system



// Include files

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"


////////////////////////////////////////////////////////////
// Random Number Generator
////////////////////////////////////////////////////////////

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

double PlaneIntersection(R3Plane *plane, R3Ray *ray)
{
  double a = plane->A();
  double b = plane->B();
  double c = plane->C();
  double d = plane->D();

  double x_1 = ray->Start().X();
  double y_1 = ray->Start().Y();
  double z_1 = ray->Start().Z();
  double x_2 = ray->Vector().X();
  double y_2 = ray->Vector().Y();
  double z_2 = ray->Vector().Z();

  double t = - (a * x_1 + b * y_1 + c * z_1 + d) / (a * x_2 + b*y_2 + c*z_2);

  return t;
}

R3Intersection BoxIntersection(R3Box *box, R3Ray *ray)
{
  R3Intersection intersection;

  R3Point p_0(box->XMax(), box->YMin(), box->ZMin());
  R3Point p_1(box->XMax(), box->YMax(), box->ZMin());
  R3Point p_2(box->XMin(), box->YMax(), box->ZMin());
  R3Point p_3(box->XMin(), box->YMin(), box->ZMin());
  R3Point p_4(box->XMax(), box->YMin(), box->ZMax());
  R3Point p_5(box->XMax(), box->YMax(), box->ZMax());
  R3Point p_6(box->XMin(), box->YMax(), box->ZMax());
  R3Point p_7(box->XMin(), box->YMin(), box->ZMax());

  R3Plane f_0(p_0, p_1, p_2);
  R3Plane f_1(p_0, p_1, p_5);
  R3Plane f_2(p_1, p_2, p_6);
  R3Plane f_3(p_2, p_3, p_6);
  R3Plane f_4(p_0, p_3, p_7);
  R3Plane f_5(p_4, p_5, p_6);
  vector<R3Plane *> faces;
  faces.push_back(&f_0);
  faces.push_back(&f_1);
  faces.push_back(&f_2);
  faces.push_back(&f_3);
  faces.push_back(&f_4);
  faces.push_back(&f_5);
  vector<R3Vector *> normals;
  R3Vector n_0(0, 0, -1);
  R3Vector n_1(1, 0, 0);
  R3Vector n_2(0, 1, 0);
  R3Vector n_3(-1, 0, 0);
  R3Vector n_4(0, -1, 0);
  R3Vector n_5(0, 0, 1);
  normals.push_back(&n_0);
  normals.push_back(&n_1);
  normals.push_back(&n_2);
  normals.push_back(&n_3);
  normals.push_back(&n_4);
  normals.push_back(&n_5);

  intersection.hit = false;

  for (int i = 0; i < (int) faces.size(); i++)
  {
    double t = PlaneIntersection(faces[i], ray);
    if (t <= 0 || (intersection.hit && t >= intersection.t))
      continue;
    R3Point plane_hit = ray->Point(t);

    double rounding_error = .0000001;
    if (plane_hit.X() >= box->XMin() - rounding_error && 
        plane_hit.X() <= box->XMax() + rounding_error &&
        plane_hit.Y() >= box->YMin() - rounding_error &&
        plane_hit.Y() <= box->YMax() + rounding_error &&
        plane_hit.Z() >= box->ZMin() - rounding_error &&
        plane_hit.Z() <= box->ZMax() + rounding_error)
    {
      intersection.t = t;
      intersection.position = plane_hit;
      intersection.hit = true;
      intersection.normal = *(normals[i]);
    }
  }

  return intersection;
}

R3Intersection SphereIntersection(R3Sphere *sphere, R3Ray *ray)
{
  R3Intersection intersection;
  R3Vector L = sphere->Center() - ray->Start();
  R3Vector V = ray->Vector();

  double t_ca = L.Dot(V);
  if (t_ca < 0)
  {
    intersection.hit = false;
    return intersection;
  }

  double r_2 = sphere->Radius() * sphere->Radius();

  double d_2 = L.Dot(L) - t_ca * t_ca;
  if (d_2 > r_2)
  {
    intersection.hit = false;
    return intersection;
  }

  double t_hc = sqrt(r_2 - d_2);
  double t = t_ca - t_hc;

  R3Point p = ray->Point(t);
  R3Vector normal = p - sphere->Center();
  normal.Normalize();

  intersection.hit = true;
  intersection.position = p;
  intersection.normal = normal;
  intersection.t = t;

  return intersection;
}

R3Intersection MeshIntersection(R3Mesh *mesh, R3Ray *ray)
{
  R3Intersection intersection;
  intersection.hit = false;

  for (int i = 0; i < (int) mesh->faces.size(); i++)
  {
    if (ray->Vector().Dot(mesh->faces[i]->plane.Normal()) == 0)
      continue;
    double t = PlaneIntersection(&(mesh->faces[i]->plane), ray);
    if (t < 0 || (intersection.hit && t >= intersection.t))
      continue;

    R3Point p = ray->Point(t);
    bool inside = true;

    for (int j = 0; j < (int) mesh->faces[i]->vertices.size() - 1; j++)
    {
      R3Point t_1 = mesh->faces[i]->vertices[j]->position;
      R3Point t_2 = mesh->faces[i]->vertices[j + 1]->position;

      R3Vector v_1 = t_1 - p;
      R3Vector v_2 = t_2 - p;

      R3Vector n_1 = v_2;
      n_1.Cross(v_1);
      n_1.Normalize();
      if (ray->Vector().Dot(n_1) < 0)
        inside = false;
    }

    R3Point t_1 = mesh->faces[i]->vertices[mesh->faces[i]->vertices.size() - 1]->position;
    R3Point t_2 = mesh->faces[i]->vertices[0]->position;

    R3Vector v_1 = t_1 - p;
    R3Vector v_2 = t_2 - p;

    R3Vector n_1 = v_2;
    n_1.Cross(v_1);
    n_1.Normalize();
    if (ray->Vector().Dot(n_1) < 0)
      inside = false;

    if (inside)
    {
      intersection.hit = true;
      intersection.position = p;
      intersection.normal = mesh->faces[i]->plane.Normal();
      intersection.t = t;

      return intersection;
    }
  }

  return intersection;;
}

R3Intersection ComputeIntersection(R3Node *node, R3Ray *ray, R3Intersection closest_intersection)
{
  if (node->shape != NULL)
  {
    switch (node->shape->type)
    {
      case R3_SPHERE_SHAPE:
      {
        R3Intersection intersection = SphereIntersection(node->shape->sphere, ray);

        if (intersection.hit && (!closest_intersection.hit || intersection.t < closest_intersection.t))
        {
          closest_intersection = intersection;
          closest_intersection.node = node;
        }
        break;
      }
      case R3_BOX_SHAPE:
      {
        R3Intersection intersection = BoxIntersection(node->shape->box, ray);

        if (intersection.hit && (!closest_intersection.hit || intersection.t < closest_intersection.t))
          {closest_intersection = intersection;
                    closest_intersection.node = node;}
        break;
      }
      case R3_MESH_SHAPE:
      {
        R3Intersection bounding_intersection = BoxIntersection(&(node->bbox), ray);
        if (!bounding_intersection.hit)
          break;

        R3Intersection intersection = MeshIntersection(node->shape->mesh, ray);

        if (intersection.hit && (!closest_intersection.hit || intersection.t < closest_intersection.t))
          {closest_intersection = intersection;
                    closest_intersection.node = node;}

        break;
      }

      default:
        break;
    }
  }

  for (int i = 0; i < (int) node->children.size(); i++)
  {
    R3Intersection closest = ComputeIntersection(node->children[i], ray, closest_intersection);
    
    if (closest.hit && (!closest_intersection.hit || closest.t < closest_intersection.t))
      closest_intersection = closest;
  }


  return closest_intersection;
}

////////////////////////////////////////////////////////////
// Generating Particles
////////////////////////////////////////////////////////////

void GenerateParticles(R3Scene *scene, double current_time, double delta_time)
{
  // Generate new particles for every source
  for (int i = 0; i < (int) scene->particle_sources.size(); i++)
  {

    R3ParticleSource *source = scene->particle_sources[i];
    int num_particles = (int) (source->rate * delta_time);
    double rem = source->rate * delta_time - (double) num_particles;
    if (RandomNumber() < rem)
      num_particles++;
 // printf("gen: %i, ", source->shape->type);
    for (int k = 0; k < num_particles; k++)
    {
    R3Vector N;
    R3Point p;
    if(source->shape->type == R3_SPHERE_SHAPE)
    {
      double z = RandomNumber();
      double r = source->shape->sphere->Radius();
      z = (z - .5) * 2 * r;
      double phi = RandomNumber() * 2 * 3.14159;
      double d = sqrt(r * r - d * d);
      p = source->shape->sphere->Center();
      N.Reset(d * cos(phi), d*sin(phi), z);
      p.Translate(N);
      N.Normalize();

    }
    if (source->shape->type == R3_CIRCLE_SHAPE)
    {
      N = source->shape->circle->Normal();
      R3Vector A(RandomNumber(), RandomNumber(), RandomNumber());
      A.Cross(N);
      A.Normalize();
      double r = RandomNumber() * source->shape->circle->Radius();
      double theta = RandomNumber() * 2 * 3.14159;
      A.Rotate(N, theta);
      A *= r;
      p = source->shape->circle->Center() + A;
    }
    // MAKE SURE THIS WORKS
      R3Vector A(RandomNumber(), RandomNumber(), RandomNumber());
      A.Cross(N);
      A.Normalize();
      double t1 = RandomNumber();
      t1 *= 2 * 3.14159;
      double t2 = RandomNumber();
      t2 *= sin(source->angle_cutoff);
      A.Rotate(N, t1);
      R3Vector VxN(A);
      VxN.Cross(N);
      A.Rotate(VxN, acos(t2));
      A.Normalize();


      R3Particle * new_particle = new R3Particle;
      new_particle->position = p;
      new_particle->velocity = A * source->velocity;
      new_particle->fixed = source->fixed;
      new_particle->drag = source->drag;
      new_particle->elasticity = source->elasticity;
      new_particle->lifetime = source->lifetime;
      new_particle->material = source->material;
      // new_particle->material = new R3Material();
      // R3Rgb source_color = source->material->kd;
      // double d1 = RandomNumber() * .2 - .1;
      // new_particle->material->kd.SetRed(source_color.Red() + d1);
      // double d2 = RandomNumber() * .2 - .1;
      // new_particle->material->kd.SetGreen(source_color.Green() + d2);
      // double d3 = RandomNumber() * .2 - .1;
      // new_particle->material->kd.SetBlue(source_color.Blue() + d3);
      // new_particle->material->kd.Clamp();
      new_particle->mass = source->mass;
      new_particle->infinite = (source->lifetime <= 0) ? true : false;

      scene->particles.push_back(new_particle);
    }
    //printf("\n");
  }
}



////////////////////////////////////////////////////////////
// Updating Particles
////////////////////////////////////////////////////////////

R3Vector f(R3Particle *particle, double delta_time, R3Scene *scene)
{
  R3Vector force(0, 0, 0);
  // Gravity
  force += particle->mass * scene->gravity;
  for (int i = 0; i < scene->NParticles(); i++)
  {
    if (scene->Particle(i) == particle)
      continue;
    R3Vector grav = scene->Particle(i)->position - particle->position;
    double r = grav.Length();
    if (r == 0)
      r = .00001;
    grav.Normalize();
    double magnitude = .0000000000667 * scene->Particle(i)->mass * particle->mass / (r * r);
    grav *= magnitude;
    force += grav;
  }
  // Drag
  force += particle->velocity * (-1 * particle->drag);
  // Sinks
  for (int j = 0; j < scene->NParticleSinks(); j++)
  {
    // Sphere Sink
    R3ParticleSink *sink = scene->ParticleSink(j);
    if (sink->shape->type == R3_SPHERE_SHAPE)
    {
      double d = (particle->position - sink->shape->sphere->Center()).Length() - sink->shape->sphere->Radius();
      
      double magnitude = sink->intensity / (sink->constant_attenuation + d * 
        sink->linear_attenuation + d * d * sink->quadratic_attenuation);
      R3Vector f_sink(sink->shape->sphere->Center() - particle->position);
      f_sink.Normalize();
      f_sink *= magnitude;
      force += f_sink;
    }
  }
  // Springs
  for (int j = 0; j < (int) particle->springs.size(); j++)
  {
    R3Particle *other_particle;
    if (particle->springs[j]->particles[0] == particle)
    {
      other_particle = particle->springs[j]->particles[1];
    }
    else
    {
      other_particle = particle->springs[j]->particles[0];
    }
    double d = (particle->position - other_particle->position).Length();
    R3Vector D = other_particle->position - particle->position;
    D.Normalize();
    double magnitude = (particle->springs[j]->ks * (d - particle->springs[j]->rest_length) + 
      (particle->springs[j]->kd * (other_particle->velocity - particle->velocity)).Dot(D));
    
    D *= magnitude;
    force += D;
  }

  return force;
}



bool UpdateParticle(R3Particle *particle, R3Scene *scene, double current_time, double delta_time, int integration_type)
{
  // Lifetime Check
  if (!particle->infinite && particle->lifetime <= 0)
  {
    return true;
  }
  if (particle->fixed)
  {
    return false;
  }
  // Collision Check
  R3Point new_position, old_position;
  old_position = particle->position;
  R3Vector old_velocity, new_velocity;
  old_velocity = particle->velocity;

  if (integration_type == MIDPOINT_INTEGRATION)
  {
    R3Point x_mid = particle->position + (delta_time / 2) * particle->velocity;
    R3Vector v_mid = particle->velocity + (delta_time / 2) * f(particle, delta_time, scene) / particle->mass;
    particle->position = x_mid;
    particle->velocity = v_mid;
    R3Vector f_mid = f(particle, delta_time, scene);
    new_position = old_position + delta_time * v_mid;
    new_velocity = old_velocity + delta_time * f_mid / particle->mass;
  }
  // default to euler integration
  else 
  {
    new_position = old_position + particle->velocity * delta_time;
    new_velocity = old_velocity + delta_time * f(particle, delta_time, scene) / particle->mass;
  }

  
  R3Vector ray_vector = new_position - old_position;

  R3Ray path(old_position, ray_vector);
  R3Intersection noHitIntersection;
  noHitIntersection.hit = false;
  // R3Intersection intersection = ComputeIntersection(scene->root, &path, noHitIntersection);

  // Sink Intersections
  for (int k = 0; k < scene->NParticleSinks(); k++)
  {
    R3ParticleSink *sink = scene->ParticleSink(k);
    if (sink->shape->type == R3_SPHERE_SHAPE)
    {
      R3Intersection intersection = SphereIntersection(sink->shape->sphere, &path);
      if (intersection.hit && intersection.t <= ray_vector.Length())
      {
        return true;
      }
    }
  }

  // if (intersection.hit && intersection.t <= ray_vector.Length())
  // {
  //   intersection.normal.Normalize();

  //   particle->position = intersection.position + (.001 * intersection.normal);
  //   R3Vector v_n(intersection.normal);
  //   v_n *= -1 * particle->velocity.Dot(intersection.normal);
  //   R3Vector v_t(particle->velocity + v_n);
  //   v_n *= particle->elasticity;
  //   particle->velocity = v_n + v_t;
  //   double dt = intersection.t / ray_vector.Length() * delta_time;
  //   particle->lifetime -= dt;
  //   return UpdateParticle(particle, scene, current_time, delta_time - dt, integration_type);
  // }



  particle->position = new_position;
  particle->velocity = new_velocity;
  
  particle->lifetime -= delta_time;
  return false;
}

void UpdateParticles(R3Scene *scene, double current_time, double delta_time, int integration_type)
{
  // Update position for every particle
  /*for (double d = 0; d < delta_time; d+=delta_time/10)
  {*/
    for (int i = 0; i < scene->NParticles(); i++)
    {
      if (UpdateParticle(scene->Particle(i), scene, current_time, delta_time, integration_type))
      {
        delete scene->Particle(i);
        scene->particles.erase(scene->particles.begin() + i);
        i--;
      }
    }
 // }
}

////////////////////////////////////////////////////////////
// Rendering Particles
////////////////////////////////////////////////////////////

void RenderParticles(R3Scene *scene, double current_time, double delta_time)
{
  double radius = .01;
  const double PI = 3.14159;
  // Draw every particle
  // REPLACE CODE HERE
  glDisable(GL_LIGHTING);
  glPointSize(5);
  glBegin(GL_POINTS);
  for (int i = 0; i < scene->NParticles(); i++) {
    R3Particle *particle = scene->Particle(i);
    glColor4f(particle->material->kd[0], particle->material->kd[1], particle->material->kd[2], .5);
    const R3Point& position = particle->position;
    // glVertex3d(position[0], position[1], position[2]);
    glPushMatrix();
    glTranslatef(position[0], position[1], position[2]);
    glBegin(GL_POLYGON);
        for(double i = 0; i < 2 * PI; i += PI / 6) //<-- Change this Value
          glVertex3f(cos(i) * radius, sin(i) * radius, 0.0);
      glEnd();
    glPopMatrix();
  }   
  glEnd();
}



