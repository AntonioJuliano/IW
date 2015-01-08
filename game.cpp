////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////

#include "R3/R3.h"
#include "R3Scene.h"
#include "particle.h"
#include "fglut/fglut.h"
#include <string>
#include <time.h>

////////////////////////////////////////////////////////////
// DEFINES
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS
////////////////////////////////////////////////////////////

// Display variables
static R3Scene *scene = NULL;
static R3Camera camera;
static double start_time = 0;
double fps;
int fps_time = 1, frame = 0, timebase = 0;

// GLUT variables
static int GLUTwindow = 0;
static int GLUTwindow_height;
static int GLUTwindow_width;
static int GLUTmodifiers = 0;

// Game varibles
static int CAST_RATE = 1 * 1000;

static bool player_cast, player_cast_success;
static int spell_num;

enum Spell {
  SPELL_NONE,
  FIREBALL,
  ICE_BLAST,
  LIGHTNING
};

static int current_spell = SPELL_NONE;
static R3ParticleSource spell_source;
static R3ParticleSink spell_sink;



// Initialize functions
R3Scene *InitializeScene();
void InitializeGame();
void castFireball();
void castIceBlast();
void castLightning();
void castEnemySpell(int val);
void DrawEnemy(R3Scene* scene);
void DrawSpell(R3Scene* scene);
void DrawParticles(R3Scene *scene);
void enemySpellHit(int val);
void playerSpellHit(int val);


////////////////////////////////////////////////////////////
// TIMER CODE
////////////////////////////////////////////////////////////

#ifdef _WIN32
#  include <windows.h>
#else
#  include <sys/time.h>
#endif

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



////////////////////////////////////////////////////////////
// RANDOM NUMBER GENERATOR
////////////////////////////////////////////////////////////

static double 
RandomNumber(void) 
{
#if defined(_WIN32)
  srand((unsigned int) clock());
  int r1 = rand();
  double r2 = ((double) rand()) / ((double) (RAND_MAX + 1));
  return (r1 + r2) / ((double) (RAND_MAX + 1));
#else
  srand48((long int) clock());
  return drand48();
#endif
}



////////////////////////////////////////////////////////////
// LOADING CODE
////////////////////////////////////////////////////////////

void LoadMatrix(R3Matrix *matrix)
{
  // Multiply matrix by top of stack
  // Take transpose of matrix because OpenGL represents vectors with 
  // column-vectors and R3 represents them with row-vectors
  R3Matrix m = matrix->Transpose();
  glMultMatrixd((double *) &m);
}

void LoadCamera(R3Camera *camera)
{
  // Set projection transformation
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(2*180.0*camera->yfov/M_PI, (GLdouble) GLUTwindow_width /(GLdouble) GLUTwindow_height, 0.01, 10000);

  // Set camera transformation
  R3Vector t = -(camera->towards);
  R3Vector& u = camera->up;
  R3Vector& r = camera->right;
  GLdouble camera_matrix[16] = { r[0], u[0], t[0], 0, r[1], u[1], t[1], 0, r[2], u[2], t[2], 0, 0, 0, 0, 1 };
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixd(camera_matrix);
  glTranslated(-(camera->eye[0]), -(camera->eye[1]), -(camera->eye[2]));
}

void LoadLights(R3Scene *scene)
{
  GLfloat buffer[4];

  // Load ambient light
  static GLfloat ambient[4];
  ambient[0] = scene->ambient[0];
  ambient[1] = scene->ambient[1];
  ambient[2] = scene->ambient[2];
  ambient[3] = 1;
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

  // Load scene lights
  for (int i = 0; i < (int) scene->lights.size(); i++) {
    R3Light *light = scene->lights[i];
    int index = GL_LIGHT0 + i;

    // Temporarily disable light
    glDisable(index);

    // Load color
    buffer[0] = light->color[0];
    buffer[1] = light->color[1];
    buffer[2] = light->color[2];
    buffer[3] = 1.0;
    glLightfv(index, GL_DIFFUSE, buffer);
    glLightfv(index, GL_SPECULAR, buffer);

    // Load attenuation with distance
    buffer[0] = light->constant_attenuation;
    buffer[1] = light->linear_attenuation;
    buffer[2] = light->quadratic_attenuation;
    glLightf(index, GL_CONSTANT_ATTENUATION, buffer[0]);
    glLightf(index, GL_LINEAR_ATTENUATION, buffer[1]);
    glLightf(index, GL_QUADRATIC_ATTENUATION, buffer[2]);

    // Load spot light behavior
    buffer[0] = 180.0 * light->angle_cutoff / M_PI;
    buffer[1] = light->angle_attenuation;
    glLightf(index, GL_SPOT_CUTOFF, buffer[0]);
    glLightf(index, GL_SPOT_EXPONENT, buffer[1]);

    // Load positions/directions
    if (light->type == R3_DIRECTIONAL_LIGHT) {
      // Load direction
      buffer[0] = -(light->direction.X());
      buffer[1] = -(light->direction.Y());
      buffer[2] = -(light->direction.Z());
      buffer[3] = 0.0;
      glLightfv(index, GL_POSITION, buffer);
    }

    // Enable light
    glEnable(index);
  }
}



//////////////////////////////////////////////////
// DRAWING FUNCTIONS
//////////////////////////////////////////////////

// Draws whole scene
void DrawScene(R3Scene *scene)
{
  DrawEnemy(scene);

  DrawSpell(scene);

  DrawParticles(scene);
}

void DrawEnemy(R3Scene* scene) {
  scene->enemy.Draw();
}

void DrawSpell(R3Scene* scene) {
  
}


void DrawParticles(R3Scene *scene)
{
  // Get current time (in seconds) since start of execution
  double current_time = GetTime();
  static double previous_time = 0;

  // printf("%i\n", scene->NParticles());

  // program just started up?
  if (previous_time == 0) previous_time = current_time;

  // time passed since starting
  double delta_time = current_time - previous_time;

  // Update particles 
  UpdateParticles(scene, current_time, delta_time, 0);

  // Generate new particles
  GenerateParticles(scene, current_time, delta_time);

  // Render particles
  RenderParticles(scene, current_time, delta_time);

  // Remember previous time
  previous_time = current_time;
}

//////////////////////////////////////////////////
// MOVEMENT
//////////////////////////////////////////////////


////////////////////////////////////////////////////////////
// DRAW HUD
////////////////////////////////////////////////////////////

void DrawHUD(R3Scene* scene)
{

  // Initialize
  glDisable(GL_DEPTH_TEST);
  glClear(GL_DEPTH_BUFFER_BIT);
  glDisable(GL_LIGHTING);
  glPushMatrix();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, GLUTwindow_width, GLUTwindow_height, 0, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Initialize text variables
  string s;
  void * font = GLUT_BITMAP_9_BY_15;

  // Draw Score
  glRasterPos2i(45,40);
  s = "Score: ";

  s += std::to_string((long double)(scene->score));
  for (string::iterator i = s.begin(); i != s.end(); ++i)
  {
    char c = *i;
    if(c == '.')
        break;
    glColor3d(1.0, 1.0, 1.0);
    glutBitmapCharacter(font, c);
  }

  // Print FPS Counter
  frame++;
  fps_time=glutGet(GLUT_ELAPSED_TIME);
  if (fps_time - timebase > 1000) {
    fps = (frame*1000.0 / (double) (fps_time-timebase));
    timebase = fps_time;
    frame = 0;
  }
  s = "FPS: ";
  s += std::to_string((long double) fps);
  glRasterPos2i(525,40);
  for (string::iterator i = s.begin(); i != s.end(); ++i)
  {
    char c = *i;
    if(c == '.')
        break;
    glColor3d(1.0, 1.0, 1.0);
    glutBitmapCharacter(font, c);
  }

  glPopMatrix();
  glEnable(GL_DEPTH_TEST);
}




////////////////////////////////////////////////////////////
// GLUT USER INTERFACE CODE
////////////////////////////////////////////////////////////

void GLUTMainLoop(void)
{
  glutMainLoop();
}

void GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Delete scene
  delete scene;

  // Exit
  exit(0);
}

void GLUTIdle(void)
{
  // Set current window
  if ( glutGetWindow() != GLUTwindow ) 
    glutSetWindow(GLUTwindow);  

  // Redraw
  glutPostRedisplay();
}

void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize camera vertical field of view to match aspect ratio of viewport
  camera.yfov = atan(tan(camera.xfov) * (double) h/ (double) w); 

  // Remember window size 
  GLUTwindow_width = w;
  GLUTwindow_height = h;

  // Redraw
  glutPostRedisplay();
}

void GLUTRedraw(void)
{
  // Initialize OpenGL drawing modes
  glEnable(GL_LIGHTING);
  glEnable(GL_BLEND);

  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(true);

  // Clear window 
  R3Rgb background = scene->background;
  
  glClearColor(background[0], background[1], background[2], background[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Load camera
  LoadCamera(&camera);

  // Load scene lights
  LoadLights(scene);

  // Draw scene surfaces
  DrawScene(scene);

  // Draw HUD
  DrawHUD(scene);

  // Swap buffers 
  glutSwapBuffers();
}

////////////////////////////////////////////////////////////
// GLUT MOUSE AND KEYBOARD
////////////////////////////////////////////////////////////

// Get key pressed
void GLUTKeyboard(unsigned char key, int, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 
  switch (key) {
    case 'w':
      if (player_cast && current_spell == FIREBALL) {
        player_cast_success = true;
        player_cast = false;
        playerSpellHit(0);
      }
      break;
    case 'z':
      if (player_cast && current_spell == ICE_BLAST) {
        player_cast_success = true;
        player_cast = false;
        playerSpellHit(0);
      }
      break;
    case 'o':
      if (player_cast && current_spell == LIGHTNING) {
        player_cast_success = true;
        player_cast = false;
        playerSpellHit(0);
      }
      break;
    case 'q':
      exit(0);
      break;
  }

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();
}

// GLUT Initialization
void GLUTInit(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  GLUTwindow = glutCreateWindow("Duel");

  // Initialize GLUT callback functions 
  glutIdleFunc(GLUTIdle);
  glutReshapeFunc(GLUTResize);
  glutDisplayFunc(GLUTRedraw);
  glutKeyboardFunc(GLUTKeyboard);

  // Initialize graphics modes 
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
}

////////////////////////////////////////////////////////////
// INITIALIZE GAME
////////////////////////////////////////////////////////////

void InitializeGame(void)
{
  start_time = GetTime();

  scene->particle_sources.push_back(&spell_source);

  scene->particle_sinks.push_back(&spell_sink);

  player_cast = false;
  spell_num = 0;
}

////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Initialize GLUT
  GLUTInit(&argc, argv);
  glutFullScreen();

  
  
  // Initialize scene
  scene = new R3Scene();
  if (!scene) exit(-1);
  camera = scene->camera;

  // Set GLUT
  GLUTwindow_height = glutGet(GLUT_WINDOW_HEIGHT);
  GLUTwindow_width = glutGet(GLUT_WINDOW_WIDTH);
  glutSetCursor(GLUT_CURSOR_NONE); 
  // glutWarpPointer(GLUTwindow_width / 2, GLUTwindow_height / 2);
  glutSetKeyRepeat(GLUT_KEY_REPEAT_OFF);

  // Initialize game settings
  InitializeGame();
  glutTimerFunc(CAST_RATE, castEnemySpell, 1000);

  // Run GLUT interface
  GLUTMainLoop();

  // Return success 
  return 0;
}

////////////////////////////////////////////////////////////
// GAME CODE
////////////////////////////////////////////////////////////

void castEnemySpell(int val) {
  scene->particles.clear();
  player_cast = true;
  player_cast_success = false;

  double r = RandomNumber();
  if (r < .33) {
    castFireball();
  }
  else if (r < .67) {
    castIceBlast();
  }
  else {
    castLightning();
  }

  glutTimerFunc(CAST_RATE * 5, enemySpellHit, ++spell_num);

  player_cast = true;
}

void enemySpellHit(int val) {
  if (player_cast_success)
    return;
  if (spell_num != val)
    return;

  player_cast = false;

  glutTimerFunc(CAST_RATE * 3, castEnemySpell, 1000);

  spell_sink.shape->sphere = new R3Sphere(R3Point(5, 0, 80), .001);

  spell_sink.intensity = 2;

}

void playerSpellHit(int val) {
  scene->particles.clear();

  glutTimerFunc(CAST_RATE * 3, castEnemySpell, 1000);

  spell_source.shape->sphere = new R3Sphere(R3Point(0, 0, 6), .01);
  spell_source.rate = 400;
  spell_source.velocity = 1;
  spell_source.angle_cutoff = 2 * 3.14;
  spell_source.mass = .1;
  spell_source.fixed = false;
  spell_source.drag = 0;
  spell_source.elasticity = 1;
  spell_source.lifetime = 2;
  spell_source.material = new R3Material();
  spell_source.material->kd.SetRed(.1);
  spell_source.material->kd.SetGreen(.9);
  spell_source.material->kd.SetBlue(.1);

  spell_sink.shape->sphere = new R3Sphere(R3Point(0, 0, 0), .001);

  spell_sink.intensity = 2;
}

void castFireball() {
  current_spell = FIREBALL;

  spell_source.shape = new R3Shape();
  spell_source.shape->type = R3_SPHERE_SHAPE;
  spell_source.shape->sphere = new R3Sphere(R3Point(-.41, .5, .1), .01);
  spell_source.rate = 400;
  spell_source.velocity = .3;
  spell_source.angle_cutoff = 2 * 3.14;
  spell_source.mass = .1;
  spell_source.fixed = false;
  spell_source.drag = 0;
  spell_source.elasticity = 1;
  spell_source.lifetime = 2;
  spell_source.material = new R3Material();
  spell_source.material->kd.SetRed(.9);
  spell_source.material->kd.SetGreen(.1);
  spell_source.material->kd.SetBlue(.1);

  spell_sink.shape = new R3Shape();
  spell_sink.shape->type = R3_SPHERE_SHAPE;
  spell_sink.shape->sphere = new R3Sphere(R3Point(-.41, .5, .1), .001);
  spell_sink.intensity = .06;
  spell_sink.constant_attenuation = 1;
  spell_sink.linear_attenuation = 0;
  spell_sink.quadratic_attenuation = 0;
}

void castIceBlast() {
  current_spell = ICE_BLAST;

  spell_source.shape = new R3Shape();
  spell_source.shape->type = R3_SPHERE_SHAPE;
  spell_source.shape->sphere = new R3Sphere(R3Point(-.41, .5, .1), .01);
  spell_source.rate = 400;
  spell_source.velocity = .3;
  spell_source.angle_cutoff = 2 * 3.14;
  spell_source.mass = .1;
  spell_source.fixed = false;
  spell_source.drag = 0;
  spell_source.elasticity = 1;
  spell_source.lifetime = 2;
  spell_source.material = new R3Material();
  spell_source.material->kd.SetRed(.1);
  spell_source.material->kd.SetGreen(.1);
  spell_source.material->kd.SetBlue(.9);

  spell_sink.shape = new R3Shape();
  spell_sink.shape->type = R3_SPHERE_SHAPE;
  spell_sink.shape->sphere = new R3Sphere(R3Point(-.41, .5, .1), .001);
  spell_sink.intensity = .06;
  spell_sink.constant_attenuation = 1;
  spell_sink.linear_attenuation = 0;
  spell_sink.quadratic_attenuation = 0;
}

void castLightning() {
  current_spell = LIGHTNING;

  spell_source.shape = new R3Shape();
  spell_source.shape->type = R3_SPHERE_SHAPE;
  spell_source.shape->sphere = new R3Sphere(R3Point(-.41, .5, .1), .01);
  spell_source.rate = 400;
  spell_source.velocity = .3;
  spell_source.angle_cutoff = 2 * 3.14;
  spell_source.mass = .1;
  spell_source.fixed = false;
  spell_source.drag = 0;
  spell_source.elasticity = 1;
  spell_source.lifetime = 2;
  spell_source.material = new R3Material();
  spell_source.material->kd.SetRed(.9);
  spell_source.material->kd.SetGreen(.9);
  spell_source.material->kd.SetBlue(.1);

  spell_sink.shape = new R3Shape();
  spell_sink.shape->type = R3_SPHERE_SHAPE;
  spell_sink.shape->sphere = new R3Sphere(R3Point(-.41, .5, .1), .001);
  spell_sink.intensity = .06;
  spell_sink.constant_attenuation = 1;
  spell_sink.linear_attenuation = 0;
  spell_sink.quadratic_attenuation = 0;
}




