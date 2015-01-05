////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////

#include "R3/R3.h"
#include "R3Scene.h"
#include "fglut/fglut.h"
#include <string>
#include <time.h>

////////////////////////////////////////////////////////////
// DEFINES
////////////////////////////////////////////////////////////

#define SHOOT_RATE 300
#define GL_GENERATE_MIPMAP 0x8191
#define CHARACTER_GEN_RATE 7
#define PLAYER_SPEED 2.7

////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS
////////////////////////////////////////////////////////////

// Video Frame Delay
static const double VIDEO_FRAME_DELAY = 1./100.;

// Program arguments
static char *output_image_name = NULL;
static const char *video_prefix = "./video-frames/";

// Materials
TerrainType currentTerrainType = TERRAIN_NONE;
bool current_blood = false;;

// Display variables
static R3Scene *scene = NULL;
static R3Camera camera;
static int save_image = 0;
static int save_video = 0;
static int num_frames_to_record = -1; 
static int quit = 0;
static double draw_distance = 75;
static double start_time = 0;
static double last_time = 0;
int fps_time = 1, frame = 0, timebase = 0;

// GLUT variables
static int GLUTwindow = 0;
static int GLUTwindow_height;
static int GLUTwindow_width;
static int GLUTmouse[4] = { 0, 0, 0, 0};
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;
static bool mouse_reset = false;

// Key states (true is down)
bool* keyStates = (bool*) calloc(256, sizeof(bool)); 
bool* movedAlready = (bool*) calloc(256, sizeof(bool)); 
double fps = 0;
double lastGen = 0;

// GLUT command list
enum {
  DISPLAY_FACE_TOGGLE_COMMAND,
  DISPLAY_EDGE_TOGGLE_COMMAND,
  DISPLAY_BBOXES_TOGGLE_COMMAND,
  DISPLAY_LIGHTS_TOGGLE_COMMAND,
  DISPLAY_CAMERA_TOGGLE_COMMAND,
  SAVE_IMAGE_COMMAND,
  SAVE_VIDEO_COMMAND,
  QUIT_COMMAND,
};

// Initialize functions
R3Scene *InitializeScene();
void InitializeGame();
void moveMainCharacter(double delta_time);


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

// Set the drawing material to that of the enemy material
void SetEnemyColor()
{
  R2Pixel pixel = R2Pixel(0.8f,0.2f,0.2f,1.0f);
  GLfloat c[4];  GLfloat d[4];  GLfloat e[4];
  c[0] = pixel.Red(); c[1] = pixel.Green(); c[2] = pixel.Blue(); c[3] = pixel.Alpha();
  glColor3f(pixel.Red(),pixel.Green(),pixel.Blue());
  d[0]=c[0]/4; d[1]=c[1]/4; d[2]=c[2]/4; d[3]=c[3];
  e[0]=0; e[1]=0; e[2]=0; e[3]=0;

  // Load Color
  glMaterialfv(GL_FRONT, GL_AMBIENT, e);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, c);
  glMaterialfv(GL_FRONT, GL_SPECULAR, d);
  glMaterialfv(GL_FRONT, GL_EMISSION, e);
  glMaterialf(GL_FRONT, GL_SHININESS, 3);
}

// Draw all characters in the scene
void DrawCharacters(R3Scene* scene)
{
  // Get current time (in seconds) since start of execution
  double current_time = GetTime();
  static double previous_time = 0;

  static double time_lost_taking_videos = 0; // for switching back and forth
					     // between recording and not
					     // recording smoothly

  // program just started up?
  if (previous_time == 0) previous_time = current_time;

  // time passed since starting
  double delta_time = current_time - previous_time;


  if (save_video) { // in video mode, the time that passes only depends on the frame rate ...
    delta_time = VIDEO_FRAME_DELAY;    
    // ... but we need to keep track how much time we gained and lost so that we can arbitrarily switch back and forth ...
    time_lost_taking_videos += (current_time - previous_time) - VIDEO_FRAME_DELAY;
  } else { // real time simulation
    delta_time = current_time - previous_time;
  }

  if(scene->mainChar && scene->characters.size() > 0)
  {
    // Set the color of the enemies
    SetEnemyColor();

    // Move Main Character
    moveMainCharacter(delta_time);
    movedAlready['w'] = false;
    movedAlready['a'] = false;
    movedAlready['s'] = false;
    movedAlready['d'] = false;

    //handle main character
    scene->mainChar->MainLoop(current_time - time_lost_taking_videos, delta_time);

    //handle enemy characters
    int n = scene->characters.size();
    for(int i = 0; i < n; i++)
      scene->characters[i]->MainLoop(current_time - time_lost_taking_videos, delta_time);
    for(int i = 0; i < n; i++)
    {

      R3Vector towards = scene->characters[i]->GetPosition() - camera.eye;
      if (towards.Length() > draw_distance)
        continue;
      R2Vector a(towards.X(), towards.Y());
      R2Vector b(camera.towards.X(), camera.towards.Y());
      if (a.Dot(b) < 0)
        continue;

      a.Normalize();
      b.Normalize();
      double theta = acos(a.Dot(b));
      if (abs(theta) > camera.xfov + .1)
        continue;

      scene->characters[i]->Draw();
    }

    scene->mainChar->Draw();

    //delete dead enemies
    vector<int> toDelete;
    for(int i = 0; i < n; i++)
    {
      if(scene->characters[i]->hp <= 0)
      {
        toDelete.push_back(i);
        scene->score += 10;
        scene->sound.playEnemyDie(scene->characters[i]->GetPosition());
      }
    }
    //actually perform deletion
    if(toDelete.size() > 0)
    {
      sort(toDelete.begin(), toDelete.end(), std::greater<int>());
      for(unsigned int i = 0; i < toDelete.size(); ++i)
      {
        // BLOOD
        int x = (int)(scene->characters[toDelete[i]]->GetPosition().X() / TERRAIN_DENSITY);
        int y = (int)(scene->characters[toDelete[i]]->GetPosition().Y() / TERRAIN_DENSITY);
        scene->has_blood[x][y][0] = true;
        scene->has_blood[x][y][1] = true;
        if (x > 0 && y > 0 && y < TERRAIN_Y_SIZE - 1 && x < TERRAIN_X_SIZE-1)
        {
          scene->has_blood[x - 1][y][1] = true;
          scene->has_blood[x][y + 1][0] = true;
          scene->has_blood[x - 1][y + 1][0] = true;
          scene->has_blood[x - 1][y + 1][1] = true;
        }
        Character* temp = scene->characters[toDelete[i]];
        scene->characters[toDelete[i]] = scene->characters[scene->characters.size()-1];
        scene->characters[scene->characters.size()-1] = temp;
        scene->characters.pop_back();
      }
    }
  }

  // Remember previous time
  previous_time = current_time;
}

// Get the color of a triangle
R2Pixel GetTriangleColor(int x, int y, int z)
{
  R2Pixel pixel;
  TerrainType type = scene->terrain_types[x][y][z];

  switch (type)
  {
    case TERRAIN_OCEAN:        pixel = R2Pixel(.33f,.49f,.67f,1.0f);      break;
    case TERRAIN_LAKE:         pixel = R2Pixel(.33f,.60f,.80f,1.0f);      break;

    case TERRAIN_DRY_COA:      pixel = R2Pixel(.95f,.95f,.60f,1.0f);      break;
    case TERRAIN_DRY_LOW:      pixel = R2Pixel(.77f,.72f,.65f,1.0f);      break;
    case TERRAIN_DRY_MID:      pixel = R2Pixel(.57f,.9f,.47f,1.0f);       break;
    case TERRAIN_DRY_HIH:      pixel = R2Pixel(.57f,.65f,.50f,1.0f);      break;
    case TERRAIN_DRY_PEA:      pixel = R2Pixel(.30f,.30f,.30f,1.0f);      break;

    case TERRAIN_COM_COA:      pixel = R2Pixel(.80f,.84f,.70f,1.0f);      break;
    case TERRAIN_COM_LOW:      pixel = R2Pixel(.69f,.72f,.65f,1.0f);      break;
    case TERRAIN_COM_MID:      pixel = R2Pixel(.605f,.9f,.57f,1.0f);      break;
    case TERRAIN_COM_HIH:      pixel = R2Pixel(.67f,.82f,.65f,1.0f);      break;
    case TERRAIN_COM_PEA:      pixel = R2Pixel(.60f,.60f,.60f,1.0f);      break;

    case TERRAIN_WET_COA:      pixel = R2Pixel(.71f,.76f,.65f,1.0f);      break;
    case TERRAIN_WET_LOW:      pixel = R2Pixel(.61f,.73f,.66f,1.0f);      break;
    case TERRAIN_WET_MID:      pixel = R2Pixel(.64f,.9f,.66f,1.0f);       break;
    case TERRAIN_WET_HIH:      pixel = R2Pixel(.77f,.90f,.80f,1.0f);      break;
    case TERRAIN_WET_PEA:      pixel = R2Pixel(.90f,.90f,.90f,1.0f);      break;
    default:                   pixel = R2Pixel(.33f,.49f,.67f,1.0f);
  }
  return pixel;
}

// Set the drawing material to the material of a certain triangle
void SetTriangleColor(int x, int y, int z)
{
  TerrainType type = scene->terrain_types[x][y][z];
  if(scene->has_blood[x][y][z] == current_blood && type == currentTerrainType)
    return;
  currentTerrainType = type;
  current_blood = scene->has_blood[x][y][z];
  R2Pixel pixel;
  if (scene->has_blood[x][y][z])
  {
    pixel = R2Pixel(.6, 0, 0, 1);
  }
  else
    pixel = GetTriangleColor(x, y, z);
  
  GLfloat c[4];
  GLfloat d[4];
  GLfloat e[4];
  c[0] = pixel.Red(); c[1] = pixel.Green(); c[2] = pixel.Blue(); c[3] = pixel.Alpha();
  glColor3f(pixel.Red(),pixel.Green(),pixel.Blue());
  d[0]=c[0]/4; d[1]=c[1]/4; d[2]=c[2]/4; d[3]=c[3];
  e[0]=0; e[1]=0; e[2]=0; e[3]=0;

  // Load Color
  glMaterialfv(GL_FRONT, GL_AMBIENT, e);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, c);
  glMaterialfv(GL_FRONT, GL_SPECULAR, d);
  glMaterialfv(GL_FRONT, GL_EMISSION, e);
  glMaterialf(GL_FRONT, GL_SHININESS, 3);
}

// Returns true if a vector is in front of the cameras FOV
bool isInFront(R3Vector toSquare, double cosXFOV, double cosYFOV)
{
  if(toSquare.Dot(camera.towards) < 0)
    return false;

  R3Vector toSquareUp = toSquare - camera.right.Dot(toSquare) * camera.right;
  toSquareUp.Normalize();
  if(!(toSquareUp.Dot(camera.towards) > cosYFOV))
    return false;
  
  R3Vector toSquareRight = toSquare - camera.up.Dot(toSquare) * camera.up;
  toSquareRight.Normalize();
  if(!(toSquareRight.Dot(camera.towards) > cosXFOV))
    return false;

  return true;
}

// Draw all the terrain triangles
void DrawTerrain(R3Scene *scene)
{
  glEnable(GL_LIGHTING);
  currentTerrainType = TERRAIN_NONE;
  current_blood = false;

  // Get Camera Details
  R3Vector XFOV = (camera.towards) - tan(camera.xfov + .01) * camera.right;
  R3Vector YFOV = (camera.towards) - tan(camera.yfov + .01) * camera.up;
  XFOV.Normalize();
  YFOV.Normalize();
  double cosXFOV = camera.towards.Dot(XFOV);
  double cosYFOV = camera.towards.Dot(YFOV);
  double xpos = camera.eye.X() / TERRAIN_DENSITY;
  double ypos = camera.eye.Y() / TERRAIN_DENSITY;

  // Draw Terrain
  double draw_dist = draw_distance / TERRAIN_DENSITY;
  int xLimit = min((double)TERRAIN_X_SIZE, xpos + draw_dist);
  int yLimit = min((double)TERRAIN_Y_SIZE, ypos + draw_dist);
  for (int i = max(0.0,xpos-draw_dist); i < xLimit - 1; i++)
  {
    for (int j = max(0.0,ypos-draw_dist); j < yLimit - 1; j++)
    {
      // Rule out far triangles
      double yDiff = ypos - j;
      double xDiff = xpos - i;
      if(xDiff * xDiff + yDiff * yDiff > draw_dist * draw_dist)
        continue;

      // Get Points on square
      R3Point none, xone, yone, both;
      none = R3Point(i * TERRAIN_DENSITY, j * TERRAIN_DENSITY, scene->terrain_heights[i][j]);
      xone = R3Point((i+1) * TERRAIN_DENSITY, j * TERRAIN_DENSITY, scene->terrain_heights[i+1][j]);
      yone = R3Point(i * TERRAIN_DENSITY, (j+1) * TERRAIN_DENSITY, scene->terrain_heights[i][j+1]);
      both = R3Point((i+1) * TERRAIN_DENSITY, (j+1) * TERRAIN_DENSITY, scene->terrain_heights[i+1][j+1]);
      //R3Point avg = (none + xone + yone + both) / 4;

      // Check close
      bool close = abs(xpos - i) + abs(ypos - j) < 15;
      
      // Check if in view triangle
      if(!close)
      {
        if(!(isInFront(none - camera.eye, cosXFOV, cosYFOV) ||
             isInFront(xone - camera.eye, cosXFOV, cosYFOV) ||
             isInFront(yone - camera.eye, cosXFOV, cosYFOV) ||
             isInFront(both - camera.eye, cosXFOV, cosYFOV)))
          continue;
      }

      // One triangle of in a square
      R3Vector normal1 = scene->terrain_normals[i][j][1];
      if(normal1.Dot(camera.towards) < cosXFOV)
      {
        glBegin(GL_TRIANGLES);
        SetTriangleColor(i,j,1);
        glNormal3f(normal1.X(), normal1.Y(), normal1.Z());
        glVertex3f(both.X(), both.Y(), both.Z());
        glVertex3f(xone.X(), xone.Y(), xone.Z());
        glVertex3f(yone.X(), yone.Y(), yone.Z());
        glEnd();
      }

      // The other triangle in a square
      R3Vector normal2 = scene->terrain_normals[i][j][0];
      if(normal2.Dot(camera.towards) < cosXFOV)
      {
        glBegin(GL_TRIANGLES);
        SetTriangleColor(i,j,0);
        glNormal3f(normal2.X(), normal2.Y(), normal2.Z());
        glVertex3f(none.X(), none.Y(), none.Z());
        glVertex3f(xone.X(), xone.Y(), xone.Z());
        glVertex3f(yone.X(), yone.Y(), yone.Z());
        glEnd();
      }
    }
  }

  // Draw Ocean
  float oceanHeight = -0.05f;
  double xStart = camera.eye.X() - draw_distance;
  double xEnd = xStart + 2 * draw_distance;
  double yStart = camera.eye.Y() - draw_distance;
  double yEnd = yStart + 2 * draw_distance;
  R3Vector normal = scene->terrain_normals[0][0][0];
  SetTriangleColor(0,0,0);
  glBegin(GL_TRIANGLES);
  glNormal3f(normal.X(), normal.Y(), normal.Z());
  glVertex3f(xStart,yStart,oceanHeight);
  glVertex3f(xStart,yEnd,oceanHeight);
  glVertex3f(xEnd,yEnd,oceanHeight);
  glEnd();
  glBegin(GL_TRIANGLES);
  glNormal3f(normal.X(), normal.Y(), normal.Z());
  glVertex3f(xStart,yStart,oceanHeight);
  glVertex3f(xEnd,yEnd,oceanHeight);
  glVertex3f(xEnd,yStart,oceanHeight);
  glEnd();
}

// Draws whole scene
void DrawScene(R3Scene *scene)
{
  // Draw terrain
  DrawTerrain(scene);
    
  // Draw characters and projectiles in the scene
  DrawCharacters(scene);
}



//////////////////////////////////////////////////
// MOVEMENT
//////////////////////////////////////////////////

// Move your character
void moveMainCharacter(double delta_time)
{
  if(scene->mainChar->hp <= 0)
    return;

  R3Vector facing;
  double speed = PLAYER_SPEED * delta_time;
  double strafe_factor = .5;

  if (keyStates['w'] && !movedAlready['w'])
  {
    movedAlready['w'] = true;
    facing = camera.towards;
    facing.SetZ(0);
    facing.Normalize();
    facing *= speed;
    camera.eye[0] += facing.X();
    if (camera.eye[0] < 0)
      camera.eye[0] = .001;
    if (camera.eye[0] > (TERRAIN_X_SIZE-1) * TERRAIN_DENSITY)
      camera.eye[0] = (TERRAIN_X_SIZE-1) * TERRAIN_DENSITY - .001;
    camera.eye[1] += facing.Y();
    if (camera.eye[1] < 0)
      camera.eye[1] = .001;
    if (camera.eye[1] > (TERRAIN_Y_SIZE-1) * TERRAIN_DENSITY)
      camera.eye[1] = (TERRAIN_Y_SIZE-1) * TERRAIN_DENSITY - .001;
    
    if (scene->LocationValid(camera.eye[0], camera.eye[1]))
    {
      camera.eye[2] = scene->getRealZValue(camera.eye[0], camera.eye[1]) + scene->mainChar->model->sphere->Radius();
      scene->mainChar->MoveTo(camera.eye);
    }
    else
    {
      camera.eye[0] = scene->mainChar->GetPosition().X();
      camera.eye[1] = scene->mainChar->GetPosition().Y();
    }

  }  
  if (keyStates['a'] && !movedAlready['a'])
  {
    movedAlready['a'] = true;

    facing = camera.right;
    facing.SetZ(0);
    facing.Normalize();
    facing *= -speed * strafe_factor;
    camera.eye[0] += facing.X();
    if (camera.eye[0] < 0)
      camera.eye[0] = .001;
    if (camera.eye[0] > (TERRAIN_X_SIZE-1) * TERRAIN_DENSITY)
      camera.eye[0] = (TERRAIN_X_SIZE-1) * TERRAIN_DENSITY - .001;
    camera.eye[1] += facing.Y();
    if (camera.eye[1] < 0)
      camera.eye[1] = .001;
    if (camera.eye[1] > (TERRAIN_Y_SIZE-1) * TERRAIN_DENSITY)
      camera.eye[1] = (TERRAIN_Y_SIZE-1) * TERRAIN_DENSITY - .001;
    if (scene->LocationValid(camera.eye[0], camera.eye[1]))
    {
      camera.eye[2] = scene->getRealZValue(camera.eye[0], camera.eye[1]) + scene->mainChar->model->sphere->Radius();
      scene->mainChar->MoveTo(camera.eye);
    }
    else
    {
      camera.eye[0] = scene->mainChar->GetPosition().X();
      camera.eye[1] = scene->mainChar->GetPosition().Y();
    }
  }
  if (keyStates['s'] && !movedAlready['s'])
  {
    movedAlready['s'] = true;

    facing = camera.towards;
    facing.SetZ(0);
    facing.Normalize();
    facing *= -speed;
    camera.eye[0] += facing.X();
    if (camera.eye[0] < 0)
      camera.eye[0] = .001;
    if (camera.eye[0] > (TERRAIN_X_SIZE-1) * TERRAIN_DENSITY)
      camera.eye[0] = (TERRAIN_X_SIZE-1) * TERRAIN_DENSITY - .001;
    camera.eye[1] += facing.Y();
    if (camera.eye[1] < 0)
      camera.eye[1] = .001;
    if (camera.eye[1] > (TERRAIN_Y_SIZE-1) * TERRAIN_DENSITY)
      camera.eye[1] = (TERRAIN_Y_SIZE-1) * TERRAIN_DENSITY - .001;
    
    if (scene->LocationValid(camera.eye[0], camera.eye[1]))
    {
      camera.eye[2] = scene->getRealZValue(camera.eye[0], camera.eye[1]) + scene->mainChar->model->sphere->Radius();
      scene->mainChar->MoveTo(camera.eye);
    }
    else
    {
      camera.eye[0] = scene->mainChar->GetPosition().X();
      camera.eye[1] = scene->mainChar->GetPosition().Y();
    }
  }
  if (keyStates['d'] && !movedAlready['d'])
  {
    movedAlready['d'] = true;

    facing = camera.right;
    facing.SetZ(0);
    facing.Normalize();
    facing *= speed * strafe_factor;
    camera.eye[0] += facing.X();
    if (camera.eye[0] < 0)
      camera.eye[0] = .001;
    if (camera.eye[0] > (TERRAIN_X_SIZE-1) * TERRAIN_DENSITY)
      camera.eye[0] = (TERRAIN_X_SIZE-1) * TERRAIN_DENSITY - .001;
    camera.eye[1] += facing.Y();
    if (camera.eye[1] < 0)
      camera.eye[1] = .001;
    if (camera.eye[1] > (TERRAIN_Y_SIZE-1) * TERRAIN_DENSITY)
      camera.eye[1] = (TERRAIN_Y_SIZE-1) * TERRAIN_DENSITY - .001;
    
    if (scene->LocationValid(camera.eye[0], camera.eye[1]))
    {
      camera.eye[2] = scene->getRealZValue(camera.eye[0], camera.eye[1]) + scene->mainChar->model->sphere->Radius();
      scene->mainChar->MoveTo(camera.eye);
    }
    else
    {
      camera.eye[0] = scene->mainChar->GetPosition().X();
      camera.eye[1] = scene->mainChar->GetPosition().Y();
    }
  }
}


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
  double hp_bar_width = 500;
  double hp_bar_height = 25;
  R2Point botleft = R2Point(50,50);
  double hp_width = hp_bar_width * (scene->mainChar->hp/scene->mainChar->max_hp);

  // Draw Game Over if dead
  if(scene->mainChar->hp <= 0)
  {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glColor4f(1.0f, 0.0f, 0.0f, 0.5f);
    glBegin(GL_QUADS);
      glVertex2f(0.0f,0.0f);
      glVertex2f(0.0f, GLUTwindow_height);
      glVertex2f(GLUTwindow_width, GLUTwindow_height);
      glVertex2f(GLUTwindow_width, 0.0f);
    glEnd();

    glDisable(GL_BLEND);
    
    void * bigfont = GLUT_BITMAP_TIMES_ROMAN_24;
    void * midfont = GLUT_BITMAP_9_BY_15;
    
    glColor3d(1.0, 0.0, 0.0);
    string s = "GAME OVER";
    glRasterPos2i(GLUTwindow_width / 2 - 100,GLUTwindow_height / 2 - 40);
    for (string::iterator i = s.begin(); i != s.end(); ++i)
      glutBitmapCharacter(bigfont, *i);
    
    glColor3d(0.8, 0.8, 0.8);
    s = "Score: ";
    s += std::to_string((long double)(scene->score));
    glRasterPos2i(GLUTwindow_width / 2 - 70,GLUTwindow_height / 2 + 20);
    for (string::iterator i = s.begin(); i != s.end(); ++i)
    {
      char c = *i;
      if(c == '.')
        break;
      glutBitmapCharacter(midfont, c);
    }
    
    glColor3d(0.8, 0.8, 0.8);
    s = "Press R to try again";
    glRasterPos2i(GLUTwindow_width / 2 - 120,GLUTwindow_height / 2 + 40);
    for (string::iterator i = s.begin(); i != s.end(); ++i)
      glutBitmapCharacter(midfont, * i);

    glEnable(GL_DEPTH_TEST);
    return;
  }

  // Draw Health bar
  glBegin(GL_TRIANGLES);
      glColor3f(1.0f, 0.0f, 0.0f);
      glVertex2f(botleft.X(), botleft.Y());
      glVertex2f(botleft.X()+hp_width, botleft.Y());
      glVertex2f(botleft.X()+hp_width, botleft.Y()+hp_bar_height);
  glEnd();
  glBegin(GL_TRIANGLES);
      glColor3f(1.0f, 0.0f, 0.0f);
      glVertex2f(botleft.X(), botleft.Y());
      glVertex2f(botleft.X(), botleft.Y()+hp_bar_height);
      glVertex2f(botleft.X()+hp_width, botleft.Y()+hp_bar_height);
  glEnd();
  glBegin(GL_TRIANGLES);
      glColor3f(1.0f, 1.0f, 1.0f);
      glVertex2f(botleft.X()+hp_width, botleft.Y());
      glVertex2f(botleft.X()+hp_bar_width, botleft.Y());
      glVertex2f(botleft.X()+hp_bar_width, botleft.Y()+hp_bar_height);
  glEnd();
  glBegin(GL_TRIANGLES);
      glColor3f(1.0f, 1.0f, 1.0f);
      glVertex2f(botleft.X()+hp_width, botleft.Y());
      glVertex2f(botleft.X()+hp_width, botleft.Y()+hp_bar_height);
      glVertex2f(botleft.X()+hp_bar_width, botleft.Y()+hp_bar_height);
  glEnd();

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

  // Print number of enemies remaining
  glRasterPos2i(45,20);
  s = "Enemies: ";
  s += std::to_string((long double)(scene->characters.size()));
  for (string::iterator i = s.begin(); i != s.end(); ++i)
  {
    char c = *i;
    if(c == '.')
        break;
    glColor3d(1.0, 1.0, 1.0);
    glutBitmapCharacter(font, c);
  }

  // Print timer
  glRasterPos2i(525,20);
  s = "Time Elapsed: ";
  if(scene->characters.size() > 0 && scene->mainChar->hp > 0)
    last_time = GetTime();
  s += std::to_string((long double)(last_time - start_time));
  for (string::iterator i = s.begin(); i != s.end(); ++i)
  {
    char c = *i;
    if(c == '.')
      break;
    glColor3d(1.0, 1.0, 1.0);
    glutBitmapCharacter(font, c);
  }

  // Draw Crosshairs
  double screen_center_x = GLUTwindow_width / 2, 
  screen_center_y = GLUTwindow_height / 2;
  float crosshair_width = 3.5, crosshair_spacing = 10, crosshair_length = 30;

  glBegin(GL_TRIANGLES);
      glColor3f(0.1f, 0.1f, 0.1f);
      glVertex2f(screen_center_x + crosshair_spacing, 
        screen_center_y);
      glVertex2f(screen_center_x + crosshair_spacing + crosshair_length
        , screen_center_y + crosshair_width);
      glVertex2f(screen_center_x + crosshair_spacing + crosshair_length
        , screen_center_y - crosshair_width);

      glVertex2f(screen_center_x - crosshair_spacing, 
        screen_center_y);
      glVertex2f(screen_center_x - crosshair_spacing - crosshair_length
        , screen_center_y + crosshair_width);
      glVertex2f(screen_center_x - crosshair_spacing - crosshair_length
        , screen_center_y - crosshair_width);

      glVertex2f(screen_center_x, 
        screen_center_y + crosshair_spacing);
      glVertex2f(screen_center_x - crosshair_width
        , screen_center_y + crosshair_spacing + crosshair_length);
      glVertex2f(screen_center_x + crosshair_width
        , screen_center_y + crosshair_spacing + crosshair_length);

      glVertex2f(screen_center_x, 
        screen_center_y - crosshair_spacing);
      glVertex2f(screen_center_x - crosshair_width
        , screen_center_y - crosshair_spacing - crosshair_length);
      glVertex2f(screen_center_x + crosshair_width
        , screen_center_y - crosshair_spacing - crosshair_length);
  glEnd();

  // Draw Minimap
  int minimap_size = 100;
  int current_x = (int) (scene->mainChar->GetPosition().X() / TERRAIN_DENSITY);
  int current_y = (int) (scene->mainChar->GetPosition().Y() / TERRAIN_DENSITY);
  int num_pixels = 3;
  int x = 50, y = GLUTwindow_height - 50 ;
  int center_x = x + minimap_size * num_pixels / 2;
  int center_y = y - minimap_size * num_pixels / 2;
  for (int i = current_y-minimap_size / 2; i < current_y+minimap_size / 2; i++)
  {
    x = 50;
    y -= num_pixels; 
    for (int j = current_x-minimap_size / 2; j < current_x+minimap_size / 2; j++)
    {
      R2Pixel pix;
      if (i < 0 || i >= TERRAIN_Y_SIZE || j < 0 || j >= TERRAIN_X_SIZE)
      {
        pix = R2Pixel(0.33f,0.49f,.67f, 1.);
      }
      else
      {
        pix = GetTriangleColor(j, i, 0);
      }

      glBegin(GL_TRIANGLES);
      glColor3f(pix.Red(), pix.Green(), pix.Blue());
      glVertex2f(x, y);
      glVertex2f(x + num_pixels, y);
      glVertex2f(x, y + num_pixels);

      glVertex2f(x + num_pixels, y + num_pixels);
      glVertex2f(x + num_pixels, y);
      glVertex2f(x, y + num_pixels);
      glEnd();

      x += num_pixels;
    }
  }

  // Draw enemies on minimap
  int nc = scene->characters.size();
  for(int i = 0; i < nc; i++)
  {
    R3Point p = scene->characters[i]->GetPosition();
    int p_x = (int) (p.X() / TERRAIN_DENSITY);
    int p_y = (int) (p.Y() / TERRAIN_DENSITY);
    int d_x = num_pixels * (p_x - (current_x - minimap_size / 2)) + 50;
    int d_y = GLUTwindow_height - num_pixels * (p_y - (current_y - minimap_size / 2)) - 50;
    if(p_x >= current_x-minimap_size / 2 && p_x < current_x+minimap_size / 2 && p_y >= current_y-minimap_size / 2 && p_y < current_y+minimap_size / 2)
    {
      glBegin(GL_POLYGON);
      glColor3f(1,0,0);
		  for(double j = 0; j < 2 * 3.14159; j +=  3.14159 / 6) 
 			  glVertex2f(d_x + cos(j) * 3, d_y + sin(j) * 3);
		  glEnd();
    }
  }

  // Draw view triangle on minimap
  glEnable( GL_BLEND );
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
  glBegin(GL_TRIANGLES);
  glColor4f(1.0f, 1.0f, 1.0f, 0.8f);
  R3Vector towards = camera.towards;
  towards.SetZ(0);
  towards.Normalize();
  towards.Rotate(R3Vector(0, 0, 1), camera.xfov);
  glVertex2f(center_x, center_y);
  glVertex2f(center_x + 100 * towards.X(), center_y - 100 * towards.Y());
  towards.Rotate(R3Vector(0, 0, 1), -2 * camera.xfov);
  glVertex2f(center_x + 100 * towards.X(), center_y - 100 * towards.Y());
  glEnd();

  // Reset GL Settings
  glDisable( GL_BLEND );


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

void GLUTSaveImage(const char *filename)
{ 
  // Create image
  R2Image image(GLUTwindow_width, GLUTwindow_height);

  // Read screen into buffer
  GLfloat *pixels = new GLfloat [ 3 * GLUTwindow_width * GLUTwindow_height ];
  glReadPixels(0, 0, GLUTwindow_width, GLUTwindow_height, GL_RGB, GL_FLOAT, pixels);

  // Load pixels from frame buffer
  GLfloat *pixelsp = pixels;
  for (int j = 0; j < GLUTwindow_height; j++) {
    for (int i = 0; i < GLUTwindow_width; i++) {
      double r = (double) *(pixelsp++);
      double g = (double) *(pixelsp++);
      double b = (double) *(pixelsp++);
      R2Pixel pixel(r, g, b, 1);
      image.SetPixel(i, j, pixel);
    }
  }

  // Write image to file
  image.Write(filename);

  // Delete buffer
  delete [] pixels;
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

  // Generate New Enemy
  if (GetTime() - lastGen > CHARACTER_GEN_RATE)
  {
    lastGen = GetTime();
    if (scene->characters.size() < MAX_ENEMIES)
      scene->SpawnCharacter();
  }

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
  GLUTmouse[1] = GLUTwindow_height / 2;
  GLUTmouse[0] = GLUTwindow_width / 2;

  // Redraw
  glutPostRedisplay();
}

void GLUTRedraw(void)
{
  // Initialize OpenGL drawing modes
  glEnable(GL_LIGHTING);
  glDisable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ZERO);
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

  // Save image
  if (save_image) {
    char image_name[256];
    static int image_number = 1;
    for (;;) {
      sprintf(image_name, "image%d.jpg", image_number++);
      FILE *fp = fopen(image_name, "r");
      if (!fp) break; 
      else fclose(fp);
    }
    GLUTSaveImage(image_name);
    printf("Saved %s\n", image_name);
    save_image = 0;
  }

  // Save video
  if (save_video) {
    char frame_name[512];
    static int next_frame = 0;
    static int num_frames_recorded = 0;
    for (;;) {
      sprintf(frame_name, "%sframe%04d.jpg", video_prefix, next_frame++);
      FILE *fp = fopen(frame_name, "r");
      if (!fp) break; 
      else fclose(fp);
    }
    GLUTSaveImage(frame_name);
    if (next_frame % 100 == 1) {
      printf("Saved %s\n", frame_name);
    }
    if (num_frames_to_record == ++num_frames_recorded) {
      save_video = 0;
      printf("Recorded %d frames, stopping as instructed.\n", num_frames_recorded);
      quit = 1;
    }
  }

  // Quit here so that can save image before exit
  if (quit) {
    if (output_image_name) GLUTSaveImage(output_image_name);
    GLUTStop();
  }

  // Draw HUD
  DrawHUD(scene);

  // Swap buffers 
  glutSwapBuffers();
}

////////////////////////////////////////////////////////////
// GLUT MOUSE AND KEYBOARD
////////////////////////////////////////////////////////////

// Mouse movement
void GLUTMouseMove(int x, int y)
{
  if(scene->mainChar->hp <= 0)
    return;

  int reset_length = 200;
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  int dx, dy;
  if (mouse_reset && (fabs((double)(x - GLUTwindow_width / 2)) > reset_length /2 || 
    fabs((double)(y - GLUTwindow_height / 2)) > reset_length /2))
  {
    // Compute mouse movement
    dx = x - GLUTmouse[2];
    dy = y - GLUTmouse[3];
  }
  else
  {
    mouse_reset = false;
    // Compute mouse movement
    dx = x - GLUTmouse[0];
    dy = y - GLUTmouse[1];
  }

  
  // Process mouse motion event
  if ((dx != 0) || (dy != 0)) {
    double vx = (double) dx / (double) GLUTwindow_width;
    double vy = (double) dy / (double) GLUTwindow_height;
    double theta_lr = -8.0 * vx;
    double theta_ud = 8.0 * vy;
    R3Vector lr_rotation_axis(0, 0, 1);
    camera.towards.Rotate(lr_rotation_axis, theta_lr);
    camera.up.Rotate(lr_rotation_axis, theta_lr);
    camera.right = camera.towards % camera.up;
    camera.towards.Rotate(camera.right, theta_ud);
    camera.up.Rotate(camera.right, theta_ud);
    camera.towards.Normalize();
    camera.up.Normalize();
    camera.right.Normalize();

    if (camera.up.Z() < .2)
    {
      camera.towards.Rotate(camera.right, -theta_ud);
      camera.up.Rotate(camera.right, -theta_ud);
      camera.towards.Normalize();
      camera.up.Normalize();
      camera.right.Normalize();
    }
    glutPostRedisplay();
  }
  if (fabs((double)(x - GLUTwindow_width / 2)) > reset_length || fabs((double)(y - GLUTwindow_height / 2)) > reset_length)
  {
    mouse_reset = true;
    glutWarpPointer(GLUTwindow_width / 2, GLUTwindow_height / 2);
    // Remember mouse position 
    GLUTmouse[0] = GLUTwindow_width / 2;
    GLUTmouse[1] = GLUTwindow_height / 2;
    GLUTmouse[2] = x;
    GLUTmouse[3] = y;
  }
  else
  {
    // Remember mouse position 
    GLUTmouse[0] = x;
    GLUTmouse[1] = y;
  }

  // Change sound position
  scene->sound.setListeningPosition(camera.eye, camera.towards, camera.up);
}

// Shoot a projectile
static void shootProjectile(int)
{

  if (GLUTbutton[0] == 1 && scene->characters.size() > 0 && scene->mainChar->hp > 0)
  {
    R3Vector v = camera.towards;
    scene->mainChar->ShootProjectile(v);
    scene->sound.playPlayerShoot();
    glutPostRedisplay(); 
  }

  glutTimerFunc(SHOOT_RATE, shootProjectile, 1); 
}

// Get mouse button pressed
void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Remember button state 
  int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
  GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();
  
   // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Redraw
  glutPostRedisplay();
}

// Get special key presses
void GLUTSpecial(int key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 
  switch (key) {
  case GLUT_KEY_F1:
    save_image = 1;
    break;
  case GLUT_KEY_F2:
    save_video = save_video ^ 1;
    break;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();
}

// Get key pressed
void GLUTKeyboard(unsigned char key, int, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  if(scene->characters.size() <= 0 || scene->mainChar->hp <= 0)
  {
    if(key == 'r' || key == 'R')
    {
      delete scene->mainChar;
      for(int i = scene->characters.size()-1; !scene->characters.empty(); i--)
      {
        delete scene->characters[i];
        scene->characters.pop_back();
      }
      scene->mainChar = NULL;
      scene->initializeTerrain();
      scene->InitCharacters();
      InitializeGame();
      scene->score = 0;
      scene->camera.eye = camera.eye;

    }
    else if(key ==  0x1B)
      exit(0);
    return;
  }
  
  // Process keyboard button event 
  switch (key) {
    case 'w':
    case 'W':
        scene->sound.playWalking();
        keyStates['w'] = true;
        
      break;

    case 'a':
    case 'A':
        scene->sound.playWalking();
        keyStates['a'] = true;
      break;

    case 's':
    case 'S':
        scene->sound.playWalking();
        keyStates['s'] = true;
      break;

    case 'd':
    case 'D':
        scene->sound.playWalking();
        keyStates['d'] = true;
      break;

    case 0x1B:
      exit(0);

    case 'f':
    case 'F':
  printf("%i x %i\n", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
      break;

  }

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();
}

// Get key no longer pressed
void keyUp (unsigned char key, int, int) {  
  switch (key) {
    case 'w':
    case 'W':
      keyStates['w'] = false;
      if (!keyStates['a'] && !keyStates['s'] &&!keyStates['d'])
        scene->sound.stopWalking();
      break;

    case 'a':
    case 'A':
      keyStates['a'] = false;
      if (!keyStates['w'] && !keyStates['s'] &&!keyStates['d'])
        scene->sound.stopWalking();
      break;

    case 's':
    case 'S':
      keyStates['s'] = false;
      if (!keyStates['a'] && !keyStates['w'] &&!keyStates['d'])
        scene->sound.stopWalking();
      break;

    case 'd':
    case 'D':
      keyStates['d'] = false;
      if (!keyStates['a'] && !keyStates['s'] &&!keyStates['w'])
        scene->sound.stopWalking();
      break;
  }
}

// GLUT Initialization
void GLUTInit(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  GLUTwindow = glutCreateWindow("Island Paintball");

  // Initialize GLUT callback functions 
  glutIdleFunc(GLUTIdle);
  glutReshapeFunc(GLUTResize);
  glutDisplayFunc(GLUTRedraw);
  glutKeyboardFunc(GLUTKeyboard);
  glutKeyboardUpFunc(keyUp);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMouseMove);
  glutPassiveMotionFunc(GLUTMouseMove);

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
  // Move main character to random position
  for (;;)
  {
    camera.eye[0] = RandomNumber() * TERRAIN_X_SIZE * TERRAIN_DENSITY;
    camera.eye[1] = RandomNumber() * TERRAIN_Y_SIZE * TERRAIN_DENSITY;
    if (scene->LocationValid(camera.eye[0], camera.eye[1]))
      break;
  }
  camera.eye[2] = scene->getRealZValue(camera.eye[0], camera.eye[1]) + scene->mainChar->model->sphere->Radius();

  scene->mainChar->MoveTo(camera.eye);

  start_time = GetTime();
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
  glutWarpPointer(GLUTwindow_width / 2, GLUTwindow_height / 2);
  glutSetKeyRepeat(GLUT_KEY_REPEAT_OFF);

  // Initialize game settings
  InitializeGame();
  glutTimerFunc(SHOOT_RATE, shootProjectile, 1);

  // Run GLUT interface
  GLUTMainLoop();

  // Return success 
  return 0;
}






