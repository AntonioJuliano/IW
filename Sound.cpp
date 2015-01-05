#include <iostream>
#include <cstdlib>
#include <string>
#include "R3/R3.h"
#include "include/irrKlang.h"
#include "Sound.h"

using namespace irrklang;
using namespace std;

#pragma comment(lib, "irrKlang.lib")

// Contructor 
Sound::Sound()
{
  // Create engine
	soundEngine = createIrrKlangDevice();
	if(!soundEngine) exit(-1);
	soundEngine->setSoundVolume(1);
  soundEngine->setRolloffFactor((irrklang::ik_f32)0.4);
}

// Destructor
Sound::~Sound() {}

// Change the 3D listening position of the viewer
void Sound::setListeningPosition(R3Point position, R3Vector towards, R3Vector up)
{
  soundEngine->setListenerPosition(
    vec3df(position.X(),position.Y(),position.Z()),
    vec3df(towards.X(),towards.Y(),towards.Z()),
    vec3df(0,0,0),
    vec3df(up.X(),up.Y(),up.Z()));
}

// // Play the player shoot sound
// void Sound::playPlayerShoot()
// {
//   soundEngine->play2D("laser.wav", false, false, true);
// }