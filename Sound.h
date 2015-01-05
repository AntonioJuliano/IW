#include <iostream>
#include <string>
#include "include/irrKlang.h"
#include "R3/R3.h"

using namespace std;
using namespace irrklang;

class Sound
{
  private:
	  ISoundEngine* soundEngine;
	  ISound* walkingSound;

  public:
	  Sound();
	  ~Sound();
  
	  void setListeningPosition(R3Point position, R3Vector towards, R3Vector up);
	  void playSound();
	  void playWalking();
	  void stopWalking();
    void playPlayerShoot();
    void playPlayerDie();
    void playPlayerOw();
    void playEnemyShoot(R3Point position);
    void playEnemyDie(R3Point position);
    void playEnemyOw(R3Point position);
    void playGameOver();
		
};