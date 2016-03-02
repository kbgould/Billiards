#include <iostream>
using namespace std;

struct VECTOR {
	double x;
	double y;
	double z;
};

// Ball Class
class Ball {

private:

	bool cueBall;
	int number;

public:

	VECTOR velocity;
	VECTOR position;
	Ball();
	Ball(double, double, bool);	
	~Ball();
		
};

// Blank constructor.  Ball will start at 0,0,0 with velocity 0,0,0.
// Ball will be initialized as a Cue Ball.  Dont use this constructor.
Ball::Ball() {

	this->velocity.x = 0;
	this->velocity.y = 0;
	this->velocity.z = 0;
	
	this->position.x = 0;
	this->position.y = 0;
	this->position.z = 0;

	this->cueBall = true;
	number = 0;

}

// Constructor in which we specify an x and z position for the ball
Ball::Ball(double x, double z, bool cue) {

	this->velocity.x = 0;
	this->velocity.y = 0;
	this->velocity.z = 0;

	this->position.x = x;
	this->position.y = 0;
	this->position.z = 0;

	this->cueBall = cue;

}

Ball::~Ball()
{

	// uhhhh... yep.

}