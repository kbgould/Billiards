/*

		Billiards

by Keith Gould and Brian Martinson

*/

#include <direct.h>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <ctime>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <GL/GL.h>
#include <assert.h>
#include <IL/il.h>
#include <IL/ilut.h>
#include <irrKlang.h>
#include "sounds.h"
#include "Billiards.h"
#include "BilliardsTexturing.h"

using namespace irrklang;

#pragma comment(lib, "C:/Program Files/irrKlang-1.3.0/lib/Win32-visualStudio/irrKlang.lib")

#define DEGREES 57.2957795
#define PI 3.14159265
#define BALL_RADIUS 0.088
#define BALL_CIRCUM 0.55292
#define POCKET_RADIUS 0.176
#define KE_FRICTION 0.00008
#define X_SIZE 5.0
#define Z_SIZE 10.0
enum {LIVE, IN_POCKET, DEAD};										// Ball statuses.
enum {TABLE_VIEW, CUE_BALL_VIEW, FROZEN};							// View modes.
enum {WAITING, TAKING_SHOT, BALLS_MOVING, DRAWING_CUE, SCRATCH, PLAYER_1_WINS, PLAYER_2_WINS};	// Game States.
enum {PLAYER_ONE, PLAYER_TWO};										// Players.

// Ball Class
class Ball {

private:

	bool cueBall;
	int number;

public:

	bool JustHitWall;
	double Rotation;
	VECTOR RotationAxis;
	VECTOR velocity;
	VECTOR position;
	VECTOR lastPosition;
	int status;
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

	this->lastPosition.x = 0;
	this->lastPosition.y = 0;
	this->lastPosition.z = 0;

	this->RotationAxis.x = 0;
	this->RotationAxis.y = 0;
	this->RotationAxis.z = 0;

	this->Rotation = 0.0;

	this->JustHitWall = false;
	this->cueBall = true;
	this->number = 0;
	this->status = LIVE;

}

// Constructor in which we specify an x and z position for the ball
Ball::Ball(double x, double z, bool cue) {

	this->velocity.x = 0;
	this->velocity.y = 0;
	this->velocity.z = 0;

	this->position.x = x;
	this->position.y = 0;
	this->position.z = z;

	this->lastPosition.x = 0;
	this->lastPosition.y = 0;
	this->lastPosition.z = 0;

	this->RotationAxis.x = 0;
	this->RotationAxis.y = 0;
	this->RotationAxis.z = 0;

	this->Rotation = 0.0;

	this->JustHitWall = false;
	this->cueBall = cue;
	this->status = LIVE;

}

Ball::~Ball()
{

	// uhhhh... yep.

}

// Global Variables

ISoundEngine* engine;  // My sound engine
char cwd[250];
char textureloc[250];
char soundloc[250];
char* ShotString;
int GameState = SCRATCH;
bool LeftButtonDown = false;
bool RightButtonDown = false;
double RotateX = 0.0;
double RotateY = 0.0;
double TableViewRotateX = 0.0;
double TableViewRotateY = 0.0;
double ShotDirection;
double CueRotateY = 0.0;
double CueOffsetZ = 0.088;
double ShotStartZ, ShotMaxZ;
double ShotSpeedTotal, ShotSpeedFrames;
long ShotStartClock, ShotFinishClock;
bool NineBallPocketed = false;
bool Breaking = true;
int TargetBall = 1;
int FirstCollision, BallCollisions;
int RailsHit, BallsPocketed;
int CurrentPlayer = PLAYER_TWO;
int MouseX, MouseY, LastMouseX, LastMouseY;
int ViewMode = TABLE_VIEW;
int PocketList, CueTipList, CueShaft1List, CueShaft2List;
Ball* cueBall;
Ball* Balls[16];
int numBalls;
double FrozenX, FrozenZ;
double EyeX, EyeZ;
double xBumpMax, xBumpMin, zBumpMax, zBumpMin;

// .7071067811
// .0622253967
// const double POCKET_CENTER_X[6] = {-2.437775, 2.437775, -2.5, 2.5, -2.437775, 2.437775};
// const double POCKET_CENTER_Z[6] = {-4.937775, 4.937775, 0.0, 0.0, 4.937775, -4.937775};

const double POCKET_CENTER_X[6] = {-2.562225, 2.562225, -2.622225, 2.622225, -2.562225,  2.562225};
const double POCKET_CENTER_Z[6] = {-5.062225, 5.062225,  0.0, 0.0, 5.062225, -5.062225};

// Function Stubs

void MouseFunction(void);
void CueBallView(double*, double*);
void StartShot(double);
void InShot(void);
void FinishShot(double);
void CancelShot(void);
void DrawCueStick(void);
void DrawBalls(void);
void DrawTable(void);
void DrawPockets(void);
void LoadPockets(void);
static void buildCueShaft1(int);
static void buildCueShaft2(int);
static void buildCueTip(int);
static void recallList(int);
void MoveBalls(void);
void Friction(VECTOR*);
void LoadWalls(void);
void CheckBallCollisions(void);
bool CheckSingleCollision(int);
void CheckWallCollisions(void);
void CheckBallsStopped(void);
void CheckPockets(void);
void Place9Ball(void);
double getDirection(VECTOR*);
double distance(VECTOR*, VECTOR*);
void midpoint(VECTOR*, VECTOR*, VECTOR*);
double dotproduct(VECTOR*, VECTOR*);
void crossproduct(VECTOR*, VECTOR*, VECTOR*);
double vectorlength(VECTOR*);
void normalize(VECTOR*, VECTOR*);
void unitvector(VECTOR*, VECTOR*);
double Magnitude(VECTOR*, VECTOR*);
int DistancePointLine(VECTOR*, VECTOR*, VECTOR*, double*);
void NextShot(void);
void NewGame(void);
void display(void);
void ChangeView(int);
void reshapeFunc(void);
void timerFunction(int);
void MouseFunction(int, int, int, int);
void MotionFunction(int, int);

////////////////////////////////
// Utility Functions ///////////
////////////////////////////////

// Get the xz direction of a vector.
double getDirection(VECTOR* v) {

	double angle;
	double tan_a = 0.0;

	if (v->x == 0.0 && v->z == 0.0) return 0;

	if (v->x == 0) {
		if (v->z < 0) {
			angle = 1.5 * PI;
			return angle * DEGREES;
		} else { 
			angle = 0.5 * PI;
			return angle * DEGREES;
		}
	}

	angle = atan(v->z / v->x);
	if (v->x < 0) angle += PI;
	if (angle < 0) angle += (2*PI);
	angle = angle * DEGREES;

	return angle;

}

// A vector distance function
double distance(VECTOR* a, VECTOR* b) {

	double x = a->x - b->x;
	double z = a->z - b->z;

	double distance = sqrt((x*x) + (z*z));

	return distance;

}

// midpoint between two vectors/vertices
void midpoint(VECTOR* v, VECTOR* a, VECTOR* b) {

	v->x = (a->x + b->x) / 2;
	v->y = 0.0;
	v->z = (a->z + b->z) / 2;

}

// vector dot product
double dotproduct(VECTOR* a, VECTOR* b) {

	double dp = (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
	return dp;

}

// vector cross product
void crossproduct(VECTOR* a, VECTOR* b, VECTOR* result) {

	double a1 = a->x;
	double a2 = a->y;
	double a3 = a->z;

	double b1 = b->x;
	double b2 = b->y;
	double b3 = b->z;

	result->x = (a2*b3) - (a3*b2);
	result->y = (a3*b1) - (a1*b3);
	result->z = (a1*b2) - (a2*b1);

}

// returns the length of vector v
double vectorlength(VECTOR* v) {

	double l = sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
	return l;

}

// returns a normalized version of vector v in result
void normalize(VECTOR* v, VECTOR* result) {

	double length = vectorlength(v);

	result->x = v->x / length;
	result->y = v->y / length;
	result->z = v->z / length;

}

// Returns a unit vector of a in v.
void unitvector(VECTOR* v, VECTOR* a) {

	double length = vectorlength(a);
	if (length == 0) {
		v->x = 0;
		v->y = 0;
		v->z = 0;
		return;
	}
	v->x = a->x / length;
	v->y = a->y / length;
	v->z = a->z / length;

}

// projection of a vector onto another.
void projection(VECTOR* v, VECTOR* a, VECTOR* b) {

	//  c = ( (a dot b) / (b dot b) ) * b

	double a_dot_b = dotproduct(a, b);
	double b_dot_b = dotproduct(b, b);
	double top = a_dot_b / b_dot_b;

	v->x = top * b->x;
	v->y = top * b->y;
	v->z = top * b->z;

}

// Calculates line intersection, but this function isnt used as it is built in to line_intersection.
void LineIntersection (VECTOR* L1Start, VECTOR* L1Fin, VECTOR* L2Start, VECTOR* L2Fin, VECTOR* Intersection) {

double x1 = L1Start->x;
double y1 = L1Start->z;
double x2 = L1Fin->x;
double y2 = L1Fin->z;
double x3 = L2Start->x;
double y3 = L2Start->z;
double x4 = L2Fin->x;
double y4 = L2Fin->z;

Intersection->x =	( (((x1*y2) - (y1*x2)) * (x3-x4) ) - ((x1-x2) * ((x3*y4) - (y3*x4))) )
								/
						(((x1 - x2) * (y3 - y4)) - ((y1-y2)*(x3-x4)));

Intersection->z =	( (((x1*y2) - (y1*x2)) * (y3-y4) ) - ((y1-y2) * ((x3*y4) - (y3*x4))) )
								/
						(((x1 - x2) * (y3 - y4)) - ((y1-y2)*(x3-x4)));


}

/* Faster Line Segment Intersection   */
/* Franklin Antonio                   */

/* return values */
#define DONT_INTERSECT 0
#define DO_INTERSECT   1
#define COLLINEAR      2

/* The use of some short working variables allows this code to run   */
/* faster on 16-bit computers, but is not essential.  It should not  */
/* affect operation on 32-bit computers.  The short working variables*/
/* to not restrict the range of valid input values, as these were    */
/* constrained in any case, due to algorithm restrictions.           */

int lines_intersect(VECTOR* L1Start, VECTOR* L1Fin, VECTOR* L2Start, VECTOR* L2Fin, VECTOR* Intersection) {

double x1 = L1Start->x;
double y1 = L1Start->z;
double x2 = L1Fin->x;
double y2 = L1Fin->z;
double x3 = L2Start->x;
double y3 = L2Start->z;
double x4 = L2Fin->x;
double y4 = L2Fin->z;

double Ax,Bx,Cx,Ay,By,Cy,d,e,f;
double x1lo,x1hi,y1lo,y1hi;

Ax = x2-x1;
Bx = x3-x4;

if(Ax<0) {						/* X bound box test*/
	x1lo=x2; x1hi=x1;
} else {
	x1hi=x2; x1lo=x1;
}

if(Bx>0) {
	if(x1hi < x4 || x3 < x1lo) return DONT_INTERSECT;
} else {
	if(x1hi < x3 || x4 < x1lo) return DONT_INTERSECT;
}

Ay = y2-y1;
By = y3-y4;

if(Ay<0) {						/* Y bound box test*/
	y1lo=y2; y1hi=y1;
} else {
	y1hi=y2; y1lo=y1;
}

if(By>0) {
	  if(y1hi < y4 || y3 < y1lo) return DONT_INTERSECT;
} else {
	  if(y1hi < y3 || y4 < y1lo) return DONT_INTERSECT;
}





Cx = x1-x3;
Cy = y1-y3;

d = By*Cx - Bx*Cy;					/* alpha numerator*/
f = Ay*Bx - Ax*By;					/* both denominator*/

if(f>0) {						/* alpha tests*/
	if(d<0 || d>f) return DONT_INTERSECT;
} else {
	if(d>0 || d<f) return DONT_INTERSECT;
}

e = Ax*Cy - Ay*Cx;					/* beta numerator*/

if(f>0) {						/* beta tests*/
	if(e<0 || e>f) return DONT_INTERSECT;
} else {
	if(e>0 || e<f) return DONT_INTERSECT;
}



/*compute intersection coordinates*/
if(f==0) return COLLINEAR;

Intersection->x =	( (((x1*y2) - (y1*x2)) * (x3-x4) ) - ((x1-x2) * ((x3*y4) - (y3*x4))) )
								/
						(((x1 - x2) * (y3 - y4)) - ((y1-y2)*(x3-x4)));

Intersection->z =	( (((x1*y2) - (y1*x2)) * (y3-y4) ) - ((y1-y2) * ((x3*y4) - (y3*x4))) )
								/
						(((x1 - x2) * (y3 - y4)) - ((y1-y2)*(x3-x4)));

return DO_INTERSECT;

}

// This function is called to initialize a new shot.
// THis is not to be confused with Drawing back the cue
// for a new shot which is StartShot() and is below this.
void InitShot() {
	CueRotateY = 0.0;
	RotateX = 0.0;
	ShotDirection = 180.0 - getDirection(&Balls[0]->position);

}

// This function is called when you start to draw back your shot.
// (aka when you first click the left mouse button.
void StartShot(double StartZ) {

	GameState = DRAWING_CUE;
	ShotStartClock = (int)clock();
	ShotStartZ = StartZ;

}

// This function is called again and again during a shot.
void InShot(double ThisZ, double LastZ) {

	if (ThisZ > ShotMaxZ) {
		ShotMaxZ = ThisZ;
	}

	if (ThisZ > LastZ) {

		ShotSpeedFrames = ShotSpeedFrames + 1.0;
			
	} else {

		ShotSpeedFrames = 0;

	}
		

}

// This function is called when you follow through with your shot.
void FinishShot(double ThisZ) {

	ShotFinishClock = clock();

	double dir = ShotDirection - RotateX + CueRotateY + 90.0;
	//double direction = ShotDirection + RotateX - CueRotateY + 90.0;

	//VECTOR* facing = new VECTOR;
	//facing->x = EyeX - Balls[0]->position.x;
	//facing->y = 0.0;
	//facing->z = EyeZ - Balls[0]->position.z;
	//double XDirOffset = getDirection(facing);

	//double direction = CueRotateY - XDirOffset - 90.0;
	double ShotSize = abs(ShotMaxZ - ShotStartZ);
	double ShotTime = (double)ShotFinishClock - (double)ShotStartClock;
	if (ShotTime > 2500) ShotTime = 2500.0;
	double PowerFraction = (2500.0 - ShotTime) / 2500;
	double ShotPower = (ShotSize * 0.0005) * PowerFraction;

	double ShotSpeedX = (ShotPower * sin(dir / DEGREES));
	double ShotSpeedZ = (ShotPower * cos(dir / DEGREES));

	Balls[0]->velocity.x = Balls[0]->velocity.x + ShotSpeedX;
	Balls[0]->velocity.z = Balls[0]->velocity.z + ShotSpeedZ;

	vec3df pos(0.0, 0.0, 0.0);
	engine->play3D("cuehit.ogg",pos);
	
	ViewMode = FROZEN;
	GameState = BALLS_MOVING;
	CueRotateY = 0.0;
	CueOffsetZ = 0.0;

}

// This function cancels the current shot.
void CancelShot(void) {

	CueRotateY = 0.0;
	CueOffsetZ = 0.0;
	GameState = TAKING_SHOT;

}

// This function builds and vertex arrays a cylinder for the tip of the cue stick.
static void buildCueTip(int numc) {

	GLfloat j, x, y, z, twopi;

	twopi = (GLfloat) (2*PI);

	glNewList(CueTipList, GL_COMPILE);

	for (int i = 0; i <= 10; i++) {

		int zpos = i-20;

		glBegin(GL_QUAD_STRIP);

		for (int k = 0; k <= 100; k++) {

			j = k * (twopi/numc);

			x = (0.01f + ((140-i)*0.0001f)) * cos(j);
			y = (0.01f + ((140-i)*0.0001f)) * sin(j);
			z = 0.01f * zpos;

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

			x = (0.01f + ((140-i)*0.0001f)) * cos(j+(twopi/numc));
			y = (0.01f + ((140-i)*0.0001f)) * sin(j+(twopi/numc));
			z = 0.01f * zpos;

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

			x = (0.01f + ((140-i)*0.0001f)) * cos(j);
			y = (0.01f + ((140-i)*0.0001f)) * sin(j);
			z = 0.01f * (zpos+1);

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

			x = (0.01f + ((140-i)*0.0001f)) * cos(j+(twopi/numc));
			y = (0.01f + ((140-i)*0.0001f)) * sin(j+(twopi/numc));
			z = 0.01f * (zpos+1);

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

		}

		glEnd();

	}
	glEndList();

}

// This function builds and vertex arrays the first part of
// the shaft of the cue stick using a cylinder.
static void buildCueShaft1(int numc) {

	GLfloat j, x, y, z, twopi;

	twopi = (GLfloat) (2*PI);

	glNewList(CueShaft1List, GL_COMPILE);

	for (int i = 0; i <= 100; i++) {

		int zpos = i-120;

		glBegin(GL_QUAD_STRIP);

		for (int k = 0; k <= 100; k++) {

			j = k * (twopi/numc);

			x = (0.02f + ((100-i)*0.0001f)) * cos(j);
			y = (0.02f + ((100-i)*0.0001f)) * sin(j);
			z = 0.01f * zpos;

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

			x = (0.02f + ((100-i)*0.0001f)) * cos(j+(twopi/numc));
			y = (0.02f + ((100-i)*0.0001f)) * sin(j+(twopi/numc));
			z = 0.01f * zpos;

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

			x = (0.02f + ((100-i)*0.0001f)) * cos(j);
			y = (0.02f + ((100-i)*0.0001f)) * sin(j);
			z = 0.01f * (zpos+1);

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

			x = (0.02f + ((100-i)*0.0001f)) * cos(j+(twopi/numc));
			y = (0.02f + ((100-i)*0.0001f)) * sin(j+(twopi/numc));
			z = 0.01f * (zpos+1);

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

		}

		glEnd();

	}
	glEndList();

}

// This function builds and vertex arrays the second part of
// the shaft of the cue stick using a cylinder.
static void buildCueShaft2(int numc) {

	GLfloat j, x, y, z, twopi;

	twopi = (GLfloat) (2*PI);

	glNewList(CueShaft2List, GL_COMPILE);

	for (int i = 0; i <= 100; i++) {

		int zpos = i-220;

		glBegin(GL_QUAD_STRIP);

		for (int k = 0; k <= 100; k++) {

			j = k * (twopi/numc);

			x = (0.02f + ((i)*0.0001f)) * cos(j);
			y = (0.02f + ((i)*0.0001f)) * sin(j);
			z = 0.01f * zpos;

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

			x = (0.02f + ((i)*0.0001f)) * cos(j+(twopi/numc));
			y = (0.02f + ((i)*0.0001f)) * sin(j+(twopi/numc));
			z = 0.01f * zpos;

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

			x = (0.02f + ((i)*0.0001f)) * cos(j);
			y = (0.02f + ((i)*0.0001f)) * sin(j);
			z = 0.01f * (zpos+1);

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

			x = (0.02f + ((i)*0.0001f)) * cos(j+(twopi/numc));
			y = (0.02f + ((i)*0.0001f)) * sin(j+(twopi/numc));
			z = 0.01f * (zpos+1);

			glNormal3f(-x, -y, 0);
			glVertex3f(x, y, z);

		}

		glEnd();

	}
	glEndList();

}

// Recalls a given display list.  This is actually totally pointless.
static void recallList(int type)
{

	glCallList(type);

}

// This function will draw the cue stick... It does some calculation
// on where the cue stick should be and uses some other functions that
// will draw pieces of the cue stick.
void DrawCueStick(void) {

	double direction = getDirection(&Balls[0]->position);
	direction = direction + RotateX - CueRotateY;

	double TmpX = (1.5 * cos(direction / DEGREES));
	double TmpZ = (1.5 * sin(direction / DEGREES));
	double TmpOffsetX = (CueOffsetZ * cos(direction / DEGREES));
	double TmpOffsetZ = (CueOffsetZ * sin(direction / DEGREES));

	glPushMatrix();

	glTranslated(Balls[0]->position.x + TmpOffsetX, 0.0, Balls[0]->position.z + TmpOffsetZ);
	// Rotate cue stick for ball / camera position
	glRotated(-direction-90.0, 0.0, 1.0, 0.0);
	// Rotate cue stick for player mouse input

	// White

	GLfloat	specular3[4] = { 2.55f,  2.55f, 2.55f, 1.0f };
	GLfloat	ambient3[4] = { 2.55f,  2.55f, 2.55f, 1.0f };
	GLfloat	diffuse3[4] = { 2.55f, 2.55f, 2.55f, 1.0f };

	// Sandy Brown: 244-164-96

	GLfloat	specular1[4]	= { 2.44f,  1.64f, 0.96f, 1.0f };
	GLfloat	ambient1[4] = { 2.44f,  1.64f, 0.96f, 1.0f };
	GLfloat	diffuse1[4] = { 0.1f, 0.1f, 0.1f, 1.0f };

	// Peru: 205-133-63

	GLfloat	specular2[4]	= { 2.05f,  0.69f, 0.19f, 1.0f };
	GLfloat	ambient2[4] = { 2.05f,  0.69f, 0.19f, 1.0f };
	GLfloat	diffuse2[4] = { 0.1f, 0.1f, 0.1f, 1.0f };

	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10.0f);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient3);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse3);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular3);

	recallList(CueTipList);

	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10.0f);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular1);

	recallList(CueShaft1List);

	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10.0f);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular2);

	recallList(CueShaft2List);
	
	glPopMatrix();

}

// This function initializes all the positions of my billiard balls.
// Read the big comment below for more info.
void InitBalls() {

	numBalls = 10;
	
	/*
	 
	  Many chunks of code get repeated for each ball, so rather than explaining the purpose of each chunk. We will
	  discuss how each ball is intitialized here. (Saves whitespace :P )

	  First we must create our ball in the Ball[] array and give it a location and status.
	  Then we must do some DevIL inits and load the image stated and bind it to a handle GLuint.
	  If the first load is successful then we must get some image information like WxH.
	  We then check our pixel map in memory of the image, if it contains data, clear it.
	  Then make a new pixel map and map each 4 bytes (RGBa) into the map.
	  Copy the pixels loaded into the map.
	  Then bind the texture to our handle for the given ball texture.
	  We then call the create sphere function and create a vector list for the sphere.
	  Following this we apply our pixel mapped image to that sphere. And BOOM, we have our textured ball.

	  Each time we draw it we just recall the vector array and rebind the texture we want for that recall.

	 */
	
	Balls[0] = new Ball(0.0, -3.0, true);
	ILuint handle;
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/cue.bmp");
	ILboolean loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		TexCue = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TDCue.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TDCue)
			delete [] TDCue;
	
		TDCue = new cRGB_Byte_Pixel[w * h * 4];
		if (!TDCue)
		{
			TexCue = -1;
		}

		//------ Copy TDCue from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TDCue);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		TexCue = handle;
	}
	ListCue = glGenLists(1);
	CreateSphere(1.0, stacks, slices, ListCue);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &TexCue);
	glBindTexture(GL_TEXTURE_2D, TexCue);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TDCue);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	
	///////////////
	// Row 1 //////
	///////////////

	Balls[1] = new Ball(0.0, 3.0, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/1.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex1 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD1.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD1)
			delete [] TD1;
	
		TD1 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD1)
		{
			Tex1 = -1;
		}

		//------ Copy TD1 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD1);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex1 = handle;
	}
	List1 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List1);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex1);
	glBindTexture(GL_TEXTURE_2D, Tex1);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD1);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	///////////////
	// Row 2 //////
	///////////////

	Balls[2] = new Ball(-0.1, 3.2, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/2.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex2 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD2.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD2)
			delete [] TD2;
	
		TD2 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD2)
		{
			Tex2 = -1;
		}

		//------ Copy TD2 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD2);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex2 = handle;
	}
	List2 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List2);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex2);
	glBindTexture(GL_TEXTURE_2D, Tex2);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD2);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	Balls[3] = new Ball(0.1, 3.2, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/3.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex3 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD3.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD3)
			delete [] TD3;
	
		TD3 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD3)
		{
			Tex3 = -1;
		}

		//------ Copy TD3 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD3);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex3 = handle;
	}
	List3 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List3);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex3);
	glBindTexture(GL_TEXTURE_2D, Tex3);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD3);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	///////////////
	// Row 3 //////
	///////////////

	Balls[4] = new Ball(-0.2, 3.4, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/4.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex4 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD4.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD4)
			delete [] TD4;
	
		TD4 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD4)
		{
			Tex4 = -1;
		}

		//------ Copy TD4 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD4);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex4 = handle;
	}
	List4 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List4);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex4);
	glBindTexture(GL_TEXTURE_2D, Tex4);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD4);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	Balls[5] = new Ball(0.2, 3.4, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/5.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex5 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD5.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD5)
			delete [] TD5;
	
		TD5 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD5)
		{
			Tex5 = -1;
		}

		//------ Copy TD5 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD5);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex5 = handle;
	}
	List5 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List5);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex5);
	glBindTexture(GL_TEXTURE_2D, Tex5);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD5);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	Balls[6] = new Ball(-0.1, 3.6, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/6.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex6 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD6.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD6)
			delete [] TD6;
	
		TD6 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD6)
		{
			Tex6 = -1;
		}

		//------ Copy TD6 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD6);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex6 = handle;
	}
	List6 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List6);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex6);
	glBindTexture(GL_TEXTURE_2D, Tex6);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD6);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	// Row 4
	Balls[7] = new Ball(0.1, 3.6, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/7.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex7 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD7.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD7)
			delete [] TD7;
	
		TD7 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD7)
		{
			Tex7 = -1;
		}

		//------ Copy TD7 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD7);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex7 = handle;
	}
	List7 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List7);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex7);
	glBindTexture(GL_TEXTURE_2D, Tex7);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD7);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	Balls[8] = new Ball(0.0, 3.8, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/8.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex8 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD8.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD8)
			delete [] TD8;
	
		TD8 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD8)
		{
			Tex8 = -1;
		}

		//------ Copy TD8 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD8);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex8 = handle;
	}
	List8 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List8);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex8);
	glBindTexture(GL_TEXTURE_2D, Tex8);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD8);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	Balls[9] = new Ball(0.0, 3.4, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/9.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex9 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD9.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD9)
			delete [] TD9;
	
		TD9 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD9)
		{
			Tex9 = -1;
		}

		//------ Copy TD9 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD9);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex9 = handle;
	}
	List9 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List9);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex9);
	glBindTexture(GL_TEXTURE_2D, Tex9);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD9);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	Balls[10] = new Ball(0.3, 3.6, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/10.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex10 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD10.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD10)
			delete [] TD10;
	
		TD10 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD10)
		{
			Tex10 = -1;
		}

		//------ Copy TD10 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD10);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex10 = handle;
	}
	List10 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List10);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex10);
	glBindTexture(GL_TEXTURE_2D, Tex10);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD10);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	// Row 5
	Balls[11] = new Ball(-0.4, 3.8, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/11.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex11 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD11.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD11)
			delete [] TD11;
	
		TD11 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD11)
		{
			Tex11 = -1;
		}

		//------ Copy TD11 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD11);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex11 = handle;
	}
	List11 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List11);
	//ApplyTexture(Tex11);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex11);
	glBindTexture(GL_TEXTURE_2D, Tex11);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD11);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	
	Balls[12] = new Ball(-0.2, 3.8, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/12.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex12 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD12.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD12)
			delete [] TD12;
	
		TD12 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD12)
		{
			Tex12 = -1;
		}

		//------ Copy TD12 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD12);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex12 = handle;
	}
	List12 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List12);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex12);
	glBindTexture(GL_TEXTURE_2D, Tex12);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD12);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	Balls[13] = new Ball(0.0, 3.8, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/13.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex13 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD13.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD13)
			delete [] TD13;
	
		TD13 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD13)
		{
			Tex13 = -1;
		}

		//------ Copy TD13 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD13);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex13 = handle;
	}
	List13 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List13);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex13);
	glBindTexture(GL_TEXTURE_2D, Tex13);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD13);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	Balls[14] = new Ball(0.2, 3.8, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/14.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex14 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD14.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD14)
			delete [] TD14;
	
		TD14 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD14)
		{
			Tex14 = -1;
		}

		//------ Copy TD14 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD14);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex14 = handle;
	}
	List14 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List14);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex14);
	glBindTexture(GL_TEXTURE_2D, Tex14);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD14);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	Balls[15] = new Ball(0.4, 3.8, true);
	ilGenImages(1, &handle);
	ilBindImage(handle);
	strcpy_s(textureloc, cwd);
	strcat_s(textureloc, "/textures/15.bmp");
	loaded = ilLoadImage(textureloc);
	if (loaded == IL_FALSE)
	{
		Tex15 = -1;
	}
	else
	{
		// THIS CODE WILL LOAD THE IMAGE'S DATA INFO INTO TD15.
		ILuint w = ilGetInteger(IL_IMAGE_WIDTH);
		ILuint h = ilGetInteger(IL_IMAGE_HEIGHT);
		TexWidth = w;
		TexHeight = h;

		if (TD15)
			delete [] TD15;
	
		TD15 = new cRGB_Byte_Pixel[w * h * 4];
		if (!TD15)
		{
			Tex15 = -1;
		}

		//------ Copy TD15 from DevIL to our pixmap.
		ilCopyPixels(0, 0, 0, w, h, 1, IL_RGBA, IL_UNSIGNED_BYTE, TD15);

		//------ Unbound image name and frees DevIL image memory.
		ilBindImage(0);
		ilDeleteImage(handle);

		Tex15 = handle;
	}
	List15 = glGenLists(1);
	CreateSphere(1.0, stacks, slices, List15);
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Tex15);
	glBindTexture(GL_TEXTURE_2D, Tex15);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, TexWidth, TexHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, TD15);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	cueBall = Balls[0];

}

// DrawBalls will draw and bind texture to each of the billiards
// balls.  It still has some remnant pieces in it from when the
// game had 16 balls.
void DrawBalls() {

	GLfloat	specular[4]	= { 0.7f,  0.7f, 0.7f, 1.0f };
	GLfloat	ambient[4] = { 0.6f, 0.6f, 0.6f, 1.0f };
	GLfloat	diffuse[4] = { 1.0f, 1.0f, 1.0f, 1.0f };

	GLfloat	specular2[4] = { 0.3f,  0.3f, 1.0f, 1.0f };
	GLfloat	ambient2[4] = { 0.3f,  0.3f, 1.0f, 1.0f };
	GLfloat	diffuse2[4] = { 0.5f, 0.5f, 1.0f, 1.0f };

	double scale = 0.088f;

	for (int i = 0; i < numBalls; i++) {

		if (Balls[i]->status == LIVE || i == 9) {

			glPushMatrix();

			if (i == 0) {
				glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10.0f);
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
				glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);

				glTranslated(Balls[i]->position.x, 0.0, Balls[i]->position.z);
				glPushMatrix();
				if (vectorlength(&Balls[i]->velocity) > 0.0) {
					crossproduct(UP_AXIS, &Balls[i]->velocity, &Balls[i]->RotationAxis);
					normalize(&Balls[i]->RotationAxis, &Balls[i]->RotationAxis);
					double speed = vectorlength(&Balls[i]->velocity);
					Balls[i]->Rotation -= (speed / BALL_CIRCUM) * 360.0;
				}
				glRotated(Balls[i]->Rotation, Balls[i]->RotationAxis.x, Balls[i]->RotationAxis.y, Balls[i]->RotationAxis.z);
				glScaled(scale, scale, scale);
				glBindTexture(GL_TEXTURE_2D, TexCue);
				DrawSphere(ListCue);
				glPopMatrix();
			} else {
				glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10.0f);
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient2);
				glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse2);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular2);

				glTranslated(Balls[i]->position.x, 0.0, Balls[i]->position.z);

				glPushMatrix();
				if (vectorlength(&Balls[i]->velocity) > 0.0) {
					crossproduct(UP_AXIS, &Balls[i]->velocity, &Balls[i]->RotationAxis);
					normalize(&Balls[i]->RotationAxis, &Balls[i]->RotationAxis);
					double speed = vectorlength(&Balls[i]->velocity);
					Balls[i]->Rotation -= (speed / BALL_CIRCUM) * 360.0;
				}
				glRotated(Balls[i]->Rotation, Balls[i]->RotationAxis.x, Balls[i]->RotationAxis.y, Balls[i]->RotationAxis.z);
				/*
				 * Here we draw the appropriate ball with the appropriate texture.
				 */
				switch(i){
					case 1:
						glBindTexture(GL_TEXTURE_2D, Tex1);
						glScaled(scale, scale, scale);
						DrawSphere(List1);
						break;
					case 2:
						glBindTexture(GL_TEXTURE_2D, Tex2);
						glScaled(scale, scale, scale);
						DrawSphere(List2);
						break;
					case 3:
						glBindTexture(GL_TEXTURE_2D, Tex3);
						glScaled(scale, scale, scale);
						DrawSphere(List3);
						break;
					case 4:
						glBindTexture(GL_TEXTURE_2D, Tex4);
						glScaled(scale, scale, scale);
						DrawSphere(List4);
						break;
					case 5:
						glBindTexture(GL_TEXTURE_2D, Tex5);
						glScaled(scale, scale, scale);
						DrawSphere(List5);
						break;
					case 6:
						glBindTexture(GL_TEXTURE_2D, Tex6);
						glScaled(scale, scale, scale);
						DrawSphere(List6);
						break;
					case 7:
						glBindTexture(GL_TEXTURE_2D, Tex7);
						glScaled(scale, scale, scale);
						DrawSphere(List7);
						break;
					case 8:
						glBindTexture(GL_TEXTURE_2D, Tex8);
						glScaled(scale, scale, scale);
						DrawSphere(List8);
						break;
					case 9:
						glBindTexture(GL_TEXTURE_2D, Tex9);
						glScaled(scale, scale, scale);
						DrawSphere(List9);
						break;
					case 10:
						glBindTexture(GL_TEXTURE_2D, Tex10);
						glScaled(scale, scale, scale);
						DrawSphere(List10);
						break;
					case 11:
						glBindTexture(GL_TEXTURE_2D, Tex11);
						glScaled(scale, scale, scale);
						DrawSphere(List11);
						break;
					case 12:
						glBindTexture(GL_TEXTURE_2D, Tex12);
						glScaled(scale, scale, scale);
						DrawSphere(List12);
						break;
					case 13:
						glBindTexture(GL_TEXTURE_2D, Tex13);
						glScaled(scale, scale, scale);
						DrawSphere(List13);
						break;
					case 14:
						glBindTexture(GL_TEXTURE_2D, Tex14);
						glScaled(scale, scale, scale);
						DrawSphere(List14);
						break;
					case 15:
						glBindTexture(GL_TEXTURE_2D, Tex15);
						glScaled(scale, scale, scale);
						DrawSphere(List15);
						break;
					default:
						break;
				}

				glPopMatrix();

			}

			glPopMatrix();

		}

	}

}

// DrawTable will draw all the parts of our table.  It's pretty
// straight forward.  Just quads and all.
void DrawTable() {

	//GLfloat	specular[4]	= { 0.34f,  1.39f, 0.34f, 1.0f };
	//GLfloat	ambient[4] = { 0.34f,  1.39f, 0.34f, 1.0f };
	//GLfloat	diffuse[4] = { 0.4f, 1.0f, 0.4f, 1.0f };
	
	GLfloat	specular[4]	= { 0.34f,  0.0f, 1.34f, 1.0f };
	GLfloat	ambient[4] = { 0.34f,  0.0f, 1.34f, 1.0f };
	GLfloat	diffuse[4] = { 0.04f, 0.1f, 0.04f, 1.0f };

	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10.0f);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);

	glPushMatrix();

	glTranslated(0.0, 0.088, 0.0);

	glBegin(GL_QUADS);

	glNormal3d(0.0, 1.0, 0.0);
	glVertex3d(X_SIZE/2 + 5*BALL_RADIUS, 0.0, Z_SIZE/2 + 5*BALL_RADIUS);
	glNormal3d(0.0, 1.0, 0.0);
	glVertex3d(X_SIZE/2 + 5*BALL_RADIUS, 0.0, -Z_SIZE/2 - 5*BALL_RADIUS);
	glNormal3d(0.0, 1.0, 0.0);
	glVertex3d(-X_SIZE/2 - 5*BALL_RADIUS, 0.0, -Z_SIZE/2 - 5*BALL_RADIUS);
	glNormal3d(0.0, 1.0, 0.0);
	glVertex3d(-X_SIZE/2 - 5*BALL_RADIUS, 0.0, Z_SIZE/2 + 5*BALL_RADIUS);

	glEnd();

	glPopMatrix();

	// Now we draw the wood panels
	// 139-69-19
	// Or in neon
	// 131, 245, 44
	//GLfloat	specular2[4] = { 1.39f,  0.69f, 0.19f, 1.0f };
	//GLfloat	ambient2[4] = { 1.39f,  0.69f, 0.19f, 1.0f };
	//GLfloat	diffuse2[4] = { 1.0f, 0.8f, 0.2f, 1.0f };


	GLfloat	specular2[4] = { 1.31f,  2.45f, 0.44f, 1.0f };
	GLfloat	ambient2[4] = { 1.31f,  2.45f, 0.44f, 1.0f };
	GLfloat	diffuse2[4] = { 0.1f, 0.08f, 0.02f, 1.0f };

	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0f);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular2);

	glPushMatrix();

	GLdouble PanelXMax = 3.102225;
	GLdouble PanelXMin = 2.602225;
	GLdouble PanelZMax = 5.602225;
	GLdouble PanelZMin = 5.102225;
	GLdouble height = 0.0;

	// Panel 1
	glBegin(GL_QUADS);

	glVertex3d(-PanelXMin, height, -PanelZMin);
	glVertex3d(-PanelXMin, height, PanelZMin);
	glVertex3d(-PanelXMax, height, PanelZMax);
	glVertex3d(-PanelXMax, height, -PanelZMax);

	glEnd();

	// Panel 2
	glBegin(GL_QUADS);

	glVertex3d(-PanelXMin, height, -PanelZMin);
	glVertex3d(PanelXMin, height, -PanelZMin);
	glVertex3d(PanelXMax, height, -PanelZMax);
	glVertex3d(-PanelXMax, height, -PanelZMax);

	glEnd();

	// Panel 3
	glBegin(GL_QUADS);

	glVertex3d(PanelXMin, height, -PanelZMin);
	glVertex3d(PanelXMin, height, PanelZMin);
	glVertex3d(PanelXMax, height, PanelZMax);
	glVertex3d(PanelXMax, height, -PanelZMax);

	glEnd();

	// Panel 4
	glBegin(GL_QUADS);

	glVertex3d(-PanelXMin, height, PanelZMin);
	glVertex3d(PanelXMin, height, PanelZMin);
	glVertex3d(PanelXMax, height, PanelZMax);
	glVertex3d(-PanelXMax, height, PanelZMax);

	glEnd();

	// Top Side Wall
	glBegin(GL_QUADS);

	glNormal3d(0.0, 0.0, 1.0);
	glVertex3d(PanelXMax, height, PanelZMax);
	glNormal3d(0.0, 0.0, 1.0);
	glVertex3d(PanelXMax, 0.8, PanelZMax);
	glNormal3d(0.0, 0.0, 1.0);
	glVertex3d(-PanelXMax, 0.8, PanelZMax);
	glNormal3d(0.0, 0.0, 1.0);
	glVertex3d(-PanelXMax, height, PanelZMax);

	glEnd();

	// Bottom Side Wall
	glBegin(GL_QUADS);

	glNormal3d(0.0, 0.0, -1.0);
	glVertex3d(PanelXMax, height, -PanelZMax);
	glNormal3d(0.0, 0.0, -1.0);
	glVertex3d(PanelXMax, 0.8, -PanelZMax);
	glNormal3d(0.0, 0.0, -1.0);
	glVertex3d(-PanelXMax, 0.8, -PanelZMax);
	glNormal3d(0.0, 0.0, -1.0);
	glVertex3d(-PanelXMax, height, -PanelZMax);

	glEnd();

	// Right Side Wall
	glBegin(GL_QUADS);

	glNormal3d(1.0, 0.0, 0.0);
	glVertex3d(PanelXMax, height, PanelZMax);
	glNormal3d(1.0, 0.0, 0.0);
	glVertex3d(PanelXMax, 0.8, PanelZMax);
	glNormal3d(1.0, 0.0, 0.0);
	glVertex3d(PanelXMax, 0.8, -PanelZMax);
	glNormal3d(1.0, 0.0, 0.0);
	glVertex3d(PanelXMax, height, -PanelZMax);

	glEnd();

	// Right Side Wall
	glBegin(GL_QUADS);

	glNormal3d(-1.0, 0.0, 0.0);
	glVertex3d(-PanelXMax, height, PanelZMax);
	glNormal3d(-1.0, 0.0, 0.0);
	glVertex3d(-PanelXMax, 0.8, PanelZMax);
	glNormal3d(-1.0, 0.0, 0.0);
	glVertex3d(-PanelXMax, 0.8, -PanelZMax);
	glNormal3d(-1.0, 0.0, 0.0);
	glVertex3d(-PanelXMax, height, -PanelZMax);

	glEnd();

	glPopMatrix();


}

// LoadPockets will vertex array the pockets.  They are just
// black circles and they are drawn with triangles.
void LoadPockets() {

	glNewList(PocketList, GL_COMPILE);

	for (int i=0; i < 6; i++) {

		glPushMatrix();
		glTranslated(POCKET_CENTER_X[i], 0.08, POCKET_CENTER_Z[i]);
		glBegin(GL_TRIANGLES);
		for (int i = 0; i <=360; i+=15) {
			GLdouble PointX = (POCKET_RADIUS * cos(i / DEGREES));
			GLdouble PointZ = (POCKET_RADIUS * sin(i / DEGREES));
			GLdouble PointX_2 = (POCKET_RADIUS * cos((i+15) / DEGREES));
			GLdouble PointZ_2 = (POCKET_RADIUS * sin((i+15) / DEGREES));

			glVertex3d(0.0, 0.0, 0.0);
			glVertex3d(PointX, 0.0, PointZ);
			glVertex3d(PointX_2, 0.0, PointZ_2);

		}
		glEnd();

		glPopMatrix();

	}

	glEndList();

}

// DrawPockets will call a vertex array list to draw the pockets
void DrawPockets() {
	
	GLfloat	specular[4]	= { 0.0f,  0.0f, 0.0f, 1.0f };
	GLfloat	ambient[4] = { 0.0f,  0.0f, 0.0f, 1.0f };
	GLfloat	diffuse[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 10.0f);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);

	glCallList(PocketList);

}

// This function will update the positions of the balls each turn.
// It will also account for friction and ball rotation.
void MoveBalls() {


	for (int i = 0; i < numBalls; i++) {

		if (Balls[i]->status == LIVE) {
		
			Friction(&Balls[i]->velocity);

			Balls[i]->lastPosition.x = Balls[i]->position.x;
			Balls[i]->lastPosition.z = Balls[i]->position.z;

			Balls[i]->position.x += Balls[i]->velocity.x;
			Balls[i]->position.z += Balls[i]->velocity.z;

		}

	}

}

// This function finds a good position for the camera behind the
// cue ball.
void CueBallView(double* EyeX, double* EyeZ) {

	double direction = getDirection(&Balls[0]->position);
	direction += RotateX;

	*EyeX = Balls[0]->position.x + (3.0 * cos(direction / DEGREES));
	*EyeZ = Balls[0]->position.z + (3.0 * sin(direction / DEGREES));

}

// This function will check each ball against all the other balls (but each pair
// only once) for a collision.  If there is a collision, the correct new velocities
// will be determined.
// Comments on how it works can be found inside this section of code as you read.
void CheckBallCollisions() {

	for (int i=0; i < numBalls; i++) {

		for (int j=i+1; j < numBalls; j++) {

			if (Balls[i]->status == LIVE && Balls[j]->status == LIVE) {

				double bdistance = distance(&Balls[i]->position, &Balls[j]->position);
				double zonesize = (2 * BALL_RADIUS);

				if (bdistance < zonesize) {

					vec3df pos((float)(Balls[i]->position.x - EyeX), 0.0, (float)(Balls[i]->position.z - EyeZ));
					int sound = (int)(rand() % 3);
					//printf("Sound: %d\n", sound);
					switch (sound) {
					case 0:
						engine->play3D("ballhithard.ogg",pos);
						break;
					case 1:
						engine->play3D("ballhitmedium.ogg",pos);
						break;
					case 2:
						engine->play3D("ballhitlight.ogg",pos);
						break;
					}
					BallCollisions++;
					if (BallCollisions == 1) FirstCollision = j;

					// printf("COLLISION DETECTED!!\n");

					// mid is the midpoint between the two balls.
					midpoint(mid, &Balls[i]->position, &Balls[j]->position);

					// Vector mid_to_B1 is the vector from the collision point to the center of ball 1
					mid_to_B1->x = Balls[i]->position.x - mid->x;
					mid_to_B1->y = 0.0;
					mid_to_B1->z = Balls[i]->position.z - mid->z;

					// Vector mid_to_B2 is the vector from the collision point to the center of ball 2
					mid_to_B2->x = Balls[j]->position.x - mid->x;
					mid_to_B2->y = 0.0;
					mid_to_B2->z = Balls[j]->position.z - mid->z;

					// Vector vc is the vector between the center of the two balls
					vc->x = Balls[j]->position.x - Balls[i]->position.x;
					vc->y = 0.0;
					vc->z = Balls[j]->position.z - Balls[i]->position.z;

					// Vector vcn is the normal of vector vc.
					vcn->x = vc->z;
					vcn->y = 0.0;
					vcn->z = -vc->x;

					// Here we are calculating unit vectors from the midpoint to
					// the center of each of the balls.  We are going to use these
					// to offset the balls by a minimum of the collision distance
					// (2 radii) each time they collide.  We need to do this to avoid
					// balls clinging to eachother when they get stuck in a "double
					// collision" (two collisions of the same ball pair back to back).
					unitvector(u_to_B1, mid_to_B1);
					unitvector(u_to_B2, mid_to_B2);

					double dir_of_m_B1 = getDirection(u_to_B1);
					double dir_of_m_B2 = getDirection(u_to_B2);
					//double dir_of_vc = getDirection(vc);
					//double dir_of_vcn = getDirection(vcn);
					//double dir_v1 = getDirection(&Balls[i]->velocity);

					// Offset the balls by 2 balls radii to ensure the tangling effect
					// does not occur
					double factor = ((zonesize - bdistance) / 2) + 0.0001;
					Balls[i]->position.x = Balls[i]->position.x + (factor * u_to_B1->x);
					Balls[i]->position.z = Balls[i]->position.z + (factor * u_to_B1->z);
					Balls[j]->position.x = Balls[j]->position.x + (factor * u_to_B2->x);
					Balls[j]->position.z = Balls[j]->position.z + (factor * u_to_B2->z);

					// Here we calculate some angles that i dont even use... ????
					// [[this is useless]]
					// double theta = abs(dir_of_vc - dir_v1);
					// double phi = abs(dir_of_vcn - dir_v1);

					// I am preparing vectors for the projections of the velocities of each ball
					// onto either vc or its normal, vcn.

					// Now I take the projections.  As an example:
					// prj_v1vc mean projection of velocity of ball 1 onto vc.
					projection(prj_v1vc, &Balls[i]->velocity, vc);
					projection(prj_v1vcn, &Balls[i]->velocity, vcn);
					projection(prj_v2vc, &Balls[j]->velocity, vc);
					projection(prj_v2vcn, &Balls[j]->velocity, vcn);

					// Ball 1's new velocity vector is the projection of its velocity onto vcn
					// plus the projection of ball 2's velocity onto vc.
					Balls[i]->velocity.x = prj_v1vcn->x + prj_v2vc->x;
					Balls[i]->velocity.y = prj_v1vcn->y + prj_v2vc->y;
					Balls[i]->velocity.z = prj_v1vcn->z + prj_v2vc->z;

					// Ball 2's new velocity vector is the projection of its velocity onto vcn
					// plus the projection of ball 1's velocity onto vc.
					Balls[j]->velocity.x = prj_v1vc->x + prj_v2vcn->x;
					Balls[j]->velocity.y = prj_v1vc->y + prj_v2vcn->y;
					Balls[j]->velocity.z = prj_v1vc->z + prj_v2vcn->z;

				}

			}

		}

	}

}

// CheckSingleCollision will check if a given BallNum is in contact
// with any other ball on the board.  It is useful for when the user
// is placing a scratched ball or if I am checking to see if the 9
// ball can be placed somewhere.
bool CheckSingleCollision(int BallNum) {

	for (int j=0; j < numBalls; j++) {

		if (j == BallNum) j++;

		if (Balls[j]->status == LIVE) {

			double bdistance = distance(&Balls[BallNum]->position, &Balls[j]->position);
			double zonesize = (2 * BALL_RADIUS);

			if (bdistance < zonesize) {

				return true;

			}

		}

	}

	return false;

}

// Place9Ball is a little algorithm that scans in concentric circles
// centering on where the 9 ball was originally placed (0.0x, 3.4z).
// When it finds an area that is not occupied by a ball it will
// exit the function, placing the 9 ball in that spot.
void Place9Ball() {

	for (int i = 0; i < 10; i++) {

		for (int j = 0; j <=360; j+=15) {

			GLdouble PointX = ((BALL_RADIUS * i) * cos(j / DEGREES));
			GLdouble PointZ = ((BALL_RADIUS * i) * sin(j / DEGREES));

			PointZ += 3.4;

			Balls[9]->position.x = PointX;
			Balls[9]->position.z = PointZ;

			bool tmp = CheckSingleCollision(9);
			if (tmp == false) {
				goto place9;
			}

		}

	}

	place9:
	return;

}

// LoadWalls loads all the positions of the bumpers.  The values were once
// calculated but I have since plugged them in manually to speed things up.
void LoadWalls() {

	////////////////////
	//// TABLE WALLS ///
	////////////////////

	BL_to_BR1->x = -2.298225;
	BL_to_BR1->y = 0.0;
	BL_to_BR1->z = -4.886225;

	BL_to_BR2->x = 2.298225;
	BL_to_BR2->y = 0.0;
	BL_to_BR2->z = -4.886225;


	BR_to_MR1->x = 2.386225;
	BR_to_MR1->y = 0.0;
	BR_to_MR1->z = -4.710225;
	
	BR_to_MR2->x = 2.386225;
	BR_to_MR2->y = 0.0;
	BR_to_MR2->z = -0.176;


	MR_to_TR1->x = 2.386225;
	MR_to_TR1->y = 0.0;
	MR_to_TR1->z = 0.176;

	MR_to_TR2->x = 2.386225;
	MR_to_TR2->y = 0.0;
	MR_to_TR2->z = 4.710225;


	TR_to_TL1->x = 2.298225;
	TR_to_TL1->y = 0.0;
	TR_to_TL1->z = 4.886225;

	TR_to_TL2->x = -2.298225;
	TR_to_TL2->y = 0.0;
	TR_to_TL2->z = 4.886225;


	TL_to_ML1->x = -2.386225;
	TL_to_ML1->y = 0.0;
	TL_to_ML1->z = 4.710225;

	TL_to_ML2->x = -2.386225;
	TL_to_ML2->y = 0.0;
	TL_to_ML2->z = 0.176;


	ML_to_BL1->x = -2.386225;
	ML_to_BL1->y = 0.0;
	ML_to_BL1->z = -0.176;

	ML_to_BL2->x = -2.386225;
	ML_to_BL2->y = 0.0;
	ML_to_BL2->z = -4.710225;

	////////////////////////
	///// POCKET WALLS /////
	////////////////////////

	// BOTTOM LEFT
	BL_Lower1->x = 2.298225;
	BL_Lower1->y = 0.0;
	BL_Lower1->z = -4.886225;

	BL_Lower2->x = -2.386225;
	BL_Lower2->y = 0.0;
	BL_Lower2->z = -5.150225;


	BL_Upper1->x = -2.386225;
	BL_Upper1->y = 0.0;
	BL_Upper1->z = -4.710225;

	BL_Upper2->x = -2.650225;
	BL_Upper2->y = 0.0;
	BL_Upper2->z = -4.886225;
	
	// BOTTOM RIGHT CORNER
	BR_Lower1->x = -2.298225;
	BR_Lower1->y = 0.0;
	BR_Lower1->z = -4.886225;

	BR_Lower2->x = 2.386225;
	BR_Lower2->y = 0.0;
	BR_Lower2->z = -5.150225;


	BR_Upper1->x = 2.386225;
	BR_Upper1->y = 0.0;
	BR_Upper1->z = -4.710225;

	BR_Upper2->x = 2.650225;
	BR_Upper2->y = 0.0;
	BR_Upper2->z = -4.886225;

	// TOP RIGHT CORNER
	TR_Lower1->x = 2.386225;
	TR_Lower1->y = 0.0;
	TR_Lower1->z = 4.710225;

	TR_Lower2->x = 2.650225;
	TR_Lower2->y = 0.0;
	TR_Lower2->z = 4.886225;


	TR_Upper1->x = 2.210225;
	TR_Upper1->y = 0.0;
	TR_Upper1->z = 4.974225;

	TR_Upper2->x = 2.386225;
	TR_Upper2->y = 0.0;
	TR_Upper2->z = 5.150225;

	// TOP LEFT CORNER
	TL_Lower1->x = -2.386225;
	TL_Lower1->y = 0.0;
	TL_Lower1->z = 4.710225;

	TL_Lower2->x = -2.650225;
	TL_Lower2->y = 0.0;
	TL_Lower2->z = 4.886225;


	TL_Upper1->x = -2.298225;
	TL_Upper1->y = 0.0;
	TL_Upper1->z = 4.886225;

	TL_Upper2->x = -2.386225;
	TL_Upper2->y = 0.0;
	TL_Upper2->z = 5.150225;

	// MID LEFT POCKET
	ML_Lower1->x = -2.386225;
	ML_Lower1->y = 0.0;
	ML_Lower1->z = -0.176;

	ML_Lower2->x = -3.0;
	ML_Lower2->y = 0.0;
	ML_Lower2->z = -0.176;


	ML_Upper1->x = -2.386225;
	ML_Upper1->y = 0.0;
	ML_Upper1->z = 0.176;

	ML_Upper2->x = -3.0;
	ML_Upper2->y = 0.0;
	ML_Upper2->z = 0.176;

	// MID RIGHT POCKET
	MR_Lower1->x = 2.386225;
	MR_Lower1->y = 0.0;
	MR_Lower1->z = -0.176;

	MR_Lower2->x = 3.0;
	MR_Lower2->y = 0.0;
	MR_Lower2->z = -0.176;


	MR_Upper1->x = 2.386225;
	MR_Upper1->y = 0.0;
	MR_Upper1->z = 0.176;

	MR_Upper2->x = 3.0;
	MR_Upper2->y = 0.0;
	MR_Upper2->z = 0.176;


}

// CheckWallCollisions is a brute force wall checking algorithm.  It uses line
// intersection and the current and last position of the ball to check whether
// the balls path has crossed any piece of the bumper.  There are six large bumpers
// and 12 little bumpers associated with each of the 6 pockets.
//
// Lots of the variables in this function such as TR_to_TL1 etc etc are declared in
// the header file and set in LoadWalls();
void CheckWallCollisions() {

	bool hit_occured = false;

	for (int i=0; i < numBalls; i++) {

		if (Balls[i]->status == LIVE && Balls[i]->JustHitWall == false) {

			int TopCollide = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, TR_to_TL1, TR_to_TL2, TopInt);
			int BottomCollide = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, BL_to_BR1, BL_to_BR2, BottomInt);

			int RightTopCollide = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, MR_to_TR1, MR_to_TR2, RightInt);
			int RightBottomCollide = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, BR_to_MR1, BR_to_MR2, RightInt);
			
			int LeftTopCollide = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, TL_to_ML1, TL_to_ML2, LeftInt);
			int LeftBottomCollide = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, ML_to_BL1, ML_to_BL2, LeftInt);

			int BLLowCorner = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, BL_Lower1, BL_Lower2, BL_Low_Int);
			int BLHighCorner = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, BL_Upper1, BL_Upper2, BL_High_Int);

			int BRLowCorner = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, BR_Lower1, BR_Lower2, BR_Low_Int);
			int BRHighCorner = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, BR_Upper1, BR_Upper2, BR_High_Int);

			int TRLowCorner = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, TR_Lower1, TR_Lower2, TR_Low_Int);
			int TRHighCorner = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, TR_Upper1, TR_Upper2, TR_High_Int);

			int TLLowCorner = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, TL_Lower1, TL_Lower2, TL_Low_Int);
			int TLHighCorner = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, TL_Upper1, TL_Upper2, TL_High_Int);

			int MLLowBumper = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, ML_Lower1, ML_Lower2, ML_Low_Int);
			int MLHighBumper = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, ML_Upper1, ML_Upper2, ML_High_Int);

			int MRLowBumper = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, MR_Lower1, MR_Lower2, MR_Low_Int);
			int MRHighBumper = lines_intersect(&Balls[i]->position, &Balls[i]->lastPosition, MR_Upper1, MR_Upper2, MR_High_Int);

			// Colliding into the high x wall
			if (RightTopCollide == DO_INTERSECT) {

				Balls[i]->position.x = MR_to_TR1->x - 0.00001;
				Balls[i]->velocity.x = -(0.8) * Balls[i]->velocity.x;
				RailsHit++;
				hit_occured = true;

			}

			if (RightBottomCollide == DO_INTERSECT) {

				Balls[i]->position.x = BR_to_MR1->x - 0.00001;
				Balls[i]->velocity.x = -(0.8) * Balls[i]->velocity.x;
				RailsHit++;
				hit_occured = true;

			}

			// Colliding into low high x wall
			if (LeftTopCollide == DO_INTERSECT) {

				Balls[i]->position.x = TL_to_ML1->x + 0.00001;
				Balls[i]->velocity.x = -(0.8) * Balls[i]->velocity.x;
				RailsHit++;
				hit_occured = true;
			}

			if (LeftBottomCollide == DO_INTERSECT) {

				Balls[i]->position.x = ML_to_BL1->x + 0.00001;
				Balls[i]->velocity.x = -(0.8) * Balls[i]->velocity.x;
				RailsHit++;
				hit_occured = true;
			}

			// Colliding into the high z wall
			if (TopCollide == DO_INTERSECT) {

				Balls[i]->position.z = TR_to_TL1->z - 0.00001;
				Balls[i]->velocity.z = -(0.8) * Balls[i]->velocity.z;
				RailsHit++;
				hit_occured = true;
			}

			// Colliding into the low z wall
			if (BottomCollide == DO_INTERSECT) {

				Balls[i]->position.z = BL_to_BR1->z + 0.00001;
				Balls[i]->velocity.z = -(0.8) * Balls[i]->velocity.z;
				RailsHit++;
				hit_occured = true;
			}

			// Collisions for the bottom left corner bumpers

			if (BLLowCorner == DO_INTERSECT) {

				double XVelocity = Balls[i]->velocity.x;

				Balls[i]->position.x = BL_Low_Int->x - 0.00001;
				Balls[i]->velocity.x = (0.8) * Balls[i]->velocity.z;

				Balls[i]->position.z = BL_Low_Int->z + 0.00001;
				Balls[i]->velocity.z = (0.8) * XVelocity;

				RailsHit++;
				hit_occured = true;
			}

			if (BLHighCorner == DO_INTERSECT) {

				double XVelocity = Balls[i]->velocity.x;

				Balls[i]->position.x = BL_High_Int->x + 0.00001; 
				Balls[i]->velocity.x = (0.8) * Balls[i]->velocity.z;

				Balls[i]->position.z = BL_High_Int->z - 0.00001;
				Balls[i]->velocity.z = (0.8) * XVelocity;

				RailsHit++;
				hit_occured = true;
			}

			// Collisions for the bottom right corner bumpers

			if (BRLowCorner == DO_INTERSECT) {

				double XVelocity = Balls[i]->velocity.x;

				Balls[i]->position.x = BR_Low_Int->x + 0.00001;
				Balls[i]->velocity.x = -(0.8) * Balls[i]->velocity.z;

				Balls[i]->position.z = BR_Low_Int->z + 0.00001;
				Balls[i]->velocity.z = -(0.8) * XVelocity;

				RailsHit++;
				hit_occured = true;
			}

			if (BRHighCorner == DO_INTERSECT) {

				double XVelocity = Balls[i]->velocity.x;

				Balls[i]->position.x = BR_High_Int->x - 0.00001; 
				Balls[i]->velocity.x = -(0.8) * Balls[i]->velocity.z;

				Balls[i]->position.z = BR_High_Int->z - 0.00001;
				Balls[i]->velocity.z = -(0.8) * XVelocity;

				RailsHit++;
				hit_occured = true;
			}

			// Collisions for the top right corner bumpers

			if (TRLowCorner == DO_INTERSECT) {

				double XVelocity = Balls[i]->velocity.x;

				Balls[i]->position.x = TR_Low_Int->x + 0.00001;
				Balls[i]->velocity.x = (0.8) * Balls[i]->velocity.z;

				Balls[i]->position.z = TR_Low_Int->z + 0.00001;
				Balls[i]->velocity.z = (0.8) * XVelocity;

				RailsHit++;
				hit_occured = true;
			}

			if (TRHighCorner == DO_INTERSECT) {

				double XVelocity = Balls[i]->velocity.x;

				Balls[i]->position.x = TR_High_Int->x - 0.00001; 
				Balls[i]->velocity.x = (0.8) * Balls[i]->velocity.z;

				Balls[i]->position.z = TR_High_Int->z - 0.00001;
				Balls[i]->velocity.z = (0.8) * XVelocity;

				RailsHit++;
				hit_occured = true;
			}

			// Collisions for the top left corner bumpers

			if (TLLowCorner == DO_INTERSECT) {

				double XVelocity = Balls[i]->velocity.x;

				Balls[i]->position.x = TL_Low_Int->x + 0.00001;
				Balls[i]->velocity.x = -(0.8) * Balls[i]->velocity.z;

				Balls[i]->position.z = TL_Low_Int->z + 0.00001;
				Balls[i]->velocity.z = -(0.8) * XVelocity;

				RailsHit++;
				hit_occured = true;
			}

			if (TLHighCorner == DO_INTERSECT) {

				double XVelocity = Balls[i]->velocity.x;

				Balls[i]->position.x = TL_High_Int->x - 0.00001; 
				Balls[i]->velocity.x = -(0.8) * Balls[i]->velocity.z;

				Balls[i]->position.z = TL_High_Int->z - 0.00001;
				Balls[i]->velocity.z = -(0.8) * XVelocity;

				RailsHit++;
				hit_occured = true;
			}

			// Collisions for mid left pocket bumpers

			if (MLHighBumper == DO_INTERSECT) {

				Balls[i]->position.z = ML_High_Int->z - 0.00001;
				Balls[i]->velocity.z = -(0.8) * Balls[i]->velocity.z;
				RailsHit++;
				hit_occured = true;
			}

			if (MLLowBumper == DO_INTERSECT) {

				Balls[i]->position.z = ML_Low_Int->z + 0.00001;
				Balls[i]->velocity.z = -(0.8) * Balls[i]->velocity.z;
				RailsHit++;
				hit_occured = true;
			}

			// Collisions for mid right pocket bumpers

			if (MRHighBumper == DO_INTERSECT) {

				Balls[i]->position.z = MR_High_Int->z - 0.00001;
				Balls[i]->velocity.z = -(0.8) * Balls[i]->velocity.z;
				RailsHit++;
				hit_occured = true;
			}

			if (MRLowBumper == DO_INTERSECT) {

				Balls[i]->position.z = MR_Low_Int->z + 0.00001;
				Balls[i]->velocity.z = -(0.8) * Balls[i]->velocity.z;
				RailsHit++;
				hit_occured = true;

			}

			if (hit_occured == true) {
				vec3df pos((float)(Balls[i]->position.x - EyeX), 0.0, (float)(Balls[i]->position.z - EyeZ));
				int sound = (int)(rand() % 3);
				//printf("Sound: %d\n", sound);
				switch (sound) {
				case 0:
					engine->play3D("bank1.ogg",pos);
					break;
				case 1:
					engine->play3D("bank2.ogg",pos);
					break;
				case 2:
					engine->play3D("bank1.ogg",pos);
					break;
				}
			}


		}

		if (hit_occured == true) {
			Balls[i]->JustHitWall = true;
		} else {
			Balls[i]->JustHitWall = false;
		}

	}

}

// Check if any of the balls have entered any of the pockets.
void CheckPockets(void) {

	for (int i = 0; i < numBalls; i++) {

		for (int j = 0; j < 6; j++) {

			thisPocket->x = POCKET_CENTER_X[j];
			thisPocket->y = 0.0;
			thisPocket->z = POCKET_CENTER_Z[j];

			double bdistance = distance(&Balls[i]->position, thisPocket);

			if (bdistance < (POCKET_RADIUS + BALL_RADIUS)) {

				if (Balls[i]->status == LIVE) {
					BallsPocketed++;
					vec3df pos((float)(Balls[i]->position.x - EyeX), 0.0, (float)(Balls[i]->position.z - EyeZ));
					engine->play3D("pocket.ogg",pos);

				}
				Balls[i]->status = IN_POCKET;
				if (i == 9) NineBallPocketed = true;

			}
		
		}

	}

}

// This function checks if all the balls have stopped moving
// if they have we start a new shot.
void CheckBallsStopped(void) {

	int liveCount = 0;
	int stopped = 0;

	for (int i=0; i < numBalls; i++) {

		if (Balls[i]->status == LIVE) liveCount++;

		if (Balls[i]->velocity.x == 0.0 && Balls[i]->velocity.z == 0.0) stopped++;

	}


	if (liveCount == stopped) {

		if (Balls[0]->status == IN_POCKET) {
			GameState = SCRATCH;
			NextShot();
		} else {
			GameState = WAITING;
			NextShot();
		}

	}

}

// NextShot() will do some calculations and prepare for the next player shot
// Depending on what happened during the last shot, it will decide who shoots
// next or if the game is over.
void NextShot() {

	// CHECK FOR ANY FOULS AND IF SO CHANGE PLAYERS AND SCRATCH THE BALL
	if ((RailsHit == 0 && BallsPocketed == 0) || GameState == SCRATCH || FirstCollision != TargetBall) {

		int Player = CurrentPlayer;

		if (Player == PLAYER_ONE) CurrentPlayer = PLAYER_TWO;
		if (Player == PLAYER_TWO) CurrentPlayer = PLAYER_ONE;

		if (Balls[9]->status == IN_POCKET) {
			Place9Ball();
			NineBallPocketed = false;
			Balls[9]->status = LIVE;
			Balls[9]->velocity.x = 0.0;
			Balls[9]->velocity.z = 0.0;
		}

		Balls[0]->position.x = 0.0;
		Balls[0]->position.z = -3.0;
		Balls[0]->velocity.x = 0.0;
		Balls[0]->velocity.z = 0.0;
		Balls[0]->status = LIVE;
		GameState = SCRATCH;
		ChangeView(TABLE_VIEW);

	} else if (BallsPocketed == 0 && NineBallPocketed == false) {

		int Player = CurrentPlayer;

		if (Player == PLAYER_ONE) CurrentPlayer = PLAYER_TWO;
		if (Player == PLAYER_TWO) CurrentPlayer = PLAYER_ONE;

	} else if (Balls[9]->status == IN_POCKET) {
		
		printf("\n\nPLAYER %d WINS!!!!\n\nPress Space to start a new game.\n\n", CurrentPlayer+1);
		
		if (CurrentPlayer == PLAYER_ONE) {
			GameState = PLAYER_1_WINS;
		} else {
			GameState = PLAYER_2_WINS;
		}

	}

	// PRINT A STRING:
	if (CurrentPlayer == PLAYER_ONE) ShotString = "PLAYER ONE SHOOTS\n";
	if (CurrentPlayer == PLAYER_TWO) ShotString = "PLAYER TWO SHOOTS\n";
	printf("%s", ShotString);

	// FIND A NEW TARGET BALL
	int LowestBall = 1;
	for (int i = 1; i < 10; i++) {
		if (Balls[i]->status == LIVE) {
			LowestBall = i;
			break;
		}
	}

	TargetBall = LowestBall;

	RailsHit = 0;
	BallsPocketed = 0;
	BallCollisions = 0;

}

//  NewGame just refreshes some variables such as ball position
// and velocity and status etc to prepare for a new game.
void NewGame() {

	Balls[0]->position.x = 0.0;
	Balls[0]->position.z = -3.0;
	Balls[0]->lastPosition.x = 0.0;
	Balls[0]->lastPosition.z = -3.0;

	Balls[1]->position.x = 0.0;
	Balls[1]->position.z = 3.0;
	Balls[1]->lastPosition.x = 0.0;
	Balls[1]->lastPosition.z = 3.0;

	Balls[2]->position.x = -0.1;
	Balls[2]->position.z = 3.2;
	Balls[2]->lastPosition.x = -0.1;
	Balls[2]->lastPosition.z = 3.2;

	Balls[3]->position.x = 0.1;
	Balls[3]->position.z = 3.2;
	Balls[3]->lastPosition.x = 0.1;
	Balls[3]->lastPosition.z = 3.2;

	Balls[4]->position.x = -0.2;
	Balls[4]->position.z = 3.4;
	Balls[4]->lastPosition.x = -0.2;
	Balls[4]->lastPosition.z = 3.4;

	Balls[5]->position.x = 0.2;
	Balls[5]->position.z = 3.4;
	Balls[5]->lastPosition.x = 0.2;
	Balls[5]->lastPosition.z = 3.4;

	Balls[6]->position.x = -0.1;
	Balls[6]->position.z = 3.6;
	Balls[6]->lastPosition.x = -0.1;
	Balls[6]->lastPosition.z = 3.6;

	Balls[7]->position.x = 0.1;
	Balls[7]->position.z = 3.6;
	Balls[7]->lastPosition.x = 0.1;
	Balls[7]->lastPosition.z = 3.6;

	Balls[8]->position.x = 0.0;
	Balls[8]->position.z = 3.8;
	Balls[8]->lastPosition.x = 0.0;
	Balls[8]->lastPosition.z = 3.8;

	Balls[9]->lastPosition.x = 0.0;
	Balls[9]->lastPosition.z = 3.4;

	for (int i = 0; i < numBalls; i++) {

		Balls[i]->velocity.x = 0;
		Balls[i]->velocity.z = 0;

		Balls[i]->Rotation = 0.0;
		Balls[i]->status = LIVE;

	}
	
	NineBallPocketed = false;
	Breaking = true;
	TargetBall = 1;

	GameState = SCRATCH;

}

// This function evaluates the friction on a velocity vector
// It's not trivial and must be done using a vector in the opposite
// direction that the ball is traveling.  Thus sins and cosines are
// involved.
void Friction(VECTOR* v) {

	double dir = getDirection(v);

	double xFactor = abs(cos(dir / DEGREES));
	double zFactor = abs(sin(dir / DEGREES));

	// Perform X friction
	if (abs(v->x) <= KE_FRICTION) {
		v->x = 0;
	} else {
		if (v->x < 0) {
			v->x += (KE_FRICTION * xFactor);
		} else {
			v->x -= (KE_FRICTION * xFactor);
		}
	}
	// Perform Z friction
	if (abs(v->z) <= KE_FRICTION) {
		v->z = 0;
	} else {
		if (v->z < 0) {
			v->z += (KE_FRICTION * zFactor);
		} else {
			v->z -= (KE_FRICTION * zFactor);
		}
	}

}

///////////////////////////////////
// TIMER, DISPLAY & GL FUNCTIONS //
///////////////////////////////////

// THis is the timer function.  Depending on GameState,
// different function will be called.
// Notice that if GameState == BALLS_MOVING that's where
// all the beefy functions get called and the 'physics engine'
// is at work.
void timerFunction(int value) {

	glutTimerFunc(16, timerFunction, 1);
	//MouseFunction();

	char* myStr = (char*)malloc(sizeof(char) * 20);
	switch (GameState) {
	case WAITING:
		break;
	case DRAWING_CUE:
		InShot((double)MouseX, (double)LastMouseX);
		break;
	case TAKING_SHOT:
		break;
	case SCRATCH:
		break;
	case BALLS_MOVING:
		MoveBalls();
		CheckPockets();
		CheckBallCollisions();
		CheckWallCollisions();
		CheckBallsStopped();
		break;
	case PLAYER_1_WINS:
		break;
	case PLAYER_2_WINS:
		break;
	}

	glutPostRedisplay();

}

// This is the display function.
// First it establishes the view based on ViewMode
// Second it draws the scene depending on GameState
// and the Game's actual state.
void display(void)
{
	

	static	GLfloat	light_position[] = { 0.0 , -1.0 , -1.0 , 1.0 };

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	double CueBallX = Balls[0]->position.x;
	double CueBallZ = Balls[0]->position.z;

	if (ViewMode == TABLE_VIEW) {
		
		gluLookAt(0.0, -15.0, -0.1 , 0.0 , 0.0 , 0.0 , 0.0 , -1.0 , 0.0);
		glRotated(TableViewRotateY, 1.0, 0.0, 0.0); // We rotate our mouse y on the x axis (odd, i know)
		glRotated(TableViewRotateX, 0.0, 1.0, 0.0); // We rotate our mouse x on the y axis (odd, i know)

	} else if (ViewMode == CUE_BALL_VIEW) {

		CueBallView(&EyeX, &EyeZ);
		gluLookAt(EyeX, -1.5, EyeZ, CueBallX, 0.0, CueBallZ , 0.0 , -1.0 , 0.0);
		FrozenX = CueBallX;
		FrozenZ = CueBallZ;
		//glRotated((RotateX % 360), CueBallX, 0.0, 0.0); // We rotate our mouse y on the x axis (odd, i know)
		//glRotated((RotateY % 360), 0.0, 1.0, 0.0); // We rotate our mouse x on the y axis (odd, i know)

	} else if (ViewMode == FROZEN) {

		gluLookAt(EyeX, -1.5, EyeZ, FrozenX, 0.0, FrozenZ , 0.0 , -1.0 , 0.0);

	}

	glLightfv(GL_LIGHT0 , GL_POSITION , light_position);

	glPushMatrix();

	if (GameState == TAKING_SHOT || GameState == DRAWING_CUE) DrawCueStick();
	DrawTable();
	DrawBalls();
	DrawPockets();

	GLdouble height = -0.02;

	glDisable(GL_LIGHTING);
	glColor3d(1.0, 0.0, 0.0);
	glPointSize(2.0f);
	glLineWidth(2.0f);

	//glBegin(GL_LINES);

	//glVertex3d(0.0, 0.0, 0.0);
	//glVertex3d(1.0, 0.0, 0.0);

	//glVertex3d(0.0, 0.0, 0.0);
	//glVertex3d(0.0, 0.0, 1.0);

	//glEnd();

	glBegin(GL_LINES);

	// BOTTOM LEFT
	glVertex3d(POCKET_CENTER_X[0] + 4*BALL_RADIUS, -height, POCKET_CENTER_Z[0] + BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[0] + 2*BALL_RADIUS, -height, POCKET_CENTER_Z[0] - BALL_RADIUS);
	
	glVertex3d(POCKET_CENTER_X[0] + BALL_RADIUS, -height, POCKET_CENTER_Z[0] + 4*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[0] - BALL_RADIUS, -height, POCKET_CENTER_Z[0] + 2*BALL_RADIUS);

	// InBETWEEN
	glVertex3d(POCKET_CENTER_X[0] + 4*BALL_RADIUS, -height, POCKET_CENTER_Z[0] + BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[5] - 4*BALL_RADIUS, -height, POCKET_CENTER_Z[5] + BALL_RADIUS);

	// BOTTOM RIGHT
	////////////////
	glVertex3d(POCKET_CENTER_X[5] - 4*BALL_RADIUS, -height, POCKET_CENTER_Z[5] + BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[5] - 2*BALL_RADIUS, -height, POCKET_CENTER_Z[5] - BALL_RADIUS);
	
	glVertex3d(POCKET_CENTER_X[5] - BALL_RADIUS, -height, POCKET_CENTER_Z[5] + 4*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[5] + BALL_RADIUS, -height, POCKET_CENTER_Z[5] + 2*BALL_RADIUS);

	// inbetween
	glVertex3d(POCKET_CENTER_X[5] - BALL_RADIUS, -height, POCKET_CENTER_Z[5] + 4*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[3] - BALL_RADIUS - .06, -height, POCKET_CENTER_Z[3] - 2*BALL_RADIUS);

	// MID RIGHT
	glVertex3d(POCKET_CENTER_X[3] - BALL_RADIUS - .06, -height, POCKET_CENTER_Z[3] + 2*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[3] + BALL_RADIUS, -height, POCKET_CENTER_Z[3] + 2*BALL_RADIUS);
	
	glVertex3d(POCKET_CENTER_X[3] - BALL_RADIUS - .06, -height, POCKET_CENTER_Z[3] - 2*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[3] + BALL_RADIUS, -height, POCKET_CENTER_Z[3] - 2*BALL_RADIUS);

	// inbetween
	glVertex3d(POCKET_CENTER_X[3] - BALL_RADIUS - .06, -height, POCKET_CENTER_Z[3] + 2*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[1] - BALL_RADIUS, -height, POCKET_CENTER_Z[1] - 4*BALL_RADIUS);

	// TOP RIGHT
	glVertex3d(POCKET_CENTER_X[1] - 4*BALL_RADIUS, -height, POCKET_CENTER_Z[1] - BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[1] - 2*BALL_RADIUS, -height, POCKET_CENTER_Z[1] + BALL_RADIUS);
	
	glVertex3d(POCKET_CENTER_X[1] - BALL_RADIUS, -height, POCKET_CENTER_Z[1] - 4*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[1] + BALL_RADIUS, -height, POCKET_CENTER_Z[1] - 2*BALL_RADIUS);

	// inbetween
	glVertex3d(POCKET_CENTER_X[1] - 4*BALL_RADIUS, -height, POCKET_CENTER_Z[1] - BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[4] + 4*BALL_RADIUS, -height, POCKET_CENTER_Z[4] - BALL_RADIUS);

	// TOP LEFT
	glVertex3d(POCKET_CENTER_X[4] + 4*BALL_RADIUS, -height, POCKET_CENTER_Z[4] - BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[4] + 2*BALL_RADIUS, -height, POCKET_CENTER_Z[4] + BALL_RADIUS);
	
	glVertex3d(POCKET_CENTER_X[4] + BALL_RADIUS, -height, POCKET_CENTER_Z[4] - 4*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[4] - BALL_RADIUS, -height, POCKET_CENTER_Z[4] - 2*BALL_RADIUS);

	// inbetween
	glVertex3d(POCKET_CENTER_X[4] + BALL_RADIUS, -height, POCKET_CENTER_Z[4] - 4*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[2] + BALL_RADIUS + .06, -height, POCKET_CENTER_Z[2] + 2*BALL_RADIUS);

	// MID LEFT
	glVertex3d(POCKET_CENTER_X[2] + BALL_RADIUS + .06, -height, POCKET_CENTER_Z[2] + 2*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[2] - BALL_RADIUS, -height, POCKET_CENTER_Z[2] + 2*BALL_RADIUS);
	
	glVertex3d(POCKET_CENTER_X[2] + BALL_RADIUS + .06, -height, POCKET_CENTER_Z[2] - 2*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[2] - BALL_RADIUS, -height, POCKET_CENTER_Z[2] - 2*BALL_RADIUS);

	// final inbetween
	glVertex3d(POCKET_CENTER_X[2] + BALL_RADIUS + .06, -height, POCKET_CENTER_Z[2] - 2*BALL_RADIUS);
	glVertex3d(POCKET_CENTER_X[0] + BALL_RADIUS, -height, POCKET_CENTER_Z[0] + 4*BALL_RADIUS);

	glEnd();
	glEnable(GL_LIGHTING);

	glutSwapBuffers();

}

// ChangeView will change the view mode to the given int mode.
void ChangeView(int mode) {

	if (mode == TABLE_VIEW) {
		TableViewRotateX = 0.0;
		TableViewRotateY = 0.0;
	}
	ViewMode = mode;

}

// Our reshape function
void reshapeFunc(int w , int h)
{

	if (h <= 0) return;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(40.0 , ((double) w) / ((double) h) , 0.01 , 30.0);
	glViewport(0 , 0 , w , h);

}

// Call this function to initialize the state of openGL.
void InitGL()
{

	glEnable(GL_LIGHTING);		// Lighting defaults to off. It must be enabled.
	glEnable(GL_LIGHT0);		// Light0 is predefined to white. It must be enabled.
	glEnable(GL_DEPTH_TEST);	// Enable depth testing.

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);
	glDepthFunc(GL_LEQUAL);
	glCullFace(GL_FRONT);
	glFrontFace(GL_CCW);
	//glEnable(GL_CULL_FACE);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GLUT_MULTISAMPLE);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.06f);

	// Initialize DevIL
	ilInit();
	iluInit();
	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_SET);
	ilEnable(IL_TYPE_SET);
	ilTypeFunc(IL_UNSIGNED_BYTE);
	ilEnable(IL_FORMAT_SET);
	ilFormatFunc(IL_RGBA);
	ilutRenderer(ILUT_OPENGL);
}

// This is our mouse function. Used when clicking a mouse
// button.
void MouseFunction(int button, int state, int x, int y) {

	switch(button) {

	case GLUT_LEFT_BUTTON:

		if (state == GLUT_DOWN) {
			LeftButtonDown = true;
			MouseX = x;
			MouseY = y;
			if (ViewMode == CUE_BALL_VIEW) {
				StartShot((double)y);
			}

		} else if (state == GLUT_UP) {
			LeftButtonDown = false;
			if (ViewMode == CUE_BALL_VIEW) {
				CancelShot();
			}
		}
		break;

	case GLUT_RIGHT_BUTTON:

		if (state == GLUT_DOWN) {
			RightButtonDown = true;
			MouseX = x;
			MouseY = y;
		} else if (state == GLUT_UP) {
			if (GameState == DRAWING_CUE) CancelShot();
			RightButtonDown = false;
		}
		break;

	}

	return;

}

// This is our mouse motion function.  Used when the mouse is moved.
void MotionFunction(int x, int y) {

	//printf("MouseX: %d || ShotStartZ: %f\n", y, ShotStartZ);

	if (ViewMode == TABLE_VIEW) {

		TableViewRotateX += MouseX - x;
		TableViewRotateY += MouseY - y;

	} else if (ViewMode == CUE_BALL_VIEW) {

		if (RightButtonDown == true) {
			if (GameState != DRAWING_CUE) {
				RotateX += ( ((double)MouseX - (double)x) / 20.0);
			}
		}

		if (LeftButtonDown == true) {
			CueRotateY += ((double) MouseX - (double) x) / 20;
			CueOffsetZ += -0.01 * (double)(MouseY - y);
			if (y < (int)(ShotStartZ - 5) && GameState == DRAWING_CUE) {
				FinishShot((double)y);
			}
		}

	} else if (ViewMode == FROZEN) {

		RotateX += MouseX - x;

	}

	MouseX = x;
	MouseY = y;
	LastMouseX = MouseX;
	LastMouseY = MouseY;


}

// This is our special function.  Basically just the arrow
// keys are used for placing scratched balls.
void specialFunc(int key, int x, int y) {

	switch (key) {

	case GLUT_KEY_UP:

		if (GameState == SCRATCH) {
			cueBall->position.z += 0.05;
			if (cueBall->position.z > (4.5 + 4*BALL_RADIUS)) cueBall->position.z = (4.5 + 4*BALL_RADIUS);
			if (Breaking == true) {
				if (cueBall->position.z > -3.0) cueBall->position.z = -3.0;
			}
		}
		break;

	case GLUT_KEY_DOWN:

		if (GameState == SCRATCH) {
			cueBall->position.z -= 0.05;
			if (cueBall->position.z < (-4.5 - 4*BALL_RADIUS)) cueBall->position.z = (-4.5 - 4*BALL_RADIUS);
		}
		break;

	case GLUT_KEY_RIGHT:

		if (GameState == SCRATCH) {
			cueBall->position.x += 0.05;
			if (cueBall->position.x > (2.0 + 4*BALL_RADIUS)) cueBall->position.x = (2.0 + 4*BALL_RADIUS);
		}
		break;

	case GLUT_KEY_LEFT:

		if (GameState == SCRATCH) {
			cueBall->position.x -= 0.05;
			if (cueBall->position.x < (-2.0 - 4*BALL_RADIUS)) cueBall->position.x = (-2.0 - 4*BALL_RADIUS);
		}
		break;

	}

}

// This is my keyboard function.
// Hitting the 'Q' key will change your view to the overhead view
// Hitting the SPACE BAR will bring you in for your next shot if 
// all the balls have stopped moving.
void keyboardFunc(unsigned char key, int x, int y) {

	switch (key) {

	case 'q':

		if (GameState != TAKING_SHOT && GameState != DRAWING_CUE) ChangeView(TABLE_VIEW);
		break;

	case ' ':

		if (GameState == WAITING || GameState == SCRATCH) {
			if (CheckSingleCollision(0) == false) {
				InitShot();
				ChangeView(CUE_BALL_VIEW);
				GameState = TAKING_SHOT;
				Breaking = false;
			}
		}

		if (GameState == PLAYER_1_WINS || GameState == PLAYER_2_WINS) {

			NewGame();


		}

		break;

	}

}

// DrawSphere calls a display list for the given sphere.
void DrawSphere(GLuint list)
{
	glCallList(list);
}

// CreateSphere loads and vertex arrays a sphere with triangle strips.
// It takes input radius, lats and longs for stacks and slices.
// It also takes a vertex array list number.
void CreateSphere(double r, int lats, int longs, GLuint list)
{
	int i, j;
	glNewList(list, GL_COMPILE);
	for(i = 0; i <= lats; i++)
	{
		double lat0 = PI * (-0.5 + (double) (i - 1) / lats);
		double z0  = sin(lat0) * r;
		double zr0 =  cos(lat0) * r;

		double lat1 = PI * (-0.5 + (double) i / lats);
		double z1 = sin(lat1) * r;
		double zr1 = cos(lat1) * r;
		
		glBegin(GL_TRIANGLE_STRIP);
		for(j = 0; j <= longs; j++)
		{
			double lng = 2 * PI * (double) (j - 1) / longs;
			double x = cos(lng);
			double y = sin(lng);

			glNormal3f((GLfloat) (x * zr0), (GLfloat) (y * zr0), (GLfloat) (z0));
			glVertex3f((GLfloat) (x * zr0), (GLfloat) (y * zr0), (GLfloat) (z0));

			GLfloat tx1 = (GLfloat)(atan2((GLfloat)(x * zr0), (GLfloat) (z0)) / (2. * PI) + 0.5);
			GLfloat ty1 = (GLfloat)(asin((GLfloat)(y * zr0)) / PI + 0.5);
			glTexCoord2f(tx1, ty1);

			glNormal3f((GLfloat) (x * zr1), (GLfloat) (y * zr1), (GLfloat) (z1));
			glVertex3f((GLfloat) (x * zr1), (GLfloat) (y * zr1), (GLfloat) (z1));

			GLfloat tx = (GLfloat) (atan2((GLfloat)(x * zr1), (GLfloat) (z1)) / (2.0 * PI) + 0.5);
			GLfloat ty = (GLfloat) (asin((GLfloat)(y * zr1)) / PI + 0.5);
			if(tx < 0.75 && tx1 > 0.75)
			tx += 1.0;
			else if(tx > 0.75 && tx1 < 0.75)
			tx -= 1.0;
			glTexCoord2f(tx, ty);
		}
		glEnd();
	}
	glEndList();
}

int	main(int argc , char * argv[])
{

	// Get the current working dir.
	_getcwd(cwd, 250);

	// Start my sound engine...
	engine = createIrrKlangDevice();
	if (!engine)
      return 0;

	// Store sounds to memory
	engine->addSoundSourceFromMemory(bank1, bank1_size, "bank1.ogg");
	engine->addSoundSourceFromMemory(bank2, bank2_size, "bank2.ogg");
	engine->addSoundSourceFromMemory(bank3, bank3_size, "bank3.ogg");
	engine->addSoundSourceFromMemory(ballhitlight, ballhitlight_size, "ballhitlight.ogg");
	engine->addSoundSourceFromMemory(ballhitmedium, ballhitmedium_size, "ballhitmedium.ogg");
	engine->addSoundSourceFromMemory(ballhithard, ballhithard_size, "ballhithard.ogg");
	engine->addSoundSourceFromMemory(cuehit, cuehit_size, "cuehit.ogg");
	engine->addSoundSourceFromMemory(pocket, pocket_size, "pocket.ogg");

	glutInit(&argc , argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(512 , 512);
	glutInitWindowPosition(0 , 0);
	glutCreateWindow("Billiards");
	glutMouseFunc(MouseFunction);
	glutMotionFunc(MotionFunction);
	glutDisplayFunc(display);
	glutReshapeFunc(reshapeFunc);
	glutSpecialFunc(specialFunc);
	glutKeyboardFunc(keyboardFunc);

	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE , 0.0);
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER , 1.0);

	// Some calculations of the X and Z bumpers
	xBumpMax = (X_SIZE / 2) - (BALL_RADIUS / 2);
	xBumpMin = -(X_SIZE / 2) + (BALL_RADIUS / 2);
	zBumpMax = (Z_SIZE / 2) - (BALL_RADIUS / 2);
	zBumpMin = -(Z_SIZE / 2) + (BALL_RADIUS / 2);

	InitGL();
	InitBalls();

	UP_AXIS->x = 0.0;
	UP_AXIS->y = 1.0;
	UP_AXIS->z = 0.0;

	PocketList = glGenLists(1);

	LoadPockets();
	LoadWalls();

	CueTipList = glGenLists(1);
	CueShaft1List = glGenLists(1);
	CueShaft2List = glGenLists(1);

	buildCueTip(100);
	buildCueShaft1(100);
	buildCueShaft2(100);

	ShotString = (char*)malloc(sizeof(char) * 50);
	NextShot();
	glutTimerFunc(16, timerFunction, 1);
	glutMainLoop();

	return 0;
}
