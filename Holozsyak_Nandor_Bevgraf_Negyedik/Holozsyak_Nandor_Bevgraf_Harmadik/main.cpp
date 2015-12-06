#include <GLFW\glfw3.h>                           // (or others, depending on the system in use)
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <Torusz.h>
#include <vector>
#include <thread>
#include <chrono>
/*======================================*/

/**
* 3D vektor
*/
typedef struct { GLdouble x, y, z; } VECTOR3;

/**
* 4D vektor
*/
typedef struct { GLdouble x, y, z, w; } VECTOR4;

/**
* 4x4 mátrix, sorfolytonosan tárolva
*/
typedef GLdouble MATRIX4[4][4];

/*======================================*/

/**
* visszadja az (x,y,z) vektort
*/
VECTOR3 initVector3(GLdouble x, GLdouble y, GLdouble z) {
	VECTOR3 P;

	P.x = x;
	P.y = y;
	P.z = z;

	return P;
}

/**
* visszadja az (x,y,z,w) vektort
*/
VECTOR4 initVector4(GLdouble x, GLdouble y, GLdouble z, GLdouble w) {
	VECTOR4 P;

	P.x = x;
	P.y = y;
	P.z = z;
	P.w = w;

	return P;
}

VECTOR3 convertToInhomogen(VECTOR4 vector) {
	return initVector3(vector.x / vector.w, vector.y / vector.w, vector.z / vector.w);
}

VECTOR4 convertToHomogen(VECTOR3 vector) {
	return initVector4(vector.x, vector.y, vector.z, 1.0);
}

GLdouble getVectorLength(VECTOR3 vec) {
	return sqrt(pow(vec.x, 2) + pow(vec.y, 2) + pow(vec.z, 2));
}

VECTOR3 normalizeVector(VECTOR3 vector) {
	GLdouble length = getVectorLength(vector);
	return initVector3(vector.x / length, vector.y / length, vector.z / length);
}

GLdouble dotProduct(VECTOR3 a, VECTOR3 b) {
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

VECTOR3 crossProduct(VECTOR3 a, VECTOR3 b) {

	return initVector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

/**
* feltölti az A mátrixot az egységmátrixszal
*/
void initIdentityMatrix(MATRIX4 A)
{
	int i, j;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			A[i][j] = 0.0f;

	for (i = 0; i < 4; i++)
		A[i][i] = 1.0f;
}

/**
* feltölti az A mátrixot a (0,0,s) középpontú centrális vetítés mátrixszával
*/
void initPersProjMatrix(MATRIX4 A, GLdouble s)
{
	initIdentityMatrix(A);

	A[2][2] = 0.0f;
	A[3][2] = -1.0f / s;
}

/**
* feltölti az A mátrixot a (wx,wy):(wx+ww,wy+wh) ablakból (vx,vy):(vx+vw,vy+vh) nézetbe
* történõ transzformáció mátrixszával
*/
void initWtvMatrix(MATRIX4 A, GLdouble wx, GLdouble wy, GLdouble ww, GLdouble wh,
	GLdouble vx, GLdouble vy, GLdouble vw, GLdouble vh)
{
	initIdentityMatrix(A);

	A[0][0] = vw / ww;
	A[1][1] = vh / wh;
	A[0][3] = -wx*(vw / ww) + vx;
	A[1][3] = -wy*(vh / wh) + vy;
}

void initViewMatrix(MATRIX4 A, VECTOR3 eye, VECTOR3 center, VECTOR3 up) {
	initIdentityMatrix(A);

	// a kamerabol az iranyt meghatarozo pontba mutato vektor
	// center az a pont,amely fele a kamerat tartjuk, eye a kamera helyét adja
	VECTOR3 centerMinusEye = initVector3(center.x - eye.x, center.y - eye.y, center.z - eye.z);

	// a fenti vektor -1 szeresének egyseg hosszra normáltja, ebbol lesz a kamera rendszerének z-tengelye
	VECTOR3 f = normalizeVector(initVector3(-centerMinusEye.x, -centerMinusEye.y, -centerMinusEye.z));

	// az up vektor es a leendo z tengelyirany vektorialis szorzata adja a kamera x-tengelyiranyat
	VECTOR3 s = normalizeVector(crossProduct(up, f));

	// a kamera y tengelyiranya a mar elkeszult z irany es x irany vektorialis szorzataként jön ki
	VECTOR3 u = crossProduct(f, s);

	A[0][0] = s.x;
	A[0][1] = s.y;
	A[0][2] = s.z;
	A[1][0] = u.x;
	A[1][1] = u.y;
	A[1][2] = u.z;
	A[2][0] = f.x;
	A[2][1] = f.y;
	A[2][2] = f.z;
	A[0][3] = -dotProduct(s, eye);
	A[1][3] = -dotProduct(u, eye);
	A[2][3] = -dotProduct(f, eye);
}


/**
* visszaadja az A mátrix és a v vektor szorzatát, A*v-t
*/
VECTOR4 mulMatrixVector(MATRIX4 A, VECTOR4 v) {
	return initVector4(
		A[0][0] * v.x + A[0][1] * v.y + A[0][2] * v.z + A[0][3] * v.w,
		A[1][0] * v.x + A[1][1] * v.y + A[1][2] * v.z + A[1][3] * v.w,
		A[2][0] * v.x + A[2][1] * v.y + A[2][2] * v.z + A[2][3] * v.w,
		A[3][0] * v.x + A[3][1] * v.y + A[3][2] * v.z + A[3][3] * v.w);
}

/**
* feltölti a C mátrixot az A és B mátrixok szorzatával, A*B-vel
*/
void mulMatrices(MATRIX4 A, MATRIX4 B, MATRIX4 C) {
	int i, j, k;

	GLdouble sum;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
			sum = 0;
			for (k = 0; k < 4; k++)
				sum = sum + A[i][k] * B[k][j];
			C[i][j] = sum;
		}
}

/*======================================*/

/**
* ablak mérete
*/
GLdouble winWidth = 800.0f, winHeight = 600.0f;


MATRIX4 view;
double **cam;

/**
* merõleges és centrális vetítések mátrixai
*/
MATRIX4 Vc;

/**
* Wtv mátrixok
*/
MATRIX4 Wc;

/**
* a fenti mátrixokból elõállított két transzformációs mátrix
*/
MATRIX4 Tcx, Tcy, Tcz, TC[4];

/**
* segédmátrix
*/
MATRIX4 Tmp, TempR, TempE;

/*
* forgatási mátrixok
*/
MATRIX4 Rx, Ry, Rz;
/*
* eltolási mátrix
*/
MATRIX4 E;
/*
* skálázás mátrixa
*/
MATRIX4 S;

MATRIX4 Sik;

GLdouble alpha = 0.0f, deltaAlpha = 3.14f / 180.0f;

/**
* nézet koordinátái
*/
//GLdouble cX = (winWidth - winHeight) / 2.0f + 150, cY = 0, cW = winHeight / 2, cH = winHeight / 2 - 100;
GLdouble cX = 100, cY = 0, cW = 600, cH = 500;
/**
* centrális vetítés középpontjának Z koordinátája
*/
GLdouble center = 6.0f;


void initEltolasMatrix(double x, double y, double z) {
	initIdentityMatrix(E);

	/*
	<<1<< 0 << 0 << x
	<< 0 << 1 << 0 << y
	<< 0 << 0 << 1 << z
	<< 0 << 0 << 0 << 1;
	*/

	E[0][3] = x;
	E[1][3] = y;
	E[2][3] = z;
}
void initSkalazasMatrix(double x, double y, double z) {
	//initIdentityMatrix(S);
	S[0][0] = x;
	S[1][1] = y;
	S[2][2] = z;
	S[3][3] = 1;
}
void initRotationMatrixes(double alpha) {
	initIdentityMatrix(Rx);
	initIdentityMatrix(Ry);
	initIdentityMatrix(Rz);

	/*X:
	<< 1 << 0 << 0 << 0
	<< 0 << std::cos(alpha) << -std::sin(alpha) << 0
	<< 0 << std::sin(alpha) << std::cos(alpha) << 0
	<< 0 << 0 << 0 << 1;
	*/
	Rx[1][1] = std::cos(alpha);
	Rx[1][2] = -std::sin(alpha);
	Rx[2][1] = std::sin(alpha);
	Rx[2][2] = std::cos(alpha);

	/*
	Y:
	<< std::cos(alpha) << 0 << std::sin(alpha) << 0
	<< 0 << 1 << 0 << 0
	<< -std::sin(alpha) << 0 << std::cos(alpha) << 0
	<< 0 << 0 << 0 << 1;
	*/

	Ry[0][0] = std::cos(alpha);
	Ry[0][2] = std::sin(alpha);
	Ry[2][0] = -std::sin(alpha);
	Ry[2][2] = std::cos(alpha);

	/*
	Z:
	<< std::cos(alpha) << -std::sin(alpha) << 0 << 0
	<< std::sin(alpha) << std::cos(alpha) << 0 << 0
	<< 0 << 0 << 1 << 0
	<< 0 << 0 << 0 << 1;
	*/
	Rz[0][0] = std::cos(alpha);
	Rz[0][1] = -std::sin(alpha);
	Rz[1][0] = std::sin(alpha);
	Rz[1][1] = std::cos(alpha);

}

double deltaA = 0.01;
/**
* elõállítja a szükséges mátrixokat
*/
double tA = 5, tB = 18, tC = 5;
double sarok = 12;
void initTransformations()
{

	// vetítési mátrixok
	initPersProjMatrix(Vc, center);

	// Wtv mátrixok
	initWtvMatrix(Wc, -0.5f, -0.5f, 1.0f, 1.0f, cX, cY, cW, cH);
	//Sík mátrixa
	//T = W×C×K
	mulMatrices(Vc, view, Tmp);
	mulMatrices(Wc, Tmp, Sik);
	//M = W×C×K×L×R, tórusz mátrixa elõbb forgatom, eltolom, kamera, centrális vetítés és w2v
	//Tóruszok
	initRotationMatrixes(deltaA);
	initEltolasMatrix(0, 5, 0);

	// egyik tórusz - X tengely
	mulMatrices(E, Rx, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, Tcx);

}

/*======================================*/


/**
* egységkocka csúcsai
*/
VECTOR4 identityCube[8] =
{
	initVector4(0.5f, 0.5f,-0.5f, 1.0f),
	initVector4(-0.5f, 0.5f,-0.5f, 1.0f),
	initVector4(-0.5f,-0.5f,-0.5f, 1.0f),
	initVector4(0.5f,-0.5f,-0.5f, 1.0f),
	initVector4(0.5f, 0.5f, 0.5f, 1.0f),
	initVector4(-0.5f, 0.5f, 0.5f, 1.0f),
	initVector4(-0.5f,-0.5f, 0.5f, 1.0f),
	initVector4(0.5f,-0.5f, 0.5f, 1.0f),
};

/**
* egységkocka lapjainak indexei
*/
GLuint faces[24] =
{
	0,1,2,3, //alsó
	0,1,5,4, //jobb
	1,2,6,5, //hátsó
	2,3,7,6, //bal
	3,0,4,7, //els?
	4,5,6,7, //fels?
};

VECTOR3 eye, up, centerVec;
double helyzet = 0.0;
Torusz A, B, C;
double R = 37;

void drawCubes(VECTOR3 color, MATRIX4 T)
{
	int i, j, id = 0;

	VECTOR4 ph, pt;
	VECTOR3 pih;

	// beállítja a kocka éleinek vastagságát
	glLineWidth(2.0f);

	// beállítja a kocka színét
	glColor3f(color.x, color.y, color.z);

	// végigmegyünk a lapokon
	for (i = 0; i < 6; i++)
	{
		glBegin(GL_LINE_LOOP);

		// végigmegyünk a csúcsokon
		for (j = 0; j < 4; j++)
		{
			// homogén koordinátás alakra hozzuk a csúcsot
			ph = initVector4(identityCube[faces[id]].x, identityCube[faces[id]].y, identityCube[faces[id]].z, 1.0f);

			// alkalmazzuk a transzformációt
			pt = mulMatrixVector(T, ph);

			// visszatérünk inhomogén alakra
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

			// kirajzoljuk a csúcsot
			glVertex2f(pih.x, pih.y);

			id++;
		}

		glEnd();
	}
}



void drawGrid(VECTOR3 color, MATRIX4 T)
{
	int i, j, id = 0;

	VECTOR4 ph, pt;
	VECTOR3 pih;


	glLineWidth(1.0f);


	glColor3f(color.x, color.y, color.z);

	glBegin(GL_LINES);
	for (double i = -sarok; i <= sarok; i += sarok / 12) {
		ph = initVector4(-sarok, 0, i, 1.0f);

		// alkalmazzuk a transzformációt
		pt = mulMatrixVector(T, ph);

		// visszatérünk inhomogén alakra
		pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

		// kirajzoljuk a csúcsot
		glVertex2f(pih.x, pih.y);

		ph = initVector4(sarok, 0, i, 1.0f);

		// alkalmazzuk a transzformációt
		pt = mulMatrixVector(T, ph);

		// visszatérünk inhomogén alakra
		pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

		// kirajzoljuk a csúcsot
		glVertex2f(pih.x, pih.y);

	}
	glEnd();
	glBegin(GL_LINES);
	for (double i = -sarok; i <= sarok; i += sarok / 12) {
		ph = initVector4(i, 0, -sarok, 1.0f);

		// alkalmazzuk a transzformációt
		pt = mulMatrixVector(T, ph);

		// visszatérünk inhomogén alakra
		pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

		// kirajzoljuk a csúcsot
		glVertex2f(pih.x, pih.y);

		ph = initVector4(i, 0, sarok, 1.0f);

		// alkalmazzuk a transzformációt
		pt = mulMatrixVector(T, ph);

		// visszatérünk inhomogén alakra
		pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

		// kirajzoljuk a csúcsot
		glVertex2f(pih.x, pih.y);
	}
	glEnd();
}



void init()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0f, winWidth, 0.0f, winHeight, 0.0f, 1.0f);

	eye = initVector3(0.0, 0.0, R); //megadja a kamera pozícióját (Ez legyen most a z tengely pozitív felén)
	eye.x = R*cos(alpha);
	eye.z = R*sin(alpha);
	centerVec = initVector3(0.0, 0.0, 0.0); //megadja, hogy merre néz a kamera (Ez legyen most az origó)
	up = initVector3(0.0, 1.0, 0.0); //megdja, hogy merre van a felfele irány (Ez legyen most az y tengely)

	initViewMatrix(view, eye, centerVec, up);

	initTransformations();
}
double PI = 3.14159265358979323846;
void drawTorusz(VECTOR3 color, MATRIX4 T) {

	double B[4][1];
	double theta, fi;
	float c = 8, a = 1;
	VECTOR4 ph, pt;
	VECTOR3 pih;
	glLineWidth(1.0);
	glPointSize(4.0);
	glColor3f(color.x, color.y, color.z);
	double delta = PI / 12;
	for (theta = 0; theta <= 2 * PI; theta += delta) {
		glBegin(GL_LINE_LOOP);
		for (fi = 0; fi <= 2 * PI; fi += delta) {
			B[0][0] = (c + a*cos(theta))*cos(fi);
			B[1][0] = a*sin(theta);
			B[2][0] = (c + a*cos(theta))*sin(fi);
			B[3][0] = 1.0;
			// homogén koordinátás alakra hozzuk a csúcsot
			ph = initVector4(B[0][0], B[1][0], B[2][0], 1.0);

			// alkalmazzuk a transzformációt
			pt = mulMatrixVector(T, ph);

			// visszatérünk inhomogén alakra
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);
			
			// kirajzoljuk a csúcsot
			glVertex2f(pih.x, pih.y);
			
		}
		glEnd();
	}


	glColor3f(1.0, 0, 0);
	for (fi = 0; fi <= 2 * PI; fi += delta) {
		glBegin(GL_LINE_LOOP);
		for (theta = 0; theta <= 2 * PI; theta += delta) {
			B[0][0] = (c + a*cos(theta))*cos(fi);
			B[1][0] = a*sin(theta);
			B[2][0] = (c + a*cos(theta))*sin(fi);
			B[3][0] = 1.0;
			// homogén koordinátás alakra hozzuk a csúcsot
			ph = initVector4(B[0][0], B[1][0], B[2][0], 1.0);

			// alkalmazzuk a transzformációt
			pt = mulMatrixVector(T, ph);

			// visszatérünk inhomogén alakra
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

			// kirajzoljuk a csúcsot
			glVertex2f(pih.x, pih.y);

		}
		glEnd();
	}
}
/**
* egy laphoz tartozó csúcsok indexei
*/
typedef struct { GLint v[4]; GLdouble cAVG; } FACE;

/**
* visszadja a (v0,v1,v2,v3) indexekhez tartozó lapot
*/
FACE initFace(GLuint v0, GLuint v1, GLuint v2, GLuint v3, GLdouble avg) {
	FACE f;
	f.v[0] = v0;
	f.v[1] = v1;
	f.v[2] = v2;
	f.v[3] = v3;
	f.cAVG = avg;
	return f;
}
//FACE faceTorusz[13 * 13];

int comparePointsbyAtlag(const void *a, const void *b) {

	if (((*(FACE*)a).cAVG < ((*(FACE*)b).cAVG))) return -1;
	if (((*(FACE*)a).cAVG == ((*(FACE*)b).cAVG))) return  0;
	if (((*(FACE*)a).cAVG > ((*(FACE*)b).cAVG))) return  1;


}
int comparePointsbyZ(const void *a, const void *b) {

	/*if (((*(FACE*)a).cAVG < ((*(FACE*)b).cAVG))) return -1;
	if (((*(FACE*)a).cAVG == ((*(FACE*)b).cAVG))) return  0;
	if (((*(FACE*)a).cAVG > ((*(FACE*)b).cAVG))) return  1;*/

	if (((*(VECTOR3*)a).z < ((*(VECTOR3*)b).z))) return -1;
	if (((*(VECTOR3*)a).z == ((*(VECTOR3*)b).z))) return  0;
	if (((*(VECTOR3*)a).z > ((*(VECTOR3*)b).z))) return  1;
}

int comparePointsbyZ4D(const void *a, const void *b) {

	/*if (((*(FACE*)a).cAVG < ((*(FACE*)b).cAVG))) return -1;
	if (((*(FACE*)a).cAVG == ((*(FACE*)b).cAVG))) return  0;
	if (((*(FACE*)a).cAVG > ((*(FACE*)b).cAVG))) return  1;*/

	if (((*(VECTOR4*)a).z < ((*(VECTOR4*)b).z))) return -1;
	if (((*(VECTOR4*)a).z == ((*(VECTOR4*)b).z))) return  0;
	if (((*(VECTOR4*)a).z >((*(VECTOR4*)b).z))) return  1;
}
/**
* visszaadja az 'a' vektorból a 'b' vektorba mutató vektort,
* azaz a 'b-a' vektort
*/
VECTOR3 vecSub(VECTOR3 a, VECTOR3 b) {
	return initVector3(
		b.x - a.x,
		b.y - a.y,
		b.z - a.z);
}


/**
* visszaadja az 'a' vektor hosszát
*/
float length(VECTOR3 a) {
	return sqrt(dotProduct(a, a));
}
/**
* visszaadja az 'a' vektor normalizáltját
*/
VECTOR3 normalize(VECTOR3 a) {
	float len = length(a);

	return initVector3(a.x / len, a.y / len, a.z / len);
}


std::vector<VECTOR4> pontok4d;
std::vector<VECTOR3> pontok;
std::vector<VECTOR3> transformedPontok;
std::vector<FACE> faceTorusz;
void drawToruszWithPlanes(VECTOR3 color, MATRIX4 T, MATRIX4 CAM) {
	Matrix<double> vectorPl(12, 12);
	double B[4][1];
	double theta = PI / 6, fi;
	float c = 3,a = 2;
	VECTOR4 ph, pt;
	VECTOR3 pih;
	int i = 0, j = 0;
	glColor3f(color.x, color.y, color.z);
	
	//2PI / (PI/6) = 12
	//2PI / (Pi/12) = 24;
	//2PI / (Pi/8) = 16;
	//2PI / (Pi/4) = 8;
	double delta = PI / 8;

	int dI = (int)((2*PI)/(delta)), dJ = (int)((2 * PI) / (delta));
	std::cout <<(double) (2 * PI) / (delta) <<"~" << dI << std::endl;
	int K = 1;
	
	for (fi = 0; fi < 2 * PI + delta; fi += delta) {
		for (theta = 0; theta < 2 * PI; theta += delta) {
			B[0][0] = (c + a*cos(theta))*cos(fi);
			B[1][0] = a*sin(theta);
			B[2][0] = (c + a*cos(theta))*sin(fi);
			B[3][0] = 1.0;
			// homogén koordinátás alakra hozzuk a csúcsot
			ph = initVector4(B[0][0], B[1][0], B[2][0], 1.0);
			//std::cout << "ph:" << ph.x << " - " << ph.y << " - " << ph.z << std::endl;
			
			pt = mulMatrixVector(CAM, ph);
			//std::cout << "Cam pt:" << pt.x << " - " << pt.y << " - " << pt.z << std::endl;
			pontok4d.push_back(pt);
			pt = mulMatrixVector(T, ph);
			//std::cout << "pt:" << pt.x << " - " << pt.y << " - " << pt.z << std::endl;
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);
			pontok.push_back(pih);
		}

	}

	/*for (fi = 0; fi < 2 * PI + delta; fi += delta) {
		//std::cout << K++ << ". 12 pont" << std::endl;
		for (theta = 0; theta < 2 * PI; theta += delta) {
			B[0][0] = (c + a*cos(theta))*cos(fi);
			B[1][0] = a*sin(theta);
			B[2][0] = (c + a*cos(theta))*sin(fi);
			B[3][0] = 1.0;
			// homogén koordinátás alakra hozzuk a csúcsot
			ph = initVector4(B[0][0], B[1][0], B[2][0], 1.0);
			//std::cout << ph.x << " - " << ph.y << " - " << ph.z << " - " << ph.w << std::endl;
			// alkalmazzuk a transzformációt
			pt = mulMatrixVector(T, ph);
			//std::cout << "z/w:" << pt.z / pt.w << std::endl;
			// visszatérünk inhomogén alakra
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);
			pontok.push_back(pih);
		//	std::cout << pih.x << " - " << pih.y << " - " << pih.z << std::endl;
			glColor3f(0, theta, 0.5);
			glBegin(GL_POINTS);
			glVertex2d(pih.x, pih.y);
			glEnd();
		}

	}*/

	//qsort(&pontok4d[0], pontok4d.size(), sizeof(VECTOR4), comparePointsbyZ4D);
	
	/*for (int l = 0; l < pontok4d.size(); l++) {
		ph = initVector4(pontok4d[l].x, pontok4d[l].y, pontok4d[l].z, 1.0);
		pt = mulMatrixVector(T, ph);
		pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);
		pontok.push_back(pih);
	}*/

	double cAVG = 0, cAVG1 = 0, cAVG2 = 0, cAVG3 = 0, cAVG4 = 0;

	for (i = 0; i < dI; i++) {
		for (int j = 0; j < dJ; j++) {
			FACE f;
			if (j == dJ-1) {
				cAVG = 
					(pontok4d[i * dI + j].z +
						pontok4d[i * dI + j + dJ].z +
						pontok4d[i * dI + j + dJ + 1 - dI].z +
						pontok4d[i * dI + j + 1 - dI].z) / 4;
				
				f = initFace((GLuint)(i * dI + j), (GLuint)(i * dI + j + dJ), (GLuint)(i * dI + j + dJ + 1 - dI), (GLuint)(i * dI + j + 1 - dI), cAVG);
			}
			else {
				cAVG  =
					(pontok4d[i * dI + j].z +
						pontok4d[i * dI + j + dJ].z +
						pontok4d[i * dI + j + dJ + 1].z +
						pontok4d[i * dI + j + 1].z) / 4;
				f = initFace((GLint)(i * dI + j), (GLint)(i * dI + j + dJ), (GLint)(i * dI + j + dJ + 1), (GLint)(i * dI + j + 1), cAVG);
			}

			faceTorusz.push_back(f);
		//	std::cout << "actual F:" << f.v[0] << " - " << f.v[1] << " - " << f.v[2] << " - " << f.v[3] << " cAVG:" << f.cAVG << std::endl;
		}
	}
	qsort(&faceTorusz[0], faceTorusz.size(), sizeof(FACE), comparePointsbyAtlag);
	
	for (int i = 0; i < faceTorusz.size(); i++) {
		//std::cout << "after sort :" << faceTorusz[i].v[0] << " - " << faceTorusz[i].v[1] << " - " << faceTorusz[i].v[2] << " - " << faceTorusz[i].v[3] << " cAVG:" << faceTorusz[i].cAVG << std::endl;
		// meghatározzuk a lap megfelel?i éleit
		VECTOR3 edges[2] =
		{
			vecSub(convertToInhomogen(pontok4d[faceTorusz[i].v[0]]),
			convertToInhomogen(pontok4d[faceTorusz[i].v[1]])),

			vecSub(convertToInhomogen(pontok4d[faceTorusz[i].v[0]]),
			convertToInhomogen(pontok4d[faceTorusz[i].v[2]]))

			/*vecSub(pontok[faceTorusz[i].v[0]],
				pontok[faceTorusz[i].v[1]]),

			vecSub(pontok[faceTorusz[i].v[0]],
				pontok[faceTorusz[i].v[2]])*/
		};

		// kiszámítjuk ezekb?l a a lap normálvektorát
		VECTOR3 normal = normalize(crossProduct(edges[0], edges[1]));

		// camerába mutató vektort
		VECTOR3 toCamera;
		//toCamera = normalize(vecSub(pontok[faceTorusz[i].v[0]], initVector3(0.0f, 0.0f, center)));
		toCamera = normalize(vecSub(convertToInhomogen(pontok4d[faceTorusz[i].v[0]]), initVector3(0.0f, 0.0f, center)));

		// láthatóság eldöntése skaláris szorzat alapján
		//std::cout << "toCamera:" << toCamera.x << " - " << toCamera.y << " - " << toCamera.z << std::endl;
		//std::cout << "normal:" << normal.x << " - " << normal.y << " - " << normal.z << std::endl;
		//std::cout << "dotProduct:" << dotProduct(normal, toCamera) << std::endl;
		if (dotProduct(normal, toCamera) < 0) {
		
			//glColor3f((double)i / faceTorusz.size(), (double)i / faceTorusz.size(), 0);
			glColor3f(0, 0, 0);
			glLineWidth(2.0);

			glBegin(GL_LINES);
			glVertex2d(pontok[faceTorusz[i].v[0]].x, pontok[faceTorusz[i].v[0]].y);
			glVertex2d(pontok[faceTorusz[i].v[1]].x, pontok[faceTorusz[i].v[1]].y);
			glEnd();
			glBegin(GL_LINES);
			glVertex2d(pontok[faceTorusz[i].v[1]].x, pontok[faceTorusz[i].v[1]].y);
			glVertex2d(pontok[faceTorusz[i].v[2]].x, pontok[faceTorusz[i].v[2]].y);
			glEnd();
			glBegin(GL_LINES);
			glVertex2d(pontok[faceTorusz[i].v[2]].x, pontok[faceTorusz[i].v[2]].y);
			glVertex2d(pontok[faceTorusz[i].v[3]].x, pontok[faceTorusz[i].v[3]].y);
			glEnd();
			glBegin(GL_LINES);
			glVertex2d(pontok[faceTorusz[i].v[3]].x, pontok[faceTorusz[i].v[3]].y);
			glVertex2d(pontok[faceTorusz[i].v[0]].x, pontok[faceTorusz[i].v[0]].y);
			glEnd();

			
			//std::cout <<(double) i / faceTorusz.size() << std::endl;
			//glColor3f((double)i / faceTorusz.size(), (double)i / faceTorusz.size(), (double)i / faceTorusz.size());
			glColor3f(color.x, color.y, color.z);
			glBegin(GL_POLYGON);
			glVertex2d(pontok[faceTorusz[i].v[0]].x, pontok[faceTorusz[i].v[0]].y);
			glVertex2d(pontok[faceTorusz[i].v[1]].x, pontok[faceTorusz[i].v[1]].y);
			glVertex2d(pontok[faceTorusz[i].v[2]].x, pontok[faceTorusz[i].v[2]].y);
			glVertex2d(pontok[faceTorusz[i].v[3]].x, pontok[faceTorusz[i].v[3]].y);
			glEnd();
			
		}

	}
	pontok4d.clear();
	pontok.clear();
	faceTorusz.clear();


}


void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);
	drawGrid(initVector3(0, 0, 0), Sik);
	drawToruszWithPlanes(initVector3(0.0f, 1.0f, 0.0f), Tcx,view);
	//drawTorusz(initVector3(0.0f, 0.0f, 1.0f), Tcx);
	deltaA += 0.002;
	initTransformations();
	glFlush();
}


void keyPressed(GLFWwindow * windows, GLint key, GLint scanCode, GLint action, GLint mods) {
	if (action == GLFW_PRESS || GLFW_REPEAT) {
		switch (key) {
		case GLFW_KEY_LEFT:
			alpha += deltaAlpha;
			eye.x = R*cos(alpha);
			eye.z = R*sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_UP:
			eye.x = R * cos(alpha);
			eye.y += 2;
			eye.z = R * sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_DOWN:
			eye.x = R* cos(alpha);
			eye.y -= 2;
			eye.z = R * sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_RIGHT:
			alpha -= deltaAlpha;
			eye.x = R*cos(alpha);
			eye.z = R*sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_KP_ADD:
			R -= 2;
			eye.x = R* cos(alpha);
			//eye.y 
			eye.z = R * sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_KP_SUBTRACT:
			R += 2;
			eye.x = R* cos(alpha);
			//eye.y 
			eye.z = R * sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		}
	}


	glfwPollEvents();
}

int main(int argc, char ** argv)
{
	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(winWidth, winHeight, "Hazi feladat", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	glfwSetKeyCallback(window, keyPressed);

	init();
	draw();
	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		draw();

		glfwSwapBuffers(window);

		glfwPollEvents();
	}

	glfwTerminate();

	return 0;
}