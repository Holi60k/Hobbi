#include <GLFW\glfw3.h>                           // (or others, depending on the system in use)
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <Torusz.h>
#include <vector>

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

void initViewMatrix(MATRIX4 A, VECTOR3 cameye, VECTOR3 center, VECTOR3 up) {
	initIdentityMatrix(A);

	// a kamerabol az iranyt meghatarozo pontba mutato vektor
	// center az a pont,amely fele a kamerat tartjuk, cameye a kamera helyét adja
	VECTOR3 centerMinuscameye = initVector3(center.x - cameye.x, center.y - cameye.y, center.z - cameye.z);

	// a fenti vektor -1 szeresének egyseg hosszra normáltja, ebbol lesz a kamera rendszerének z-tengelye
	VECTOR3 f = normalizeVector(initVector3(-centerMinuscameye.x, -centerMinuscameye.y, -centerMinuscameye.z));

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
	A[0][3] = -dotProduct(s, cameye);
	A[1][3] = -dotProduct(u, cameye);
	A[2][3] = -dotProduct(f, cameye);
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
double tA = 5,tB = 18,tC = 5;
double sarok = 12;
void initTransformations()
{

	// vetítési mátrixok
	initPersProjMatrix(Vc, center);

	// Wtv mátrixok
	initWtvMatrix(Wc, -0.5f, -0.5f, 1.0f, 1.0f, cX, cY, cW, cH);
	//initWtvMatrix(Wc, -4, 0, 4, 0, cX, cY, cW, cH);
	//Sík mátrixa
	//T = W×C×K
	mulMatrices(Vc, view, Tmp);
	mulMatrices(Wc, Tmp, Sik);
	
	//T = W×C×K×E×S
	//Azaz ha S a skálázás, E az eltolás mátrixát jelöli egy hasáb esetén, akkor az alábbi transzformációs lánc adja a végsõ mátrixot.
	//5×5 egység alapterületû, 18 egységnyi magasságú hasáb
	//Téglalapok
	initSkalazasMatrix(tA, tB, tC);
	// elsõ téglatest
	initEltolasMatrix(-sarok+tA/2, 9, sarok - tA / 2);
	mulMatrices(E, S, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, TC[0]);



	// mádosik téglatest
	initEltolasMatrix(sarok - tA / 2, 9, sarok - tA / 2);
	mulMatrices(E, S, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, TC[1]);


	// harmadik
	initEltolasMatrix(-sarok + tA / 2, 9, -sarok + tA / 2);
	mulMatrices(E, S, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, TC[2]);


	// negyedik
	initEltolasMatrix(sarok - tA / 2, 9, -sarok + tA / 2);
	mulMatrices(E, S, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, TC[3]);

	//M = W×C×K×L×R, tórusz mátrixa elõbb forgatom, eltolom, kamera, centrális vetítés és w2v
	//Tóruszok
	initRotationMatrixes(deltaA);
	initEltolasMatrix(0, 0, 0);

	// egyik tórusz - X tengely
	mulMatrices(E, Rx, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, Tcx);


	// második - Y tengely
	mulMatrices(E, Ry, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, Tcy);

	//harmadik - Z tengely
	mulMatrices(E, Rz, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, Tcz);
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

VECTOR3 up, centerVec;
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
	for (double i = -sarok; i <= sarok; i+=sarok/12) {
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
	for (double i = -sarok ; i <= sarok; i+=sarok/12) {
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


/**
* kamera paraméterek
*/
VECTOR3 cameye = initVector3(0.0, 0.0f, R);
VECTOR3 cameraCenter = initVector3(0.0f, 0.0f, 0.0f);
VECTOR3 cameraUp = initVector3(0.0f, 1.0f, 0.0f);

void init()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0f, winWidth, 0.0f, winHeight, 0.0f, 1.0f);

	//cameye = initVector3(0.0, 0.0,R); //megadja a kamera pozícióját (Ez legyen most a z tengely pozitív felén)
	//cameye.x = R*cos(alpha);
	//cameye.z = R*sin(alpha);
	centerVec = initVector3(0.0, 0.0, 0.0); //megadja, hogy merre néz a kamera (Ez legyen most az origó)
	up = initVector3(0.0, 1.0, 0.0); //megdja, hogy merre van a felfele irány (Ez legyen most az y tengely)

	//initViewMatrix(view, cameye, centerVec, up);
	initViewMatrix(view, cameye, cameraCenter, cameraUp);

	initTransformations();
}

double PI = 3.14159265358979323846;


void drawTorusz(VECTOR3 color, MATRIX4 T) {

	double B[4][1];
	float theta, fi;
	float c = 8, a = 1;
	VECTOR4 ph, pt;
	VECTOR3 pih;
	glLineWidth(1.0);
	glColor3f(color.x, color.y, color.z);
	for (theta = 0; theta <= 2 * PI + 0.00001; theta += PI/6) {
		glBegin(GL_LINE_LOOP);
		for (fi = 0; fi <= 2 * PI + 0.00001; fi += PI/6 ) {
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



	for (fi = 0; fi <= 2 * PI + 0.000001; fi += PI / 6) {
		glBegin(GL_LINE_LOOP);
		for (theta = 0; theta <= 2 * PI + 0.0000001; theta += PI / 6) {
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
typedef struct { GLint v[4]; } FACE;

/**
* visszadja a (v0,v1,v2,v3) indexekhez tartozó lapot
*/
FACE initFace(GLint v0, GLint v1, GLint v2, GLint v3) {
	FACE f;

	f.v[0] = v0;
	f.v[1] = v1;
	f.v[2] = v2;
	f.v[3] = v3;

	return f;
}
FACE faceTorusz[13 * 13];
int comparePointsZ(const void *a, const void *b) {
	if ( (*(VECTOR3*)a).z < (*(VECTOR3*)b).z) return -1;
	if ( (*(VECTOR3*)a).z == (*(VECTOR3*)b).z) return  0;
	if ( (*(VECTOR3*)a).z >(*(VECTOR3*)b).z) return  1;
}

void drawToruszWithPlanes(VECTOR3 color, MATRIX4 T) {
	Matrix<double> vectorPl(12,12);
	double B[4][1];
	double theta = PI/6, fi;
	float c = 8, a = 1;
	VECTOR4 ph, pt;
	VECTOR3 pih;
	glLineWidth(1.0);
	int i = 0;
	glColor3f(color.x, color.y, color.z);
	double delta = PI / 6;
	std::vector<VECTOR3> pontok;
	for (theta = 0; theta <= 2 * PI + 0.00001; theta += delta) {
		glBegin(GL_LINE_LOOP);
		for (fi = 0; fi <= 2 * PI + 0.00001; fi += delta) {
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
			pontok.push_back(pih);
			i++;
		}
		glEnd();
	}

	/*for (theta = 0; theta <= 2 * PI; theta += delta) {
		
		glPointSize(4);
		glBegin(GL_POINTS);
		for (fi = 0; fi <= 2 * PI + 0.00001; fi += delta) {
			//std::cout << i++ << " " << std::endl;
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
		}
		glEnd();
	}*/

	qsort(pontok.data(), pontok.size()-1, sizeof(VECTOR3), comparePointsZ);
	
	for (i = 0; i < 13; i++) {
		for (int j = 0; j < 13; j++) {
			faceTorusz[i * 13 + j] = initFace( (GLint) (i*13+j), (GLint)(i * 13 + j+13), (GLint)(i * 13 + j + 13+1), (GLint)(i*13 + j+1));
			/*std::cout << i * 12 + j << " lap pontjainak indexe: " << faceTorusz[i * 12 + j].v[0] << " - "
				<< faceTorusz[i * 12 + j].v[1] << " - "
				<< faceTorusz[i * 12 + j].v[2] << " - "
				<< faceTorusz[i * 12 + j].v[3] << std::endl;*/
		}
	}

	for (i =0; i < 12; i++) {
		for (int j = 0; j <12; j++) {
			glColor3f(0.0, 0, 0);
			glLineWidth(3.0);
			glBegin(GL_LINES);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[0]].x, pontok[faceTorusz[i * 13 + j].v[0]].y);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[1]].x, pontok[faceTorusz[i * 13 + j].v[1]].y);
			glEnd();
			glBegin(GL_LINES);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[1]].x, pontok[faceTorusz[i * 13 + j].v[1]].y);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[2]].x, pontok[faceTorusz[i * 13 + j].v[2]].y);
			glEnd();
			glBegin(GL_LINES);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[2]].x, pontok[faceTorusz[i * 13 + j].v[2]].y);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[3]].x, pontok[faceTorusz[i * 13 + j].v[3]].y);
			glEnd();
			glBegin(GL_LINES);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[3]].x, pontok[faceTorusz[i * 13 + j].v[3]].y);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[0]].x, pontok[faceTorusz[i * 13 + j].v[0]].y);
			glEnd();
			glColor3f(1.0, 0, 0);
			/*switch (j)
			{
			case 0:glColor3f(0.1, 0.5, 1.0);
				break;
			case 1:glColor3f(0.2, 0.5, 1.0);
				break;
			case 2:glColor3f(0.3, 0.5, 1.0);
				break;
			case 3:glColor3f(0.4, 0.5, 1.0);
				break;
			case 4:glColor3f(0.5, 0.5, 1.0);
				break;
			case 5:glColor3f(0.6, 0.5, 1.0);
				break;
			case 6:glColor3f(0.7, 0.5, 1.0);
				break;
			case 7:glColor3f(0.8, 0.5, 1.0);
				break;
			case 8:glColor3f(0.9, 0.5, 1.0);
				break;
			case 9:glColor3f(0.9, 0.7, 1.0);
				break;
			case 10:glColor3f(0.9, 0.8, 1.0);
				break;
			case 11:glColor3f(0.9, 0.9, 1.0);
				break;
			default:
				break;
			}*/
			glBegin(GL_POLYGON);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[0]].x, pontok[faceTorusz[i * 13 + j].v[0]].y);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[1]].x, pontok[faceTorusz[i * 13 + j].v[1]].y);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[2]].x, pontok[faceTorusz[i * 13 + j].v[2]].y);
			glVertex2d(pontok[faceTorusz[i * 13 + j].v[3]].x, pontok[faceTorusz[i * 13 + j].v[3]].y);
			glEnd();


		}
	}
	pontok.clear();
	glColor3f(color.x, color.y, color.z);
	for (fi = 0; fi <= 2*PI + 0.000001; fi += delta) {
		glBegin(GL_LINE_LOOP);
		for (theta = 0; theta <= 2 * PI + 0.0000001; theta += delta) {
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
			//glVertex2f(pih.x, pih.y);
			pontok.push_back(pih);

		}
		glEnd();
	}
}


void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);

	//drawGrid(initVector3(0.2f, 0.2f, 0.2f), Sik);
	/*for (int i = 0; i < 4; i++) {
		drawCubes(initVector3(0.0f, 0.0f, 0.0f), TC[i]);
	}*/
	drawToruszWithPlanes(initVector3(0.0f, 1.0f, 0.0f), Tcx);
	//drawToruszWithPlanes(initVector3(0.0f, 1.0f, 0.0f), Tcy);
	//drawToruszWithPlanes(initVector3(0.0f, 1.0f, 0.0f), Tcz);
	//drawTorusz(initVector3(1.0f, 0.0f, 0.0f), Tcy);
	//drawTorusz(initVector3(0.0f, 0.0f, 1.0f), Tcz);
	deltaA += 0.002;
	initTransformations();
	glFlush();
}

void keyPressed(GLFWwindow * windows, GLint key, GLint scanCode, GLint action, GLint mods) {
	if (action == GLFW_PRESS || GLFW_REPEAT) {
		switch (key) {
		case GLFW_KEY_LEFT:
			alpha += deltaAlpha;
			cameye.x = R*cos(alpha);
			cameye.z = R*sin(alpha);
			initViewMatrix(view, cameye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_UP:
			cameye.x = R * cos(alpha);
			cameye.y += 2;
			cameye.z = R * sin(alpha);
			initViewMatrix(view, cameye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_DOWN:
			cameye.x = R* cos(alpha);
			cameye.y -= 2;
			cameye.z = R * sin(alpha);
			initViewMatrix(view, cameye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_RIGHT:
			alpha -= deltaAlpha;
			cameye.x = R*cos(alpha);
			cameye.z = R*sin(alpha);
			initViewMatrix(view, cameye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_KP_ADD:
			R -= 2;
			cameye.x = R* cos(alpha);
			//cameye.y 
			cameye.z = R * sin(alpha);
			initViewMatrix(view, cameye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_KP_SUBTRACT:
			R += 2;
			cameye.x = R* cos(alpha);
			//cameye.y 
			cameye.z = R * sin(alpha);
			initViewMatrix(view, cameye, centerVec, up);
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
	//draw();
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
