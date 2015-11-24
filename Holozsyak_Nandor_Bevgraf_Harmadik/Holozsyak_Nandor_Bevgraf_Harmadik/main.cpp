#include <GLFW\glfw3.h>                           // (or others, depending on the system in use)
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <Torusz.h>
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
* 4x4 m�trix, sorfolytonosan t�rolva
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
* felt�lti az A m�trixot az egys�gm�trixszal
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
* felt�lti az A m�trixot a (0,0,s) k�z�ppont� centr�lis vet�t�s m�trixsz�val
*/
void initPersProjMatrix(MATRIX4 A, GLdouble s)
{
	initIdentityMatrix(A);

	A[2][2] = 0.0f;
	A[3][2] = -1.0f / s;
}

/**
* felt�lti az A m�trixot a (wx,wy):(wx+ww,wy+wh) ablakb�l (vx,vy):(vx+vw,vy+vh) n�zetbe
* t�rt�n� transzform�ci� m�trixsz�val
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
	// center az a pont,amely fele a kamerat tartjuk, eye a kamera hely�t adja
	VECTOR3 centerMinusEye = initVector3(center.x - eye.x, center.y - eye.y, center.z - eye.z);

	// a fenti vektor -1 szeres�nek egyseg hosszra norm�ltja, ebbol lesz a kamera rendszer�nek z-tengelye
	VECTOR3 f = normalizeVector(initVector3(-centerMinusEye.x, -centerMinusEye.y, -centerMinusEye.z));

	// az up vektor es a leendo z tengelyirany vektorialis szorzata adja a kamera x-tengelyiranyat
	VECTOR3 s = normalizeVector(crossProduct(up, f));

	// a kamera y tengelyiranya a mar elkeszult z irany es x irany vektorialis szorzatak�nt j�n ki
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
* visszaadja az A m�trix �s a v vektor szorzat�t, A*v-t
*/
VECTOR4 mulMatrixVector(MATRIX4 A, VECTOR4 v) {
	return initVector4(
		A[0][0] * v.x + A[0][1] * v.y + A[0][2] * v.z + A[0][3] * v.w,
		A[1][0] * v.x + A[1][1] * v.y + A[1][2] * v.z + A[1][3] * v.w,
		A[2][0] * v.x + A[2][1] * v.y + A[2][2] * v.z + A[2][3] * v.w,
		A[3][0] * v.x + A[3][1] * v.y + A[3][2] * v.z + A[3][3] * v.w);
}

/**
* felt�lti a C m�trixot az A �s B m�trixok szorzat�val, A*B-vel
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
* ablak m�rete
*/
GLdouble winWidth = 800.0f, winHeight = 600.0f;


MATRIX4 view;
double **cam;

/**
* mer�leges �s centr�lis vet�t�sek m�trixai
*/
MATRIX4 Vc;

/**
* Wtv m�trixok
*/
MATRIX4 Wc;

/**
* a fenti m�trixokb�l el��ll�tott k�t transzform�ci�s m�trix
*/
MATRIX4 Tcx, Tcy, Tcz, TC[4];

/**
* seg�dm�trix
*/
MATRIX4 Tmp, TempR, TempE;

/*
* forgat�si m�trixok
*/
MATRIX4 Rx, Ry, Rz;
/*
* eltol�si m�trix
*/
MATRIX4 E;
/*
* sk�l�z�s m�trixa
*/
MATRIX4 S;

MATRIX4 Sik;

GLdouble alpha = 0.0f, deltaAlpha = 3.14f / 180.0f;

/**
* n�zet koordin�t�i
*/
GLdouble cX = (winWidth - winHeight) / 2.0f + 150, cY = 0, cW = winHeight / 2, cH = winHeight / 2 - 100;
/**
* centr�lis vet�t�s k�z�ppontj�nak Z koordin�t�ja
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
* el��ll�tja a sz�ks�ges m�trixokat
*/
double tA = 5,tB = 18,tC = 5;
double sarok = 12;
void initTransformations()
{

	// vet�t�si m�trixok
	initPersProjMatrix(Vc, center);

	// Wtv m�trixok
	initWtvMatrix(Wc, -0.5f, -0.5f, 1.0f, 1.0f, cX, cY, cW, cH);
	//S�k m�trixa
	//T = W�C�K
	mulMatrices(Vc, view, Tmp);
	mulMatrices(Wc, Tmp, Sik);
	
	//T = W�C�K�E�S
	//Azaz ha S a sk�l�z�s, E az eltol�s m�trix�t jel�li egy has�b eset�n, akkor az al�bbi transzform�ci�s l�nc adja a v�gs� m�trixot.
	//5�5 egys�g alapter�let�, 18 egys�gnyi magass�g� has�b
	//T�glalapok
	initSkalazasMatrix(tA, tB, tC);
	// els� t�glatest
	initEltolasMatrix(-sarok+tA/2, 9, sarok - tA / 2);
	mulMatrices(E, S, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, TC[0]);



	// m�dosik t�glatest
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

	//M = W�C�K�L�R, t�rusz m�trixa el�bb forgatom, eltolom, kamera, centr�lis vet�t�s �s w2v
	//T�ruszok
	initRotationMatrixes(deltaA);
	initEltolasMatrix(0, 9, 0);

	// egyik t�rusz - X tengely
	mulMatrices(E, Rx, TempR);
	mulMatrices(view, TempR, TempE);
	mulMatrices(Vc, TempE, Tmp);
	mulMatrices(Wc, Tmp, Tcx);


	// m�sodik - Y tengely
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
* egys�gkocka cs�csai
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
* egys�gkocka lapjainak indexei
*/
GLuint faces[24] =
{
	0,1,2,3, //als�
	0,1,5,4, //jobb
	1,2,6,5, //h�ts�
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

	// be�ll�tja a kocka �leinek vastags�g�t
	glLineWidth(2.0f);

	// be�ll�tja a kocka sz�n�t
	glColor3f(color.x, color.y, color.z);

	// v�gigmegy�nk a lapokon
	for (i = 0; i < 6; i++)
	{
		glBegin(GL_LINE_LOOP);

		// v�gigmegy�nk a cs�csokon
		for (j = 0; j < 4; j++)
		{
			// homog�n koordin�t�s alakra hozzuk a cs�csot
			ph = initVector4(identityCube[faces[id]].x, identityCube[faces[id]].y, identityCube[faces[id]].z, 1.0f);

			// alkalmazzuk a transzform�ci�t
			pt = mulMatrixVector(T, ph);

			// visszat�r�nk inhomog�n alakra
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

			// kirajzoljuk a cs�csot
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

	// be�ll�tja a kocka �leinek vastags�g�t
	glLineWidth(2.0f);

	// be�ll�tja a kocka sz�n�t
	glColor3f(color.x, color.y, color.z);

	glBegin(GL_LINES);
	for (double i = -sarok; i <= sarok; i+=sarok/12) {
			ph = initVector4(-sarok, 0, i, 1.0f);

			// alkalmazzuk a transzform�ci�t
			pt = mulMatrixVector(T, ph);

			// visszat�r�nk inhomog�n alakra
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

			// kirajzoljuk a cs�csot
			glVertex2f(pih.x, pih.y);

			ph = initVector4(sarok, 0, i, 1.0f);

			// alkalmazzuk a transzform�ci�t
			pt = mulMatrixVector(T, ph);

			// visszat�r�nk inhomog�n alakra
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

			// kirajzoljuk a cs�csot
			glVertex2f(pih.x, pih.y);
			//glVertex2f(-12,i);
			//glVertex2f(12, -i);

		}
	glEnd();
	glBegin(GL_LINES);
	for (double i = -sarok ; i <= sarok; i+=sarok/12) {
		ph = initVector4(i, 0, -sarok, 1.0f);

		// alkalmazzuk a transzform�ci�t
		pt = mulMatrixVector(T, ph);

		// visszat�r�nk inhomog�n alakra
		pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

		// kirajzoljuk a cs�csot
		glVertex2f(pih.x, pih.y);

		ph = initVector4(i, 0, sarok, 1.0f);

		// alkalmazzuk a transzform�ci�t
		pt = mulMatrixVector(T, ph);

		// visszat�r�nk inhomog�n alakra
		pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

		// kirajzoljuk a cs�csot
		glVertex2f(pih.x, pih.y);
		//glVertex2f(-12,i);
		//glVertex2f(12, -i);

	}
	glEnd();
}



void init()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0f, winWidth, 0.0f, winHeight, 0.0f, 1.0f);

	eye = initVector3(0.0, 0.0, R); //megadja a kamera poz�ci�j�t (Ez legyen most a z tengely pozit�v fel�n)
	eye.x = R*cos(alpha);
	eye.z = R*sin(alpha);
	centerVec = initVector3(0.0, 0.0, 0.0); //megadja, hogy merre n�z a kamera (Ez legyen most az orig�)
	up = initVector3(0.0, 1.0, 0.0); //megdja, hogy merre van a felfele ir�ny (Ez legyen most az y tengely)

	initViewMatrix(view, eye, centerVec, up);

	initTransformations();
}
double PI = 3.14159265358979323846;
void drawTorusz(VECTOR3 color, MATRIX4 T) {

	double B[4][1];
	float theta, fi;
	float c = 8, a = 1;
	VECTOR4 ph, pt;
	VECTOR3 pih;
	glColor3f(color.x, color.y, color.z);
	for (theta = 0; theta <= 2 * PI + 0.00001; theta += PI/6) {
		glBegin(GL_LINE_LOOP);
		for (fi = 0; fi <= 2 * PI + 0.00001; fi += PI / 6) {
			B[0][0] = (c + a*cos(theta))*cos(fi);
			B[1][0] = a*sin(theta);
			B[2][0] = (c + a*cos(theta))*sin(fi);
			B[3][0] = 1.0;
			// homog�n koordin�t�s alakra hozzuk a cs�csot
			ph = initVector4(B[0][0], B[1][0], B[2][0], 1.0);

			// alkalmazzuk a transzform�ci�t
			pt = mulMatrixVector(T, ph);

			// visszat�r�nk inhomog�n alakra
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

			// kirajzoljuk a cs�csot
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
			// homog�n koordin�t�s alakra hozzuk a cs�csot
			ph = initVector4(B[0][0], B[1][0], B[2][0], 1.0);

			// alkalmazzuk a transzform�ci�t
			pt = mulMatrixVector(T, ph);

			// visszat�r�nk inhomog�n alakra
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

			// kirajzoljuk a cs�csot
			glVertex2f(pih.x, pih.y);

		}
		glEnd();
	}
}

void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);

	drawGrid(initVector3(0.0f, 0.0f, 0.0f), Sik);
	for (int i = 0; i < 4; i++) {
		drawCubes(initVector3(0.0f, 0.0f, 0.0f), TC[i]);
	}
	drawTorusz(initVector3(0.0f, 1.0f, 0.0f), Tcx);
	drawTorusz(initVector3(1.0f, 0.0f, 0.0f), Tcy);
	drawTorusz(initVector3(0.0f, 0.0f, 1.0f), Tcz);
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
