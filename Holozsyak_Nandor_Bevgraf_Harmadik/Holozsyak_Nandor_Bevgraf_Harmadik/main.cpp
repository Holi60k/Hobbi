#include <Torusz.h>
#include <math.h>
#include <time.h> 
#define PI 3.141592653589793

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

	for (i = 0; i<4; i++)
		for (j = 0; j<4; j++)
			A[i][j] = 0.0f;

	for (i = 0; i<4; i++)
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
* történ? transzformáció mátrixszával
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

/**
* mer?leges és centrális vetítések mátrixai
*/
MATRIX4 Vc;

/**
* Wtv mátrixok
*/
MATRIX4 Wc;

/**
* a fenti mátrixokból el?állított két transzformációs mátrix
*/
MATRIX4 Tc;

/**
* segédmátrix
*/
MATRIX4 Tmp;

GLdouble alpha = 0.0f, deltaAlpha = 3.14f / 180.0f;

/**
* nézet koordinátái
*/
GLdouble cX = (winWidth - winHeight) / 2.0f, cY = 0.0f, cW = winHeight, cH = winHeight;
/**
* centrális vetítés középpontjának Z koordinátája
*/
GLdouble center = 2.0f;

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

/**
* el?állítja a szükséges mátrixokat
*/
void initTransformations()
{
	// vetítési mátrixok
	initPersProjMatrix(Vc, center);

	// Wtv mátrixok
	initWtvMatrix(Wc, -0.5f, -0.5f, 1.0f, 1.0f, cX, cY, cW, cH);

	// centrális mátrixok
	mulMatrices(Vc, view, Tmp);
	mulMatrices(Wc, Tmp, Tc);
}

/*======================================*/

VECTOR3 eye, up, centerVec;

double R = 1;

void init()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0f, winWidth, 0.0f, winHeight, 0.0f, 1.0f);

	eye = initVector3(0.0, 0.0, R); //megadja a kamera pozícióját (Ez legyen most a z tengely pozitív felén)
	centerVec = initVector3(0.0, 0.0, 0.0); //megadja, hogy merre néz a kamera (Ez legyen most az origó)
	up = initVector3(0.0, 1.0, 0.0); //megdja, hogy merre van a felfele irány (Ez legyen most az y tengely)

	initViewMatrix(view, eye, centerVec, up);

	initTransformations();
}

/*
void init()
{
	glClearColor(1.0, 1.0, 1.0, 0.0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(0.f, winWidth, 0.f, winHeight, 0.f, 1.f);
	//Origó központúak vagyunk most már
	glOrtho(-(GLdouble)winWidth / 2, (GLdouble)winWidth / 2, (GLdouble)winHeight / 2, -(GLdouble)winHeight / 2, 0, 1);
}*/

void MatrixMul(Matrix<double> *A, Matrix<double> *B, Matrix<double> *C) {
	double Var = 0;
	if (A->GetY() == B->GetX()) {

		for (int i = 0; i < C->GetY();i++) {
			for (int j = 0; j < C->GetX(); j++) {

				for (int l = 0; l < A->GetY();l++) {
					Var += A->GetValue(j, l) * B->GetValue(l, i);
				}
				C->FillMatrix(Var, j, i);
				Var = 0;

			}

		}

	}
	else {
		std::cout << "Sorry I can not multiplicate these two Matrixes..." << std::endl;

	}

}

GLdouble updateFrequency = 0.01, lastUpdate;
Torusz A(PI/2,PI/6,0), B(0, 0,0), C(0, 0,0);
void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0f, 0.4f, 0.2f);
	double now = glfwGetTime();
	//if (now - lastUpdate > updateFrequency) {
		glColor3f(1.0f, 0.0f, 0.0f);
		A.draw(1);
		glColor3f(0.0f, 1.0f, 0.0f);
		B.draw(2);
		glColor3f(0.0f, 0.0f, 1.0f);
		C.draw(3);
		//lastUpdate = now;
	//}

	glFlush();
}

void keyPressed(GLFWwindow * windows, GLint key, GLint scanCode, GLint action, GLint mods) {
	if (action == GLFW_PRESS || GLFW_REPEAT) {
		switch (key) {
		case GLFW_KEY_LEFT:
			alpha += deltaAlpha;
			eye.x = R*cos(alpha);
			eye.y = 0.0;
			eye.z = R*sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_UP:
			R -= 0.1f;
			eye.x = R * cos(alpha);
			eye.y = 0.0;
			eye.z = R * sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_DOWN:
			R += 0.1f;
			eye.x = R* cos(alpha);
			eye.y = 0.0;
			eye.z = R * sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_RIGHT:
			alpha -= deltaAlpha;
			eye.x = R*cos(alpha);
			eye.y = 0.0;
			eye.z = R*sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		}
	}

	glfwPollEvents();
}

int main(int argc, char** argv)
{

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(winWidth, winHeight, "Circle", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	init();
	glfwSetKeyCallback(window, keyPressed);
	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		draw(); /* Render here */

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	glfwTerminate();

	return 0;
}