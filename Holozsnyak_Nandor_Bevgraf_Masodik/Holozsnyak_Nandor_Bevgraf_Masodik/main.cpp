#include <GLFW\glfw3.h>
#include <matrix.h>
#include <math.h>
#include <time.h> 

#define PI 3.141592653589793

/*

A második házi feladat a következõ.

Vegyünk fel az egérrel hét pontot a képernyõn.
Ez a hét pont egy 5 kontrollpontos Bézier-görbe és egy "3 pont 1 érintõvektor" típusú Hermite-ív alappontja lesz,
melyek C1 folytonosan csatlakoznak egymáshoz.
A Bézier görbe végpontja és az Hermite-ív kezdõpontja ugyanaz lesz tehát,
ezért van csak hét pont.
Meg kell rajzolni az elsõ öt pont által meghatározott Bézier görbét, de Casteljau algoritmussal,
az Hermite-ívet pedig a GMT formula segítségével CT mátrixszorzással úgy,
hogy a kezdõponthoz 0,
a középsõ ponthoz 1,
a végponthoz 2 paramétert rendelünk.
Az Hermite-ív kezdõpontbeli érintõvektorát úgy kell beállítani, hogy teljesüljön a C1 folytonos csatlakozás a tanultak szerint.

A kontrollpontokat és az Hermite-ív pontjait és érintõvektorát is lehessen egérrel módosítani pontok mozgatása által.
Mindenféle változtatás során meg kell maradnia a C1 folytonos csatlakozásnak.
Ha a közös pontot mozgatjuk, maradhat fixen az utolsó elõtti kontrollpont, vizonyíthatunk ahhoz, így az Hermite-ív érintõvektorát kell majd mindig ilyen esetben újraszámolni.

Az 1-7 billentyûk megnyomására a görbe a megadott sorszámú pont körül kezdjen el forogni, illetve álljon meg ismételt megnyomásra.

*/


typedef struct vector2d { GLdouble x, y; } VECTOR2D;

typedef struct point2d { GLdouble x, y; } POINT2D;

GLsizei winWidth = 800, winHeight = 600;
Matrix<GLdouble> T(4, 4);
VECTOR2D initVector2D(GLdouble x, GLdouble y) {
	VECTOR2D P;
	P.x = x;
	P.y = y;
	return P;
}

POINT2D initPoint2D(GLdouble x, GLdouble y) {
	POINT2D P;
	P.x = x;
	P.y = y;
	return P;
}

void init()
{
    glClearColor (1.0, 1.0, 1.0, 0.0);	
    glMatrixMode (GL_PROJECTION);	
	glLoadIdentity();         
	glPointSize(10.0);
	glOrtho(-(GLdouble)winWidth / 2, (GLdouble)winWidth / 2, (GLdouble)winHeight / 2, -(GLdouble)winHeight / 2, 0, 1);
}

#define MAX_DEGREE 20

int numControls = 5, numControls2;                   /* control pontok számai */
double controls[MAX_DEGREE + 1][2];
double controls2[MAX_DEGREE + 1][3];


GLdouble param = 0.5;
#define NUM_CONTROLS 7
GLdouble controlPoints[NUM_CONTROLS][2];
int givenPoints = 0;

void allBernstein(int n, double u, double B[]) {
	int j, k;
	double u1 = 1.0 - u;
	double saved;
	B[0] = 1.0;
	for (j = 1; j <= n; j++) {
		saved = 0.0;
		for (k = 0; k < j; k++) {
			double temp = B[k];
			B[k] = saved + u1*temp;
			saved = u*temp;
		}
		B[j] = saved;
	}
}

void pointOnBezierCurve(double P[][2], int n, double u, double C[2]) {
	double B[20];
	int k;
	allBernstein(n, u, B);
	C[0] = C[1] = 0.0;
	for (k = 0; k <= n; k++) {
		C[0] += B[k] * P[k][0];
		C[1] += B[k] * P[k][1];
		//C[2] += B[k] * P[k][2];
	}
}
double*  findfixpoint(double P[][3], int n, int m) {
	int i, j;
	double retpoint[1][3];
	retpoint[0][0] = (P[n - 1][0] - P[n - 2][0]) + P[n - 1][0];
	retpoint[0][1] = (P[n - 1][1] - P[n - 2][1]) + P[n - 1][1];
	retpoint[0][2] = 0;

	return retpoint[0];

}

void Casteljau(int n, double P[][2], double t) {
	double tomb[7][2];
	int i, j;
	for (i = 0; i<n; i++) {
		tomb[i][0] = P[i][0];
		tomb[i][1] = P[i][1];
	}

	glColor3f(0.0, 0.0, 1.0);
	glPointSize(3.0);
	glBegin(GL_POINTS);
	
	for (i = n-1;i>0;i--) {
	
		for (j = 0;j <= i ;j++) {
			tomb[j][0] = (1 - t)*tomb[j][0] + t * tomb[j + 1][0];
			tomb[j][1] = (1 - t)*tomb[j][1] + t * tomb[j + 1][1];
		}		
	}
	glVertex2dv(tomb[0]);
	glEnd();
	//glBegin(GL_POINTS);
	//glVertex3dv(tomb[0]);
	//glEnd();
}

void drawLineBetweenControlPoints() {
	glBegin(GL_LINES);
	for (int i = 0; i < givenPoints - 1; i++) {
		glColor3f(1.0, 0.0, 1.0);
		glVertex2f(controlPoints[i][0], controlPoints[i][1]);
		glVertex2f(controlPoints[i + 1][0], controlPoints[i + 1][1]);
	}
	glEnd();

}

void HermiteIv() {
	int t1 = 0, t2 = 1, t3 = 2;
	Matrix<float> T(4, 4),O(4,4),G(4,4);
	Matrix<float> M(4,4);
	double _4th;
	for (int i = 3; i >=0;i--) {
		if (isinf(pow(t1, i - 1))) {
			_4th = 0;
		}
		else {
			_4th = pow(t1, i - 1);
		}
		
		O << pow(t1, i) << pow(t2, i) << pow(t3, i) <<_4th;
	}
	T = O;
	M = O.Inverz_Matrix();
	M = M*O;


	std::cout << M;
}

void draw()
{
    glClear (GL_COLOR_BUFFER_BIT);
	//glColor3f (0.0f, 0.4f, 0.2f);
	//controls[0][1] = 10;
	//controls[0][2] = 15;
	//controls[0][3] = 10.0;
	//controls[1][1] = 30;
	//controls[1][2] = 70;
	//controls[1][3] = 34.0;
	//controls[2][1] = 120;
	//controls[2][2] = 150;
	//controls[2][3] = 76.0;
	//controls[3][1] = 200;
	//controls[3][2] = 180;
	//controls[3][3] = 145.0;
	//controls[4][1] = 450;
	//controls[4][2] = 600;
	//controls[4][3] = 600.0;
	//glBegin(GL_LINE_STRIP);
	/*int i, u;
	for (i = 0, u = 0.0; i < SAMPLES; i++, u += DU) {
		double C[2];
		pointOnBezierCurve(controlPoints, givenPoints - 1, u, C);
		//glVertex3dv(C);
	}*/
	//glEnd();
	for (double j = 0; j <= 1; j += 0.001) {
		Casteljau(givenPoints <6?givenPoints:5, controlPoints, j);
	}
	glPointSize(8.0);
	
	glBegin(GL_POINTS);
	
	for (int i = 0; i < givenPoints; i++) {
		glColor3f(1.0, 0.0, 0.0);
		glVertex2d(controlPoints[i][0], controlPoints[i][1]); //display a point
															//glVertex2i((it++)->GetValue(0, 0), (it++)->GetValue(1, 0)); //display a poin
	}
	HermiteIv();
	glEnd();
	drawLineBetweenControlPoints();
	glFlush ( );
}
double xpos, ypos;

GLint dragged = -1;

GLint Round(GLfloat n) { return (GLint)(n + 0.5); }

bool movePoint = false;
GLint selectedPointIndex = -1;

GLfloat dist2(POINT2D P1, POINT2D P2) {
	GLfloat t1 = P1.x - P2.x;
	GLfloat t2 = P1.y - P2.y;
	return t1 * t1 + t2 * t2;
}


GLint getActivePoint1(GLdouble p[][2], GLint size, GLint sens, GLint x, GLint y) {
	GLint i, s = sens * sens;
	POINT2D P = initPoint2D(x, y),pP;

	for (i = 0; i < size; i++) {
		pP = initPoint2D(p[i][0], p[i][1]);
		if (dist2(pP,P) < s) {
			std::cout << i << std::endl;
			return i;
		}
	}
	return -1;
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (movePoint) {
		if (selectedPointIndex != -1) {
			controlPoints[selectedPointIndex][0] = xpos - winWidth * 0.5f; //set the selected point new coordinates
			controlPoints[selectedPointIndex][1] = ypos - winHeight * 0.5f;
			
		}
		else {
			selectedPointIndex = getActivePoint1(controlPoints, 10, 8, xpos - winWidth * 0.5f, ypos - winHeight * 0.5f); //which point was clicked by user 
		}
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && givenPoints < 7) {
		glfwGetCursorPos(window, &xpos, &ypos);
		GLint convertedX = xpos - winWidth * 0.5f;
		GLint convertedY = ypos- winHeight * 0.5f;
		controlPoints[givenPoints][0] = convertedX;
		controlPoints[givenPoints][1] = convertedY;		
		std::cout << "X:" << controlPoints[givenPoints][0] << " Y:" << controlPoints[givenPoints][1] << std::endl;
		givenPoints++;
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		movePoint = true; //left button was presseb, so we need to do some interact inside the mouseMove function
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		movePoint = false; //left mouse button was released so we can take a rest
		selectedPointIndex = -1; //no point is selected currently
	}

	
		//std::cout << "kattintott" << std::endl;
}


int main (int argc, char** argv)
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

	glfwSetCursorPosCallback(window, cursor_position_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	init();

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