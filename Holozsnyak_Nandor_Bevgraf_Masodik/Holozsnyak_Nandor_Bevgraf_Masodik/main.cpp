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
int keyPressed = -1;

GLdouble param = 0.5;
#define NUM_CONTROLS 7+1
GLdouble controlPoints[NUM_CONTROLS][2];
int givenPoints = 0;
bool manualErinto = false;

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
}

void drawLineBetweenControlPoints() {
	glBegin(GL_LINES);
	for (int i = 0; i < givenPoints -2 ; i++) {
		glColor3f(1.0, 0.0, 1.0);
		glVertex2f(controlPoints[i][0], controlPoints[i][1]);
		glVertex2f(controlPoints[i + 1][0], controlPoints[i + 1][1]);
	}
	glEnd();

}

void HermiteIv() {
	int t1 = 0, t2 = 1, t3 = 2;
	Matrix<float>O(4,4),G(2,4);
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
	//T = O;
	M = O.Inverz_Matrix();
	//M = M*T;
	givenPoints = 8;
	
	if (manualErinto == false) {
		controlPoints[7][0] = (controlPoints[4][0] + 2 * (controlPoints[4][0] - controlPoints[3][0]));
		controlPoints[7][1] = (controlPoints[4][1] + 2 * (controlPoints[4][1] - controlPoints[3][1]));
	}
	G << controlPoints[4][0] << controlPoints[5][0] << controlPoints[6][0] << controlPoints[7][0];
	G << controlPoints[4][1] << controlPoints[5][1] << controlPoints[6][1] << controlPoints[7][1];
	
	glPointSize(8.0);
	glBegin(GL_POINTS);
	glColor3f(1.0, 0.0, 0.0);
	glVertex2d(controlPoints[7][0], controlPoints[7][1]); //display a point
	glEnd();

	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex2d(controlPoints[4][0], controlPoints[4][1]);
	glVertex2d(controlPoints[7][0], controlPoints[7][1]); //display a point
	glEnd();
	Matrix<float> C = G*M,Z(C.GetY(),T.GetX());
	Matrix<float> T(4, 1);
	//std::cout << "C:" << C;
	glColor3f(0.0, 0.0, 1.0);
	glPointSize(3.0);
	glBegin(GL_POINTS);
	float x, y;

	for (float i = t1; i <= t2; i += 0.001)
	{
		
		T << 3*i*i;
		T << 2*i;
		T << 1;
		T << 0;
		Z = C*T;
		x = Z.GetValue(0, 0);
		y = Z.GetValue(1, 0);
		if (x == controlPoints[4][0] && y == controlPoints[4][1]) {
			
			//std::cout << Z;
			glVertex2f(x, y);
		}
	}

	for (float i = t2; i <= t3; i += 0.001)
	{

		T << i*i*i;
		T << i*i;
		T << i;
		T << 1;
		Z = C*T;
		x = Z.GetValue(0, 0);
		y = Z.GetValue(1, 0);
		//std::cout << Z;
		glVertex2f(x, y);
	}
	glEnd();
	
	
	///std::cout << fin;
}
void turnActualPoint(int actual) {
	double alpha = -3.14159265 / (10 * 360);
	Matrix<double> m(3, 1);
	for (int i = 0; i < givenPoints; i++) {
		if (i != actual) {
			m << controlPoints[i][0] - controlPoints[actual][0] << controlPoints[i][1] - controlPoints[actual][1] << 1;
			m.RotateX(alpha);
			controlPoints[i][0] = m.GetValue(0, 0) + controlPoints[actual][0];
			controlPoints[i][1] = m.GetValue(1, 0) + controlPoints[actual][1];
		}
	}
}

void draw()
{
    glClear (GL_COLOR_BUFFER_BIT);
	
	if (givenPoints == 7 || givenPoints == 8) {
	for (double j = 0; j <= 1; j += 0.001) {
			Casteljau(givenPoints <6 ? givenPoints : 5, controlPoints, j);
	}
	
	HermiteIv();
	
	drawLineBetweenControlPoints();
		if (keyPressed != -1) {
			turnActualPoint(keyPressed);
		}
	}

	glPointSize(8.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < givenPoints; i++) {
		glColor3f(1.0, 0.0, 0.0);
		glVertex2d(controlPoints[i][0], controlPoints[i][1]); //display a point
															  //glVertex2i((it++)->GetValue(0, 0), (it++)->GetValue(1, 0)); //display a poin
	}
	glEnd();

	
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
			
			
			if (selectedPointIndex == 7) {
				manualErinto = true;
				//controlPoints[3][0] = controlPoints[7][0] - controlPoints[5][0];
				//controlPoints[3][1] = controlPoints[5][1] -controlPoints[7][1];
				controlPoints[selectedPointIndex][0] = xpos - winWidth * 0.5f; //set the selected point new coordinates
				controlPoints[selectedPointIndex][1] = ypos - winHeight * 0.5f;
			}
			else {
				manualErinto = false;
				controlPoints[selectedPointIndex][0] = xpos - winWidth * 0.5f; //set the selected point new coordinates
				controlPoints[selectedPointIndex][1] = ypos - winHeight * 0.5f;
			}

		}
		else {
			selectedPointIndex = getActivePoint1(controlPoints, givenPoints, 8, xpos - winWidth * 0.5f, ypos - winHeight * 0.5f); //which point was clicked by user 
		}
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && givenPoints < 8) {
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


void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	
	
		if (key == GLFW_KEY_1 && action==GLFW_PRESS)
		{
			if (keyPressed == 0) {
				keyPressed = -1;
			}
			else {
				keyPressed = 0;
			}

		}
		if (key == GLFW_KEY_2 && action == GLFW_PRESS)
		{
			if (keyPressed == 1) {
				keyPressed = -1;
			}
			else {
				keyPressed = 1;
			}

		}
		if (key == GLFW_KEY_3 && action == GLFW_PRESS)
		{
			if (keyPressed == 2) {
				keyPressed = -1;
			}
			else {
				keyPressed = 2;
			}

		}
		if (key == GLFW_KEY_4 && action == GLFW_PRESS)
		{
			if (keyPressed == 3) {
				keyPressed = -1;
			}
			else {
				keyPressed = 3;
			}

		}
		if (key == GLFW_KEY_5 && action == GLFW_PRESS)
		{
			if (keyPressed == 4) {
				keyPressed = -1;
			}
			else {
				keyPressed = 4;
			}

		}
		if (key == GLFW_KEY_6 && action == GLFW_PRESS)
		{
			if (keyPressed == 5) {
				keyPressed = -1;
			}
			else {
				keyPressed = 5;
			}

		}

		
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
	glfwSetKeyCallback(window, key_callback);
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