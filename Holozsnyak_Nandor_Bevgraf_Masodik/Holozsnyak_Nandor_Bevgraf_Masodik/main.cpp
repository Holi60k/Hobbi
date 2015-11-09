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
//Matrix<GLdouble> T(4, 4);
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
	glEnable(GL_POINT_SMOOTH);
	
	
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
	for (i = n-1;i>0;i--) {
	
		for (j = 0;j <= i ;j++) {
			tomb[j][0] = (1 - t)*tomb[j][0] + t * tomb[j + 1][0];
			tomb[j][1] = (1 - t)*tomb[j][1] + t * tomb[j + 1][1];
		}		
	}
	glColor3f(0.0, 0.0, 1.0);
	glLineWidth(3.0);
	glVertex2dv(tomb[0]);
}

void drawLineBetweenControlPoints() {
	glBegin(GL_LINES);
	for (int i = 0; i < 4 ; i++) {
		glColor3f(1.0, 0.0, 1.0);
		glVertex2f(controlPoints[i][0], controlPoints[i][1]);
		glVertex2f(controlPoints[i + 1][0], controlPoints[i + 1][1]);
	}
	glEnd();

}
Matrix<GLdouble>O(4, 4), G(2, 4);
Matrix<GLdouble> M(4, 4);
Matrix<GLdouble> T(4, 1);
Matrix<GLdouble> C(G.GetX(), M.GetY());
Matrix<GLdouble> Z(C.GetY(), T.GetX());
//Csoba István - Piazza
GLdouble erintoX;
GLdouble erintoY;
void HermiteIv() {
	int t1 = 0, t2 = 1, t3 = 2;
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
	M = O.Inverz_Matrix();
	givenPoints = 8;
	//Csoba István - Piazza
	erintoX = 4 * (controlPoints[4][0] - controlPoints[3][0]);
	erintoY = 4 * (controlPoints[4][1] - controlPoints[3][1]);
	G << controlPoints[4][0] << controlPoints[5][0] << controlPoints[6][0] << erintoX;
	G << controlPoints[4][1] << controlPoints[5][1] << controlPoints[6][1] << erintoY;
	if (manualErinto == false) {
		//controlPoints[7][0] = (controlPoints[4][0] + 2 * (controlPoints[4][0] - controlPoints[3][0]));
		//controlPoints[7][0] = ;
		//controlPoints[7][1] =  4 * ( controlPoints[4][1] - controlPoints[3][1]);
		//controlPoints[7][1] = (controlPoints[4][1] + 2 * (controlPoints[4][1] - controlPoints[3][1]));
	}
	controlPoints[7][0] = erintoX + controlPoints[4][0];
	controlPoints[7][1] = erintoY + controlPoints[4][1];
	glPointSize(8.0);
	glBegin(GL_POINTS);
		glColor3f(1.0, 0.0, 0.0);
	//	glVertex2d(controlPoints[7][0], controlPoints[7][1]); //display a point
	glEnd();

	glBegin(GL_LINE_STRIP);
		glColor3f(1.0, 0.0, 0.0);
		glVertex2d(controlPoints[4][0], controlPoints[4][1]);
		glVertex2d(controlPoints[7][0], controlPoints[7][1]);
		//controlPoints[7][0] += controlPoints[4][0];
		//controlPoints[7][1] += controlPoints[4][1];
	glEnd();
	
	C = G*M;
	
	glColor3f(0.0, 0.0, 1.0);
	glPointSize(1.0);
	glBegin(GL_LINE_STRIP);
	GLdouble x,y;
	for (float i = t1; i <= t3; i += 0.05)
	{
		T << i*i*i;
		T << i*i;
		T << i;
		T << 1;
		Z = C*T;
		x = Z.GetValue(0, 0);
		y = Z.GetValue(1, 0);
		glVertex2d(x,y);		
	}
	
	glEnd();
	
}
Matrix<double> m(3, 1);
void turnActualPoint(int actual) {
	double alpha = -3.14159265 / (10 * 360);
	
	for (int i = 0; i < givenPoints; i++) {
		if (i != actual) {
			m << controlPoints[i][0] - controlPoints[actual][0] << controlPoints[i][1] - controlPoints[actual][1] << 1;
			m.RotateX(alpha);
			controlPoints[i][0] = m.GetValue(0, 0) + controlPoints[actual][0];
			controlPoints[i][1] = m.GetValue(1, 0) + controlPoints[actual][1];
		}
	}
}
void Bezier() {
	glBegin(GL_LINE_STRIP);
	for (double j = 0; j <= 1; j += 0.001) {
		Casteljau(givenPoints <6 ? givenPoints : 5, controlPoints, j);
	}
	glEnd();
}
void Points() {
	glPointSize(8.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < givenPoints; i++) {
		glColor3f(1.0, 0.0, 0.0);
		//if (i == 7) {
			//controlPoints[7][0] += controlPoints[4][0];
			//controlPoints[7][1] += controlPoints[4][1];
		//}
		//else {
			glVertex2d(controlPoints[i][0], controlPoints[i][1]);
		//}
	}
	glEnd();
}
void draw()
{
    glClear (GL_COLOR_BUFFER_BIT);
	Points();
	if (givenPoints == 7 || givenPoints == 8) {
		Bezier();
		HermiteIv();
		drawLineBetweenControlPoints();
		if (keyPressed != -1) {
			turnActualPoint(keyPressed);
		}
	}

	glFlush ( );
}


GLint dragged = -1;
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
double xpos, ypos;
static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (movePoint) {
		if (selectedPointIndex != -1) {
			
			
			if (selectedPointIndex == 7) {
				manualErinto = true;
				erintoX = 4 * (controlPoints[4][0] - controlPoints[3][0]);
				erintoY = 4 * (controlPoints[4][1] - controlPoints[3][1]);
				controlPoints[selectedPointIndex][0] = xpos - winWidth * 0.5f; //set the selected point new coordinates
				controlPoints[selectedPointIndex][1] = ypos - winHeight * 0.5f;
				controlPoints[3][0] = (4 * controlPoints[4][0] - controlPoints[7][0]) / 4;
				controlPoints[3][1] = (4 * controlPoints[4][1] - controlPoints[7][1]) / 4;				
				
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
GLint convertedX, convertedY;
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && givenPoints < 8) {
		glfwGetCursorPos(window, &xpos, &ypos);
		convertedX = xpos - winWidth * 0.5f;
		convertedY = ypos- winHeight * 0.5f;
		controlPoints[givenPoints][0] = convertedX;
		controlPoints[givenPoints][1] = convertedY;		
		std::cout << "X:" << controlPoints[givenPoints][0] << " Y:" << controlPoints[givenPoints][1] << std::endl;
		givenPoints++;
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		movePoint = true; //left button was presseb, so we need to do some interact inside the mouseMove function
		glfwGetCursorPos(window, &xpos, &ypos);
		convertedX = xpos - winWidth * 0.5f;
		convertedY = ypos - winHeight * 0.5f;
		std::cout << "X:" << xpos << " Y:" << ypos << std::endl;
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
		if (key == GLFW_KEY_7 && action == GLFW_PRESS)
		{
			if (keyPressed == 6) {
				keyPressed = -1;
			}
			else {
				keyPressed = 6;
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