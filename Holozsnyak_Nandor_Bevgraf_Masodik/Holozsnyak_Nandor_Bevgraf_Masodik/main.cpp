#include <GLFW\glfw3.h>
#include <matrix.h>
#include <math.h>
#include <time.h> 

#define PI 3.141592653589793

typedef struct vector2d { GLdouble x, y; } VECTOR2D;

typedef struct point2d { GLdouble x, y; } POINT2D;

GLsizei winWidth = 800, winHeight = 600;

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
	glOrtho(0.f, winWidth, 0.f, winHeight, 0.f, 1.f);
}

void semicircle(POINT2D O, GLdouble r) {

	glBegin(GL_LINE_STRIP);
	for (GLdouble t = 0; t <= PI; t += 0.01)
		glVertex2d(O.x + r * cos(t), O.y + r * sin(t));
	glEnd();
}

void circle(POINT2D O, GLdouble r) {

	glBegin(GL_LINE_LOOP);
	for (GLdouble t = 0; t <= 2 * PI; t += 0.01)
		glVertex2d(O.x + r * cos(t), O.y + r * sin(t));
	glEnd();
}

void DrawXCircles(int xx) {
	
	srand(time(NULL));
	int x;
	int y;
	for (int i = 0; i < xx; i++) {
		 x = rand() % winWidth + 1;
		 y = rand() % winHeight + 1;
		circle(initPoint2D(x, y), 30);
	}	
}

void draw()
{
    glClear (GL_COLOR_BUFFER_BIT);
	glColor3f (0.0f, 0.4f, 0.2f);
	
	glFlush ( );
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