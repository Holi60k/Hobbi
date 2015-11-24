#include <Torusz.h>
#include <math.h>
#include <time.h> 
#define PI 3.141592653589793

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
GLsizei winWidth = 800, winHeight = 600;
/*void init()
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho((GLdouble)winWidth / 2, -(GLdouble)winWidth / 2, -(GLdouble)winHeight / 2, (GLdouble)winHeight / 2, 0, 1);
	//glOrtho(0.0f, winWidth, 0.0f, winHeight, 0.0f, 1.0f);
}
Torusz A, B, C;

void drawGrid() {
	int k = 10;
	int egyseg = 24*10;
	Matrix<double> Pont(4, 1);

	glColor3f(0.0, 0.0, 0.0);
//	Pont << 100 << 0 << 100 << 1;
	glBegin(GL_LINES);
	for (int i = -egyseg; i < egyseg; i+=2) {
		Pont << -i << 0 << egyseg << 1;
		Pont.CVetites(2);
		glVertex3d(Pont.GetValue(0, 0) / Pont.GetValue(3, 0), Pont.GetValue(1, 0) / Pont.GetValue(3, 0), Pont.GetValue(2, 0) / Pont.GetValue(3, 0));
		//glVertex3d(-egyseg,0,i);
		Pont << i << 0 << egyseg << 1;
		glVertex3d(egyseg,0,i);
	}
	glEnd();
	for (int i = -egyseg; i < egyseg; i+=2) {
		glVertex3d(i, 0, -egyseg);
		glVertex3d(i, 0,egyseg);
	}
	glEnd();
	
}
void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0f, 0.0f, 0.0f);
	//A.draw(1);
	glColor3f(0.0f, 1.0f, 0.0f);
	//B.draw(2);
	glColor3f(0.0f, 0.0f, 1.0f);
	//C.draw(3);
	
	drawGrid();
	glFlush();
}

void keyPressed(GLFWwindow * windows, GLint key, GLint scanCode, GLint action, GLint mods) {
	glfwPollEvents();
}*/

