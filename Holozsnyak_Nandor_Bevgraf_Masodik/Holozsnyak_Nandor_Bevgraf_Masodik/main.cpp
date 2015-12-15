#include <GLFW\glfw3.h>
#include <matrix.h>
#include <math.h>
#include <time.h> 

#define PI 3.141592653589793
#define NUM_CONTROLS 7+1
#define MAX_DEGREE 20
/*
A m�sodik h�zi feladat a k�vetkez�.

Vegy�nk fel az eg�rrel h�t pontot a k�perny�n.
Ez a h�t pont egy 5 kontrollpontos B�zier-g�rbe �s egy "3 pont 1 �rint�vektor" t�pus� Hermite-�v alappontja lesz,
melyek C1 folytonosan csatlakoznak egym�shoz.
A B�zier g�rbe v�gpontja �s az Hermite-�v kezd�pontja ugyanaz lesz teh�t,
ez�rt van csak h�t pont.
Meg kell rajzolni az els� �t pont �ltal meghat�rozott B�zier g�rb�t, de Casteljau algoritmussal,
az Hermite-�vet pedig a GMT formula seg�ts�g�vel CT m�trixszorz�ssal �gy,
hogy a kezd�ponthoz 0,
a k�z�ps� ponthoz 1,
a v�gponthoz 2 param�tert rendel�nk.
Az Hermite-�v kezd�pontbeli �rint�vektor�t �gy kell be�ll�tani, hogy teljes�lj�n a C1 folytonos csatlakoz�s a tanultak szerint.

A kontrollpontokat �s az Hermite-�v pontjait �s �rint�vektor�t is lehessen eg�rrel m�dos�tani pontok mozgat�sa �ltal.
Mindenf�le v�ltoztat�s sor�n meg kell maradnia a C1 folytonos csatlakoz�snak.
Ha a k�z�s pontot mozgatjuk, maradhat fixen az utols� el�tti kontrollpont, vizony�thatunk ahhoz, �gy az Hermite-�v �rint�vektor�t kell majd mindig ilyen esetben �jrasz�molni.

Az 1-7 billenty�k megnyom�s�ra a g�rbe a megadott sorsz�m� pont k�r�l kezdjen el forogni, illetve �lljon meg ism�telt megnyom�sra.

*/

//strukt�r�k majd a haszn�land� dolgokhoz
typedef struct vector2d { GLdouble x, y; } VECTOR2D;

typedef struct point2d { GLdouble x, y; } POINT2D;
void initHermiteIvPont();
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
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glPointSize(10.0);
	glEnable(GL_POINT_SMOOTH);
	glOrtho(-(GLdouble)winWidth / 2, (GLdouble)winWidth / 2, (GLdouble)winHeight / 2, -(GLdouble)winHeight / 2, 0, 1);
	glLineStipple(1, 0x00FF);
	initHermiteIvPont();
}

//melyik gomb van lenyomva, kezdetben egy sem
int keyPressed = -1;
//melyik pont van kiv�lasztva, kezdetben egyik sem.
GLint selectedPointIndex = -1;

GLdouble param = 0.5;
//kontrollpontok sz�ma
GLdouble controlPoints[NUM_CONTROLS][2];
//az eddigi megadott pontok sz�ma, ezt minden egyes kattint�sra n�velj�k.
int givenPoints = 0;
bool manualErinto = false;
//de casteljau algoritmus, a B�zier-g�rbe el��ll�t�s�hoz.
//3 param�terrel, a t�mb m�ret�vel, mag�val a t�mbbel �s az aktu�lis param-mal ami a bezi�r g�rbe aktu�lis param�tere, a pont l�trehoz�s�hoz.
void Casteljau(int n, double P[][2], double t) {
	//ez egy seg�d t�mb lesz, el�sz�r ide �tpakoljuk a P pontjainkat �s ezekkel sz�molunk tov�bb.
	double segedTomb[7][2];

	for (int i = 0; i < n; i++) {
		segedTomb[i][0] = P[i][0];
		segedTomb[i][1] = P[i][1];
	}

	//egy dupla for ciklussal folyamatosan sz�molgatjuk a pontokat a szakaszokon a megadott k�plet alapj�n.
	for (int i = n - 1;i > 0;i--) {
		for (int j = 0;j <= i;j++) {
			//(1-t)*P[j][0] + t*P[j+1][0]
			segedTomb[j][0] = (1 - t)*segedTomb[j][0] + t * segedTomb[j + 1][0];
			//(1-t)*P[j][0] + t*P[j+1][0]
			segedTomb[j][1] = (1 - t)*segedTomb[j][1] + t * segedTomb[j + 1][1];
		}
	}
	glColor3f(0.0, 0.0, 1.0);
	glLineWidth(3.0);
	//amint megvan a k�v�nt pont azt kirajzoljuk.
	glVertex2dv(segedTomb[0]);
}

//vonalak megrajzol�sa a kontrollpontok k�z�tt
void drawLineBetweenControlPoints() {
	glBegin(GL_LINES);
	for (int i = 0; i < 4; i++) {
		glColor3f(1.0, 0.0, 1.0);
		glVertex2f(controlPoints[i][0], controlPoints[i][1]);
		glVertex2f(controlPoints[i + 1][0], controlPoints[i + 1][1]);
	}
	glEnd();

}

void MatrixMul(Matrix<double> *A, Matrix<double> *B, Matrix<double> *C) {
	double Var = 0;

	if (A->GetY() == B->GetX()) {
		for (int i = 0; i < A->GetY();i++) {
			for (int j = 0; j < A->GetX(); j++) {

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

//Megfelel� m�trixok l�trehoz�sa majd a GMT szorz�shoz
//Ahol a G egy 2x4 es m�trix, az M egy 4x4 -es
//illetve a T egy 4x1
//C = G*M lesz (2x4)
//Z = C*T (2*1)
//O m�trix egy seg�d m�trix amibe kezdetben a param�terek k�l�b�z� hatv�nyait tessz�k bele, �s ut�na pedig invert�ljuk majd �s azt az M m�trixba rakjuk
Matrix<GLdouble> O(4, 4), G(2, 4);
Matrix<GLdouble> M(4, 4);
Matrix<GLdouble> T(4, 1);
Matrix<GLdouble> C(G.GetX(), M.GetY());
Matrix<GLdouble> Z(C.GetY(), T.GetX());

//Csoba Istv�n - Piazza ---- az �rint�nek az x,y koordin�t�i
GLdouble erintoX;
GLdouble erintoY;
//a feladat ki�r�sban szerepl� param�terek.
int t1 = 0, t2 = 1, t3 = 2;
void initHermiteIvPont() {
	
	//egy seg�dv�ltoz�t
	double _4th;
	//szimpla t�mb felt�lt�s, a megadott m�trixforma alapj�n.
	// O << t_i^3
	//	 << t_i^2
	//	 << t_i^1
	//	 << 1
	//illetve az utols� oszlopa az �rint� adatai
	//O << 3 * t_i ^ 2
	//  << 2 * t_i
	//	<< 1
	//  << 0;
	//V�gs� alak, ami majd a for ciklus ut�n lesz
	//O << t_1 ^ 3 << t_2 ^ 3 << t_3 ^ 3 << 3 * t_1 ^ 2
	//     t_1 ^ 2 << t_2 ^ 2 << t_3 ^ 2 << 2 * t_1 ^ 1
	//	   t_1 ^ 1 << t_2 ^ 1 << t_3 ^ 1 << 1
	//	   1 << 1 << 1 << 0
	std::cout << "------" << std::endl;
	for (int i = 3; i >= 0;i--) {
		//egy kissebb hiba kik�sz�b�l�se, hogy a pow az ne lehessen inf ez�rt kell egy kis felt�tel.

		if (isinf(pow(t1, i - 1))) {
			std::cout << "-inf-" << std::endl;
			_4th = 0;
		}
		else {
			_4th = i*pow(t1, i - 1);
		}
		std::cout << _4th << std::endl;
		O << pow(t1, i) << pow(t2, i) << pow(t3, i) << _4th;
	}
	std::cout << "------" << std::endl;
	//std::cout << O;
	//m�trix invert�l�s.
	M = O.Inverz_Matrix();
}
//Hermite�v rajzol�s�hoz 
void HermiteIv() {
	
	//itt m�r a givenPoints 8 mivel az �rint�nek is fel kell venni a pontj�t a megjelen�t�shez.
	givenPoints = 8;
	//Csoba Istv�n - Piazza
	erintoX = 4 * (controlPoints[4][0] - controlPoints[3][0]);
	erintoY = 4 * (controlPoints[4][1] - controlPoints[3][1]);
	//ha a megv�laszott pont az eg�rrel val� mozgat�sn�l nem a 7-es akkor folyamatosan sz�moljuk a hetes pontot, ha pl a 4-es pontot mozgatjuk akkor kell is.
	if (selectedPointIndex != 7) {
		controlPoints[7][0] = erintoX + controlPoints[4][0];
		controlPoints[7][1] = erintoY + controlPoints[4][1];
	}

	//G m�trix felt�lt�se, amiben a 4,5,6 pont van benne, a Hermite �v pontjai illetve az �rint�vektor
	G << controlPoints[4][0] << controlPoints[5][0] << controlPoints[6][0] << erintoX;
	G << controlPoints[4][1] << controlPoints[5][1] << controlPoints[6][1] << erintoY;

	/*glPointSize(8.0);
	glBegin(GL_POINTS);
		glColor3f(1.0, 0.0, 0.0);
	glEnd();*/


	//Itt rajzoljuk ki az 5. pontb�l az �rint�be az egyenest, ami szagatott lesz.
	glEnable(GL_LINE_STIPPLE);
	glBegin(GL_LINE_STRIP);
	glColor3f(0.0, 0.0, 0.0);
	glVertex2d(controlPoints[4][0], controlPoints[4][1]);
	glVertex2d(controlPoints[7][0], controlPoints[7][1]);
	glEnd();
	glDisable(GL_LINE_STIPPLE);
	//C = G*M (G a gemoteria adatok, M a param�terek helyes behelyettes�t�se ut�n inver�lt m�trix)
	//C = G*M;
	//sajnos memory leak miatt ezt az alakot kell haszn�lni.
	MatrixMul(&G, &M, &C);

	glColor3f(0.0, 0.0, 1.0);
	glPointSize(1.0);
	glBegin(GL_LINE_STRIP);
	GLdouble x, y;
	//itt fogjuk kirajzolni az �vet, felt�lj�k egyes�vel a T m�trixot ami ugyeb�r egy 4x1-es m�trix,
	//ezzel pedig folyamatosan fogjuk szorozni a m�r l�trej�tt C m�trixot Z = C*T �s a t�rt vonalakat folyamatosan fogjuk kirajzolni.
	//T m�trix (i^3 << i^2 << i << 1) transzpon�ltja a C m�trix pedig egy m�r el�re �sszeszorzott


	for (float i = t1; i <= t3; i += 0.05)
	{
		T << i*i*i;
		T << i*i;
		T << i;
		T << 1;
		//Z = C*T;
		MatrixMul(&C, &T, &Z);
		x = Z.GetValue(0, 0);
		y = Z.GetValue(1, 0);
		glVertex2d(x, y);
	}
	glEnd();
}
//forgat�s aktu�lis pont k�r�l.
Matrix<double> m(3, 1),Forgatas(3,3),EltolasOrigoba(3,3),EltolasPontba(3,3),Temp1(3,3),MozgatForgat(3,3);

void initEltolasForgatasMatrix(double alpha, double x, double y) {


	//
	EltolasPontba << 1 << 0 << x;
	EltolasPontba << 0 << 1 << y;
	EltolasPontba << 0 << 0 << 1;
	
	Forgatas << cos(alpha) << -sin(alpha) << 0;
	Forgatas << sin(alpha) << cos(alpha) << 0;
	Forgatas << 0 << 0 << 1;

	MatrixMul(&EltolasPontba, &Forgatas, &Temp1);

	EltolasOrigoba << 1 << 0 << -x;
	EltolasOrigoba << 0 << 1 << -y;
	EltolasOrigoba << 0 << 0 << 1;
	MatrixMul(&Temp1, &EltolasOrigoba, &MozgatForgat);
	
}
void turnActualPoint(int actual) {
	double alpha = -3.14159265 / (10 * 360);

	for (int i = 0; i < givenPoints; i++) {
		if (i != actual) {
			//az eltol�s �s forgat�si m�trixokat inicializ�ljuk �s v�g�l pedig ezt benyomjuk egy szorz�sba, ami egy nagy szorz�sk�nt van meg.
			initEltolasForgatasMatrix(alpha, controlPoints[actual][0], controlPoints[actual][1]);
			m << controlPoints[i][0] << controlPoints[i][1] << 1;
			MatrixMul(&MozgatForgat, &m, &m);
			controlPoints[i][0] = m.GetValue(0, 0);
			controlPoints[i][1] = m.GetValue(1, 0);
		}
	}
}

//bezi�r g�rbe kirajzol�sa,ezen bel�l fut a for ciklus a v�ltoz� param�ter sz�mmal amiket �sszek�t�get.
void Bezier() {
	glBegin(GL_LINE_STRIP);
	for (double j = 0; j <= 1; j += 0.001) {
		Casteljau(givenPoints < 6 ? givenPoints : 5, controlPoints, j);
	}
	glEnd();
}

//pontok kirajzol�sa
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
	glClear(GL_COLOR_BUFFER_BIT);
	Points();
	if (givenPoints == 7 || givenPoints == 8) {
		Bezier();
		HermiteIv();
		drawLineBetweenControlPoints();
		if (keyPressed != -1) {
			turnActualPoint(keyPressed);
		}
	}

	glFlush();
}


GLint dragged = -1;
bool movePoint = false;

//k�t pont k�z�tti t�vols�gn�gyzet
GLfloat dist2(POINT2D P1, POINT2D P2) {
	GLfloat t1 = P1.x - P2.x;
	GLfloat t2 = P1.y - P2.y;
	return t1 * t1 + t2 * t2;
}

//megfogunk egy pontot,
GLint getActivePoint1(GLdouble p[][2], GLint size, GLint sens, GLint x, GLint y) {
	GLint i, s = sens * sens;
	POINT2D P = initPoint2D(x, y), pP;

	for (i = 0; i < size; i++) {
		pP = initPoint2D(p[i][0], p[i][1]);
		if (dist2(pP, P) < s) {
			std::cout << i << std::endl;
			return i;
		}
	}
	return -1;
}

double xpos, ypos;
double seged[2][2];
static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (movePoint) {
		if (selectedPointIndex != -1) {
			//ha a hetes pontot fogtuk meg az azt jelenti hogy az az �rint�nek a pontja, ilyenkor pedig a manu�lis sz�mol�s kell a 4. pont eset�n.
			if (selectedPointIndex == 7) {
				//azaz be�ll�tjuk az �j �rint� pontot.
				controlPoints[7][0] = xpos - winWidth * 0.5f; //set the selected point new coordinates
				controlPoints[7][1] = ypos - winHeight * 0.5f;
				//egy seg�dv�ltoz�ba let�roljuk az �rint� pontj�nak �s az 5. pontnak a t�vols�g�tt
				seged[0][0] = controlPoints[7][0] - controlPoints[4][0];
				seged[0][1] = controlPoints[7][1] - controlPoints[4][1];
				//itt pedig kisz�moljuk a 4. pontot ami a k�vetkez�k alapj�n alakul
				controlPoints[3][0] = (seged[0][0] / -4) + controlPoints[4][0];
				controlPoints[3][1] = (seged[0][1] / -4) + controlPoints[4][1];
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
		convertedY = ypos - winHeight * 0.5f;
		controlPoints[givenPoints][0] = convertedX;
		controlPoints[givenPoints][1] = convertedY;
		std::cout << "X:" << controlPoints[givenPoints][0] << " Y:" << controlPoints[givenPoints][1] << std::endl;
		givenPoints++;
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		movePoint = true; //a bal eg�rgomb lett lenyomva �s ekkor most ak�r mozoghat is egy pont
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


	if (key == GLFW_KEY_1 && action == GLFW_PRESS)
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



int main(int argc, char** argv)
{

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(winWidth, winHeight, "2. hazi feladat", NULL, NULL);
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