#include <GLFW\glfw3.h>
#include <matrix.h>
#include <math.h>
#include <time.h> 

#define PI 3.141592653589793
#define NUM_CONTROLS 7+1
#define MAX_DEGREE 20
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

//struktúrák majd a használandó dolgokhoz
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
//melyik pont van kiválasztva, kezdetben egyik sem.
GLint selectedPointIndex = -1;

GLdouble param = 0.5;
//kontrollpontok száma
GLdouble controlPoints[NUM_CONTROLS][2];
//az eddigi megadott pontok száma, ezt minden egyes kattintásra növeljük.
int givenPoints = 0;
bool manualErinto = false;
//de casteljau algoritmus, a Bézier-görbe elõállításához.
//3 paraméterrel, a tömb méretével, magával a tömbbel és az aktuális param-mal ami a beziér görbe aktuális paramétere, a pont létrehozásához.
void Casteljau(int n, double P[][2], double t) {
	//ez egy segéd tömb lesz, elõször ide átpakoljuk a P pontjainkat és ezekkel számolunk tovább.
	double segedTomb[7][2];

	for (int i = 0; i < n; i++) {
		segedTomb[i][0] = P[i][0];
		segedTomb[i][1] = P[i][1];
	}

	//egy dupla for ciklussal folyamatosan számolgatjuk a pontokat a szakaszokon a megadott képlet alapján.
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
	//amint megvan a kívánt pont azt kirajzoljuk.
	glVertex2dv(segedTomb[0]);
}

//vonalak megrajzolása a kontrollpontok között
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

//Megfelelõ mátrixok létrehozása majd a GMT szorzáshoz
//Ahol a G egy 2x4 es mátrix, az M egy 4x4 -es
//illetve a T egy 4x1
//C = G*M lesz (2x4)
//Z = C*T (2*1)
//O mátrix egy segéd mátrix amibe kezdetben a paraméterek külöbözõ hatványait tesszük bele, és utána pedig invertáljuk majd és azt az M mátrixba rakjuk
Matrix<GLdouble> O(4, 4), G(2, 4);
Matrix<GLdouble> M(4, 4);
Matrix<GLdouble> T(4, 1);
Matrix<GLdouble> C(G.GetX(), M.GetY());
Matrix<GLdouble> Z(C.GetY(), T.GetX());

//Csoba István - Piazza ---- az érintõnek az x,y koordinátái
GLdouble erintoX;
GLdouble erintoY;
//a feladat kiírásban szereplõ paraméterek.
int t1 = 0, t2 = 1, t3 = 2;
void initHermiteIvPont() {
	
	//egy segédváltozót
	double _4th;
	//szimpla tömb feltöltés, a megadott mátrixforma alapján.
	// O << t_i^3
	//	 << t_i^2
	//	 << t_i^1
	//	 << 1
	//illetve az utolsó oszlopa az érintõ adatai
	//O << 3 * t_i ^ 2
	//  << 2 * t_i
	//	<< 1
	//  << 0;
	//Végsõ alak, ami majd a for ciklus után lesz
	//O << t_1 ^ 3 << t_2 ^ 3 << t_3 ^ 3 << 3 * t_1 ^ 2
	//     t_1 ^ 2 << t_2 ^ 2 << t_3 ^ 2 << 2 * t_1 ^ 1
	//	   t_1 ^ 1 << t_2 ^ 1 << t_3 ^ 1 << 1
	//	   1 << 1 << 1 << 0
	std::cout << "------" << std::endl;
	for (int i = 3; i >= 0;i--) {
		//egy kissebb hiba kiküszöbölése, hogy a pow az ne lehessen inf ezért kell egy kis feltétel.

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
	//mátrix invertálás.
	M = O.Inverz_Matrix();
}
//HermiteÍv rajzolásához 
void HermiteIv() {
	
	//itt már a givenPoints 8 mivel az érintõnek is fel kell venni a pontját a megjelenítéshez.
	givenPoints = 8;
	//Csoba István - Piazza
	erintoX = 4 * (controlPoints[4][0] - controlPoints[3][0]);
	erintoY = 4 * (controlPoints[4][1] - controlPoints[3][1]);
	//ha a megválaszott pont az egérrel való mozgatásnál nem a 7-es akkor folyamatosan számoljuk a hetes pontot, ha pl a 4-es pontot mozgatjuk akkor kell is.
	if (selectedPointIndex != 7) {
		controlPoints[7][0] = erintoX + controlPoints[4][0];
		controlPoints[7][1] = erintoY + controlPoints[4][1];
	}

	//G mátrix feltöltése, amiben a 4,5,6 pont van benne, a Hermite ív pontjai illetve az érintõvektor
	G << controlPoints[4][0] << controlPoints[5][0] << controlPoints[6][0] << erintoX;
	G << controlPoints[4][1] << controlPoints[5][1] << controlPoints[6][1] << erintoY;

	/*glPointSize(8.0);
	glBegin(GL_POINTS);
		glColor3f(1.0, 0.0, 0.0);
	glEnd();*/


	//Itt rajzoljuk ki az 5. pontból az érintõbe az egyenest, ami szagatott lesz.
	glEnable(GL_LINE_STIPPLE);
	glBegin(GL_LINE_STRIP);
	glColor3f(0.0, 0.0, 0.0);
	glVertex2d(controlPoints[4][0], controlPoints[4][1]);
	glVertex2d(controlPoints[7][0], controlPoints[7][1]);
	glEnd();
	glDisable(GL_LINE_STIPPLE);
	//C = G*M (G a gemoteria adatok, M a paraméterek helyes behelyettesítése után inverált mátrix)
	//C = G*M;
	//sajnos memory leak miatt ezt az alakot kell használni.
	MatrixMul(&G, &M, &C);

	glColor3f(0.0, 0.0, 1.0);
	glPointSize(1.0);
	glBegin(GL_LINE_STRIP);
	GLdouble x, y;
	//itt fogjuk kirajzolni az ívet, feltöljük egyesével a T mátrixot ami ugyebár egy 4x1-es mátrix,
	//ezzel pedig folyamatosan fogjuk szorozni a már létrejött C mátrixot Z = C*T és a tört vonalakat folyamatosan fogjuk kirajzolni.
	//T mátrix (i^3 << i^2 << i << 1) transzponáltja a C mátrix pedig egy már elõre összeszorzott


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
//forgatás aktuális pont körül.
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
			//az eltolás és forgatási mátrixokat inicializáljuk és végül pedig ezt benyomjuk egy szorzásba, ami egy nagy szorzásként van meg.
			initEltolasForgatasMatrix(alpha, controlPoints[actual][0], controlPoints[actual][1]);
			m << controlPoints[i][0] << controlPoints[i][1] << 1;
			MatrixMul(&MozgatForgat, &m, &m);
			controlPoints[i][0] = m.GetValue(0, 0);
			controlPoints[i][1] = m.GetValue(1, 0);
		}
	}
}

//beziér görbe kirajzolása,ezen belül fut a for ciklus a változó paraméter számmal amiket összekötöget.
void Bezier() {
	glBegin(GL_LINE_STRIP);
	for (double j = 0; j <= 1; j += 0.001) {
		Casteljau(givenPoints < 6 ? givenPoints : 5, controlPoints, j);
	}
	glEnd();
}

//pontok kirajzolása
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

//két pont közötti távolságnégyzet
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
			//ha a hetes pontot fogtuk meg az azt jelenti hogy az az érintõnek a pontja, ilyenkor pedig a manuális számolás kell a 4. pont esetén.
			if (selectedPointIndex == 7) {
				//azaz beállítjuk az új érintõ pontot.
				controlPoints[7][0] = xpos - winWidth * 0.5f; //set the selected point new coordinates
				controlPoints[7][1] = ypos - winHeight * 0.5f;
				//egy segédváltozóba letároljuk az érintõ pontjának és az 5. pontnak a távolságátt
				seged[0][0] = controlPoints[7][0] - controlPoints[4][0];
				seged[0][1] = controlPoints[7][1] - controlPoints[4][1];
				//itt pedig kiszámoljuk a 4. pontot ami a következõk alapján alakul
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
		movePoint = true; //a bal egérgomb lett lenyomva és ekkor most akár mozoghat is egy pont
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