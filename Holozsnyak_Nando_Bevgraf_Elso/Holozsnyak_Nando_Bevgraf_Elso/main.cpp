#include <GLFW\glfw3.h>
#include "Header.h"
#include <math.h>
#include <time.h> 
#include <thread>       
#include <chrono> 
#define PI 3.141592653589793

typedef struct vector2d { GLdouble x, y; } VECTOR2D;

typedef struct point2d { GLdouble x, y; } POINT2D;

typedef struct circle2d { GLdouble x, y, r; } CIRCLE2D;

GLFWwindow* window;


//definiálunk egy kör struktúrák, amiben tároljuk a kör x,y koordinátáját, a tömegét (mass) illetve a sebességvektorát ami M,
//illetve egy cS-el azt mondjuk meg hogy milyen állapota lesz azaz inkább a színét, szürke ha 0, fehér ha 1 és fekete ha 2
typedef struct circle {
	Vector<GLdouble> *M = new Vector<GLdouble>(2);
	GLdouble x, y, r, mass;
	int cS; //circle State - ez mutatja hogy a kör szürke 0 , fehér 1 vagy fekete 2 
} CIRCLE;


//itt adjuk meg hogy fut-e a játék vagy sem.
bool gameRunning = true;
//a gyõztes játékos száma
int winner = 0;
//hány darab körünk legyen
int circleNum = 15;
//hány darab kör van hátra, ez minden egyes körlap a körlappal való ütközésre csökkeni fog.
int circleRemaining = circleNum - 1;
//mekkorák legyen a körök sugara
int circleRadius = 15;
//milyen legyen a körök kezdeti súlya
GLdouble circleMass = 1;
//ablakbeállítások
GLsizei winWidth = 800, winHeight = 600;
//az updatefrequency megmondja hogy milyen gyakran kell majd az ablakot frissíteni, a lastupdate pedig hogy mikor volt utoljára frissítve
GLdouble updateFrequency = 0.01, lastUpdate;
//ebben lesznek a körök cArr(circle Arrays)
CIRCLE *cArr = new CIRCLE[circleNum];
//ebben lesznek a játékokos
CIRCLE *Players = new CIRCLE[2];

//játékosok pontjaik amiket majd meccs közben szerezhetnek
int PlayerPoints[2] = { 0 };

//A1 = (winWidth,winHeight/2);
//A2 = (0,winHeight/2);
//B1 = (winWidth/2,0);
//B1 = (winWidth/2,winHeight);

//e = (a,b,c) = (a2 ? b2, b1 ? a1, a1b2 ? b1a2)


Vector<GLdouble> e1(2), e2(2), e3(2), e4(2);
//Egy vonal struktúra amiben egy vonalnak tudjuk megadni a kezdõ és végponjainak a koordinátáit késõbb ebbõl tudunk majd számolni fal és körlap közötti távolságot.
typedef struct line {
	int x1, x2, y1, y2;
	//Vector<GLdouble> *N;
}LINE2D;
//4 vonalunk lesz amelyek ebben a lines nevû tömmben lesznek tárolva
LINE2D lines[4];
//Egy új lap inicializálásra való függvény
LINE2D newLine(int x1, int y1, int x2, int y2) {
	LINE2D line;
	line.x1 = x1;
	line.y1 = y1;
	line.x2 = x2;
	line.y2 = y2;
	return line;
}


//egy vektor és egy pont távolsága
double lineDistance(Vector<GLdouble> normal, POINT2D P) {
	//Vector <GLdouble> V(2);
	return std::abs((normal.GetValue(0)*P.x + normal.GetValue(1)*P.y)) / (normal.Distance());

}
//random pontokat generál számunkra.
double Randomizer(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}
//2d vector inicializálása
VECTOR2D initVector2D(GLdouble x, GLdouble y) {
	VECTOR2D P;
	P.x = x;
	P.y = y;
	return P;
}
//2d pont inicializálása
POINT2D initPoint2D(GLdouble x, GLdouble y) {
	POINT2D P;
	P.x = x;
	P.y = y;
	return P;
}
//2d kör inicializálása
CIRCLE2D initCircle2D(GLdouble x, GLdouble y, GLdouble r) {
	CIRCLE2D P;
	P.x = x;
	P.y = y;
	P.r = r;
	return P;
}

bool checkCircleSpawn() {
	bool spawnedTogether = false;
	return spawnedTogether;
}
//a négy egyenes inicializálása
void initLines() {
	LINE2D line1, line2, line3, line4;
	//felsõ
	line1.x1 = 0;
	line1.y1 = 0;
	line1.x2 = winWidth;
	line1.y2 = 0;
	//jobb oldali
	line2.x1 = winWidth;
	line2.y1 = 0;
	line2.x2 = winWidth;
	line2.y2 = winHeight;
	//alsó
	line3.x1 = winWidth;
	line3.y1 = winHeight;
	line3.x2 = 0;
	line3.y2 = winHeight;
	//bal oldali
	line4.x1 = 0;
	line4.y1 = winHeight;
	line4.x2 = 0;
	line4.y2 = 0;

	lines[0] = line1;
	lines[1] = line2;
	lines[2] = line3;
	lines[3] = line4;
}

void init()
{
	glClearColor(0.3, 0.3, 0.3, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.f, winWidth, 0.f, winHeight, 0.f, 1.f);
	initLines();
}
//egy kör inicializálása
void circle(CIRCLE O) {

	glBegin(GL_TRIANGLE_FAN);
	for (GLdouble t = 0; t <= 2 * PI; t += 0.01)
		glVertex2f(O.x + O.r * cos(t), O.y + O.r * sin(t));
	glEnd();
}

int szor = 15;

//ez csak debug jelleggel van itt, berajzolja számunkra a körök sebességvektorát
void circleVector(CIRCLE O) {
	/*Vector<GLdouble> a(2);
	a << O.r << O.r;
	double C = (Vector<GLdouble>::innerMultiply(a, *O.M))/(a.Distance()*O.M->Distance());*/

	glBegin(GL_LINES);
	glVertex2f(O.x, O.y);
	glVertex2f(O.x + O.M->GetValue(0)*O.r, O.y + O.M->GetValue(1)*O.r);
	glEnd();

	/*glBegin(GL_TRIANGLE_FAN);
	for (GLdouble t = 0; t <= 2 * PI; t += 0.01)
		glVertex2f(O.x + O.M->GetValue(0) * cos(t), O.y + O.M->GetValue(1) * sin(t));
	glEnd();*/
}
//x tengelyre való tükrözés
void Mirror2X(Vector<GLdouble> *T, Vector<GLdouble> *norm) {
	Vector<double> normal(2);
	normal.Reset();
	normal << 1 << 0;
	if (norm != NULL) {
		*T = T->Mirror(*norm);
	}
	else {
		*T = T->Mirror(normal);
	}
}
//y tengelyre való tükrözés
void Mirror2Y(Vector<GLdouble> *T, Vector<GLdouble> *norm) {
	Vector<double> normal(2);
	normal.Reset();
	normal << 0 << 1;
	if (norm != NULL) {
		*T = T->Mirror(*norm);
	}
	else {
		*T = T->Mirror(normal);
	}
}

//két pont közötti távolságot adja vissza
double pointDistance(POINT2D p1, POINT2D p2) {
	double sum = 0;
	sum += (p2.x - p1.x)*(p2.x - p1.x);
	sum += (p2.y - p1.y)*(p2.y - p1.y);
	return (double)std::sqrt(sum);
}

//körök mozgatása, végigmegynük az összes körön
void MoveBalls() {

	for (int id = 0; id < circleNum; id++) {
		//e = (a,b,c) = (a2 ? b2 b1 ? a1 a1b2 ? b1a2)
		//ezek majd a 4 vonalból jönnek.

		double A;
		float d;
		for (int i = 0; i < 4; i++) {

			//Vector<GLdouble> Vector2D(2), vonal(2), vonal90(2);
			//Vector2D << (cArr[id].x+cArr[id].M->GetValue(0)) - lines[i].x1 << (cArr[id].y+ cArr[id].M->GetValue(1)) - lines[i].y1;
			//vonal << lines[i].x2 - lines[i].x1 << lines[i].y2 - lines[i].y1;
			//vonal90 << vonal.GetValue(1) << vonal.GetValue(0)*-1;
			//double L = Vector2D.PROJ(vonal90);

			//egyenes paraméteres egyenlete igazából ez lesz itt :)
			Vector<GLdouble>egyenes(3), normal(2);
			//Egyenes elõállítása két vektor vektorális szorzataként,
			//itt a koordináták y2-y1,x2-x1,x1*y2-x2*y1 - ezt a képlet alapján levezetve kapjuk meg.
			//Könyebb úgy számolni ha felírjuk a pontokt egy 3x3-mas mátrixba utolsó pontjukat 1-re választjuk meg azaz homogén koordinátába tesszük
			//és ezek után pedig determinánst számolunk és kijönnek ezek az adatok.
			//ami kijön egyenes itt, ezek koordinátái az egyenesek A,B,C paraméterei lesznek.
			egyenes << lines[i].y1 - lines[i].y2 << lines[i].x2 - lines[i].x1 << lines[i].x1*lines[i].y2 - lines[i].x2 * lines[i].y1;
			//az egyenes normálvektora 
			normal << egyenes.GetValue(0) << egyenes.GetValue(1);

			//az akutális egyenesnek  és az aktuális körlapnak a távolságát adjuk meg ebben a d változóban.
			//a kiadott pdf-ben lévõ képlet alapján
			// |Ax_+By_+C|/||n_e||
			d = std::fabs(egyenes.GetValue(0)*cArr[id].x + egyenes.GetValue(1)*cArr[id].y + egyenes.GetValue(2)) / normal.Distance();
			//itt az x tengellyel párhuzamos vonalakat nézzük,
			if (std::fabs(d) <= cArr[id].r && (i % 2) == 1) {
				*cArr[id].M = *cArr[id].M * (-1);
				Mirror2X(cArr[id].M, &normal);
			}//itt pedig egyértelmûen azokat amelyek az y tengellyel párhuzamosak
			else if (std::fabs(d) <= cArr[id].r && (i % 2) == 0) {
				*cArr[id].M = *cArr[id].M * (-1);
				Mirror2Y(cArr[id].M, &normal);
			}
		}
		//itt fogjuk lényegesen mozgatni a golyókat,növeljük a x illetve y koordinátájukat a sebességvektoruk koordinátájával
		cArr[id].x += cArr[id].M->GetValue(0);
		cArr[id].y += cArr[id].M->GetValue(1);
	}

}

void MovePlayers() {
	//Itt a játékosainkat mozgatjuk, ugyan ezen az elven.
	for (int id = 0; id < 2; id++) {
		double d;
		for (int i = 0; i < 4; i++) {

			//Vector<GLdouble> Vector2D(2), vonal(2), vonal90(2);
			//Vector2D << (Players[id].x+Players[id].M->GetValue(0)) - lines[i].x1 << (Players[id].y+ Players[id].M->GetValue(1)) - lines[i].y1;
			//vonal << lines[i].x2 - lines[i].x1 << lines[i].y2 - lines[i].y1;
			//vonal90 << vonal.GetValue(1) << vonal.GetValue(0)*-1;
			//double L = Vector2D.PROJ(vonal90);

			//egyenes paraméteres egyenlete igazából ez lesz itt :)

			//Itt ugyan ebbõl az elgondolásból lesz az egyenes paraméteres elõállítása
			Vector<GLdouble>egyenes(3), normal(2);
			egyenes << lines[i].y1 - lines[i].y2 << lines[i].x2 - lines[i].x1 << lines[i].x1*lines[i].y2 - lines[i].x2 * lines[i].y1;
			normal << egyenes.GetValue(0) << egyenes.GetValue(1);

			//az aktuális kör középpontjától az egyensig való távolság
			d = std::fabs(egyenes.GetValue(0)*Players[id].x + egyenes.GetValue(1)*Players[id].y + egyenes.GetValue(2)) / normal.Distance();

			if (std::fabs(d) <= Players[id].r && (i % 2) == 1) {
				*Players[id].M = *Players[id].M * (-1);
				Mirror2X(Players[id].M, &normal);

			}
			else if (std::fabs(d) <= Players[id].r && (i % 2) == 0) {
				*Players[id].M = *Players[id].M * (-1);
				Mirror2Y(Players[id].M, &normal);

			}
			/*if (Players[id].x + Players[id].r > 800 || Players[id].x - Players[id].r < 0)
			{
				*Players[id].M = *Players[id].M * (-1);
				Mirror2X(Players[id].M);
			}
			if (Players[id].y + Players[id].r > 600 || Players[id].y - Players[id].r < 0)
			{
				*Players[id].M = *Players[id].M * (-1);
				Mirror2Y(Players[id].M);
			}*/
		}
		Players[id].x += Players[id].M->GetValue(0);
		Players[id].y += Players[id].M->GetValue(1);

	}


}


//ebben a függvényben fogjuk két kört ütköztetni.
void collidee(CIRCLE & a, CIRCLE & b) {


	//ez a w vektor ami a két körlap középpontjának különbségeként fog létrejönni. a mozgóból az állóba fog mutatni
	Vector<GLdouble> w(2);
	//a lassabban mozgó körlap sebességvektora lesz a smaller
	Vector<GLdouble> smaller(2);
	//megnézzük melyik a lassabb körlap
	if (a.M->Distance() > b.M->Distance())
	{
		//ha az a körlap a gyorsabb akkor a b körlapot vesszük fel mint lasabb körlap
		smaller = *b.M;
		//a w vektor koordinátáinak meghatározása
		GLdouble x = b.x - a.x, y = b.y - a.y;
		w << x << y;
	}
	else {
		smaller = *a.M;
		GLdouble x = a.x - b.x, y = a.y - b.y;
		w << x << y;
	}
	//itt kivonjuk a már létezõ sebességvektorolból a smaller vektort hogy megtudjuk ekkor melyik fog megállni
	*a.M = *a.M - smaller;
	*b.M = *b.M - smaller;

	Vector<GLdouble> m_1(2);//m1'
	Vector<GLdouble> u1(2);//u1
	Vector<GLdouble> m2(2);//m2 az új sebességvektor az alapból "mozgó" körlapnak
	Vector<GLdouble> a2(2);//a2 az új sebességvektor az alapból "álló" körlapnak
	GLdouble r;

	if (a.M->Distance() > 0) //ha nem nulla akkor A mozog 
	{
		
		//ekkor b az álló kör :)

		//a megadott képletek alapján kiszámoljuk a vektorokat
		//r = az energiamedmardáshoz kell
		r = a.mass / b.mass;
		//m1' = a-nak a w-re vett merõleges vetülete
		m_1 = a.M->projectTo(w);
		//u1 = a ból ki kell vonni m1'-t
		u1 = *a.M - m_1;
		//a megadott képlet alapján újra számoljuk a sebességvektorokat
		m2 = (m_1 + u1) * ((r - 1) / (r + 1)); //mozgó kör új sebességvektora
		a2 = (m_1 * (2 * r)) / (r + 1); // az álló kör új sebességvektora
		//minkettõhöz hozzáadjuk az álló körlap sebességvektorát
		m2 = m2 + smaller;
		a2 = a2 + smaller;
		//ezeket pedig visszadobjuk a körlapokba
		*a.M = m2;
		*b.M = a2;


	}
	else //fordítva
	{
		r = b.mass / a.mass;
		m_1 = b.M->projectTo(w);
		u1 = *b.M - m_1;
		m2 = (m_1 + u1) * ((r - 1) / (r + 1)); // the moving circle's new speedvector
		a2 = (m_1 * (2 * r)) / (r + 1); // the standing circle's new speedvector
		m2 = m2 + smaller;
		a2 = a2 + smaller;
		*a.M = a2;
		*b.M = m2;

	}
}

void PlayerCollideWithCircle(int id) {
	//egy for ciklussal végig megyünk a körökön és megnézzük hogy ütköztethetõek-e vagy sem, illetve itt kapják meg a pontokat a játékosok
	for (int i = 0; i < circleNum; i++) {
		//két pont közötti hosszt vizsgáljuk az aktuális játékos és egy körlap középpontjainak a távolságát, ha ez kisebb mint a játékos és a szürke körlap sugarának az összege akkor ütköztetünk illetve pontokat adunk a játékosnak a kör aktuális cS értéke szerint.
		if (pointDistance(initPoint2D(cArr[i].x, cArr[i].y), initPoint2D(Players[id].x, Players[id].y)) < Players[id].r + cArr[i].r && cArr[i].cS == 0) {
			collidee(cArr[i], Players[id]);
			cArr[i].r = 0;
			cArr[i].x = -100;
			cArr[i].y = -100;
			cArr[i].M->Reset();
			PlayerPoints[id]++;
			circleRemaining--;
			std::cout << "Piros: " << PlayerPoints[0] << "\t Kék:" << PlayerPoints[1] << std::endl;

		}
		//fehér golyó
		if (pointDistance(initPoint2D(cArr[i].x, cArr[i].y), initPoint2D(Players[id].x, Players[id].y)) < Players[id].r + cArr[i].r && cArr[i].cS == 1) {
			collidee(cArr[i], Players[id]);
			cArr[i].r = 0;
			cArr[i].x = -100;
			cArr[i].y = -100;
			cArr[i].M->Reset();
			PlayerPoints[id] += 5;
			circleRemaining--;
			std::cout << "feher golyó kilõve" << std::endl;
			std::cout << "Piros: " << PlayerPoints[0] << "\t Kék:" << PlayerPoints[1] << std::endl;

		}
		//fekete golyó
		if (pointDistance(initPoint2D(cArr[i].x, cArr[i].y), initPoint2D(Players[id].x, Players[id].y)) < Players[id].r + cArr[i].r && cArr[i].cS == 2) {
			collidee(cArr[i], Players[id]);
			cArr[i].r = 0;
			cArr[i].x = -100;
			cArr[i].y = -100;
			cArr[i].M->Reset();
			std::cout << "fekete golyó kilõve" << std::endl;
			winner = id == 1 ? 0 : 1;
			gameRunning = false;
		}


	}
	//std::cout << id << "-edik játkos pontjai:" << PlayerPoints[id] << std::endl;
}

void checkCollide() {
	//a = m1
	//b = w
	//ap = m_1
	//am = u1
	//Az a vektor b-vel párhuzamos komponense : a_p=((a*b)/|b|^2)*b
	//Az a vektor b - re merõleges komponense : a_m = a-a_p
	//Ellenõrzizzük a játék közbeni szürke köröknek az ütköztetését egymással 
	for (int i = 0; i < circleNum; i++) {
		for (int j = i; j < circleNum; j++) {
			if (i != j)
				if (pointDistance(initPoint2D(cArr[i].x + cArr[i].M->GetValue(0), cArr[i].y + cArr[i].M->GetValue(1)), initPoint2D(cArr[j].x + cArr[j].M->GetValue(0), cArr[j].y + cArr[j].M->GetValue(1))) < cArr[i].r + cArr[j].r) {
					collidee(cArr[i], cArr[j]);
				}
		}
	}
	//illetve itt a játékos játékossal való ütközését és a szürke körlap játékossal való ütköztetését
	for (int i = 0; i < 2; i++) {
		PlayerCollideWithCircle(i);
		for (int j = i; j < 2; j++) {
			if (i != j) {
				if (pointDistance(initPoint2D(Players[i].x + Players[i].M->GetValue(0), Players[i].y + Players[i].M->GetValue(1)), initPoint2D(Players[j].x + Players[j].M->GetValue(0), Players[j].y + Players[j].M->GetValue(1))) < Players[i].r + Players[j].r) {
					collidee(Players[i], Players[j]);
				}
			}
		}
	}
}



//mozgások, billentyûvel való mozágsok
void simpleKeyCheck(GLFWwindow * window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_A && action == GLFW_PRESS) {
		Mirror2X(Players[0].M, NULL);
	}
	if (key == GLFW_KEY_S && action == GLFW_PRESS) {
		Mirror2Y(Players[0].M, NULL);
	}
	if (key == GLFW_KEY_K && action == GLFW_PRESS) {
		Mirror2X(Players[1].M, NULL);
	}
	if (key == GLFW_KEY_L && action == GLFW_PRESS) {
		Mirror2Y(Players[1].M, NULL);
	}
}

//Ki nyert?! ha a winner = 0, akkor a piros játékos nyert, ha 1 akkor kék ha 3 akkor pedig döntetlen :O
void winnerColor() {
	glClear(GL_COLOR_BUFFER_BIT);
	if (winner == 0) {
		glClearColor(1, 0, 0, 0);
	}
	else if (winner == 1) {
		glClearColor(0, 0, 1, 0);
	}
	else if (winner == 3) {

		if (time(NULL) % 2) {
			glClearColor(1, 0, 0, 0);
		}
		else {
			glClearColor(0, 0, 1, 0);
		}
	}
}

void draw()
{
	glfwSetKeyCallback(window, simpleKeyCheck);
	double now = glfwGetTime();
	if (now - lastUpdate > updateFrequency) {
		MoveBalls();
		MovePlayers();
		checkCollide();

		if (circleRemaining == 0) {
			winner = PlayerPoints[0] > PlayerPoints[1] ? 0 : 1;
			if (PlayerPoints[0] == PlayerPoints[1]) winner = 3;
			gameRunning = false;
		}
		lastUpdate = now;
	}

	glClear(GL_COLOR_BUFFER_BIT);

	for (int i = 0; i < circleNum - 2; i++) {
		glColor3f(0.6f, 0.6f, 0.6f);
		circle(cArr[i]);
		//glColor3f(1.0, 0.4, 0.2);
		//circleVector(cArr[i]);
	}
	//fekete
	glColor3f(0, 0, 0);
	circle(cArr[circleNum - 2]);

	//glColor3f(0.0, 0.4, 0.2);
	//circleVector(cArr[circleNum - 2]);

	//fehér
	glColor3f(1, 1, 1);
	circle(cArr[circleNum - 1]);

	//glColor3f(0.0, 0.4, 0.2);
	//circleVector(cArr[circleNum - 1]);

	//red - 0
	glColor3f(1, 0, 0);
	circle(Players[0]);

	//glColor3f(0.0, 0.4, 0.2);
	//circleVector(Players[0]);

	//blue - 1
	glColor3f(0, 0, 1);
	circle(Players[1]);

	//glColor3f(0.0, 0.4, 0.2);
	//circleVector(Players[1]);

	glFlush();
}

//itt teszteljük a körök x,y koordinátáinak "valódiságá" hogy nem e szerepel már valakinél ez a koordináta úgymond,
//de valójában a pontok közötti távolságot vizsgáljuk
bool isItFineXY(GLdouble & x, GLdouble & y, int i) {

	for (int j = 0; j < i; j++) {
		//cArr[j].x == x || cArr[j].x + circleRadius == x || cArr[j].x - circleRadius == x || cArr[j].y == y || cArr[j].y + circleRadius == y || cArr[j].y - circleRadius == y
		if (pointDistance(initPoint2D(cArr[j].x, cArr[j].y), initPoint2D(x, y)) <= circleRadius * 2) {
			std::cout << "x,y ujrageneralasa a " << i << " kornel" << std::endl;
			//std::cout << i << " - " << "X: " << x << " Y:" << y << std::endl;
			//j = 0;
			return false;
		}
	}
	return true;
}

void initCircle(int i, int Code) {

	GLdouble x, y, k, l;
	//x,y generálása, x esetén egy 15 és 785 közötti számot generálunk, ne spawonljunk a falbe se kört, y esetén ugyan ez a móka csak más számmal.
	x = rand() % (winWidth - circleRadius * 2) + circleRadius;
	y = rand() % (winHeight - circleRadius * 2) + circleRadius;
	//std::cout << i << " - " << "X: " << x << " Y:" << y << std::endl;
	//ha jó az x,y koordináta akkor spawonljuk a kört, azaz ez a koordináta nem fog belelõgni egyetlen egy körbe sem.
	if (isItFineXY(x, y, i)) {
		//ha a kód 0, azaz szürke golyó akkor ezeket az adatokat kapja meg
		k = Randomizer(2, -2);
		l = Randomizer(2, -2);
		
		if (Code == 0 || Code == 1 || Code == 2) {
			cArr[i].r = circleRadius;
			cArr[i].x = x;
			cArr[i].y = y;
			cArr[i].M->AddItem(k);
			cArr[i].M->AddItem(l);
			cArr[i].cS = Code;
			cArr[i].mass = circleMass;
		}
		else if (Code == 3) {//játékos inicializálás
			Players[i].r = circleRadius;
			Players[i].x = x;
			Players[i].y = y;
			Players[i].M->AddItem(k);
			Players[i].M->AddItem(l);
			Players[i].cS = 2;
			Players[i].mass = circleMass / 2;
		}
	}
	else {//ha nem jó akkor hívjuk meg újra hogy generáljon nekünk egy újabb x,y koordinátát.
		initCircle(i, Code);
	}
}

//A1 = (winWidth,winHeight/2);
//A2 = (0,winHeight/2);
//B1 = (winWidth/2,0);
//B1 = (winWidth/2,winHeight);
//a = (a1,a2) 
//b = (b1,b2)
//e = (a,b,c) = (a2 ? b2, b1 ? a1, a1b2 ? b1a2)

int main(int argc, char** argv)
{

	srand(time(NULL));
	//13 kör inicializálása
	for (int i = 0; i < circleNum - 2; i++) {
		initCircle(i, 0);
	}
	//egy fehér és egy fekete
	initCircle(circleNum - 2, 2);//fekete
	initCircle(circleNum - 1, 1);//fehér 
	initCircle(0, 3);//piros
	initCircle(1, 3);//kék

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
	char windowTitle[50];
	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{

		if (gameRunning) {
			draw(); /* Render here */
			snprintf(windowTitle, sizeof(windowTitle), "Piros: %d - kek: %d", PlayerPoints[0], PlayerPoints[1]);
		}
		else {
			snprintf(windowTitle, sizeof(windowTitle), "Gyoztes jatekos: %s", (winner == 0 ? "Piros" : (winner == 1 ? "Kek" : "Dontetlen")));
			winnerColor();
		}
		glfwSetWindowTitle(window, windowTitle);
		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}
	std::cout << "Pontok:" << std::endl;
	glfwTerminate();

	return 0;
}
