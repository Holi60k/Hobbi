#include <includes.h>
#define _USE_MATH_DEFINES

class Torusz
{
public:
	Torusz(double alphax, double alphay,double alphaz) {
		this->alphax = alphax;
		this->alphay = alphay;
		this->alphaz = alphaz;
	};
	~Torusz() {
	};
	void draw(int rotation) {

		Vector<double> B(4);
		float theta, fi;
		float c = 8, a = 1;
		for (theta = 0; theta <= 2 * PI + 0.00001; theta += 0.6) {
			glBegin(GL_LINE_LOOP);
			for (fi = 0; fi <= 2 * PI+0.00001; fi += 0.3) {
				B << (c + a*cos(theta))*cos(fi)
					<< a*sin(theta)
					<< (c + a*cos(theta))*sin(fi)
					<< 1.0;
				
				B.Rotate(rotation, alpha);
				
				
				B.MoveTo(0, 9, 0);
				B.CVetites(B, 6);
				B.Window2Viewport();
				glVertex3d(B.GetValue(0, 0) / B.GetValue(3, 0), B.GetValue(1, 0) / B.GetValue(3, 0), B.GetValue(2, 0) / B.GetValue(3, 0));
			}
			glEnd();
		}
		


		for (fi = 0; fi <= 2 * PI + 0.000001; fi += 0.6) {
			glBegin(GL_LINE_LOOP);
			for (theta = 0; theta <= 2 * PI + 0.0000001; theta += 0.3) {
				B << (c + a*cos(theta))*cos(fi)
				<< a*sin(theta)
				<< (c + a*cos(theta))*sin(fi)
				<< 1.0;

				B.Rotate(rotation, alpha);
				B.MoveTo(0, 9, 0);
				B.CVetites(B, 6);
				B.Window2Viewport();
				glVertex3d(B.GetValue(0, 0) / B.GetValue(3, 0), B.GetValue(1, 0) / B.GetValue(3, 0), B.GetValue(2, 0) / B.GetValue(3, 0));

			}
			glEnd();
		}
		alpha += 0.1;


	}
	void rotateX()
	{
		alphax += 0.1;
		alpha += 0.1;
	}
	void rotateY()
	{
		alphay += 0.1;
		alpha += 0.1;
	}
	void rotateZ()
	{
		alphaz += 0.1;
		alpha += 0.1;
	}
private:
	double PI = 3.14159265358979323846;
	double alphax, alphay, alphaz;
	double alpha;

};