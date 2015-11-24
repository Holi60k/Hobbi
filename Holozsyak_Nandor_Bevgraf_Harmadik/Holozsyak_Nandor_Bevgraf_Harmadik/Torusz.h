#include <includes.h>
#define _USE_MATH_DEFINES

class Torusz
{
public:
	Torusz(double alphax, double alphay, double alphaz) {
		this->alphax = alphax;
		this->alphay = alphay;
		this->alphaz = alphaz;
	};
	Torusz() {};
	~Torusz() {
	};
	void draw(int rotation) {
		Vector<double> B(4);
		float theta, fi;
		float c = 8, a = 1;
		for (theta = 0; theta <= 2 * PI + 0.00001; theta += 0.6) {
			glBegin(GL_LINE_LOOP);
			for (fi = 0; fi <= 2 * PI + 0.00001; fi += 0.6) {
				B << (c + a*cos(theta))*cos(fi)
					<< a*sin(theta)
					<< (c + a*cos(theta))*sin(fi)
					<< 1.0;
				//Mátrix szorzási lista!
				B.Rotate(rotation, alpha);
				B.MoveTo(0, 9, 0);
				//B.MatrixMul(Camera, &B, &B);
				B.CVetites(100);
				B.Window2Viewport();
				//std::cout << B;
				glVertex3d(B.GetValue(0, 0) / B.GetValue(3, 0), B.GetValue(1, 0) / B.GetValue(3, 0), B.GetValue(2, 0) / B.GetValue(3, 0));
			}
			glEnd();
		}



		for (fi = 0; fi <= 2 * PI + 0.000001; fi += 0.6) {
			glBegin(GL_LINE_LOOP);
			for (theta = 0; theta <= 2 * PI + 0.0000001; theta += 0.6) {
				B << (c + a*cos(theta))*cos(fi)
					<< a*sin(theta)
					<< (c + a*cos(theta))*sin(fi)
					<< 1.0;

				B.Rotate(rotation, alpha);
				B.MoveTo(0, 9, 0);
			//	B.MatrixMul(Camera, &B,&B);
				B.Vetites();
				B.Window2Viewport();
				//std::cout << B;
				glVertex3d(B.GetValue(0, 0) / B.GetValue(3, 0), B.GetValue(1, 0) / B.GetValue(3, 0), B.GetValue(2, 0) / B.GetValue(3, 0));

			}
			glEnd();
		}

		alpha += 0.1;


	}
	void rotateX()
	{
		alphax += 0.1;

	}
	void rotateY()
	{
		alphay += 0.1;

	}
	void rotateZ()
	{
		alphaz += 0.1;

	}

	void setCameraMatrix(Matrix<double> C) {
		Camera = &C;
		//std::cout << "Camera Matrix X:" << Camera->GetX() << " Y:" << Camera->GetY() << std::endl;
	}
private:
	double PI = 3.14159265358979323846;
	double alphax, alphay, alphaz;
	double alpha;
	Matrix<double> *Camera;
};