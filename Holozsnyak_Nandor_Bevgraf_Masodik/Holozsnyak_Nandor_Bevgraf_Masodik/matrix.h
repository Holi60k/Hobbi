//  Copyright 2015 Nándor Kristóf Holozsnyák
//  Matrix class made by Holi60k
//  Just for fun
//  Nevermind if it has some bugs, they may be fixed :)

#include <iostream>
#include <limits>
#include <cmath>

template<typename T>
class Matrix {
 public:
	// Konstuktor
	Matrix(int a = 0, int b = 0):siX(a), siY(b) {
		
		//std::cout << "Matrix ctor" << std::endl;
		// inicializálás itt történik, memóriában foglalunk helyet ennek a sok szép számnak
		// beállítjuk a mutatónkat hogy mutasson egy mutatóra... igen ez most így lesz
		Matx = new T*[a];
		for(int i = 0; i < a; i++) {
			// aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] = new T[b];
		}
		Reset();
		if(siX == siY) {

			Squared = true;
		} else {
			Squared = false;
		}
		
		// determináns előjele
		Det_Sign = 1;
	}
	// Destruktor
	~Matrix() {
		// destruktor törli a sok szép számunkat
		// logikus hogy előbb azokat az adatokat töröljük amiket legutoljára raktunk a memóriába
		// pointer to pointer, A to B... B(k) törlése előbb a memóriából aztán pedig az A-t
		for (int i = 0; i < this->GetX(); i++) {
			try
			{
				//if(Matx[i] != NULL)
					//delete[] Matx[i];
			}
			catch (const std::exception&e)
			{
				std::cout << e.what() << std::endl;
			}
			

		}
		if(Matx != NULL)
		delete [] Matx;

	}
	void Destroy() {
		this->~Matrix();
	}
	// mozgató szemantika
	Matrix (Matrix<T> && A) {

		//std::cout << "move ctor" << std::endl;
		Matx = new T*[A.GetX()];
		for(int i = 0; i < A.GetX(); i++) {
			// aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] = A.Matx[i];
			A.Matx[i] = nullptr;
		}
		siX = A.GetX();
		siY = A.GetY();
		Squared = A.GetSquared();
		Det_Sign = A.Get_Det_Sign();
	}
	// mozgató értékadá
	Matrix & operator= (Matrix<T> && A) {

		//std::cout << "move assingment" << std::endl;
		siX = A.GetX();
		siY = A.GetY();
		Squared = A.GetSquared();
		Det_Sign = A.Get_Det_Sign();
		// egy előre lefoglalt kibaszott 0 terület törlése :'(
		//delete [] Matx;


		Matx = new T*[siX];
		for (int i = 0; i < siX; i++) {
			// aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] =  A.Matx[i];
			A.Matx[i] = NULL;
		}
		return *this;

	}

	// értékadó operátorunk, másoló értékadás!
	Matrix & operator= (Matrix & A) {

		// std::cout << "copy assingment" << std::endl;
		siX = A.siX;
		siY = A.siY;
		// egy előre lefoglalt kibaszott 0 terület törlése :'(
		delete [] Matx;


		Matx = new T*[siX];
		for (int i = 0; i < siX; i++) {
			// aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] = new T[siY];
		}

		for (int i = 0; i < A.GetX();i++) {
			for(int j = 0; j < A.GetY();j++) {
				this->FillMatrix(A.GetValue(i,j),i,j);
			}
		}

		Squared = A.Squared;
		return *this;

	}

	Matrix<T> & operator=(const Matrix<T> & A) {
		// std::cout << "copy assingment" << std::endl;
		siX = A.siX;
		siY = A.siY;
		// egy előre lefoglalt kibaszott 0 terület törlése :'(
		//delete [] Matx;


		Matx = new T*[siX];
		for(int i = 0; i < siX; i++) {
			// aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] = new T[siY];
		}

		for(int i = 0; i < A.GetX();i++) {
			for(int j = 0; j < A.GetY();j++) {
				this->FillMatrix(A.GetValue(i,j),i,j);
			}
		}

		Squared = A.Squared;
		return *this;
	}
	// Másoló konstuktor
	Matrix (const Matrix<T> & t) {

		// std::cout << "copy ctor " << std::endl;
		siX = t.GetX();
		siY = t.GetY();

		Matx = new T*[siX];
		// std::cout << "Old: " << &t.Matx << std::endl;
		// std::cout << "New: " << &Matx << std::endl;
		for(int i = 0; i < siX; i++) {
			// aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] = new T[siY];
			// std::cout << "Old: "<< i << " - " << &t.Matx[i] << std::endl;
			// std::cout << "New: "<< i << " - " << &Matx[i] << std::endl;
		}

		for (int i = 0; i < t.GetX();i++) {
			for (int j = 0; j < t.GetY();j++) {
				FillMatrix(t.GetValue(i,j),i,j);
			}
		}
		// GetWholeMatrix();
		Squared = t.Squared;
		Det_Sign = t.Det_Sign;

	}
	
	// Matrix (const Matrix&);
	// két mátrix összeadására szolgáló operátorunk
	Matrix operator+ (const Matrix & t) {

		if(this->GetY() == t.GetX()) {

			Matrix<T> Result(this->GetX(),t.GetY());

			for(int i = 0; i < Result.GetX();i++) {
				for(int j = 0; j < Result.GetY(); j++) {				
					Result.FillMatrix(this->GetValue(i,j)+t.GetValue(i,j),i,j);
				}

			}
			return Result;
		} else {
			std::cout << "Sorry I can not add these two Matrixes..."<< std::endl;

		}
		return *this;
	}

	Matrix operator* (const Matrix & A) {
		T Var = 0;
		//Matrix<T> Result(this->GetX(), A.GetY());
		
		if (this->GetY() == A.GetX()) {
			
			Matrix<T> *Result = new Matrix<T>(this->GetX(), A.GetY());
			for (int i = 0; i < Result->GetY();i++) {
				for (int j = 0; j < Result->GetX(); j++) {

					for (int l = 0; l < this->GetY();l++) {
						Var += this->GetValue(j, l) * A.GetValue(l, i);
					}
					Result->FillMatrix(Var, j, i);
					Var = 0;

				}

			}
			return *Result;
		}
		else {
			std::cout << "Sorry I can not multiplicate these two Matrixes..." << std::endl;

		}
		return 0;
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


	Matrix & operator*= (const T & B) {
		

		for(int i = 0; i < this->siX; i++) {
			for(int j = 0; j < this->siY; j++) {
				this->Matx[i][j] = this->Matx[i][j]*B;
			}
		}

		return *this;
	}
	void ResetToNull() {
		
		for (int i = 0; i < siX; i++) {
			// aztán pedig egyesével lefoglalunk a mutatók számára dinamikus memóriát
			Matx[i] = new T[1];
			Matx[i] = NULL;
		}
		Matx = NULL;

	}
	bool GetSquared() const {
		return Squared;
	}

	// Mátrix feltöltése, v - érték, x - x koordináta, y - y koordináta
	void FillMatrix(T v, int x, int y) const {
		Matx[x][y] = v;		
	}
	// Visszaadja a mátrix sorainak számát
	int GetX() const {
		return siX;
	}
	// Visszadja a mátrix oszlopainak számát
	int GetY() const {
		return siY;
	}
	// visszaadja a mátrix x.sorának y.edik elemét
	T GetValue(int x, int y) const {
		return Matx[x][y];
	}
	// két sor cseréje, elemi sorművelet ugyebár amit majd a Gauss(-Jordan) elimináció során fogunk haszálni
	// természetesen ezek csak mutató cserék mivel a sorainak egyben mutatók ahova tároljuk az eleminket
	// ennél a függvénynél található egy Det_Sign *= -1 is amivel már a determinánsunk előjelét is tudjuk váltani,
	//  mivel definíció szerint ha sort vagy oszlopot cserélünk akkor a determináns értéke a -1-szeresére változik.
	void ChangeRows(int a, int b) {
		T *cs;
		cs = Matx[a];
		Matx[a] = Matx[b];
		Matx[b] = cs;
		Det_Sign *= -1;
		
		
	}
	// Visszaadja a mátrix i.edik sorát, igazából nem volt használva de jó hogyha van egy ilyenünk is
	T* GetRows(int i) { 
		// std::cout<< i << " - " << Matx[i] << std::endl;
		return Matx[i];
	}
	// visszadja a determinánsunk előjelét
	int Get_Det_Sign() {
		return Det_Sign;
	}
	// Igaz hogy GaussElimination néven van a függvény de ez lesz nekünk a Gauss-Jordan eliminálónk, mivel felül illetve alul is eliminálunk
	// Működése elég egyértelmű de rövid mondatokban
	// 1 - a vizsgált sort cseréljük az első sorral(ha az első persze akkor is) és beosztjuk mindig a sorunkat az aktuális i.edik elemmel ez mindig az aktuálisan feldolgozandó sor i.edik eleme
	// 2 - ezután ennek a sornak az X-szeresét (x = a következő sorok i.edik pozíciójában álló eleme) az alatta levő sorokhoz
	// 3 - visszacseréljük a két sort és újraindul az algoritmus
	// végülis Gauss-Jordan elimináció csak itt nem tűnik fel az alul-felül való elimináció mert mindig sort cserélünk és csak lefelé eliminálunk :)
	void GaussElimination() {
		//std::cout << "Gauss-Jordan" << std::endl;

		for(int i = 0; i < siX; i++) {
			ChangeRows(0,i);
			MultiplicateRow(0,Matx[0][i] != 0?1/Matx[0][i]:1);
			for(int j = 1; j < siX; j++) {

				AddRow2Row(j,0,(-1)*Matx[j][i]);
			}
			ChangeRows(0,i);
		}

	}
	// Ugyan az mint a rendes GaussElimination() fgv,
	// de ez inkább felső 3-szög alakra hozza a mátrixot :)
	void GaussElimination_2() {
		// std::cout << "Gauss_2" << std::endl;
		T pivot;
		int csI;

		for (int i = 0; i < siX; i++) {

			for(int j = i; j < siY; j++) {

				csI = Abs_Max(j);
				ChangeRows(j,csI);


			}
			//std::cout << *this;
			//if (WhereIsItsRow(i) == i){
				pivot = Matx[i][i] != 0?1/Matx[i][i]:1;
				MultiplicateRow(i,pivot);

				for (int j = i+1; j < siX; j++) {

					AddRow2Row(j,i,(-1)*Matx[j][i]);
				}

				MultiplicateRow(i,1/pivot);
			//} else {
			//	ChangeRows(i,WhereIsItsRow(i));
			//	GaussElimination_2();
			//}
		}


	}
	// Egy sor value-szeresét hozzáadjuk egy másik sor-hoz
	// a-hoz adjuk hozzá b value szeresét
	void AddRow2Row(int a, int b, T value) {

		for(int i = 0; i < siY; i++) {

			Matx[a][i] += value * Matx[b][i];
		}

	}
	// Nos igen, az eddigi eliminációk igazából csak homogén rendszerekre volt jó, ennél a függvénynél hozzádobunk egy "vektort" is
	// De a működés itt is természetesen a Gauss-Jordan...
	// Szabad paraméterek is működnek elviekben de még nem jól... :(
	void GaussEliminationWVector(Matrix & b) {
		//Matrix<T> Result = *this;
		//std::cout << "Gauss_W_Vector" << std::endl;
		T pivot;
		int Free_Param = 0;
		int csI = 0;
		for(int i = 0; i < siX; i++) {

			for(int j = i; j < siY; j++) {
				csI = Abs_Max(j);
				ChangeRows(j,csI);
				b.ChangeRows(j,csI);
			}

			ChangeRows(0,i);
			b.ChangeRows(0,i);
			pivot = Matx[0][i] != 0?1/Matx[0][i]:1;
			
			MultiplicateRow(0,pivot);

			b.MultiplicateRow(0,pivot);

			for(int j = 1; j < siX; j++) {
				
				pivot = (-1)* Matx[j][i];
				
				b.AddRow2Row(j,0,pivot);
				AddRow2Row(j,0,pivot);
			}

         	ChangeRows(0,i);
			b.ChangeRows(0,i);
		}
		//Result.GaussElimination_2();
		//Free_Param = siY - Result.Rank();
		//std::cout << "Free Parameters: " << Free_Param << std::endl;

	}

	int Abs_Max(int oIndex) {
		int max = 0;
		int index;
		for(int i = 0; i < GetY(); i++) {
			if(std::abs(Matx[i][oIndex]) > max) {
				max = Matx[i][oIndex];
				index = i;
			}
		}
		return index;

	}
	// Determináns számítás
	float Determinant() {

		// std::cout << "Determinant" << std::endl;
		double Det = 1;
		Matrix<T> Result = *this;
		if(Squared) {

			/*if(!Result.IsTriangle()) {
				if(Result.CanIReOrderToTriangle()) 
					Result.ReOrderTriangle();
			}*/
			
			Result.GaussElimination_2();
			std::cout << Result;
			for(int i = 0; i < siX; i++) {
				Det *= Result.Matx[i][i];
			  // std::cout << Result.Matx[i][i] << std::endl;
				//std::cout << Det;
			}
			//std::cout <<"Det:" << Det*Result.Get_Det_Sign() << std::endl;
			return Det*Result.Get_Det_Sign();

		} else {
			std::cout << "Sorry I can not count its determinant cause this matrix is non-squared." << std::endl;
			return -1;
		}
	}
	// Rang számítás... még nem tökéletes a koncepció
	int Rank() {

		Matrix<T> Result = *this;
		int rank = 0;
		Result.GaussElimination();

		for(int i = 0; i < siX; i++) {

			if(Result.CountZeroItems(i) == 0) rank++;
			// std::cout << "R: " << rank << std::endl;
		}

		return rank;

	}
	// megszámolja egy sorban lévő 0 elemeket, rangszámításhoz kell...
	int CountZeroItems(int a) {
		int Zero = 0;

		if(Matx[a][a] != 0 && Squared) {
			return 0;
		}

		for(int i = 0; i < siY; i++) {
			if(Matx[a][i] == 0) Zero++;
		}
		return Zero;
	}
	/*
	// Majdnem jó lett de mégsem...
	void ReOrderRows() { 
		int Zero_Num[siX] = {0};

		for(int i = 0; i < siX; i++) {
		Zero_Num[i] = 0;
			for(int j= 0; j < siY; j++) {
				if(Matx[i][j] == 0)
					Zero_Num[i]++;
			}
		}
	}
	*/
	// A következő pár függvény a felső háromszög alakra hozásban játszik szerepet
	// IsFineRow fügvénnyel meg tudjuk nézni hogy az aktuális sorunk jó helyen van-e, úgymond megfelelő a nullák száma és poziciója
	bool IsFineRow(int a) {

		int Zero_Num = 0;

		for(int i = 0; i < a; i++) {
			if(Matx[a][i] == 0)
				Zero_Num++;
		}
		if(Zero_Num == a)
			return true;
		else
			return false;

	}
	// Hol kellene hogy az aktuálisan kapott sor legyen, hanyadik sorban
	int WhereIsItsRow(int a) {

		int Zero_Num = 0;
		for(int j = 0; j < a; j++) {

			if(Matx[a][j] == 0) {
				Zero_Num++;
			}
		}
		return Zero_Num;

	}
	// Rendezzük át a kis mátrixunkat, ez akkor lényeges ha már eleve egy olyan mátrixot kapunk amit egyből felső háromszög alakra tudunk hozni elimináció nélkül
	void ReOrderTriangle() {
		//  0 0 1
		//  1 0 1
		//  0 1 1

		//  1 0 1
		//  0 1 1
		//  0 0 1

		int Places[siX];
		int Changes = 0;
		for(int i = 0; i < siX; i++) {
			Places[i] = i;
			if(!IsFineRow(i)) {
				ChangeRows(i,WhereIsItsRow(i));
				Places[i] = WhereIsItsRow(i);
			} 
			
		}
	}
	// Már felső háromszög alakú a mátrix?
	bool IsTriangle() {

		int Zero_Num[siX];
		bool Row_OK[siX];
		int Rows = 0;

		for(int i = 0; i < siX; i++) {
			Row_OK[i] = true;
			for(int j = 0; j < i; j++) {
				if(Matx[i][j] == 0) {
					Row_OK[i] = true;
				} else {
					Row_OK[i] = false;
					break;
				}
			}
		}

		for(int i = 0; i < siX; i++){
			if(Row_OK[i] == true) {
				Rows++;
				
			}
		}

		if(Rows == siX)
			return true;
		else
			return false;

	}
	// Illetve hogy tudunk-e megfelelő felső háromszög alakra rendezni... 
	bool CanIReOrderToTriangle() {
		int Zero_Num[siX];
		
		int Rows = 0;
		int a = 0;
		for(int i = 0; i < siX; i++) {
			Zero_Num[i] = 0;
			while(Matx[i][a] == 0 && a < siY) {

				Zero_Num[i]++;
				a++;
			}
			
			a = 0;
			
		}

		for(int i = 0; i < siX; i++) {

			if(Zero_Num[i] > 0){
				Rows++;
			}

		}
		if(Rows == siX-1)
			return true;
		else
			return false;

	}
	// a mátrix egy sorának szorzása egy megfelelő számmal
	void MultiplicateRow(int a, T Multi) {

		for(int i = 0; i < siY; i++) {
			// std::cout << "SubRow:(in) " << i  << " - " << Matx[a][i] << std::endl;
			Matx[a][i] *= Multi;

		}
	}
	// Kiiratja mátrixunkat
	void GetWholeMatrix() const {
		std::cout << "---" << std::endl;
		// std::cout << "Get the whole Martix..." << std::endl;
		for(int i = 0; i < siX; i++)
		{
			for(int j = 0; j < siY; j++)
			{
				// Matx[i][j] = v;
				std::cout << Matx[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "---" << std::endl;

	}
	// Kiiratja a mátrixunkat illetve a paraméterül kapott vektorunkat is 
	void GetWholeMatrixVector(Matrix & b) const {
		// std::cout << "Get the whole Martix..." << std::endl;

		for(int i = 0; i < siX; i++)
		{
			for(int j = 0; j < siY; j++)
			{
				// Matx[i][j] = v;
				std::cout.precision(15);
				std::cout << Matx[i][j] << " ";
			}
			std::cout.precision(15);
			std::cout << b.Matx[i][0];
			std::cout << std::endl;
		}

	}
	// ???? lényegtelen még...
	Matrix Zero_Vector() {

	}

	Matrix Transparent() {

		Matrix<T> Transparent(siX,siY);
		for(int i = 0; i < siX; i++) {
			for(int j = 0; j < siY; j++) {
				//fordított sorrendben feltöltjük a mátrixot
				Transparent.FillMatrix(GetValue(j,i),i,j);
			}
		}
		return Transparent;
	}

	void LU() {

	}

	void Cholesky() {
		//A = LU
		//L mátrix alsó háromszög alakú, főátlójában csak egyes van
		//azok alatt pedig a pivot elemek (?!)
		//U = felső háromszög alakú mátrix ami a gauss elimináció után marad
		//U főátlóbeli elemeit diagonális mátrixba gyűjtve D
		//A = LDU
		//A = LL^T

		Matrix<T> L(siX,siY),U(siX,siY),D(siX,siY),Lv(siX,siY);
		L.Make_Identity();
		std::cout << L;
		U.Reset();
		T pivot;
		

		for (int i = 0; i < siX; i++) {
			
			pivot = Matx[i][i] != 0?1/Matx[i][i]:1;
			for(int j = i;j < siX;j++) {
				L.Matx[j][i] = Matx[i][j]*pivot;
			}

			MultiplicateRow(i,pivot);

			for (int j = i+1; j < siX; j++) {

				AddRow2Row(j,i,(-1)*Matx[j][i]);
			}

			MultiplicateRow(i,1/pivot);
	    }

	    U = *this;
		std::cout << L;
		std::cout << U;

		D.Reset();
		for(int i = 0; i < siX; i++) {
			D.FillMatrix(std::sqrt(U.GetValue(i,i)),i,i);
		}
		Lv = L*D;
		std::cout << D;
		D = Lv*Lv.Transparent();
		std::cout << D;
	}

	void Cholesky_Scratch(Matrix & b) {

		Matrix<T> y(siX,1);
		y.Reset();
		//std::cout << b;
		std::cout << "Cholesky start.." << std::endl;
		for(int K = 0; K < siX; K++) {

			if(Matx[K][K] > 10e-15) {
				Matx[K][K] = sqrt(Matx[K][K]);

			}

			for(int i = K+1; i < siX; i++) {
				Matx[i][K] = Matx[i][K] / Matx[K][K];

				for(int j = K+1; j <=i; j++) {
					Matx[i][j] = Matx[i][j] - Matx[i][K]*Matx[j][K];
				}
			}
		}
		double sum = 0;
		for(int i = 0; i < siX; i++) {
			sum = 0;
			for(int j = 0; j < i; j++) {
				sum += Matx[i][j]*b.Matx[j][0];
			}
			
			
			b.Matx[i][0] = (b.Matx[i][0] - sum)/Matx[i][i];
	
		}


		std::cout << b;

		for(int i = 0; i < siX; i++) {
			sum = 0;
			for(int j = i+1; j < siX; j++) {
				sum += Matx[i][j]*b.Matx[j][0];
			}
			
			
			b.Matx[i][0] = (b.Matx[i][0] - sum)/Matx[i][i];
	
		}
		std::cout << b;
		std::cout << "Cholesky finnish." << std::endl;

	}

	void Reset() {
		for(int i = 0; i < siX; i++) {
			for(int j = 0; j < siY; j++) {
				FillMatrix(0,i,j);
			}
		}
	}

	void Make_Identity() {
		if(Squared) {
			Reset();
			for(int i = 0; i < siX; i++) {
				FillMatrix(1,i,i);
			}
		}
	}

	

	Matrix Inverz_Matrix() {
		Matrix<T> ID(siX, siY);

		//O = *this;
		ID.Make_Identity();
		T pivot;
		int csI = 0;
		
		for(int i = 0; i < this->GetX(); i++) {

			//for (int j = i; j < this->GetY(); j++) {
			//	csI = Abs_Max(j);
				ChangeRows(i, csI);
				ID.ChangeRows(i, csI);
			//}

			ChangeRows(0,i);
			ID.ChangeRows(0,i);
			pivot = GetValue(0,i)!= 0?1/GetValue(0,i):1;
			//pivot = 1 / O.GetValue(0, i);
			
			MultiplicateRow(0,pivot);

			ID.MultiplicateRow(0,pivot);

			for(int j = 1; j < siX; j++) {
				
				pivot = (-1)* GetValue(j,i);
				
				ID.AddRow2Row(j,0,pivot);
				AddRow2Row(j,0,pivot);
			}

			ChangeRows(0,i);
			ID.ChangeRows(0,i);

			
		}
	return ID;
	}
	// kimeneti operátor << túlterhelése
	friend std::ostream & operator<< ( std::ostream & os, Matrix & A ) {
		A.GetWholeMatrix();
		return os;
	}
	Matrix & operator<<(T value) {
		static int i = 0, j = 0;
		FillMatrix(value,i,j++);
		if (j == siY) { j = 0; i++;}
		if (i == siX) i = 0;
		return *this;
	}

	void RotateX(double alpha) {
		Matrix<T> Forgat(3, 3);
		Forgat << cos(alpha) << -sin(alpha) << 0;
		Forgat << sin(alpha) << cos(alpha) << 0;
		Forgat << 0 << 0 << 1;
		//*this = Forgat*(*this);
		MatrixMul(&Forgat, this, this);

	}

	void moveTo(double x, double y) {
		Matrix<T> Eltolas(3, 3);
		Eltolas << 1 << 0 << x;
		Eltolas << 0 << 1 << y;
		Eltolas << 0 << 0 << 1;
		MatrixMul(&Eltolas, this, this);
	}

protected:
	// egy pár privát tag.
	int siX;
	int siY;
	bool Squared;
	bool Null_Matrix;
	int Det_Sign = 1;
	T** Matx = NULL;
};


template <typename T>
class Vector:public Matrix<T> {

public:
	Vector(int x):Matrix<T>(x,1) {
		//std::cout << "Vector ctor " << std::endl;
		this->siX = x;
		this->siY = 1;
	};
	~Vector() {
	//std::cout << "Vector dtor" << std::endl;
	};

	void GetVector() const {
		for(int i = 0; i < this->siX; i++)
		{
			std::cout.precision(8);
			std::cout << std::fixed << this->Matx[i][0] << " ";
		}
		std::cout << std::endl;

	}

	friend std::ostream & operator<< ( std::ostream & os, Vector & A ) {
		A.GetVector();
		return os;
	}


	
private:

};