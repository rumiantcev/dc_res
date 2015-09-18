// ---------------------------------------------------------------------------

#include <stdlib>
#include <system.hpp>
// #include <sysutils.hpp>
#include <assert.h>

#pragma hdrstop
#include "Net.h"

using namespace std;

// ------------------------------- operator () --------------------------------//
double __fastcall TNet:: operator()(long i, int j) {
	if (v != NULL)
		if ((v->v->linkCount > 1) && (!updated))
			update();
	return getIJ(i, j);
};
/* */
// ------------------------------- operator [] --------------------------------//
///*!inline*/  const double& __fastcall TNet::operator [](long i) const {f->v[i];};
// ------------------------------- getVecAt -----------------------------------//

Vector* __fastcall TNet::getVecAt(long i) {
	int j;
	// Vector result(Dim);
	checkCacheVector();

	if (!updated) {
		if (isVirtual)
			dynUpdate();
		else
			update();
	}
	for (j = 0; j < Dim; j++)
		cache->v->v[j] = getIJ(i, j);
	return cache;
}

/* */// ------------------------------default constructor -------------------------//
__fastcall TNet::TNet() {
	initNetDefault();
}

/* */// ------------------------------Parameters Initializatin----------------------//

void __fastcall TNet::initNetParams(const int& res, const int& _res) {
	long i;

	virtDim = Dim;
	halfRes = (double)Res / 2;
	dRes = 1 / (float)Res;

	SurfaceDim = Dim - 1;
	NumOfSur = Dim * 2;

	// (сетка сокращается только по одной коодинате)для двумерных сеток сразу считаем разрешение без одной точки. Хоть и небольшая экономия, но всё равно приятно
	NumOfPoints = _res;
	for (i = 2; i <= SurfaceDim; i++)
		NumOfPoints *= (res);
	/* */
	Count = NumOfPoints * NumOfSur;
}

/* */// ------------------------------constructor ----------------------------------//
__fastcall TNet::TNet(int dim, long perf, long res, bool virt)
	: perfomance(perf), Res(res), isVirtual(virt), Dim(dim) {
	// int i;
	long _res = res; // - 1;
	// Для экономии памяти, в случае если сетка реальная, число точек на грани можно сделать на одну меньше.

	initNetDefault(); // инициаизируем по умолчанию
	initNetParams(res, _res);
	buildPowerVectors(_res);
	// кешируем значения степеней для ускорения расчёта значений сетки
	if (!isVirtual)
		create(Dim /* , perf, res */);
	else
		v = NULL;
}

/* */// ------------------------------constructor ----------------------------------//
__fastcall TNet::TNet(long mm, long nn, long perf, long res) {
	initNetDefault(); // инициаизируем по умолчанию
	create(mm, nn, perf, res);
};

/* */// ------------------------copy constructor-----------------------------------//
__fastcall TNet::TNet(const TNet& Net) {
	initNetDefault();
	copyNetFrom(Net);
}

/* */// ------------------------------- destructor ---------------------------------//
/* void TNet:: operator delete(void *p) {
 TNet *ptr = (TNet*)p;

 if (ptr->v == NULL || ptr->isVirtual)
 delete(void*) p;
 else
 p = NULL;
 }

/* */// ---------------------------- destructor-------------------------------------//
__fastcall TNet::~TNet() {
	destroy();
}

/* */// ------------------------------ destroy -------------------------------------//
void __fastcall TNet::destroy() {
	if (v != NULL)
		if (--v->v->linkCount == 0) {
			delete v;
			v = NULL;
		}

	if (u_mx != NULL) {
		delete u_mx;
		u_mx = NULL;
	}
	if (u_vec != NULL) {
		delete u_vec;
		u_vec = NULL;
	}
	if (powVec != NULL) {
		delete powVec;
		powVec = NULL;
	}
	if (powVec_1 != NULL) {
		delete powVec_1;
		powVec_1 = NULL;
	}
	if (cache != NULL) {
		delete cache;
		cache = NULL;
	}

}
/* */// ------------------------------ MVMul ---------------------------------------//

void __fastcall TNet::MVMul(double* vsrc, double* vdst) {
	long i, j;
	double sum;
	if (vdst != NULL) {
		double* vs;
		vs = new double[u_mx->v->m];
		for (i = 0; i < u_mx->v->m; i++) {
			sum = 0;
			for (j = 0; j < u_mx->v->n; j++)
				sum += vsrc[j] * u_mx->v->v[i][j];
			vs[i] = sum;
		}
		for (i = 0; i < u_mx->v->m; i++)
			vdst[i] = vs[i];
		delete[]vs;
	}
}

/* */// ------------------------------ MVMul ---------------------------------------//

void __fastcall TNet::MVMul(Vector* vsrc, Vector* vdst) {
	long i, j;
	double sum, *temp;
	if (vdst->v->v != NULL) {
		double* vs;
		vs = new double[u_mx->v->m];
		for (i = 0; i < u_mx->v->m; i++) {
			sum = 0;
			for (j = 0; j < u_mx->v->n; j++)
				sum += vsrc->v->v[j] * u_mx->v->v[i][j];
			vs[i] = sum;
		}
		temp = vdst->v->v;
		delete[]temp;
		vdst->v->v = vs;
		vdst->v->size = u_mx->v->m;
	}
}

/* */// ------------------------------------Dynamically Update----------------------------------//

void __fastcall TNet::dynUpdate() {
	int j;
	if (cache == NULL) {
		cache = new Vector(Dim);
		cacheCurrent = -1;
	}
	if (upd != 1.0)
		for (j = 0; j < Dim; j++)
			(*cache)[j] *= upd;
	if (!umx) {
		if (!u_mx->updated)
			u_mx->update();
		MVMul(cache, cache);
		// umx = true;
	};
	if (!uvec) {
		for (j = 0; j < Dim; j++)
			(*cache)[j] += u_vec->v->v[j];
		// uvec = true;
	}
}

/* */// ------------------------------------Update----------------------------------//

void __fastcall TNet::update() {
	Vector c(v->v->n);
	Vector vec(u_mx->v->m);
	Matrix *vv;
	long i, j;
	long k;

	detach();
	if (upd != 1.0) {
		for (i = 0; i < Count; i++)
			for (j = 0; j < Dim; j++)
				v->v->v[i][j] = upd * v->v->v[i][j];
		upd = 1;
	}
	if (!umx) {
		if (u_mx->v->n == Dim) {
			for (k = 0; k < Count; k++)
				MVMul(v->v->v[k], v->v->v[k]);
		}
		else {
			vv = new Matrix(Count, u_mx->v->n);
			for (k = 0; k < vv->v->m; k++)
				MVMul(v->v->v[k], vv->v->v[k]);
			delete v;
			v = vv;
		};
		umx = true;
	};
	if (!uvec) {
		// MVMul(u_vec->v->v, u_vec->v->v);
		for (i = 0; i < Count; i++)
			for (j = 0; j < Dim; j++)
				v->v->v[i][j] += u_vec->v->v[j];
		uvec = true;
	}
	updated = true;
}

/* */// ------------------------------------Get current plate---------------------//

/* TODO -cОшибка :
 По-моему логика работы нарушена
 Надо проверить */
int __fastcall TNet::getCurrentPlate(long ID) {
	int i;
	long ind = 0;
	for (i = 0; i < 2 * Dim; i++) {
		ind += pow(Res, i * 0.5);
		if (ID < ind)
			return i - 1;
	}
	return -1;
}

/* */// ------------------------------------Cet coord vector-----------------------//
/* !inline */ /* Vector __fastcall TNet::getCoordVectorByID(long ID)
 {
 return parseCoordinate(ID);
 }
/* */// ------------------------------------Create----------------------------------//
void __fastcall TNet::create(long Dim /* long perf, long Res */) {

	long i, j, k, ind, curr, mm, md; // , ff;
	// long coordNum; // номер текущей координаты по которой идёт приращение
	long normalDir;
	double val;
	double step = 2.0 / double(Res /* -1/* */); // шаг приращения сетки
	virtDim = Dim;

	Vector exclude(Dim), baseV(Dim);

	// Count = NumOfPoints * NumOfSur;  //Общее количество точек
	v = new Matrix(Count, Dim); // Создаём матрицу координат
	// cout<<*v<<endl;

	for (j = 0; j < NumOfSur; j++) {
		ind = long(j * 0.5); // № текущей координаты
		md = (ind + 1) % Dim; // номер координаты по которой происходит экономия памяти
		(ind == 0 ? mm = 1 : mm = 0);  //с какой координаты начинаются изменения - если нулевая координата фиксирована, то с первой,
		//в противном случае  с нулевой
		normalDir = (j % 2) * 2 - 1;
		// направление нормали, а также какой конец будем исключать из построения сетки - верхний или нижний
		for (k = 0; k < Dim; k++) // заполняем базовый вектор -1.0
				baseV[k] = -1.0;
		if (normalDir > 0) {	// корректируем базовый вектор с учётом фиксированной грани и отсупа для экономии памати
			baseV[ind] = 1.0;
			baseV[md] += step;
			// если нормаль направлена в положительную сторону - исключаем начальные значения,
			// если в отрицательную - исключаем конечные значения по координате с +1 по модулю размерности индексом;
			// иначе в размеры матрицы можно не влезть
		}
		curr = (j * NumOfPoints); // первая точка текущей гиперплоскости
		// coordNum = 0;
		for (i = 0; i < NumOfPoints; i++) {
			// cout<<baseV<<endl;
			memcpy(v->v->v[curr], baseV.v->v, Dim*sizeof(LDouble));   //Копируем базовый вектор
			val = v->v->v[curr][mm] + step;    //приращиваем текущую координату
			if (i == NumOfPoints - 1)//если дошли до конца грани - выходим
				break;
			if (val > 1 || (val >= 1 && normalDir < 0)) {  //если координату прирастили до res - переходим к следующей
				mm++;
				if ((mm == ind) && (mm < Dim))
					// если чего - проскакиваем текущую кооддинату
						mm = (mm + 1) % Dim;
				for (k = 0; k < mm; k++) {    // и перестраиваем базовый вектор
					baseV.v->v[k] = -1.0;
				}
				if (normalDir > 0)
					// корректируем базовый вектор с учётом фиксированной грани и отсупа для экономии памати
				{
					baseV[ind] = 1.0;
					baseV[md] += step;
					// если нормаль направлена в положительную сторону - исключаем начальные значения, если в отрицательную - исключаем конечные значения по координате с +1 по модулю размерности индексом;
					// иначе в размеры матрицы можно не влезть
				}
				baseV.v->v[mm] += step;
			}
			else
				baseV.v->v[mm] = val; //корректируем базовый вектор (текущую координату) для использования при следующем копировании
			// (ind == 0 ? mm = 1 : mm = 0);
			curr++;
		}

	}

	/* for (j = 0; j < NumOfSur; j++) {
	 ind = long(j * 0.5); //№ текущей координаты
	 md = j % 2;  //??направление нормали
	 curr = (j * NumOfPoints);  //первая точка текущей гиперплоскости
	 p = 2 * nn * md - nn;
	 for (mm = 0; mm < Dim; mm++) {
	 exclude[mm] = nn + p;
	 p = -p;
	 }
	 cout<<exclude<<endl;
	 for (i = 1; i <= NumOfPoints; i++) {
	 if (j - ind * 2 == 0)              //заполнение 1 и -1 на концах, некорректно работает для прямоугольников
	 v->v->v[curr][ind] = -1.0;
	 else
	 v->v->v[curr][ind] = 1.0;

	 ff = i - 1;
	 mm = 0;
	 for (k = 0; k < Dim; k++) {
	 s = (ff % Res);
	 ff = ff / Res;
	 val = 2 / (double)(Res - 1) * s - 1;
	 if (k != ind)
	 val += exclude[++mm];
	 if (v->v->v[curr][k] == 0) {

	 v->v->v[curr][k] = val;
	 }
	 else {
	 if ((++k) < Dim) {
	 v->v->v[curr][k] = val;
	 }
	 }
	 }
	 curr++;
	 }
	 } /* */
}

/* */// ------------------------------------Create----------------------------------//
void __fastcall TNet::create(long mm, long nn, short perf, long res) {
	perfomance = perf;
	Res = res;
	halfRes = (double)Res / 2;
	v = new Matrix(mm, nn);
	Count = mm;
	Dim = nn;
	virtDim = Dim;
	isVirtual = false;
}

/* */// ------------------------------------Get i,j-----------------------------//
/*
 * Вычисляет значение находящееся в массиве сетки
 */
double __fastcall TNet::getIJ(long current, int coordNumber) {
	double res;
	if (isVirtual) {
		checkCacheVector();
		if (current != cacheCurrent) {
			cache = parseCoordinate(current);
			cacheCurrent = current;
			cache->norm(halfRes);
		}

		if (!updated)
			dynUpdate();
		res = cache->v->v[coordNumber];
	}
	else {
		if (!updated)
			update();
		res = v->v->v[current][coordNumber];
	}
	return res;
}

/* */// ------------------------------------Parse coord-----------------------------//
// protected:
Vector* __fastcall TNet::parseCoordinate(long current) {
	if (!isVirtual)
		return getVecAt(current);
	else
		return parseCoordinateForShift(current);
}

/* */// ------------------------------------Parse coord-----------------------------//
// protected:
Vector* __fastcall TNet::parseCoordinateForShift(long current) {
	// Кешируем, при необходимости, целый вектор
	int i;
	checkCacheVector();
	if (current != cacheCurrent) {
		cacheCurrent = current;

		_mod = current / NumOfPoints;
		// Номер координаты*2 , нормалью к которой будет построенный вектор - номер стороны гиперкуба
		__mod = _mod % 2; // Определяем знак нормали   (0,1)
		coordNumber = _mod * 0.5; // Определяем  номер координаты

		current -= _mod * NumOfPoints; // Вычитаем сдвиг от начала  стороны
		for (i = 0; i < coordNumber; i++) {
			cache->v->v[i] = current % Res;
			current *= dRes;
		}
		cache->v->v[coordNumber] = __mod * (Res - 1);
		// в текущей координате либо 0, либо res-1 до нормирования
		if (coordNumber < Dim)
			for (i = coordNumber + 1; i < Dim; i++) {
				cache->v->v[i] = current % Res;
				current *= dRes;
			}
	}
	return cache;
}

// public:
/* */// ------------------------------------Shift-----------------------------------//
/*
 * Осуществляет переход на соседнюю точку по поверхности сетки в заданном направлении
 */
long __fastcall TNet::shift(long current, int coordNumber, int step,
	bool& borderChanged) throw(exInvalidMoveDirection) {
	// Vector vv= *powVec,vv_1=*powVec_1;
	int currPlateNorm = current / (NumOfPoints * 2);
	int _dim = Dim == virtDim ? Dim : virtDim;
	long ind = current;
	long tOldCoordShift;
	// зациклить при попытке перехода за нижнюю или верхнюю границу массива
	if ((ind < 0)||(ind >= Count))
		return labs(ind) % Count;

	// вернуть ошибку при попытке перехода не по поверхности куба
	if ((currPlateNorm == coordNumber) && (Dim > 1))
		// throw(exInvalidMoveDirection());
			return current;
	// вернуть ошибку если при попытке перехода на не сущ. координату
	if (coordNumber < 0 || coordNumber > _dim)
		// throw(exInvalidMoveDirection());
			return current;

	Vector* c = parseCoordinateForShift(current);
	long posWithShift = step + c->v->v[coordNumber];
	if ((posWithShift < Res) && (posWithShift >= 0)) {
		if (coordNumber > currPlateNorm)
			tOldCoordShift = step * powVec->v->v[coordNumber - 1];
		else
			tOldCoordShift = step * powVec->v->v[coordNumber];

		minBorder = current - (current % NumOfPoints);
		maxBorder = minBorder + NumOfPoints - 1;
		ind = current + tOldCoordShift;
	   //	borderChanged = false;
	}
	else {
		minBorder = -1;
		maxBorder = minBorder;
		//borderChanged = true;
	}
	// Если не выходим за пределы текущей грани, то возврашаем индекс соотв. вектора
	if ((ind >= minBorder) && (ind <= maxBorder)) {
		return ind;
	}
	else {
		borderChanged = true;
		int newNormSign = (step / labs(step) + 1) * 0.5;
		// направление нормали к грани, на которую осуществляется сдвиг (0 - если отрицательный шаг, 1 - если шаг положителен) ;
		//новой нормалью становится та координата, по направлению которой происходил сдвиг
		// вычисляем нижнюю границу грани, на которую осуществляется переход переход
		long newMin = ((2 * coordNumber + newNormSign) /*% (2 * _dim)*/)
			* NumOfPoints;
		// нормаль к грани  - индекс соотв. направления сдвига
		int newPlateNorm = coordNumber;
		// нормаль к грани, на которую осуществляется переход

		ind = newMin;  //индекс - начало отсчёта - начальный индекс новой грани
		// Вычисляем точку, на которую осуществляем переход
		int i, j = _dim - 1;
		for (i = _dim - 1; i > newPlateNorm; --i) {
			ind +=  powVec_1->v->v[j];// * c->v->v[i];
			--j;
		}
		for (i = newPlateNorm - 1; i >= 0; --i) {
			ind +=  powVec_1->v->v[j];// * c->v->v[i];
			--j;
		}
		if (ind == Count) {
			ind--;
			// cout<<"err"<<endl;
			cout << current << ", " << coordNumber << ", " << step << endl;
			// cout<<(*this)<<endl;
		}
		return ind;
	};
}

/* */// ------------------------------------Detach----------------------------------//
void __fastcall TNet::detach() {
	long i, j;
	if (v != NULL)
		if (v->v->linkCount > 1) {
			v->v->linkCount--;
			Matrix *vv;
			vv = v;
			v = NULL;
			create(vv->v->m, vv->v->n, perfomance, Res);
			for (i = 0; i < Count; i++)
				for (j = 0; j < Dim; j++)
					v->v->v[i][j] = vv->v->v[i][j];
			v->v->linkCount = 1;
		};
	if (u_vec != NULL)
		u_vec->detach();
	if (u_mx != NULL)
		u_mx->detach();
	if (powVec != NULL)
		powVec->detach();
	if (powVec_1 != NULL)
		powVec_1->detach();
	if (cache != NULL)
		cache->detach();
};

/* */// -------------------------------- << ----------------------------------------//
ostream& __fastcall operator << (ostream& out_data, TNet& C) {
	long SurfaceDim = C.Dim - 1;
	long NumOfSur = C.Dim * 2;
	long NumOfPoints = C.Res;
	long i;
	int j;

	if (C.isVirtual)
		C.dynUpdate();
	else
		C.update();
	for (i = 2; i <= SurfaceDim; i++)
		NumOfPoints *= C.Res;

	for (i = 0; i < NumOfPoints * NumOfSur; i++) {
		out_data << "[";
		for (j = 0; j < C.Dim - 1; j++)
			out_data << C.getIJ(i, j) << ",";
		out_data << C.getIJ(i, j);
		out_data << "]";
		if (i == NumOfPoints*NumOfSur - 1)
			out_data << ";";
		out_data << endl;
	}
	return out_data;
};

/* */// ----------------------------------- >> -------------------------------------//
istream& __fastcall operator >> (istream& in_data, TNet& C) {

	if (C.v == NULL)
		C.v = new Matrix(C.Count, C.Res);
	in_data >> *C.v;
	return in_data;
};
/* */// ----------------------------------- = --------------------------------------//
TNet& __fastcall TNet:: operator = (const TNet & Net) {
	if (this == &Net)
		return*this;

	if (!isVirtual)
		if (--v->v->linkCount == 0) {
			delete v;
			v = NULL;
		}
	copyNetFrom(Net);
	return *this;
};

/* */// ----------------------------------- * --------------------------------------//
TNet& __fastcall TNet:: operator *= (const Matrix& A) {
	bool isDimChanged = (Dim != A.m()) ? true : false;
	if (!uvec)
		update();
	if (u_mx != NULL)
		* u_mx = *u_mx * A;
	else
		*u_mx = A;
	Dim = A.m();
	if (isDimChanged)
		buildPowerVectors(Dim);

	if (!isVirtual) {
	}
	updated = false;
	umx = false;
	return *this;
}

/* */// ----------------------------------- * --------------------------------------//
const TNet __fastcall operator*(const Matrix &A, const TNet& B) {
	TNet result = B;
	return TNet(B) *= A;
};

/* */// ----------------------------------- * --------------------------------------//
/* TNet& __fastcall TNet:: operator *= (const double &a) {
 upd *= a;
 updated = false;
 return *this;
 }

/* */// ----------------------------------- * --------------------------------------//
const TNet __fastcall operator *(const double& a, const TNet& B) {
	return TNet(B) *= a;
};

/* */// ----------------------------------- + --------------------------------------//
TNet& __fastcall TNet:: operator += (const TNet& B) {
	long i, j;
	TNet *pB = (TNet*)&B;
	double coeff = 1 / upd;
	if (isVirtual) // если сеть была виртуальная, то создаем физическую
			create(Count, Dim, perfomance, Res);
	if (v == B.v) {
		upd += B.upd;
	}
	else {
		// if(updated) update();
		if (!B.updated)
			pB->update();
		// for(i=0;i<A.m;i++)
		// result.f->v[i]=A.f->v[i]+B.f->v[i];
		for (i = 0; i < Count; i++)
			for (j = 0; j < Dim; j++)
				v->v->v[i][j] += coeff * pB->v->v->v[i][j];
	}
	updated = false;
	return *this;
}

/* */// ----------------------------------- + --------------------------------------//
const TNet __fastcall operator +(const TNet& A, const TNet& B) {
	return TNet(A) += B;
};

/* */// ----------------------------------- + --------------------------------------//
TNet& __fastcall TNet:: operator += (const Vector& A) {
	if (umx)
		update();
	if (u_vec != NULL) {
		*u_vec = *u_vec + A;
	}
	else {
		u_vec = new Vector(A);
		// *u_vec=A;
	};
	updated = false;
	uvec = false;
	return *this;
}

/* */// ----------------------------------- + --------------------------------------//
const TNet __fastcall operator +(const Vector& A, const TNet& B) {
	return TNet(B) += A;
};

/* */// ------------------------------ Clear ---------------------------------------//

void __fastcall TNet::Clear() {
	long i, j;
	if (!isVirtual)
		for (i = 0; i < Count; i++)
			for (j = 0; j < Dim; j++)
				v->v->v[i][j] = 0;
}

/* */// ------------------------------inintialize pointers -------------------------//
void __fastcall TNet::initNetDefault() {
	updated = true;
	uvec = true;
	umx = true;
	upd = 1.0;
	u_vec = NULL;
	u_mx = NULL;
	cache = NULL;
	cacheCurrent = -1;
	// parsedCache=NULL;
	powVec = NULL;
	powVec_1 = NULL;
	v = NULL;
}

/* */// ------------------------------copy variables -------------------------//
void __fastcall TNet::copyNetFrom(const TNet& Net) {
	isVirtual = Net.isVirtual;
	if ((Net.v != NULL) && (!isVirtual))
		v = Matrix::copy(Net.v, v);
	Res = Net.Res;
	Dim = Net.Dim;
	virtDim = Net.virtDim;
	halfRes = Net.halfRes;
	dRes = Net.dRes;
	NumOfPoints = Net.NumOfPoints;
	Count = Net.Count;
	perfomance = Net.perfomance;
	updated = Net.updated;
	upd = Net.upd;
	cached = Net.cached;
	halfRes = Net.halfRes;
	cacheCurrent = Net.cacheCurrent;
	cache = Vector::copy(Net.cache, cache);

	uvec = Net.uvec;
	u_vec = Vector::copy(Net.u_vec, u_vec);

	umx = Net.umx;
	u_mx = Matrix::copy(Net.u_mx, u_mx);

	powVec = Vector::copy(Net.powVec, powVec);
	powVec_1 = Vector::copy(Net.powVec_1, powVec_1);
}

/* */// ----------------------build vectors with power values----------------------------------------------//
void __fastcall TNet::buildPowerVectors(int _res) {
	int i, j = 0;
	if (powVec != NULL)
		delete powVec;
	if (powVec_1 != NULL)
		delete powVec_1;
	powVec = new Vector(Dim);
	powVec_1 = new Vector(Dim);

	(*powVec_1)[0] = 1;
	for (i = 1; i < Dim; i++) {
		(*powVec_1)[i] = pow((LDouble)_res, j);
		(*powVec)[j++] = (*powVec_1)[i];
	}
	(*powVec)[j] = pow((float)_res, j);
}

/* */// ----------------------Checks cached Vector----------------------------------------------//
void __fastcall TNet::checkCacheVector() {
	if (cache == NULL) {
		cache = new Vector(Dim);
		cacheCurrent = -1;
	}
	else if (cache->size() != Dim) {
		delete cache;
		cache = new Vector(Dim);
		cacheCurrent = -1;
	}
}
/* */

/* */// ----------------------метка недействительности кеша----------------------------------------------//
void TNet::wasteCache() {
    cacheCurrent = -1;
}
