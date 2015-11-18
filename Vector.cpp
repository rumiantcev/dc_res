// ---------------------------------------------------------------------------

#include <assert.h>
#include <math.h>
#include <algorithm>
// #include <dstring.h>
#pragma hdrstop
#include "Vector.h"

using namespace std;

// ------------------------------- size () ------------------------------------//
long __fastcall Vector::size() const {
	return v->size;
}

// ------------------------------- destructor ---------------------------------//
/* void Vector::operator delete(void *p)
 {
 Vector *ptr = (Vector*)p;

 if (ptr->v==NULL)
 delete (void*) p;
 else
 p=NULL;
 } */
// ------------------------------- destructor ---------------------------------//
__fastcall Vector::~Vector() {
		delete v;
		v = NULL;
}

// ------------------------------- constructor --------------------------------//
__fastcall Vector::Vector(long sz) {
	create(sz);
}

// ------------------------------- constructor --------------------------------//
__fastcall Vector::Vector() {
	create(0);
}

// ------------------------------- copy constr --------------------------------//
__fastcall Vector::Vector(const Vector &C) : v(C.v), upd(C.upd),
	updated(C.updated) {
	v->linkCount++;
}

// ------------------------------- copy constr --------------------------------//
__fastcall Vector::Vector(const double* vv, long sz) : upd(1.0), updated(true) {
	v = new sVec(vv, sz);
}

// ------------------------------- string reader --------------------------------//
__fastcall Vector::Vector(const string& str, long size) {
	long k;
	string::const_iterator i, j;
	create(size);
	i = find(str.begin(), str.end(), '[');
	j = ++i;
	string _val;
	for (k = 0; k < size; k++) {
		while (*j != ',')
			j++;
		_val.assign(i, j);
		v->v[k] = atof(_val.c_str());
		i = ++j;
	}
}
// ------------------------------- = ------------------------------------------//
Vector& __fastcall Vector:: operator = (const Vector & C) {
	if (this == &C)
		return*this;
	delete v;
	v = C.v;
	upd = C.upd;
	updated = C.updated;
	v->linkCount++;
	return *this;
} /* */

// ------------------------------- = ------------------------------------------//
Vector& __fastcall Vector:: operator = (const double & c) {
	long i;

	for (i = 0; i < v->size; i++)
		v->v[i] = c;
	return *this;
}

// ---------------------------- += --------------------------------------------//
Vector& __fastcall Vector:: operator += (const Vector& A) {
	long i;
	double coeff = A.upd / upd;
	if (A.v == v)
		upd += A.upd;
	else {
		if (v->linkCount > 1)
			detach();
		for (i = 0; i < A.v->size; i++)
			v->v[i] += coeff * A.v->v[i];
	}
	updated = false;
	return *this;
}

Vector& __fastcall Vector::vSum(const Vector& A) // added for profiling
{
	long i;
	if (A.v == v) {
		upd += A.upd;
	}
	else {
		if (v->linkCount > 1)
			detach();
		for (i = 0; i < A.v->size; i++)
			v->v[i] += A.v->v[i];
	}
	updated = false;
	return *this;
}

// ---------------------------- + ---------------------------------------------//
const Vector __fastcall operator +(const Vector& A, const Vector &B) {
	return (Vector)A += B;
}

// ----------------------------- * -------------------------------------------//
const double __fastcall operator *(const Vector& A, const Vector &B) {
	long i;
	double res = 0;
	// assert(A.v->size==B.v->size);
	for (i = 0; i < A.v->size; i++)
		res += A.v->v[i] * B.v->v[i];
	res *= B.upd * A.upd;
	return res;

}

const double __fastcall scMul(const Vector& A, const Vector &B)
{
	long i;
	double res = 0;
	for (i = 0; i < A.v->size; i++)
		res += A.v->v[i] * B.v->v[i];
	return res;
}

// ------------------------------ * -------------------------------------------//
Vector& __fastcall Vector:: operator *= (const double& scalar) {
	// !if (v->linkCount>1) detach();
	/* DONE -orum -cCheck : ѕроверить корректность оптимизации. */
	updated = false;
	upd *= scalar;
	return *this;
}

// ------------------------------ * -------------------------------------------//
const Vector __fastcall operator*(const double &scalar, const Vector &A) {
	Vector result = A;
	result.updated = false;
	result.upd *= scalar;
	return result;
}

// -------------------------------- << ----------------------------------------//
ostream& __fastcall operator << (ostream& out_data, Vector& C) {
	long i;
	// ≈сли прос€т вывести пустой вектор - возвращаем NULL
	if ((&C == NULL) || (C.v == NULL)) {
		out_data << "NULL;";
		return out_data;
	}
	// выводим данные дл€ нормально заданных векторов
	if (!C.updated)
		C.update();
	out_data << "[";
	for (i = 0; i < C.v->size - 1; i++)
		out_data << C.v->v[i] << ",";
	out_data << C.v->v[C.v->size - 1];
	out_data << "]";
	out_data << ';' << endl;
	return out_data;
}

// -------------------------------- >> ----------------------------------------//
istream& __fastcall operator >> (istream& in_data, Vector& C) {
	char c = 0;
	double tempVal;
	long i;

	// »нициализируем несчитываемые значени€
	C.upd = 1.0;
	C.updated = true;

	// «аполн€ем матрицу данными из потока

	if (!C.updated)
		C.update();
	for (i = 0; i < C.v->size;) {
		in_data.get(c);
		switch (c) {
		case '[':
			break;
		case ']':
			break;
		case ',':
			break;
		case '\n':
			break;
		default: {
				in_data.putback(c);
				in_data >> tempVal;
				C.v->v[i] = tempVal;
				++i;
			}
		}
	}
	// C=vec;
	in_data.get(c);
	while (c != ';' && c != EOF)
		in_data.get(c);
	in_data.get(c);
	return in_data;
}

// ------------------------ GetSubVector --------------------------------------//
Vector __fastcall Vector::GetSubVector(long start, long end) {
	long i;
	Vector result(end - start + 1);

	update();
	for (i = start - 1; i < end; i++)
		result.v->v[i - start + 1] = v->v[i];
	result.upd = upd;
	result.updated = false;
	return result;
}

// ------------------------------- create ------------------------------------//
/* !inline */ void __fastcall Vector::create(long size) const {
	try {
		upd = 1;
		updated = true;
		v = new sVec(size);
	}
	catch (...) {
		cout << "Could not allocmte. Bye ...";
		exit(-1);
	}
}

// ------------------------------- detach -------------------------------------//
Vector& __fastcall Vector::detachT() {
	long i;

	if (v->linkCount > 1) {
		Vector *vv = new Vector(v->size);
		for (i = 0; i < v->size; i++)
			vv->v->v[i] = v->v[i];
		delete v;
		return *vv;
	}
	return *this;
}

// ------------------------------- detach -------------------------------------//
void __fastcall Vector::detach() const {
	long i;

	if (v->linkCount > 1) {
		sVec *vv = v;
		delete v;
		v = new sVec(vv->size);
		for (i = 0; i < v->size; i++)
			v->v[i] = vv->v[i];
	}
}

// ------------------------------- Update -------------------------------------//
void __fastcall Vector::update() const {
	long i;
	detach();
	if (upd != 1.0)
		for (i = 0; i < v->size; i++)
			v->v[i] *= upd;
	upd = 1;
	updated = true;
}

Vector* __fastcall Vector::copy(Vector* src, Vector* dst) {
	if ((src != NULL) && (dst != NULL) && (src->v == dst->v))
		return dst;
	if (dst != NULL) {
		delete dst;
		dst = NULL;
	}
	if (src != NULL)
		dst = new Vector(*src);
	return dst;
}

// ------------------------------------Get real Vector-------------------------//
// ¬ычисл€ет норму вектора относительно нул€, сдвига€ гиперкубическую сеть так,
// ноль был в центре, т.е. на половину разрешени€ сетки
// в принципе т.к. это расчЄт нормали относ€ценс€ только к гиперкубическим сеткам,
// то имеет смысл перенести это в соотв. модуль
/* TODO -orum -crepair :  ѕеренести функцию в TNet */
void __fastcall Vector::norm(const int& halfRes) {
	int i;
	double acc = 0;
	for (i = 0; i < v->size; ++i) {
		v->v[i] -= halfRes;
		acc += (v->v[i]*v->v[i]);
	}
	acc = sqrt(acc);
	if (acc != 0)
		for (i = 0; i < v->size; ++i)
			v->v[i] /= acc;
}

// ------------------------------------Get real Vector norm-------------------//
// ¬ычисл€ет норму  вектора
LDouble __fastcall Vector::norm() {
	int i;
	double acc = 0;
	for (i = 0; i < v->size; ++i)
		acc += (v->v[i]*v->v[i]);
	return sqrt(acc);
}
