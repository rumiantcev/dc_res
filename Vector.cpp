// ---------------------------------------------------------------------------

//#include <assert.h>
//#include <math.h>
//#include <algorithm>
// #include <dstring.h>
#pragma hdrstop
#include "Vector.h"

using namespace std;

// ------------------------------- size () ------------------------------------//
//long __fastcall Vector::size() const {
//	return v->size;
//}

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
__fastcall Vector::Vector(unsigned long sz) {
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
__fastcall Vector::Vector(const LDouble* vv, unsigned long sz) : upd(1.0), updated(true) {
	v = new sVec(vv, sz);
}

// ------------------------------- string reader --------------------------------//
__fastcall Vector::Vector(const string& str, unsigned long size) {
	unsigned long k;
	string::const_iterator i, j;
	create(size);
	i = find(str.begin(), str.end(), '[');
	j = ++i;
	string _val;
	for (k = 0; k < size; k++) {
		while (*j != ',')
			++j;
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
Vector& __fastcall Vector:: operator = (const LDouble & c) {
   unsigned	long i;

	for (i = 0; i < v->size; i++)
		v->v[i] = c;
	return *this;
}

// ---------------------------- += --------------------------------------------//
Vector& __fastcall Vector:: operator += (const Vector& A) {
	LDouble coeff = A.upd / upd;
	if (A.v == v)
		upd += A.upd;
	else {
		unsigned long i;
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
	if (A.v == v) {
		upd += A.upd;
	}
	else {
		unsigned long i;
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
const LDouble __fastcall operator *(const Vector& A, const Vector &B) {
	unsigned long i;
	LDouble res = 0;
	// assert(A.v->size==B.v->size);
	for (i = 0; i < A.v->size; i++)
		res += A.v->v[i] * B.v->v[i];
	res *= B.upd * A.upd;
	return res;

}

const LDouble __fastcall scMul(const Vector& A, const Vector &B)
{
	unsigned long i;
	LDouble res = 0;
	for (i = 0; i < A.v->size; i++)
		res += A.v->v[i] * B.v->v[i];
	return res;
}

// ------------------------------ * -------------------------------------------//
Vector& __fastcall Vector:: operator *= (const LDouble& scalar) {
	// !if (v->linkCount>1) detach();
	/* DONE -orum -cCheck : Проверить корректность оптимизации. */
	updated = false;
	upd *= scalar;
	return *this;
}

// ------------------------------ * -------------------------------------------//
const Vector __fastcall operator*(const LDouble &scalar, const Vector &A) {
	Vector result = A;
	result.updated = false;
	result.upd *= scalar;
	return result;
}

// -------------------------------- << ----------------------------------------//
ostream& __fastcall operator << (ostream& out_data, Vector& C) {
	unsigned long i;
	// Если просят вывести пустой вектор - возвращаем NULL
	if ((&C == NULL) || (C.v == NULL)) {
		out_data << "NULL;";
		return out_data;
	}
	// выводим данные для нормально заданных векторов
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
	LDouble tempVal;
	unsigned long i;

	// Инициализируем несчитываемые значения
	C.upd = 1.0;
	C.updated = true;

	// Заполняем матрицу данными из потока

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
	while ((c != ';') /*&& c != ','*/ && (!in_data.eof()))
		in_data.get(c);
	//in_data.get(c);
	return in_data;
}

// ------------------------ GetSubVector --------------------------------------//
Vector __fastcall Vector::GetSubVector(unsigned long start, unsigned long end) {
	unsigned long i;
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
	
	if (v->linkCount > 1) {
		unsigned long i;
		
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

	if (v->linkCount > 1) {
		unsigned long i;
		
		sVec *vv = v;
		delete v;
		v = new sVec(vv->size);
		for (i = 0; i < v->size; i++)
			v->v[i] = vv->v[i];
	}
}

// ------------------------------- Update -------------------------------------//
void __fastcall Vector::update() const {
	detach();
	if (upd != 1.0){
		unsigned long i;
		
		for (i = 0; i < v->size; i++)
			v->v[i] *= upd;
	}
	upd = 1;
	updated = true;
}

Vector* __fastcall Vector::copy(Vector* src, Vector* dst) {
	if ((src != NULL) &&
		(dst != NULL) &&
		(src->v == dst->v))
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
// Вычисляет норму вектора относительно нуля, сдвигая гиперкубическую сеть так,
// ноль был в центре, т.е. на половину разрешения сетки
// в принципе т.к. это расчёт нормали относяценся только к гиперкубическим сеткам,
// то имеет смысл перенести это в соотв. модуль
/* TODO -orum -crepair :  Перенести функцию в TNet */
void __fastcall Vector::norm(const LDouble& halfRes) {
	unsigned long i;
	LDouble acc = 0;
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
// Вычисляет норму  вектора
LDouble __fastcall Vector::norm() {
	unsigned long i;
	LDouble acc = 0;
	for (i = 0; i < v->size; ++i)
		acc += (v->v[i]*v->v[i]);
	return sqrt(acc);
}

//-------------евклидово расстояние между векторами---------------------------//
const LDouble eu_dist(const Vector& a, const Vector& b){
 unsigned long i;
 LDouble dist=0;

 for (i = 0; i < a.size(); i++)
	dist += (a[i]-b[i]) * (a[i]-b[i]);
 return sqrt(dist);
}
//------------------------- заполнение константами ---------------------------//
void Vector::one() {
	// 	memset(v->v,1.0,v->size*sizeof(LDouble));
	for (unsigned long i = 0; i < v->size; i++)
		v->v[i] = 1.0;
}

 void Vector::zero() {
	memset(v->v,0,v->size*sizeof(LDouble));
  //for (unsigned long i = 0; i < v->size; i++)
  //		v->v[i] = 0.0;
}
//------------------------- двоичное представление val -----------------------//
int Vector::toBinary(unsigned long val) {
	unsigned long len(log(val)/log(2));
	unsigned long i, rem, quot;

	if (size() <= len) {
		return -1; //не хватит места в массиве
	}
	quot = val;
    zero();
	for (i = v->size-1; i>=0; i--) {
		quot >>= 1;
		rem = val - (quot<<1);
		v->v[i] = rem;
		val = quot;
	}
    return 0;
}
