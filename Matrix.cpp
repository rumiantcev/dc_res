// ---------------------------------------------------------------------------

#include <assert.h>
#include <algorithm>
#pragma hdrstop
#include "matrix.h"
#include "smx.h"
#include "vector.h"

// -------------------------------------- m () --------------------------------//
long Matrix::m() const {
	return v->m;
}

// -------------------------------------- n () --------------------------------//
long Matrix::n() const {
	return v->n;
}

// ------------------------------- operator () --------------------------------//
double& Matrix:: operator()(long i, long j) {
	if ((v->linkCount > 1) && (!updated))
		update();
	return v->v[i][j];
}

// ------------------------------- destructor ---------------------------------//
/* void Matrix::operator delete(void *p)
 {
 Matrix *ptr = (Matrix*)p;

 if (ptr->v==NULL)
 delete (void*) p;
 else
 p=NULL;
 }/* */
// ------------------------------- destructor ---------------------------------//
__fastcall Matrix::~Matrix() {
	// if (--v->linkCount==0)
   //	{
		delete v;
		v = NULL;
   //	}
}

// ----------------------------- constructor ----------------------------------//
__fastcall Matrix::Matrix(long mm, long nn) {
	create(false, mm, nn);
}

// ------------------------------- copy constr --------------------------------//
__fastcall Matrix::Matrix(const Matrix &C) : v(C.v), upd(C.upd),
	updated(C.updated) {
	v->linkCount++;
}

// ----------------------------- constructor ----------------------------------//
//__fastcall Matrix::Matrix(double **vv, long mm, long nn) : upd(1), updated(true)
//{
//	v = new sMx(vv, mm, nn);
//}

// -------------------------- E constructor -----------------------------------//
__fastcall Matrix::Matrix(long mm) {
	create(true, mm, mm);
}

// ------------------------------- string reader --------------------------------//
__fastcall Matrix::Matrix(const string& str, long m, long n) {
	long k, l;
	string::const_iterator i, j;
	create(false, m, n);
	// i = &str[str.find('[', 0)];
	i = find(str.begin(), str.end(), '[');
	j = ++i;
	// j=i;
	string _val;
	for (l = 0; l < m; l++) {
		for (k = 0; k < n; k++) {
			while ((*j != ',') && (*j != ']'))
				++j;
			_val.assign(i, j);
			v->v[l][k] = atof(_val.c_str());
			i = ++j;
		}
		while ((*i != '[') && (i != str.end()))
			++i;
		j = ++i;
	}
}

// ----------------------------- += --------------------------------------------//
Matrix& __fastcall Matrix:: operator += (const Matrix& A) {
	
	// Matrix *pA = (Matrix*)(&A);
	double coeff = A.upd / upd;
	if (A.v == v) {
		upd += A.upd;
	}
	else {
		long i, j;
		if (v->linkCount > 1)
			detach();
		for (i = 0; i < A.v->m; i++)
			for (j = 0; j < A.v->n; j++)
				v->v[i][j] += coeff * A.v->v[i][j];
	}
	updated = false;
	return *this;
}

// ----------------------------- + --------------------------------------------//
const Matrix __fastcall operator +(const Matrix& A, const Matrix& B) {
	return (Matrix)A += B;
}

// ------------------------- * -----------------------------------------------//
Matrix& __fastcall Matrix:: operator *= (const Matrix& B) {
	long i, j, k;

	sMx *result = new(sMx)(v->m, B.v->n);
	double sum;

	if (v->linkCount > 1)
		detach();
	for (i = 0; i < v->m; i++)
		for (j = 0; j < B.v->n; j++) {
			sum = 0;
			for (k = 0; k < B.v->m; k++)
				sum += v->v[i][k] * B.v->v[k][j];
			result->v[i][j] = sum;
		}
	upd *= B.upd;
	updated = false;
	delete v;
	v = result;
	return *this;
}

// ------------------------- * -----------------------------------------------//
const Matrix __fastcall operator *(const Matrix& A, const Matrix &B) {
	return Matrix(A) *= B;
}

// ------------------------------- * ------------------------------------------//
Matrix& __fastcall Matrix:: operator *= (const double &scalar) {
	// if (v->linkCount>1) detach();
	/* DONE -orum -cCheck : Проверить корректность оптимизации для матриц. */
	updated = false;
	upd *= scalar;
	return *this;
}

// ------------------------------- * ------------------------------------------//
const Matrix __fastcall operator *(const double &scalar, const Matrix &A)
	// Matrix __fastcall operator *(const double &scalar, const Matrix &A)
{
	return Matrix(A) *= scalar;
}

// ----------------------------- * -------------------------------------------//
const Vector __fastcall operator *(const Matrix &A, const Vector& B) {
	long i, j;
	Vector result(A.v->m);

	for (i = 0; i < A.v->m; i++) {
		double sum = 0;
		for (j = 0; j < A.v->n; j++)
			sum += B.v->v[j] * A.v->v[i][j];
		result.v->v[i] = sum;
	}
	result.upd = A.upd * B.upd;
	result.updated = A.updated && B.updated;
	return result;
}

// ----------------------------- * -------------------------------------------//
const Vector __fastcall operator *(Vector& B, Matrix &A) {
	long i, j;
	Vector result = Vector(A.v->m);

	assert(B.v->size != A.v->m);

	for (j = 0; j < A.v->n; j++) {
		double sum = 0;
		for (i = 0; i < A.v->m; i++)
			sum += B.v->v[i] * A.v->v[i][j];
		result.v->v[j] = sum;
	}
	result.upd = A.upd * B.upd;
	result.updated = false;
	return result;
}

// ----------------------------------- << -------------------------------------//
ostream& __fastcall operator << (ostream& out_data, Matrix& C) {
	long i, j;
	// Если просят вывести пустую матрицу - возвращаем NULL
	if ((&C == NULL) || (C.v == NULL)) {
		out_data << "NULL;";
		return out_data;
	}
	// выводим данные для нормально заданных матриц

	if (!C.updated)
		C.update();
	for (i = 0; i < C.v->m; i++) {
		out_data /* <<i /* */ << "[";
		for (j = 0; j < C.v->n - 1; j++)
			out_data << C.v->v[i][j] << ",";
		out_data << C.v->v[i][C.v->n - 1];
		out_data << "]";
		if (i == C.v->m - 1)
			out_data << ';';
		out_data << '\n';
	}
	return out_data;
}

// ----------------------------------- >> -------------------------------------//
istream& __fastcall operator >> (istream& in_data, Matrix& C) {
	// Matrix mx(C.v->m,C.v->n);
	char c = 0;
	double tempVal;
	long i, j; // ,k=0,l=0;

	// Инициализируем несчитываемые значения
	C.upd = 1;
	C.updated = true;

	// Заполняем матрицу данными из потока
	for (i = 0; i < C.v->m;)
		for (j = 0; j < C.v->n;) {
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
						C.v->v[i][j] = tempVal;
						j++;
						if (j >= C.v->n)
							i++;
					}
			}
		}
	// C=mx;
	in_data.get(c);
	while (c != ';' && c != EOF)
		in_data.get(c);
	in_data.get(c);
	return in_data;
}
// -------------------------------- = -----------------------------------------//
Matrix& __fastcall Matrix:: operator = (const Matrix & C) {
	if (this == &C)
		return*this;
	delete v;
	v = C.v;
	upd = C.upd;
	updated = C.updated;
	v->linkCount++;

	return *this;
}

// ----------------------------- Transpose -----------------------------------//
Matrix __fastcall Transpose(const Matrix &A) {
	long i, j;

	Matrix result(A.v->m, A.v->n);
	for (i = 0; i < A.v->m; i++)
		for (j = 0; j < A.v->n; j++)
			result.v->v[j][i] = A.v->v[i][j];
	result.upd = A.upd;
	result.updated = false;
	return result;
} /* */

// ---------------------------- Норма Фробениуса матрицы----------------------//
double __fastcall Matrix::Norm() {
	long i, j;
	double res = 0;
	if (!updated)
		update();
	for (j = 0; j < v->n; j++) {
		for (i = 0; i < v->m; i++)
			res += v->v[i][j] * v->v[i][j];

	}
	res = sqrt(res);
	return res;
}

// ---------------------------- Exponential ----------------------------------//
Matrix __fastcall Exponential(Matrix& A, double t, double delta) {
	if (!A.updated)
		A.update();
	Vector x0(A.v->m);
	Matrix E(A.v->m), B(A.v->m, A.v->n), Bn(A.v->m, A.v->n); // ,
	long i, k, md;
	double tp;
	double norm = abs(t) * A.Norm();
	if (norm >= 1)
		md = ceil(log10(norm) / log10(2.0));
	else
		md = 0;
	tp = t / pow(2.0, double(md));
	B = E;
	Bn = E;
	k = 0;
	while ((Bn.Norm() > delta) && (k < 10)) {
		k = k + 1;
		Bn *= (tp / k)* (A);
		B = B + Bn;
	}
	for (i = 0; i < md; i++) {
		B = B * B;
	}
	return B;
}

// -------------------------- Get row -----------------------------------------//
Vector __fastcall Matrix::GetRow(long i) {
	Vector result(v->n);
	long j;

	for (j = 0; j < v->n; j++)
		result.v->v[j] = v->v[i][j] * upd;
	result.updated = false;
	return result;
}

// -------------------------- Get col -----------------------------------------//
Vector __fastcall Matrix::GetCol(long j) {
	Vector result(v->m);
	long i;

	for (i = 0; i < v->m; i++)
		result.v->v[i] = v->v[i][j] * upd;
	result.updated = false;
	return result;
}

// -------------------------- Set row -----------------------------------------//
void __fastcall Matrix::SetRow(const Vector& vec, long i) {

	if (upd != 0)
		for (long j = 0; j < v->n; j++)
			v->v[i][j] = vec.v->v[j] / upd;
	updated = false;
}

// -------------------------- Set col -----------------------------------------//
void __fastcall Matrix::SetCol(const Vector& vec, long j) {

	if (upd != 0)
		for (long i = 0; i < v->m; i++)
			v->v[i][j] = vec.v->v[i] / upd;
	updated = false;
}

// ---------------------------- Solve -----------------------------------------//
Vector __fastcall Solve(const Matrix& A, const Vector& b, LDouble epsilon) {
	Vector y(b), result(A.m());
	Matrix a(A);
	LDouble temp;
	long i,j, k, n;

	if (!a.updated) a.update();
	if (!y.updated) y.update();
	n = A.m();
	k = 0;

	while (k < n) {
		LDouble max;
        long index;
		// Поиск строки с максимальным A[i][k]
		max = abs(a.v->v[k][k]);
		index = k;
		for (i = k + 1; i < n; i++)
			if (abs(a.v->v[i][k]) > max) {
				max = abs(a.v->v[i][k]);
				index = i;
			}

		// Перестановка строк
		if (max < epsilon) {
			// нет ненулевых диагональных элементов
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return y;
		}
		for (j = 0; j < n; j++) {
			temp = a.v->v[k][j];
			a.v->v[k][j] = a.v->v[index][j];
			a.v->v[index][j] = temp;
		}
			temp = y[k];
			y[k] = y[index];
			y[index] = temp;
		// Нормализация уравнений
		for (i = k; i < n; i++) {
			temp = a.v->v[i][k];
			if (abs(temp) < epsilon) continue; // для нулевого коэффициента пропуск
			for (j = 0; j < n; j++)
				a.v->v[i][j] /=  temp;
			y[i] /= temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (j = 0; j < n; j++)
				a.v->v[i][j] -= a.v->v[k][j];
			y[i] -= y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--) {
		result[k] = y[k];
		for (i = 0; i < k; i++)
			y[i] -= a.v->v[i][k] * result[k];
	}
	return result;
}

// ---------------------------- Mirror ----------------------------------------//
Matrix __fastcall Mirror(Matrix& A) {
	if (!A.updated)
		A.update();
	long j;
	Matrix result(A.v->m), temp = A;
	Vector tmp(A.v->m);

	for (j = 0; j < A.v->n; j++) {
		tmp = result.GetCol(j);
		tmp = Solve(temp, tmp, 0.0001);
		result.SetCol(tmp, j);
	}
	return result;
}

// ------------------------------ GetSubMatrix --------------------------------//
Matrix __fastcall Matrix::GetSubMatrix(long a, long b, long c, long d) {
	long i, j;

	Matrix result(b - a + 1, d - c + 1);
	for (i = a - 1; i < b; i++)
		for (j = c - 1; j < d; j++)
			result.v->v[i - a + 1][j - c + 1] = v->v[i][j];
	result.upd = upd;
	result.update();
	return result;
}

// ------------------------------- Update -------------------------------------//
void __fastcall Matrix::update() {
	detach();
	if (upd != 1){
	   long  i, j;
		for (i = 0; i < v->m; i++)
			for (j = 0; j < v->n; j++)
				v->v[i][j] = upd * v->v[i][j];
	}
	upd = 1;
	updated = true;
}

// ------------------------------- create ------------------------------------//
/* !inline */ void __fastcall Matrix::create(bool isE, long mm, long nn) {
	try {
		upd = 1.0;
		updated = true;

		if (isE) {
			v = new sMx(mm);
		}
		else {
			v = new sMx(mm, nn);
		}
	}
	catch (...) {
		cout << "Could not allocate. Bye ...";
		exit(-1);
	}
}

// ------------------------------- detach -------------------------------------//
void __fastcall Matrix::detach() {

	if (v->linkCount > 1) {
		long  i, j;
		sMx *vv;
		v->linkCount--;
		vv = v;
		v = NULL;
		create(false, vv->m, vv->n);
		for (i = 0; i < vv->m; i++)
			for (j = 0; j < vv->n; j++)
				v->v[i][j] = vv->v[i][j];
		v->linkCount = 1;
	}
}
// ---------------------------------------------------------------------------

Matrix* __fastcall Matrix::copy(Matrix* src, Matrix* dst) {
	if ((src != NULL) && (dst != NULL) && (src->v == dst->v))
		return dst;
	if (dst != NULL) {
		delete dst;
		dst = NULL;
	}
	if (src != NULL)
		dst = new Matrix(*src);
	return dst;
}
