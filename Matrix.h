// ---------------------------------------------------------------------------
#ifndef matrixH
#define matrixH

#pragma once
// #include <stream>

#include "svec.h"
#include "vector.h"
#include "smx.h"
#include "general.h"

/**
 * ------------------------------- Matrix -------------------------------------<url>element://model:project::Project2/cpp:e_field:src:Project2:sVec.v</url>
 */
class Matrix {
public:
	// int		m, n;
	sMx *v;
	double upd;
	bool updated;
	// ,cached; Reserved

	__fastcall Matrix(long, long);
	__fastcall Matrix(long);
	__fastcall Matrix(const Matrix&);
	__fastcall Matrix(double **, long, long);
	__fastcall Matrix(const string&, long, long);
	// void operator delete(void *p);

	void __fastcall update();
	/* !inline */ void __fastcall create(bool, long, long);
	void __fastcall detach();
	virtual __fastcall ~Matrix();

	long m() const ; // число строк
	long n() const ; // число столбцов
	double& operator()(long i, long j);
	Matrix& __fastcall operator = (const Matrix&);

	friend ostream& __fastcall operator << (ostream&, Matrix&);
	friend istream& __fastcall operator >> (istream&, Matrix&);

	/* !inline */ Matrix& __fastcall operator += (const Matrix&);
	/* !inline */ Matrix& __fastcall operator *= (const double&);
	/* !inline */ Matrix& __fastcall operator *= (const Matrix&);

	friend const Matrix __fastcall operator +(const Matrix&, const Matrix&);
	friend const Matrix __fastcall operator *(const double&, const Matrix&);
	friend const Vector __fastcall operator *(const Matrix&, const Vector&);
	friend const Matrix __fastcall operator *(const Vector&, const Matrix&);
	friend const Matrix __fastcall operator *(const Matrix&, const Matrix&);
	Vector __fastcall GetRow(long);
	Vector __fastcall GetCol(long);
	void __fastcall SetRow(const Vector&, long);
	void __fastcall SetCol(const Vector&, long);
	Matrix __fastcall GetSubMatrix(long, long, long, long);
	friend Matrix __fastcall Mirror(const Matrix&);
	friend Matrix __fastcall Transpose(const Matrix&);
	friend Matrix __fastcall Exponential(Matrix&, double, double);
	friend Vector __fastcall Solve(const Matrix& A, const Vector& b, LDouble epsilon);
	static Matrix* __fastcall copy(Matrix* src, Matrix* dst);
	double __fastcall Norm();
};
// ---------------------------------------------------------------------------
#endif
