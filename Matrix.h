/*
* Copyright (c) 2001-2017, Alexey Rumyantsev.
* e-mail: rumiantcev@yandex.ru
* All rights reserved.
*
*/
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
	LDouble upd;
	bool updated;
	// ,cached; Reserved

	__fastcall Matrix(unsigned long, unsigned long);
	explicit __fastcall Matrix(unsigned long);
	__fastcall Matrix(const Matrix&);
   //	__fastcall Matrix(LDouble **, long, long);
	__fastcall Matrix(const string&, unsigned long, unsigned long);
	// void operator delete(void *p);

	void __fastcall update();
	/* !inline */ void __fastcall create(bool, unsigned long, unsigned long);
	void __fastcall detach();
	virtual __fastcall ~Matrix();

	unsigned long m() const ; // число строк
	unsigned long n() const ; // число столбцов
	LDouble& operator()(unsigned long i, unsigned long j);
	Matrix& __fastcall operator = (const Matrix&);

	friend ostream& __fastcall operator << (ostream&, Matrix&);
	friend istream& __fastcall operator >> (istream&, Matrix&);

	/* !inline */ Matrix& __fastcall operator += (const Matrix&);
	/* !inline */ Matrix& __fastcall operator *= (const LDouble&);
	/* !inline */ Matrix& __fastcall operator *= (const Matrix&);

	friend const Matrix __fastcall operator +(const Matrix&, const Matrix&);
	friend const Matrix __fastcall operator *(const LDouble&, const Matrix&);
	friend const Vector __fastcall operator *(const Matrix&, const Vector&);
	friend const Matrix __fastcall operator *(const Vector&, const Matrix&);
	friend const Matrix __fastcall operator *(const Matrix&, const Matrix&);
	Vector __fastcall GetRow(unsigned long);
	Vector __fastcall GetCol(unsigned long);
	void __fastcall SetRow(const Vector&, unsigned long);
	void __fastcall SetCol(const Vector&, unsigned long);
	Matrix __fastcall GetSubMatrix(unsigned long, unsigned long, unsigned long, unsigned long);
	friend Matrix __fastcall Mirror(const Matrix&);
	friend Matrix __fastcall Transpose(const Matrix&);
	friend Matrix __fastcall Exponential(Matrix&, LDouble, LDouble);
	friend Vector __fastcall Solve(const Matrix& A, const Vector& b, LDouble epsilon);
	static Matrix* __fastcall copy(Matrix* src, Matrix* dst);
	LDouble __fastcall Norm();
};
// ---------------------------------------------------------------------------
#endif
