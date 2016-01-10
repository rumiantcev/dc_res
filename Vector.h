// ---------------------------------------------------------------------------
#ifndef vectorH
#define vectorH
#pragma once
// ----------------------------------- Vector ---------------------------------//
//#include <iostream>
//#include <fstream>
#include "svec.h"
#include "smx.h"
 #include "general.h"
using namespace std;

/**
 * <url>element://model:project::Project2/cpp:e_field:src:Project2:sVec.v</url>
 */
class Vector {
public:
	// int	   	size;
	mutable sVec *v;
	mutable double upd;
	mutable bool updated; // ,

	// cached; Reserved
	__fastcall Vector();
	__fastcall Vector(long);
	__fastcall Vector(const Vector&);
	__fastcall Vector(const double*, long);
	__fastcall Vector(const string&, long);
	virtual __fastcall ~Vector();
	// void operator delete(void* p);

	void __fastcall update() const ;
	inline long __fastcall size() const {return v->size;};
	/* !inline */ void __fastcall create(long size) const ;
	void __fastcall detach() const ;
	Vector& __fastcall detachT();

	friend ostream& __fastcall operator << (ostream&, Vector&);
	friend istream& __fastcall operator >> (istream&, Vector&);

	Vector& __fastcall operator = (const Vector & C);
	Vector& __fastcall operator = (const double & c);

	/* !inline */ Vector& __fastcall operator += (const Vector&);
	/* !inline */ Vector& __fastcall operator *= (const double&);

	friend const Vector __fastcall operator +(const Vector&, const Vector&);
	friend const Vector __fastcall operator *(const double&, const Vector&);
	friend const double __fastcall operator *(const Vector&, const Vector&);
	// ---------------------------------------------------------------------------
	friend const double __fastcall scMul(const Vector& A, const Vector &B);
	// added for profiling
	Vector& __fastcall vSum(const Vector&); // added for profiling
	// ---------------------------------------------------------------------------

	/* !inline */ Vector __fastcall GetSubVector(long, long);

	/* версия с проверкой ошибок */
	// inline double& operator [](long i) {if((v->linkCount>1)||!updated) update(); if (i < v->size) return v->v[i]; else throw "out of bounds";}
	inline double& operator[](long i) {
		if ((v->linkCount > 1) || !updated)
			update();
		return v->v[i];
	} /* версия без проверки ошибок */

	inline double& operator()(long i) {
		if ((v->linkCount > 1) || !updated)
			update();
		return v->v[i];
	}

	inline const double& operator[](long i) const {
		if ((v->linkCount > 1) || !updated)
			update();
		return v->v[i];
	}
	// /*!inline*/ const double operator [](long i) const{return v->v[i];}
	static Vector* __fastcall copy(Vector* src, Vector* dst);
	void __fastcall norm(const int& halfRes);
	LDouble __fastcall norm();
	friend const LDouble eu_dist(const Vector&, const Vector&);
};
// ---------------------------------------------------------------------------
#endif
