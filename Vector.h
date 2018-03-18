/*
* Copyright (c) 2001-2017, Alexey Rumyantsev.
* e-mail: rumiantcev@yandex.ru
* All rights reserved.
*
*/
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
	mutable LDouble upd;
	mutable bool updated; // ,

	// cached; Reserved
	__fastcall Vector();
	explicit __fastcall Vector(unsigned long);
	__fastcall Vector(const Vector&);
	__fastcall Vector(const LDouble*, unsigned long);
	__fastcall Vector(const string&, unsigned long);
	virtual __fastcall ~Vector();
	// void operator delete(void* p);

	void __fastcall update() const ;
	inline unsigned long __fastcall size() const {return v->size;};
	/* !inline */ void __fastcall create(long size) const ;
	void __fastcall detach() const ;
	Vector& __fastcall detachT();

	friend ostream& __fastcall operator << (ostream&, Vector&);
	friend istream& __fastcall operator >> (istream&, Vector&);

	Vector& __fastcall operator = (const Vector & C);
	Vector& __fastcall operator = (const LDouble & c);

	/* !inline */ Vector& __fastcall operator += (const Vector&);
	/* !inline */ Vector& __fastcall operator -= (const Vector&);
	/* !inline */ Vector& __fastcall operator *= (const LDouble&);

	friend const Vector __fastcall operator +(const Vector&, const Vector&);
    friend const Vector __fastcall operator -(const Vector&, const Vector&);
	friend const Vector __fastcall operator *(const LDouble&, const Vector&);
	friend const LDouble __fastcall operator *(const Vector&, const Vector&);
	// ---------------------------------------------------------------------------
	friend const LDouble __fastcall scMul(const Vector& A, const Vector &B);
	// added for profiling
	Vector& __fastcall vSum(const Vector&); // added for profiling
	// ---------------------------------------------------------------------------

	/* !inline */ Vector __fastcall GetSubVector(unsigned long, unsigned long);

	/* версия с проверкой ошибок */
	// inline LDouble& operator [](long i) {if((v->linkCount>1)||!updated) update(); if (i < v->size) return v->v[i]; else throw "out of bounds";}
	 LDouble& operator[](unsigned long i) {
		if ((v->linkCount > 1) || !updated)
			update();
		return v->v[i];
	} /* версия без проверки ошибок */

	inline LDouble& operator()(unsigned long i) {
		if ((v->linkCount > 1) || !updated)
			update();
		return v->v[i];
	}

	inline const LDouble& operator[](unsigned long i) const {
		if ((v->linkCount > 1) || !updated)
			update();
		return v->v[i];
	}

	void one();  //заполнение единицами
	void zero(); //заполнение нулями
	int toBinary(unsigned long val);
	// /*!inline*/ const LDouble operator [](long i) const{return v->v[i];}
	static Vector* __fastcall copy(Vector* src, Vector* dst);
	void __fastcall norm(const LDouble& halfRes);
	LDouble __fastcall norm();
	friend const LDouble eu_dist(const Vector&, const Vector&);
	friend bool operator== ( const Vector&, const Vector&);
};
// ---------------------------------------------------------------------------
#endif
