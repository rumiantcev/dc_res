/*
* Copyright (c) 2001-2017, Alexey Rumyantsev.
* e-mail: rumiantcev@yandex.ru
* All rights reserved.
*
*/
// ---------------------------------------------------------------------------
#ifndef traectoryH
#define traectoryH

#include <algorithm>
#include<vector>
#include"Vector.h"
#include "Matrix.h"
#include"Netfunc.h"
#include"general.h"
using namespace std;


// ---------------------------------------------------------------------------
class Traectory {
public:
	Vector x0; // ,x;
	Vector psi;
	LDouble mu, T;
	long ind;
	vector<TNetF*>NetList;
	VecOfLong psiExtr; // Значения psi на которых достигается extr

  //	VecOfVec uu_i, vv_i, xx_i; //траектория и управления для НЕфиксированного времени
	Matrix *u_i, *v_i, *x_i; //траектория и управления для фиксированного времени
	VecOfVec vecBrige;
	VecOfVec w;

	Traectory();
  	Traectory(const Traectory &from);
	~Traectory();


};
// ---------------------------------------------------------------------------
#endif
