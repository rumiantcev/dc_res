/*
* Copyright (c) 2001-2017, Alexey Rumyantsev.
* e-mail: rumiantcev@yandex.ru
* All rights reserved.
*
*/
// ---------------------------------------------------------------------------

#pragma hdrstop

#include "traectory.h"

// ---------------------------------------------------------------------------

#pragma package(smart_init)

Traectory::Traectory() {
	T = 0;
	x_i = NULL;
	u_i = NULL;
	v_i = NULL;
	mu =0.0;
	ind =0;
}

Traectory::Traectory(const Traectory &from):
	x0(from.x0),psi(from.psi),mu(from.mu),T(from.T),ind(from.ind){

	//переносим "альтернированный интеграл"
	NetList.reserve(from.NetList.size());
	long i;
	TNetF *net;
	for (i = 0; i < from.NetList.size(); i++) {
		net = new TNetF(*from.NetList[i]);
		NetList.push_back(net);
	}

	psiExtr.reserve(from.psiExtr.size()); // переносим значения psi на которых достигается extr
	for (i = 0; i < from.psiExtr.size(); i++) {
		psiExtr.push_back(from.psiExtr[i]);
	}

	x_i = NULL;
	u_i = NULL;
	v_i = NULL;
  // переносим траекторию и управления для фиксированного времени
	if (from.u_i != NULL)
		u_i = new Matrix (*from.u_i);
	if (from.v_i != NULL)
		v_i = new Matrix (*from.v_i);
	if (from.x_i != NULL)
		x_i= new Matrix (*from.x_i);

	vecBrige.reserve(from.vecBrige.size());
	Vector *vec;
	for (i = 0; i < from.vecBrige.size(); i++) {
		vec = new Vector(*from.vecBrige[i]);
		vecBrige.push_back(vec);
	}

	w.reserve(from.w.size());
	for (i = 0; i < from.w.size(); i++) {
		vec = new Vector(*from.w[i]);
		w.push_back(vec);
	}

}  /**/

Traectory::~Traectory() {
	for_each(NetList.begin(), NetList.end(), DeleteObj());
	for_each(vecBrige.begin(), vecBrige.end(), DeleteObj());
	for_each(w.begin(), w.end(), DeleteObj());

	if (x_i != NULL)
		delete x_i;
	if (u_i != NULL)
		delete u_i;
	if (v_i != NULL)
		delete v_i;
}



