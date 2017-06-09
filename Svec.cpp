/*
* Copyright (c) 2001-2017, Alexey Rumyantsev.
* e-mail: rumiantcev@yandex.ru
* All rights reserved.
*
*/
// ---------------------------------------------------------------------------

#include <stddef.h>
#include <iostream>
// #include <excpt>
#pragma hdrstop

#include "svec.h"

using namespace std;

// ------------------------------- destructor ---------------------------------//
void sVec:: operator delete(void *p) {
	if(p == NULL)
		return;

	sVec *ptr = static_cast<sVec*>(p);

	if (ptr->v == NULL)
		delete p;
	else
		p = NULL;
}

// ------------------------------- destructor ---------------------------------//
__fastcall sVec::~sVec() {
	if (--linkCount == 0) {
		delete[]v;
		v = NULL;
	}
}

// ------------------------------- destroy ------------------------------------//
///*!inline*/ void __fastcall sVec::destroy(){delete[] v; size =0;linkTo=NULL;}
// ------------------------------- create ------------------------------------//
void __fastcall sVec::create() {
	try {
		v = new LDouble[size];
		linkCount = 1;
	}
	catch (...) {
		cout << "Could not allocate. Bye ...";
		exit(-1);
	}
}

// ------------------------------- constructor --------------------------------//
__fastcall sVec::sVec(unsigned long sz) : size(sz) {
	create();
}

// ------------------------------- copy constr --------------------------------//
sVec::sVec(const sVec &C): size(C.size){
	create();
 	memcpy(v, C.v, size*sizeof(LDouble));
}

// ------------------------------- copy constr --------------------------------//
__fastcall sVec::sVec(const LDouble *vv, unsigned long sz) : size(sz) {
	create();
	memcpy(v, vv, size*sizeof(LDouble));
}
// ---------------------------------------------------------------------------
