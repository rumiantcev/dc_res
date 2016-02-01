// ---------------------------------------------------------------------------
#ifndef sMxH
#define sMxH
#pragma once

#include <iostream>
// #include"sVec.h"

/**
 * ------------------------------- sMx -------------------------------------<url>element://model:project::Project2/cpp:e_field:src:Project2:sVec.v</url>
 */
 //Класс нужен собственно для того, чтобы можно было быстро, без копирования массива, перезавать охватывающий класс по значению
class sMx {
public:
	long m, n;
	double **v, *vv;
	long linkCount;

	__fastcall sMx(long, long);
	explicit __fastcall sMx(long);
	__fastcall sMx(sMx&);
   //	__fastcall sMx(double **vv, long mm, long nn);
	virtual __fastcall ~sMx();
	void operator delete(void *p);
	/* !inline */ void __fastcall create();

};
// ---------------------------------------------------------------------------
#endif
