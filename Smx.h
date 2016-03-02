// ---------------------------------------------------------------------------
#ifndef sMxH
#define sMxH
#pragma once

#include <iostream>
#include "general.h"

/**
 * ------------------------------- sMx -------------------------------------<url>element://model:project::Project2/cpp:e_field:src:Project2:sVec.v</url>
 */
 //Класс нужен собственно для того, чтобы можно было быстро, без копирования массива, передавать охватывающий класс по значению
class sMx {
public:
	unsigned long m, n;
	LDouble **v, *vv;
	long linkCount;

	__fastcall sMx(unsigned long, unsigned long);
	explicit __fastcall sMx(unsigned long);
	__fastcall sMx(sMx&);
   //	__fastcall sMx(LDouble **vv, long mm, long nn);
	virtual __fastcall ~sMx();
	void operator delete(void *p);
	/* !inline */ void __fastcall create();

};
// ---------------------------------------------------------------------------
#endif
