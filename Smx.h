/*
* Copyright (c) 2001-2017, Alexey Rumyantsev.
* e-mail: rumiantcev@yandex.ru
* All rights reserved.
*
*/
// ---------------------------------------------------------------------------
#ifndef sMxH
#define sMxH
#pragma once

#include <iostream>
#include "general.h"

/**
 * ------------------------------- sMx -------------------------------------
 */
 //Класс нужен для того, чтобы можно было быстро, без копирования массива, передавать охватывающий класс по значению
 //реализован несколько извращённый подсчёт ссылок (создавалось, когда auto_ptr ещё не было в проекте)
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
