/*
* Copyright (c) 2001-2017, Alexey Rumyantsev.
* e-mail: rumiantcev@yandex.ru
* All rights reserved.
*
*/
// ---------------------------------------------------------------------------
#ifndef sVecH
#define sVecH
#pragma once

 #include "general.h"
/* Описание простого вектора с возможностью счёта собственно сам вектор */
//Класс нужен для того, чтобы можно было быстро, без копирования массива, передавать охватывающий класс по значению
 //реализован несколько извращённый подсчёт ссылок (создавалось, когда auto_ptr ещё не было в проекте)
/**
 * -----------------------------------SimpleVec--------------------------------
 */
class sVec {
public:
	unsigned long size;
	
	LDouble *v;
	long linkCount;

	explicit __fastcall sVec(unsigned long);
	__fastcall sVec(const LDouble *vv, unsigned long sz);
	virtual __fastcall ~sVec();
	void operator delete(void* p);

private:
	void __fastcall create();
	sVec(const sVec&);
};
// ---------------------------------------------------------------------------
#endif
