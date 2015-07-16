// ---------------------------------------------------------------------------
#ifndef sVecH
#define sVecH
#pragma once

 #include "general.h"
/* Описание простого вектора с возможностью счёта собственно сам вектор */
//Класс нужен собственно для того, чтобы можно было быстро, без копирования массива, перезавать охватывающий класс по значению
/**
 * -----------------------------------SimpleVec--------------------------------<url>element://model:project::Project2/cpp:e_field:src:Project2:sVec.v</url>
 */
class sVec {
public:
	long size;
	/**
	 * <url>element://model:project::Project2/cpp:e_class:src:Project2:TNet</url>
	 * <url>element://model:project::Project2/cpp:e_class:src:Project2:Matrix</url>
	 * <url>element://model:project::Project2/cpp:e_class:src:Project2:Vector</url>
	 * <url>element://model:project::Project2/cpp:e_class:src:Project2:sMx</url>
	 * <url>element://model:project::Project2/cpp:e_class:src:Project2:sVec</url>
	 */
	double *v;
	long linkCount;

	__fastcall sVec(long);
	__fastcall sVec(const double *vv, long sz);
	virtual __fastcall ~sVec();
	void operator delete(void* p);

private:
	void __fastcall create();
	sVec(const sVec&);
};
// ---------------------------------------------------------------------------
#endif
