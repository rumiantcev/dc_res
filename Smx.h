// ---------------------------------------------------------------------------
#ifndef sMxH
#define sMxH
#pragma once

#include <iostream>
// #include"sVec.h"

/**
 * ------------------------------- sMx -------------------------------------<url>element://model:project::Project2/cpp:e_field:src:Project2:sVec.v</url>
 */
 //����� ����� ���������� ��� ����, ����� ����� ���� ������, ��� ����������� �������, ���������� ������������ ����� �� ��������
class sMx {
public:
	unsigned long m, n;
	double **v, *vv;
	long linkCount;

	__fastcall sMx(unsigned long, unsigned long);
	explicit __fastcall sMx(unsigned long);
	__fastcall sMx(sMx&);
   //	__fastcall sMx(double **vv, long mm, long nn);
	virtual __fastcall ~sMx();
	void operator delete(void *p);
	/* !inline */ void __fastcall create();

};
// ---------------------------------------------------------------------------
#endif
