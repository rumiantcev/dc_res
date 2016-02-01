// ---------------------------------------------------------------------------

#include <stddef.h>
// #include <math>
#include <fstream>
#pragma hdrstop
#include "sMx.h"


using namespace std;

// ------------------------------- destructor ---------------------------------//
void sMx:: operator delete(void *p) {
    if(p==0) return;

	sMx *ptr = static_cast<sMx*>(p);

	if (ptr->v == NULL)
		delete(void*) p;
	else
		p = NULL;
}

// ------------------------------ destructor ----------------------------------//
__fastcall sMx::~sMx() {
	if (--linkCount == 0) {
		delete[]vv;
		vv = NULL;
		delete[]v;
		v = NULL;
	}
}

// ------------------------------- create -------------------------------------//
/* !inline */
void __fastcall sMx::create()
	// ручное создание двумерного массива со страничками по 16384
{
	long quant /* , ind */ ;
	if (m > 16384)
		quant = 16384;
	else
		quant = m;

	try {
		long i, j, k, mm, count, last;
		v = new double*[m]; // массив указателей (на строки)
		vv = new double[m * n]; // массив значений
		count = long(m / quant); // число страниц по 16384
		last = m - count * quant; // mod   //остаток от числа страниц
		for (i = 0; i < count; i++) {
			k = i * quant;
			mm = k * n;
			for (j = 0; j < quant; j++) {
				v[k + j] = &vv[mm + j * n];
				// заполняем массив указателями на строки в рамках страниц
			}
		}
		if ((count > 1) && (i == count)) {
			k = i * quant;
			mm = k * n;
			for (j = 0; j < last; j++) {
				v[k + j] = &vv[mm + j * n];
				////заполняем массив указателями на строки в хвосте
			}
		}

	}
	catch (...) {
		cout << "Could not allocate. Bye ...";
		exit(-1);
	}
}

// ----------------------------- constructor ----------------------------------//
__fastcall sMx::sMx(long mm, long nn) : m(mm), n(nn), linkCount(1) {
	create();
}

// ------------------------------- copy constr --------------------------------//
__fastcall sMx::sMx(sMx &C) : m(C.m), n(C.n), linkCount(1) {
	memcpy(vv, C.vv, m*n*sizeof(double));
	memcpy(v, C.v, m*sizeof(double*));
}

// ------------------------------- copy constr --------------------------------//
//__fastcall sMx::sMx(double **vv, long mm, long nn) : m(mm), n(nn), linkCount(1)
//{
//	long i, j;
//	for (i = 0; i < m; i++)
//		for (j = 0; j < n; j++)
//			v[i][j] = vv[i][j];
//}

// -------------------------- E constructor -----------------------------------//
__fastcall sMx::sMx(long mm) : m(mm), n(mm), linkCount(1) {
	long i, j;
	create();
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			if (i == j)
				v[i][i] = 1;
			else
				v[i][j] = 0;
}
// ---------------------------------------------------------------------------
