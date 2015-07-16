// ---------------------------------------------------------------------------

#include <stddef.h>
// #include <math>
#include <fstream>
#pragma hdrstop
#include "sMx.h"


using namespace std;

// ------------------------------- destructor ---------------------------------//
void sMx:: operator delete(void *p) {
	sMx *ptr = (sMx*)p;

	if (ptr->v == NULL)
		delete(void*) p;
	else
		p = NULL;
}; /* */

// ------------------------------ destructor ----------------------------------//
__fastcall sMx::~sMx() {
	if (--linkCount == 0) {
		delete[]vv;
		vv = NULL;
		delete[]v;
		v = NULL;
	}
};

// ------------------------------- create -------------------------------------//
/* !inline */
void __fastcall sMx::create()
	// ручное создание двумерного массива со страничками по 16384
{
	long i, j, k, mm;

	long quant, count, last /* , ind */ ;
	if (m > 16384)
		quant = 16384;
	else
		quant = m;

	try {
		v = new double*[m]; // массив указателей (на строки)
		vv = new double[m * n]; // массив значений
		// cout<<"Matrix size: "<<m*n<<endl;
		count = long(m / quant); // число страниц по 16384
		last = m - count * quant; // mod   //остаток от числа страниц
		for (i = 0; i < count; i++) {
			k = i * quant;
			mm = k * n;
			for (j = 0; j < quant; j++) {
				v[k + j] = &vv[mm + j * n];
				// заполняем массив указателями на строки в рамках страниц
			}
		};
		if ((count > 1) && (i == count)) {
			k = i * quant;
			mm = k * n;
			for (j = 0; j < last; j++) {

				v[k + j] = &vv[mm + j * n];
				////заполняем массив указателями на строки в хвосте
			};
		};
		// !       for (j = 0; j < m; j++)
		// !        v[j] = &vv[j*n];
		// for (j = 0; j < m; j++)
		// v[j] = new double[n];
	}
	catch (...) {
		cout << "Could not allocate. Bye ...";
		exit(-1);
	};
}

// ----------------------------- constructor ----------------------------------//
__fastcall sMx::sMx(long mm, long nn) : m(mm), n(nn), linkCount(1) {
	create();
};

// ------------------------------- copy constr --------------------------------//
__fastcall sMx::sMx(sMx &C) : m(C.m), n(C.n), linkCount(1) {
	// ng i, j;
	memcpy(vv, C.vv, m*n*sizeof(double));
	// убить если не заработает рабочий код ниже
	memcpy(v, C.v, m*sizeof(double*));
	// vv=C.vv;
	// for (i = 0; i < m; i++)
	// for (j = 0; j < n; j++)
	// v[i][j] = C.v[i][j];
};

// ------------------------------- copy constr --------------------------------//
__fastcall sMx::sMx(double **vv, long mm, long nn) : m(mm), n(nn), linkCount(1)
{
	long i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			v[i][j] = vv[i][j];
};

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
};
// ---------------------------------------------------------------------------
