// ---------------------------------------------------------------------------
#ifndef NetH
#define NetH

#pragma once

#include<ostream>
#include<istream>
#include<math.h>

#include "vector.h"
#include "matrix.h"
#include "general.h"
// #include "Netfunc.h"

class exInvalidMoveDirection {
};

// сетка элеметнтов поверхности гиперкуба.
// формируется в виде массива координат, по умолчанию 0 в центре гиперкуба, соотв. нулевая грань с 1 в первой координате, а первая с -1 в первой координате

class TNet {
public:
	long Res;
	// Число точек на грани в физическом массиве с учётом оптимизации памяти

	short perfomance;
	Matrix *v; // массив координат
	long Count; // общее число точек
	long NumOfPoints; // число точек на грани!!!
	int Dim; // размерность сети
	long minBorder, maxBorder;
	long SurfaceDim; // размерность повернности гиперкуба
	long NumOfSur; // число поверхностей гиперкуба

	bool updated, // признак актуальности
		umx, cached, isVirtual;
	// isVirtual - признак расчёной сети т.е. если true, то сетка не хранится в памяти а рассчитывается при обращении к ней
	double upd; // кешированный коэффициент умножения
	Matrix *u_mx; // кешированная матрица

private:
	bool uvec;
	Vector *u_vec; // кешированный вектор

protected:
	int currPlateNorm, currPlateDir;
	Vector *powVec,
		// вектор степеней размерности- необходим для расчётов в виртуализированном виде
		*powVec_1,
		// вектор степеней размерности- необходим для расчётов в виртуализированном виде  (Для поверхностей?)
		*cache; // , *parsedCache;
	long cacheCurrent; // ,
	// parsedCacheCurrent;
	int virtDim;
	double halfRes, dRes;
	long coordNumber; // , i;
	int _mod, __mod;

	// protected:
public:
	__fastcall TNet(int, long, long, bool);
	__fastcall TNet(long, long, long, long);
	__fastcall TNet(const TNet&);
	__fastcall TNet();
	virtual __fastcall ~TNet();
	// void operator delete(void *p);

	virtual /* */ void __fastcall update();
	virtual /* */ void __fastcall dynUpdate();
	virtual /* */ void __fastcall detach();
	virtual /* */ void __fastcall destroy();
	Vector* __fastcall parseCoordinate(long current);
	Vector* __fastcall parseCoordinateForShift(long current);
	long __fastcall shift(long current, int coordNumber, int step,
		bool& borderChanged) throw(exInvalidMoveDirection);
	virtual void wasteCache();

private:
	virtual /* */ void __fastcall create(long dim /* , long perf, long res */);
	void __fastcall create(long mm, long nn, short perf, long res);

	int __fastcall getCurrentPlate(long ID);

	// Vector __fastcall getCoordVectorByID(long ID);
public:
	double __fastcall getIJ(long ID, int coordNumber);

	friend ostream& __fastcall operator << (ostream&,  TNet&);
	friend istream& __fastcall operator >> (istream&, TNet&);

	/* !inline */ TNet& __fastcall operator *= (const Matrix&);
	/* !inline */// TNet& __fastcall operator *= (const double&);
	/* !inline */ TNet& __fastcall operator += (const TNet&);
	/* !inline */ TNet& __fastcall operator += (const Vector&);

	friend const TNet __fastcall operator *(const Matrix&, const TNet&);
	friend const TNet __fastcall operator *(const double &, const TNet&);
	friend const TNet __fastcall operator +(const TNet&, const TNet&);
	friend const TNet __fastcall operator +(const Vector&, const TNet&);

	TNet& __fastcall operator = (const TNet&);

	virtual void __fastcall Clear();
	double __fastcall operator()(long i, int j);
	// возвращаем значение по коорд i,j
	void /* !inline */ __fastcall copyNetFrom(const TNet&);
	void /* !inline */ __fastcall buildPowerVectors(int);
	virtual /* !inline */ Vector* __fastcall getVecAt(long i);

private:
	void /* !inline */ __fastcall MVMul(Vector* vsrc, Vector* vdst);
	// перемножение кешированной матрицы на заданый вектор
	void /* !inline */ __fastcall MVMul(double* vsrc, double* vdst);
	// перемножение кешированной матрицы на заданный вектор

	void inline __fastcall initNetDefault();
	void /* !inline */ __fastcall initNetParams(const int& res,
		const int& _res);

public:
	void /* !inline */ __fastcall checkCacheVector();
};
// ---------------------------------------------------------------------------
#endif
