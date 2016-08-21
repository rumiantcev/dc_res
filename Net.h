// ---------------------------------------------------------------------------
#ifndef NetH
#define NetH

#pragma once

//#include<ostream>
//#include<istream>
//#include<math.h>

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
	unsigned long Res;
	// Число точек на грани в физическом массиве с учётом оптимизации памяти

	short perfomance;
	Matrix *v; // массив координат
	unsigned long Count; // общее число точек
	unsigned long NumOfPoints; // число точек на грани!!!
	unsigned long Dim; // размерность сети
	unsigned long minBorder, maxBorder;
	unsigned int SurfaceDim; // размерность повернности гиперкуба
	unsigned int NumOfSur; // число поверхностей гиперкуба

	bool updated, // признак актуальности
		umx, cached, isVirtual;
	// isVirtual - признак расчёной сети т.е. если true, то сетка не хранится в памяти а рассчитывается при обращении к ней
	LDouble upd; // кешированный коэффициент умножения
	Matrix *u_mx; // кешированная матрица

private:
	bool uvec;
	Vector *u_vec; // кешированный вектор

protected:
	int currPlateNorm, currPlateDir;
	Vector *powVec,
		// вектор степеней размерности- необходим для расчётов в виртуализированном виде
		*powVec_1,
		// вектор степеней размерности- необходим для расчётов в виртуализированном виде  (Для поверхностей)
		*cache; // , *parsedCache;
	unsigned long cacheCurrent; // ,
	// parsedCacheCurrent;
	unsigned long virtDim;
	LDouble halfRes, //половинка разрешения сетки
		dRes; //=1/Res
	unsigned long coordNumber; // , i;
	unsigned long _mod, __mod;

	// protected:
public:
	__fastcall TNet(unsigned long, long, unsigned long, bool);
	__fastcall TNet(unsigned long, unsigned long, long, unsigned long);
	__fastcall TNet(const TNet&);
	//__fastcall TNet();
	virtual __fastcall ~TNet();
	// void operator delete(void *p);

	virtual /* */ void __fastcall update();
	virtual /* */ void __fastcall dynUpdate();
	virtual /* */ void __fastcall detach();
	virtual /* */ void __fastcall destroy();
	Vector* __fastcall parseCoordinate(unsigned long current);
	Vector* __fastcall parseCoordinateForShift(unsigned long current);
	unsigned long __fastcall shift(unsigned long current, unsigned long coordNumber, int step,
		bool& borderChanged) throw(exInvalidMoveDirection);
	virtual void wasteCache();

private:
	virtual /* */ void __fastcall create( unsigned long /* , long perf, long res */);
	void __fastcall create(unsigned long mm, unsigned long nn, short perf, unsigned long res);

	unsigned long __fastcall getCurrentPlate(unsigned long ID);

	// Vector __fastcall getCoordVectorByID(long ID);
public:
	LDouble __fastcall getIJ(unsigned long ID, unsigned long coordNumber);

	friend ostream& __fastcall operator << (ostream&,  TNet&);
	friend istream& __fastcall operator >> (istream&, TNet&);

	/* !inline */ TNet& __fastcall operator *= (const Matrix&);
	/* !inline */// TNet& __fastcall operator *= (const LDouble&);
	/* !inline */ TNet& __fastcall operator += (const TNet&);
	/* !inline */ TNet& __fastcall operator += (const Vector&);

	friend const TNet __fastcall operator *(const Matrix&, const TNet&);
	//friend const TNet __fastcall operator *(const LDouble &, const TNet&);
	friend const TNet __fastcall operator +(const TNet&, const TNet&);
	friend const TNet __fastcall operator +(const Vector&, const TNet&);

	TNet& __fastcall operator = (const TNet&);

	virtual void __fastcall Clear();
	LDouble __fastcall operator()(unsigned long i, unsigned long j);
	// возвращаем значение по коорд i,j
	void /* !inline */ __fastcall copyNetFrom(const TNet&);
	void /* !inline */ __fastcall buildPowerVectors(unsigned long);
	virtual /* !inline */ Vector* __fastcall getVecAt(unsigned long i);

private:
	void /* !inline */ __fastcall MVMul(Vector* vsrc, Vector* vdst);
	// перемножение кешированной матрицы на заданый вектор
	void /* !inline */ __fastcall MVMul(LDouble* vsrc, LDouble* vdst);
	// перемножение кешированной матрицы на заданный вектор

	void inline __fastcall initNetDefault();
	void /* !inline */ __fastcall initNetParams(const unsigned long& res);

public:
	void /* !inline */ __fastcall checkCacheVector();
};
// ---------------------------------------------------------------------------
#endif
