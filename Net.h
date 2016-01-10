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

// ����� ���������� ����������� ���������.
// ����������� � ���� ������� ���������, �� ��������� 0 � ������ ���������, �����. ������� ����� � 1 � ������ ����������, � ������ � -1 � ������ ����������

class TNet {
public:
	long Res;
	// ����� ����� �� ����� � ���������� ������� � ������ ����������� ������

	short perfomance;
	Matrix *v; // ������ ���������
	long Count; // ����� ����� �����
	long NumOfPoints; // ����� ����� �� �����!!!
	int Dim; // ����������� ����
	long minBorder, maxBorder;
	long SurfaceDim; // ����������� ����������� ���������
	long NumOfSur; // ����� ������������ ���������

	bool updated, // ������� ������������
		umx, cached, isVirtual;
	// isVirtual - ������� �������� ���� �.�. ���� true, �� ����� �� �������� � ������ � �������������� ��� ��������� � ���
	double upd; // ������������ ����������� ���������
	Matrix *u_mx; // ������������ �������

private:
	bool uvec;
	Vector *u_vec; // ������������ ������

protected:
	int currPlateNorm, currPlateDir;
	Vector *powVec,
		// ������ �������� �����������- ��������� ��� �������� � ������������������ ����
		*powVec_1,
		// ������ �������� �����������- ��������� ��� �������� � ������������������ ����  (��� ������������?)
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
	// ���������� �������� �� ����� i,j
	void /* !inline */ __fastcall copyNetFrom(const TNet&);
	void /* !inline */ __fastcall buildPowerVectors(int);
	virtual /* !inline */ Vector* __fastcall getVecAt(long i);

private:
	void /* !inline */ __fastcall MVMul(Vector* vsrc, Vector* vdst);
	// ������������ ������������ ������� �� ������� ������
	void /* !inline */ __fastcall MVMul(double* vsrc, double* vdst);
	// ������������ ������������ ������� �� �������� ������

	void inline __fastcall initNetDefault();
	void /* !inline */ __fastcall initNetParams(const int& res,
		const int& _res);

public:
	void /* !inline */ __fastcall checkCacheVector();
};
// ---------------------------------------------------------------------------
#endif
