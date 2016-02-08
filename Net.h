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
	unsigned long Res;
	// ����� ����� �� ����� � ���������� ������� � ������ ����������� ������

	short perfomance;
	Matrix *v; // ������ ���������
	unsigned long Count; // ����� ����� �����
	unsigned long NumOfPoints; // ����� ����� �� �����!!!
	unsigned long Dim; // ����������� ����
	unsigned long minBorder, maxBorder;
	unsigned int SurfaceDim; // ����������� ����������� ���������
	unsigned int NumOfSur; // ����� ������������ ���������

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
	unsigned long cacheCurrent; // ,
	// parsedCacheCurrent;
	unsigned long virtDim;
	double halfRes, dRes;
	unsigned long coordNumber; // , i;
	unsigned int _mod, __mod;

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
	double __fastcall getIJ(unsigned long ID, unsigned long coordNumber);

	friend ostream& __fastcall operator << (ostream&,  TNet&);
	friend istream& __fastcall operator >> (istream&, TNet&);

	/* !inline */ TNet& __fastcall operator *= (const Matrix&);
	/* !inline */// TNet& __fastcall operator *= (const double&);
	/* !inline */ TNet& __fastcall operator += (const TNet&);
	/* !inline */ TNet& __fastcall operator += (const Vector&);

	friend const TNet __fastcall operator *(const Matrix&, const TNet&);
	//friend const TNet __fastcall operator *(const double &, const TNet&);
	friend const TNet __fastcall operator +(const TNet&, const TNet&);
	friend const TNet __fastcall operator +(const Vector&, const TNet&);

	TNet& __fastcall operator = (const TNet&);

	virtual void __fastcall Clear();
	double __fastcall operator()(unsigned long i, unsigned long j);
	// ���������� �������� �� ����� i,j
	void /* !inline */ __fastcall copyNetFrom(const TNet&);
	void /* !inline */ __fastcall buildPowerVectors(unsigned long);
	virtual /* !inline */ Vector* __fastcall getVecAt(unsigned long i);

private:
	void /* !inline */ __fastcall MVMul(Vector* vsrc, Vector* vdst);
	// ������������ ������������ ������� �� ������� ������
	void /* !inline */ __fastcall MVMul(double* vsrc, double* vdst);
	// ������������ ������������ ������� �� �������� ������

	void inline __fastcall initNetDefault();
	void /* !inline */ __fastcall initNetParams(const unsigned long& res,
		const unsigned long& _res);

public:
	void /* !inline */ __fastcall checkCacheVector();
};
// ---------------------------------------------------------------------------
#endif
