// ---------------------------------------------------------------------------
#ifndef NetfuncH
#define NetfuncH
//#include <string>
//#include <map>
//#include <list>
//#include <iostream>
//#include <sstream>
//#include <assert.h>
//#include <system.hpp>
//#include <time.h>
#include "general.h"

using namespace std;
//#include<math.h>
// #include <omp.h>

#include "dll.h"
using namespace Dll;
#include "SICx.h"
// #include "Vector.h"
// #include "matrix.h"
#include "net.h"
// #include "rapideval\RapidEvaluator.hpp"


class TNetF;

typedef LDouble  (*cCrit) (long, LDouble, const TNetF*);
typedef LDouble __fastcall (*scM)(long, const Vector&, TNetF*, alphType* );

class TNetF : public TNet {
public:

	OpType operType;
	bool /* toPoint, */ alphaMode, is_empty; // ,updated;
	// double upd;  //������������ ����������� ���������
	// Matrix *u_mx;   	//������������ �������
	// ,isPoint; //��������� ��� �� ��������� ����� �������� �� ������������ �����������
	LDouble zeroPrecision; // �������� ���������� �� 0;
	LDouble maxRad;
	alphType alpha; // ����� ������� ��������� �� ������������ �����������
	Vector *f;
	TNet *net;
	string fStr;
	string lib_name;

	// double updSum;

	TDll *dll; // ���������� � �������������
	TSIC_Data sic; // ��� �����������
	// TRapidEvaluator       *func;
	double *vars;
	double t;

	__fastcall TNetF(int, long, long); // Ok
	__fastcall TNetF(int, long, long, string); // Ok
	__fastcall TNetF(long, int, long, long); // Ok
	explicit __fastcall TNetF(TNet&); // Ok
	// __fastcall TNetF(TNet&,TSIC_Data&);
	__fastcall TNetF(TNet&, const string&);
	__fastcall TNetF(const TNetF&); // Ok
	// __fastcall TNetF(TNetF&,TSIC_Data&);
	__fastcall TNetF(TNetF&, const string&);
	__fastcall ~TNetF();
	void operator delete(void *p); // Ok

	virtual void __fastcall update();
	virtual void __fastcall dynUpdate();
	virtual void __fastcall create(long /* , long, long */);
	virtual void __fastcall detach();

	friend ostream& __fastcall operator << (ostream&, TNetF&); // Ok
	friend istream& __fastcall operator >> (istream&, TNetF&); // Ok

	inline  double __fastcall getIJ(long current, int coordNumber);

	TNetF& __fastcall operator = (const TNetF&); // Ok

	/* !inline */ TNetF& __fastcall operator *= (const Matrix&); // Ok
	/* !inline */ TNetF& __fastcall operator *= (const LDouble &);
	/* !inline */ TNetF& __fastcall operator += (const TNetF&);
	/* !inline */ TNetF& __fastcall operator -= (const TNetF&);
	/* !inline */
	// TNetF& __fastcall operator += (const double &); �� ����� ������

	friend const TNetF __fastcall operator *(const Matrix&, const TNetF&);
	friend const TNetF __fastcall operator *(const LDouble, const TNetF&);
	friend const TNetF __fastcall operator +(const TNetF&, const TNetF&); // Ok
	friend const TNetF __fastcall operator -(const TNetF&, const TNetF&);
	// friend const TNetF __fastcall operator +(const double, const TNetF&); //�� ����� ������
	// friend TNetF __fastcall operator *(double &,TNetF&); // Changes in update()

	LDouble __fastcall oporn(const Vector &x, LDouble t, int sign);
	LDouble __fastcall oporn(const Vector &x, LDouble t, const string &funcStr,
		int sign);
	// LDouble __fastcall oporn(const Vector &x,  LDouble t,
	// const TRapidEvaluator &func,int sign);
	void __fastcall oporn(LDouble t, int sign);
	TNetF __fastcall opornx(LDouble t, int sign);
	LDouble __fastcall oporn(const Vector &x, LDouble t, int sign,
		const Matrix &A);
	void __fastcall oporn(LDouble t, int sign, const Matrix &A);

	long __fastcall GetExtrDirection(const Vector& vec, scM scmul,
	cCrit crit, OpType extrOper, ZeroAware isZeroAware, long index, alphType* coeff,
	LDouble& extr, 	TNetF& net);
	// ����� ���������� � ����������� �������� �������� vec
	long __fastcall GetExtrGlobal(OpType extrOper, long index, LDouble& extr);
	// ����� ���������� �� ������� �������

	void __fastcall Clear(); // Ok
	// /*!inline*/ double& Vector operator[](long i);
	// /*!inline*/  const double& operator [](long i)const;

	void /* !inline */ __fastcall SetFunc(const string& fstr); // Ok
	void __fastcall Conv(bool *); // �������� ��������
    void __fastcall ConvTimS(bool *); // �������� ��������
	void __fastcall ConvSerial(bool *, const Vector & a);
	TNet __fastcall Points(bool compactPoints);
	// ����������� �� ������� ������� � ����� �� ����� �� ������� ���������
	void __fastcall saveAsVrml(string);
	void __fastcall smoothFunction(double epsilon);
	virtual /* !inline */ Vector* __fastcall getVecAt(long i);
	long /* !inline */ selectExtrX(const Vector& vec, scM scmul, cCrit crit,
		long current, long& result, LDouble &extr, OpType extrOper,
		ZeroAware isZeroAware, bool &isExtrExist, TNetF& net);

protected:
	// ����� ���������� � ����������� �������� �������� vec ���������
	long __fastcall findExtrSlowDirection(const Vector& vec, scM scmul,
	cCrit crit,  ZeroAware isZeroAware,OpType extrOper, long index, alphType* coeff, LDouble& extr, TNetF& net);
	// ����� ���������� � ����������� �������� �������� vec ������� ������
	long __fastcall findExtrAnnealingDirection(const Vector& vec, scM scmul,
	cCrit crit,  ZeroAware isZeroAware,OpType extrOper, long index, alphType* coeff, LDouble& extr, TNetF& net);


	long __fastcall findExtrSlowGlobal(OpType extrOper,LDouble &extr);
	// ����� ���������� �� ������� ������� ���������
	long __fastcall findExtrAnnealingGlobal(OpType extrOper, LDouble &extr);
	// ����� ���������� �� ������� ������� ������� ������

	void __fastcall AddVariables();
	// void __fastcall SetVariables(const DynamicArray<double> &vv);

	virtual void __fastcall makeAlpha(alphType& alpha,  bool* L,
		 TNetF& net);


	long __fastcall findExtrFastXDirection(const Vector& vec, scM scmul,
	cCrit crit, OpType extrOper, ZeroAware isZeroAware, long index, alphType* coeff,
	LDouble& extr, 	TNetF& net);
	// ����� ���������� � ����������� �������� �������� vec
	long __fastcall findExtrFastXGlobal(OpType extrOper,
		 long index, LDouble& extr);
	// ����� ���������� �� ��. ������� � ����� ��������� �����������

public:
	void /* !inline */ __fastcall initNetFDefault();
	void /* !inline */ __fastcall copyNetFFrom(const TNetF&, bool);
	Vector __fastcall getBorderPoint(long index, const Vector& psi);

};

// --convex Criterias
// public:
extern LDouble  convCriteria(long num, LDouble sk, const TNetF *v);
extern LDouble  convCriteria1(long num, LDouble sk, const TNetF *v);
extern LDouble  convCriteria2(long num, LDouble sk, const TNetF *v);
// --Scalar Multex
extern LDouble __fastcall  scm (long num, const Vector &vec, TNetF *v, alphType* coeff);
extern LDouble __fastcall scm1(long num, const Vector &vec, TNetF *v, alphType* coeff);
extern LDouble __fastcall scm2(long num, const Vector &vec, TNetF *v, alphType* coeff);

#endif
