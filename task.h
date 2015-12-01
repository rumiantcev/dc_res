// ---------------------------------------------------------------------------
#ifndef taskH
#define taskH
// ---------------------------------------------------------------------------
#include<classes.hpp>
#include"vector.h"
#include"matrix.h"
#include"net.h"
#include"netfunc.h"
// #include "RapidEvaluator.hpp"
#include"traectory.h"
#include<vector>
#include<iterator>
#include<fstream>
#include<math.h>

using namespace std;


typedef enum TaskState {
	maxTimeReached = 0
} TaskState;
// typedef vector<Vector> VecOfVec;
// typedef vector<Vector>::iterator VIter;
// typedef vector<long> VecOfLong;
// typedef vector<double> VecOfDouble;

class Task {
public:
	long Level, IndI;
	short perfomance, v_control;
	TaskState state;
	Matrix A, B, C, PP;
	LDouble currMu, currT;
	vector<Traectory>tr_s;

	TNetF *cP, *cQ, *cM;
	// ������� ������� �� ������ ��� �������� ��������� ���������������, ���������� � ������������� ���������


	long dim_x, dim_u, dim_v, dim_m1, dim_m, steps;
	// long maxWayLength;      //����� ������������� ���� �� ������� � ����� �� tau.
	double epsilon,
		// ��������, ������������ ������� ���������� �������� ����������
		// � ������������� ���������
		precision;
	// �������� �������� ��� ��������� ����������� (����. ������������)
	double tau, // ��� �� �������.
		maxTime, t0, T, tmpT;
	int priority, method; // , traectoryCount;
	string description;

	Task();
	Task(long, long, long, long, LDouble, double, double, LDouble, LDouble,
		long, short, long, int);
	Task(const Task&);
	virtual ~Task();
	Vector rungeCutt(const Vector& xn, const Vector& un, const Vector& vn);
	Vector rungeCutt(const Vector& xn, const Vector& un, const Vector& vn, LDouble tau);
	// virtual LDouble timeCriteria(const TNet&,long,LDouble);
	//virtual LDouble TimeCalc();
	virtual LDouble TimeCalc_Pontryagin(int trNum);
	virtual LDouble TimeCalc_AltInt(int trNum);
	//virtual void Control();
	virtual void Control_AltInt(int trNum);
	virtual void Control_Pontryagin(int trNum);
	virtual void Control_R1(int trNum);
	virtual void Control_R2(int trNum);
	// virtual long ImageWake()

	void saveTask(char *);
	static Task* loadTask(/*istream&*/ string);
	// friend istream& operator >>(istream& in_f, Task* res);
	
	void __fastcall setFuncToNetF(TNetF& net, string func);
	virtual long ImageUp();

protected:
	long psi0Index; // ������, ����������� ��� ���������� ����������.
};
#endif
