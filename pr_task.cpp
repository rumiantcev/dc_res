//---------------------------------------------------------------------------

#pragma hdrstop

#include "pr_task.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
// ------------------------------ desttructor --------------------------------//

PR_Task::~PR_Task() {
	for_each(PursuerList.begin(), PursuerList.end(), DeleteObj());
}

// ------------------------------ constructor --------------------------------//
PR_Task::PR_Task(long dimX, long dimU, long dimV, long dimM, LDouble ts, double prec,
	double eps, LDouble delta, LDouble tmax, long st, short perf, long stat,
	int tr_count):Task(dimX,  dimU,  dimV,  dimM,  ts,  prec,
	 eps,  delta,  tmax,  st,  perf,  stat,  tr_count) {
	}

//TNetF& __fastcall TNetF::

//TNetF& __fastcall PR_Task::Find_Ns(int trNum) {
void PR_Task::Find_Ns(int trNum) {
	vector<TNetF*> tmplist;     //'input'
	tmplist = PursuerList;

	vector<TNetF*> NList;       //'output'

	bool *L;
	LDouble t = t0;
	Matrix EtA(A.m()), PiEtA(A.m(), A.n()), PiEtAB(B.m(), B.n()), PiEtAC(C.m(), C.n()), intEtA(A.m());
	//Vector PiEtAx0(A.m()), psi(A.v->m), Fs(A.v->m), vmin(A.v->m);

	TNetF c(PP.m(), perfomance, steps), tmpPNet(*cP), tmpQNet(*cQ), tmpNet(c);

	//long i, j;

	intEtA = (0.5 * tau) * (Exponential(A, 0, precision) + Exponential(A, -tau, precision));

	while (!tmplist.empty()) {

		TNetF NiNet(*(tmplist.back()));
		NiNet.oporn(t0, 1);

        EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA * intEtA;

		while (1) {
			t += tau;

			PiEtAB = PiEtA * B;
			tmpPNet = *cP;
			tmpPNet.oporn(t, -1);
			tmpPNet *= PiEtAB;

			PiEtAC = PiEtA * C;
			tmpQNet = *(tmplist.back());
			tmpQNet.oporn(t, -1);
			tmpQNet *= PiEtAC;

			tmpNet = tmpQNet - tmpPNet;
		  /* */
			cout << tmpNet.is_empty << endl;
			if (tmpNet.is_empty) {
				break;
			}

			NiNet += (tau) * tmpNet;

		}
		L = new bool[NiNet.Count];
		NiNet.Conv(L);
		delete[]L;
		NList.push_back(&NiNet);
		tmplist.pop_back();

		cout << t << endl;
	}
	//return *this;        //?
		/**/
		//
}

