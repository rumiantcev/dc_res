/*
* Copyright (c) 2015-2017, Alexey Rumyantsev, Julia Filimonova.
* e-mail: rumiantcev@yandex.ru, jul305@gmail.com
* All rights reserved.
*
*/
//---------------------------------------------------------------------------

#ifndef pr_taskH
#define pr_taskH

#include "task.h"
#include "traectory.h"
#include "pursuerType.h"
#include "pursuer.h"
//#include "matrix.h"

#include "general.h"

using namespace std;

typedef vector<pursuerType*>VecOfPTypes;
typedef vector<Pursuer*>VecOfPursuers;

class PR_Task: public Task {
	public:
		vector<TNetF*>PursuerList;
		vector<TNetF*> NList;

		VecOfPTypes pTypes;
		VecOfPursuers Pursuers;

		PR_Task(long ,  long dimU,  long dimV,  long dimM, LDouble ts, LDouble prec,
			LDouble eps, LDouble delta, LDouble tmax,  long st, short perf,
			long stat,	int tr_count);

		PR_Task(const PR_Task&);
		virtual ~PR_Task();

		virtual void calcPursuerSets(int);
		virtual void Find_Ns(int);
		virtual LDouble TimeCalc_PR(int trNum);
		virtual void Control_PR(int trNum);
		virtual int Control_PR_fullSets(int trNum);

		virtual void plot(int trNum);
	private:
		void calcNextAltInt(int trNum, LDouble t,  /*Matrix& PiEtA,*/ Matrix& PiEtAk, Matrix& PiEtAkB, /* Matrix& PiEtAkC,*/ TNetF& c);
		void calcNextAltInt(int trNum, LDouble t);
		void buildFunnels(pursuerType* pT);
	protected:
		//рефакторенные функции
		TNetF &point_oporn(const Vector &PiEtAx0, TNetF &x0Net) const;
		void pointGeomDiff(int trNum, unsigned long k, const Vector &PiEtAx, TNetF &c, TNetF &x_Net) const;
		void storeResults(int trNum,  VecOfVec &vx_i, VecOfVec &vu_i, VecOfVec &vv_i);
		void storeLocResults(const Vector &u_i, const Vector &v_i, const Vector &x_i, VecOfVec &vx_i, VecOfVec &vu_i, VecOfVec &vv_i,
					Vector *r_i) const;

	public:
		friend  void Control_PR_fullSets_smooth(int trNum, PR_Task& mt);
        friend  void recusionStep(int trNum, PR_Task &mt, const Vector &v_i, VecOfVec &vx_i, VecOfVec &vu_i, VecOfVec &vv_i,
				  VecOfLong &vec_j, VecOfLong &vec_k,  alphType &vec_t, unsigned long &j,
				  unsigned long &k, LDouble &t, Vector &u_i, Vector &x_i);

};
//---------------------------------------------------------------------------


#endif

