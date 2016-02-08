
//---------------------------------------------------------------------------

#ifndef pr_taskH
#define pr_taskH

#include "task.h"
#include "traectory.h"
#include "pursuerType.h"
#include "pursuer.h"

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

		PR_Task(long ,  long dimU,  long dimV,  long dimM, LDouble ts, double prec,
			double eps, LDouble delta, LDouble tmax,  long st, short perf,
			long stat,	int tr_count);
		virtual ~PR_Task();

		virtual void calcPursuerSets(int);
		virtual void Find_Ns(int);
		virtual LDouble TimeCalc_PR(int trNum);
		virtual void Control_PR(int trNum);
		virtual void Control_PR_fullSets(int trNum);
	private:
		void calcNextAltInt(int trNum, LDouble t,  /*Matrix& PiEtA,*/ Matrix& PiEtAk, Matrix& PiEtAkB, /* Matrix& PiEtAkC,*/ TNetF& c);
		void calcNextAltInt(int trNum, LDouble t);
		void buildFunnels(pursuerType* pT);


};
//---------------------------------------------------------------------------
#endif

