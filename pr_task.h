
//---------------------------------------------------------------------------

#ifndef pr_taskH
#define pr_taskH

#include "task.h"
#include "traectory.h"
#include "pursuerType.h"
#include "pursuer.h"
#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>

#include "general.h"

using namespace std;

typedef vector<pursuerType*>VecOfPTypes;
typedef vector<Pursuer*>VecOfPursuers;

class PR_Task: public Task {
	public:
		vector<TNetF*>PursuerList;
		vector<TNetF*> NList;
		//vector<pair<LDouble, LDouble> > xy;

		VecOfPTypes pTypes;
		VecOfPursuers Pursuers;

		PR_Task(long ,  long dimU,  long dimV,  long dimM, LDouble ts, LDouble prec,
			LDouble eps, LDouble delta, LDouble tmax,  long st, short perf,
			long stat,	int tr_count);
		virtual ~PR_Task();

		virtual void calcPursuerSets(int);
		virtual void Find_Ns(pursuerType* pT);
		virtual LDouble TimeCalc_PR(int trNum);
		virtual void Control_PR(int trNum);
		virtual void Control_PR_fullSets(int trNum);

		virtual void plot(int trNum);
	private:
		void calcNextAltInt(int trNum, LDouble t,  /*Matrix& PiEtA,*/ Matrix& PiEtAk, Matrix& PiEtAkB, /* Matrix& PiEtAkC,*/ TNetF& c);
		void calcNextAltInt(int trNum, LDouble t);
		void buildFunnels(pursuerType* pT);

};
//---------------------------------------------------------------------------


#endif

