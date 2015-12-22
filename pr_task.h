
//---------------------------------------------------------------------------

#ifndef pr_taskH
#define pr_taskH

#include "task.h"
class PR_Task: public Task {
	public:
		vector<TNetF*>PursuerList;

        PR_Task(long , long dimU, long dimV, long dimM, LDouble ts, double prec,
	double eps, LDouble delta, LDouble tmax, long st, short perf, long stat,
	int tr_count);
		~PR_Task();
	 //	PR_Task* loadTask(string szFileName);
	 virtual void PR_Task::Find_Ns(int trNum);
	 //TNetF& __fastcall PR_Task::Find_Ns(int trNum);
};
//---------------------------------------------------------------------------
#endif

