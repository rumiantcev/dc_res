//---------------------------------------------------------------------------

#ifndef pursuerTypeH
#define pursuerTypeH
//---------------------------------------------------------------------------

//#include "pr_task.h"
#include "Netfunc.h"
#include "traectory.h"


#include "general.h"

typedef vector<TNetF*>VecOfNetF;
using namespace std;

class pursuerType{
	public:
		VecOfVec centres;

		//даннык для создания сеток
		string func_m, func_v;
		int dim;
		short perfomance;
		long steps;

		pursuerType(int d, short p, long s);
		virtual ~pursuerType();


	friend istream& __fastcall operator >> (istream&, pursuerType&);
	   //	friend ostream& __fastcall operator << (ostream&, pursuerType&);
   //	private:



  };
#endif
