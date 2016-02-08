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

		//������ ��� �������� �����
		string func_m, func_v;
		unsigned long dim;
		short perfomance;
		unsigned long steps;

		pursuerType(unsigned long d, short p, unsigned long s);
		virtual ~pursuerType();


	friend istream& __fastcall operator >> (istream&, pursuerType&);
	   //	friend ostream& __fastcall operator << (ostream&, pursuerType&);
   //	private:



  };
#endif
