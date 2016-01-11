//---------------------------------------------------------------------------

#ifndef pursuerH
#define pursuerH

#include "Vector.h"
#include "pursuerType.h"
#include "general.h"
//---------------------------------------------------------------------------
class Pursuer{
	public:
		Vector* center;
		LDouble funnelDepth;
		LDouble radarVisibility;
		pursuerType* pType;
		long currInd;
		VecOfNetF funnel;

		Pursuer();
		virtual ~Pursuer();

		LDouble checkVisibility(const Vector& from);
};
#endif
