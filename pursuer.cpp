//---------------------------------------------------------------------------

#pragma hdrstop

#include "pursuer.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

 Pursuer::~Pursuer(){
	for_each(funnel.begin(), funnel.end(), DeleteObj());
	if (center != NULL) 
		delete center;
 }

 Pursuer::Pursuer():currInd(0){
 }
