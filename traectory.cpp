// ---------------------------------------------------------------------------

#pragma hdrstop

#include "traectory.h"

// ---------------------------------------------------------------------------

#pragma package(smart_init)

Traectory::Traectory() {
	T = 0;
	x_i = NULL;
	u_i = NULL;
	v_i = NULL;
}

Traectory::~Traectory() {
	for_each(NetList.begin(), NetList.end(), DeleteObj());
   //	 for_each(uu_i.begin(),uu_i.end(),DeleteObj());
   //	 for_each(vv_i.begin(),vv_i.end(),DeleteObj());
   //	 for_each(xx_i.begin(),xx_i.end(),DeleteObj());
	for_each(vecBrige.begin(), vecBrige.end(), DeleteObj());
	for_each(w.begin(), w.end(), DeleteObj());
	if (x_i != NULL)
		delete x_i;
	if (u_i != NULL)
		delete u_i;
	if (v_i != NULL)
		delete v_i;
}



