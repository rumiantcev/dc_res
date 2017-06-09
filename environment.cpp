/*
* Copyright (c) 2001-2017, Alexey Rumyantsev.
* e-mail: rumiantcev@yandex.ru
* All rights reserved.
*
*/
//---------------------------------------------------------------------------

#pragma hdrstop

#include <iostream>
#include "environment.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

class LocEnv : public Environment {
 // public:
 //	void op1()
  //		{  }
};

Environment::Environment()	{
	_extr_e_param  = -0.8;
	_extr_t0_param  = 1.0;
	_extr_tmin_param  = 0.0001;
	_extr_e_val = exp(_extr_e_param);
	ldZeroDf = 0.0000000001f;
	localPath = "..//Data";
    _pinMark = -999999;
}

Environment::~Environment()
	{/* std:: cout << "Destroyed" << std:: endl;*/ }

Environment *Environment::instance_ = 0;

Environment &Environment::instance() {
    if( !instance_ )
		instance_ = new LocEnv;
    return *instance_;
}

// Simple mechanism for Singleton cleanup.  Often works,
// but is far from foolproof.
class Destroyer {
  public:
	~Destroyer() { delete Environment::instance_; }
};
namespace {
Destroyer d;
}
