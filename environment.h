//---------------------------------------------------------------------------

#ifndef environmentH
#define environmentH

#include <string>
#include "general.h"
//---------------------------------------------------------------------------
class Environment { // Singleton
  public:
	static Environment &instance();
	LDouble _extr_e_param;
	LDouble _extr_t0_param;
	LDouble _extr_tmin_param;
	LDouble _extr_e_val;
	LDouble ldZeroDf;
    LDouble _pinMark;
	std::string localPath;
   // virtual void op1() = 0;
    // other operations...
  protected:
    Environment();
    virtual ~Environment();
  private:
	static Environment *instance_;

	friend class Destroyer;
};


#endif
