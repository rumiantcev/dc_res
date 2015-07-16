//---------------------------------------------------------------------------

#pragma hdrstop

#include "exception.h"
#include <tchar.h>
#include <stdio.h>
#include <string.h>
#include <windows.h>
//---------------------------------------------------------------------------

 Exception::Exception(const AnsiString s = "Unknown"){what = s);      }
 //Exception::Exception(int i){what = strdup("Unknown");      }

 Exception::Exception(const Exception& e ){what = e.wha); }

Exception::~Exception()                   {         }

 AnsiString Exception::msg() const             {return what;           }


#pragma package(smart_init)
