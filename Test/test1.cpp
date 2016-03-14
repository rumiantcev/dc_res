#include <DUnitX.TestFramework.hpp>
#include <stdio.h>

#include "Vector.h"
#include "netfunc.h"
#pragma option --xrtti

class __declspec(delphirtti) dc_res_Test : public TObject
{
public:

  virtual void __fastcall SetUp();
  virtual void __fastcall TearDown();

__published:
  void __fastcall Test1();
  void __fastcall Test2();
};


void __fastcall dc_res_Test::SetUp()
{
}

void __fastcall dc_res_Test::TearDown()
{
}

void __fastcall dc_res_Test::Test1()
{
  // TODO
   Vector r(2);
   TNetF *net;

  String s("Hello");
  Dunitx::Testframework::Assert::IsTrue(s == "Hello");
}

void __fastcall dc_res_Test::Test2()
{
  // TODO
  String s("Hello");
  Dunitx::Testframework::Assert::IsTrue(s == "Bonjour"); // This fails for illustrative purposes
}

static void registerTests()
{
  TDUnitX::RegisterTestFixture(__classid(dc_res_Test));
}
#pragma startup registerTests 33