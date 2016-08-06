//#include <DUnitX.TestFramework.hpp>
//#include <stdio.h>
//
#include "Vector.h"
#include "netfunc.h"
//#pragma option --xrtti
//
//class __declspec(delphirtti) dc_res_Test : public TObject
//{
//public:
//
//  virtual void __fastcall SetUp();
//  virtual void __fastcall TearDown();
//
//__published:
//  void __fastcall TestNetTimSConv();
//  void __fastcall Test2();
//};
//
//
//void __fastcall dc_res_Test::SetUp()
//{
//}
//
//void __fastcall dc_res_Test::TearDown()
//{
//}
//
//void __fastcall dc_res_Test::TestNetTimSConv()
//{
//  // TODO
//   TNetF c(2, optTimS, 20), n(c), r(c);
//   bool *L = new bool[c.Count];
//
//   c.SetFunc("2.0*sqrt(p0*p0+p1*p1)");
//   n.SetFunc("1.8*sqrt(p0*p0+p1*p1)");
//   c.oporn(0, 1);
//   n.oporn(0, 1);
//
//   r = c;// + (-1)*n; // арифметическая разность
//   cout << r;
//   r.ConvTimS(L);//овыпукление
//   //cout << r;
//
//   r = c - n;//геометричческая разность
//   //cout << r;
//
//
//  String s("Hello");
//  Dunitx::Testframework::Assert::IsTrue(s == "Hello");
//  delete L;
//}
//
//void __fastcall dc_res_Test::Test2()
//{
//  // TODO
//  String s("Hello");
//  Dunitx::Testframework::Assert::IsTrue(s == "Bonjour"); // This fails for illustrative purposes
//}
//
//static void registerTests()
//{
//  TDUnitX::RegisterTestFixture(__classid(dc_res_Test));
//}
//#pragma startup registerTests 33

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(Multiply)
{
//BOOST_CHECK(multiply(4, 5) == 20);
	BOOST_CHECK(1 == 2);
}
