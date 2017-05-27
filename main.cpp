#pragma hdrstop
#pragma argsused

#ifdef _WIN32
#include <tchar.h>
#else
  typedef char _TCHAR;
  #define _tmain main
#endif


#include "general.h"
#include "task.h"
#include "taskloader.h"
#include "CDataFile.h"
#include "environment.h"


using namespace Dll;
using namespace::std;

int _test();

int _tmain(int argc, _TCHAR* argv[]) {

//-------------
	randomize();   //�������������� ������������
	setlocale( LC_ALL,"Russian.1572" ); //������������� ������ �������� ������ ���, ����� ���������� ������������ ���� �����, � �� �������
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	string localPath = getenv("DC_PATH");  //���������� ���������� ��� ��������� ������ � ����. ����� ... � ����� ���������� ��������
	Environment::instance().localPath=localPath;
	TaskLoader tl(localPath);
	cout<<localPath<<endl;

	cout << "Begin.." << endl;
	cout << "������: 07.01.2017" << endl;
	tl.loadGlobalParameters();
	cout << "Parameters are loaded" << endl;

	tl.load_and_calc_tasks();
   //	_test();

	return 0;
}

/*
int _test(){//��� ����������������� �������� ����������� � dc_res_test

   TNetF c(2, optTimS, 20), n(c), r(c);
   bool *L = new bool[c.Count];

   c.SetFunc("2.0*sqrt(p0*p0+p1*p1)");
   n.SetFunc("1.8*sqrt(p0*p0+p1*p1)");
   c.oporn(0, 1);
   n.oporn(0, 1);

   r = c;// + (-1)*n; // �������������� ��������
   cout << r;
   r.ConvTimS(L);//�����������
   //cout << r;

   r = c - n;//�������������� ��������
   //cout << r;

  delete[] L;
}
*/
