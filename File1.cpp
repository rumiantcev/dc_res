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



//#include <boost/program_options/config_file.hpp>

using namespace Dll;
using namespace::std;


int _tmain(int argc, _TCHAR* argv[]) {

	randomize();   //�������������� ������������
	setlocale( LC_ALL,"Russian.1572" ); //������������� ������ �������� ����� ���, ����� ���������� ������������ ���� �����, � �� �������
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	string localPath = getenv("DC_PATH");  //���������� ���������� ��� ��������� ������ � ����. ����� ... � ����� ���������� ��������
	Environment::instance().localPath=localPath;
	TaskLoader tl(localPath);
	cout<<localPath<<endl;


	cout << "Begin.." << endl;
	cout << "������: 01.02.2016" << endl;
	tl.loadGlobalParameters();
	cout << "Parameters are loaded" << endl;
	tl.load_and_calc_tasks();

	return 0;
}

