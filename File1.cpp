#pragma hdrstop
#pragma argsused

#ifdef _WIN32
#include <tchar.h>
#else
  typedef char _TCHAR;
  #define _tmain main
#endif


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <locale.h>
#include <windows.h>
#include<time.h>
#include <io.h>

#include "task.h"
#include "taskloader.h"
#include "CDataFile.h"


//#include <boost/program_options/config_file.hpp>

using namespace Dll;
using namespace::std;

int _tmain(int argc, _TCHAR* argv[]) {

	randomize();   //�������������� ������������
	setlocale( LC_ALL,"Russian.1572" ); //������������� ������ �������� ����� ���, ����� ���������� ������������ ���� �����, � �� �������
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);


	string localPath = getenv("DC_PATH");  //���������� ���������� ��� ��������� ������ � ����. ����� ... � ����� ���������� ��������
	TaskLoader tl(localPath);
	cout<<localPath<<endl;

	cout << "Begin.." << endl;
	tl.loadGlobalParameters();
	cout << "Parameters are loaded" << endl;
	tl.load_and_calc_tasks();

	return 0;
}
