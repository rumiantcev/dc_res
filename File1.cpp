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
#include <locale>
//#include <chrono>
#include<time.h>
#include <io.h>

//#include "Matrix.h"
//#include "Net.h"
//#include "Netfunc.h"
#include "task.h"
//#include "dll.h"
//#include "SICx.h"
#include "taskloader.h"

//#include <float.h>	// needed for the FLT_MIN define

#include "CDataFile.h"


//#include <boost/program_options/config_file.hpp>

using namespace Dll;
using namespace::std;
typedef vector<Task*>VecOfTask;

int _tmain(int argc, _TCHAR* argv[]) {

	randomize();
  //	  setlocale( LC_ALL,"Russian.1572" );
  //setlocale( LC_ALL,"RUS" );
  //		SetConsoleCP(1251);
  //	SetConsoleOutputCP(1251);
  //	string localPath(argv[0]),pathHelper="\\Data";
 //	cout<<localPath<<endl;
    string localPath = getenv("DC_PATH");  //переменная определяет где находятся данные и конф. файлы ... а также результаты расчётов
	TaskLoader tl(localPath);

	//localPath.append(pathHelper);
	cout<<localPath<<endl;

	cout << "Begin.." << endl;
	tl.loadGlobalParameters();
	cout << "Parameters are loaded" << endl;
	tl.load_and_calc_tasks();

	return 0;
}
