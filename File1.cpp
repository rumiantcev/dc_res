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

	randomize();   //инициализируем рандомизатор
	setlocale( LC_ALL,"Russian.1572" ); //устанавливаем локаль русского язяка так, чтобы десятичным разделителем была точка, а не запятая
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	string localPath = "C://Users/Julia/Documents/Embarcadero/Studio/Projects/Test1/";  //переменная определяет где находятся данные и конф. файлы ... а также результаты расчётов
	TaskLoader tl(localPath);
	cout<<localPath<<endl;

	cout << "Begin.." << endl;
	tl.loadGlobalParameters();
	cout << "Parameters are loaded" << endl;
	cout << endl;

	string fileName = std::string(localPath + "Task_pr_2kr2.inp");
	PR_Task *t = NULL;
	t =  (PR_Task *) tl.loadTask(fileName);
	vector<TNetF*> tmp;
	tmp = ((PR_Task *)t)->PursuerList;

	/*while (!tmp.empty()) {
		cout << "!" << endl << *(tmp.back()) << endl;
		tmp.pop_back();
	} /**/
	//cout << "!!" << endl << *(t->cM) << endl;

	/**/

	int trNum;
	trNum = 0;
  ((PR_Task *)t)->Find_Ns(trNum);

	return 0;
}
