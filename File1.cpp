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
//#include <chrono>
#include<time.h>

#include "Vector.h"
#include "Matrix.h"
#include "Net.h"
#include "Netfunc.h"
#include "task.h"
#include "dll.h"
#include "SICx.h"
// #include "TESTx.h"


#include <float.h>	// needed for the FLT_MIN define

#include "CDataFile.h"


//#include <boost/program_options/config_file.hpp>

using namespace Dll;
//using namespace::std;
typedef vector<Task*>VecOfTask;

void loadGlobalParameters()
{
	CDataFile DFile(".\\dconsole.ini");

	Report(E_INFO, "[doSomething] The file <control.ini> contains %d sections, & %d keys.",
				   DFile.SectionCount(), DFile.KeyCount());

	_extr_e_param = DFile.GetFloat("extr_e_param","Параметры метода отжига");
	_extr_t0_param  = DFile.GetFloat("extr_t0_param", "Параметры метода отжига");

	_extr_tmin_param  = DFile.GetFloat("extr_tmin_param", "Параметры метода отжига");

	_extr_e_val = exp(_extr_e_param);
	/**/
  //  DFile.Clear();
}

void load_and_calc_tasks(){
	VecOfTask taskList;
	int j = 0;
	ofstream out_f;
	Task* t = NULL;

	CDataFile DFile(".\\control.ini");
	Report(E_INFO, "[doSomething] The file <control.ini> contains %d sections, & %d keys.",
				   DFile.SectionCount(), DFile.KeyCount());
	SectionItor i = DFile.GetFirstSectionIter();
	//i++;  //пропускаем пустую секцию
	t_Section* section;
	string fileName, method, opt,v_control;
	string sectionName;
	i++;

	while (i!= NULL) {
		section = (t_Section*)&(*i);
		fileName = ".\\"+section->Keys[0].szValue;//fileName = DFile.GetValue(section->szName, "file");
		t =  Task::loadTask(fileName);
		cout << "Task "<< t->description<< "  loaded" << endl;
		taskList.push_back(t);
		fileName = ".\\"+section->Keys[4].szValue;//DFile.GetValue(section->szName, "result");
		out_f.open(fileName.c_str());
		//загрузка настроечных параметров задачи
		method =  section->Keys[1].szValue;//DFile.GetValue(section->szName, "method");
		opt =  section->Keys[2].szValue;//DFile.GetValue(section->szName, "search_optimisation");
		v_control =  section->Keys[3].szValue;
		//Настройка параметров поиска экстремума на сетке
		if (opt == "none")
			taskList[j]->perfomance = optNone;
		if (opt == "annealing")
			taskList[j]->perfomance = optAnnealing;
		if (opt == "gradient")  //только для выпуклых множеств - использовать осторожно
			taskList[j]->perfomance = optGradient;
		taskList[j]->cP->perfomance  = taskList[j]->perfomance;
		taskList[j]->cQ->perfomance  = taskList[j]->perfomance;
		taskList[j]->cM->perfomance  = taskList[j]->perfomance;

		if (method == "pontryagin")
			taskList[j]->TimeCalc_Pontryagin(0);
		if (method == "alt_int")
			taskList[j]->TimeCalc_AltInt(0);
		cout << "Time evaluated.." << endl;
		out_f << "Time evaluated.." << endl;
		cout << "time is: " << taskList[j]->tr_s[0].T << endl;
		out_f << "time is: " << taskList[j]->tr_s[0].T << endl;
		if (method == "pontryagin")
			taskList[j]->Control_Pontryagin(0);
		if (method == "alt_int")
		   //	taskList[j]->Control_AltInt(0); //Временно
		  taskList[j]->Control_R1(0);
		cout << "Control evaluated.." << endl;
		cout << "Traectory: " << endl;
		cout<< *taskList[j]->tr_s[0].x_i;
		out_f << "Control evaluated.." << endl;
		out_f << "Traectory: " <<  endl;
		out_f << *taskList[j]->tr_s[0].x_i;
		out_f.flush();
		out_f.close();
		j++;
		i = DFile.GetNextSectionIter(i);
   }

	for_each(taskList.begin(), taskList.end(), DeleteObj());
}

int _tmain(int argc, _TCHAR* argv[]) {
	//	setlocale (LC_ALL, "RUS");
	/*auto t1 = std::chrono::high_resolution_clock::now();
		int s=0;
		for (int i=0; i<10; ++i)
			s+=i;
		auto t2 = std::chrono::high_resolution_clock::now();
		long dt = ((std::chrono::nanoseconds)(t2-t1)).count();
		std::cout << dt << std::endl;
	/**/


	randomize();
	cout << "Begin.." << endl;
	loadGlobalParameters();
	cout << "Parameters are loaded" << endl;
	
	load_and_calc_tasks();
   //	char c;
   // cin >> c;

	return 0;
}
