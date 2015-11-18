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

//#include "Vector.h"
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
using namespace::std;
typedef vector<Task*>VecOfTask;

void loadGlobalParameters(string localPath)
{
	string fileName = localPath+"\\dconsole.ini";
	cout<<fileName<<endl;
  // 	CDataFile DFile("C:\\Users\\alex\\Documents\\MSU\\dc_res\\Data\\dconsole.ini");
	CDataFile DFile(fileName);
	Report(E_INFO, "[doSomething] The file <control.ini> contains %d sections, & %d keys.",
				   DFile.SectionCount(), DFile.KeyCount());

	_extr_e_param = DFile.GetFloat("extr_e_param","��������� ������ ������");
	_extr_t0_param  = DFile.GetFloat("extr_t0_param", "��������� ������ ������");

	_extr_tmin_param  = DFile.GetFloat("extr_tmin_param", "��������� ������ ������");

	_extr_e_val = exp(_extr_e_param);
	/**/
  //  DFile.Clear();
}

void load_and_calc_tasks(string localPath){
	VecOfTask taskList;
	int j = 0;
	ofstream out_f;
	Task* t = NULL;
	clock_t before;
	double elapsed;

	CDataFile DFile(std::string(localPath+"\\control.ini").c_str());
	Report(E_INFO, "[doSomething] The file <control.ini> contains %d sections, & %d keys.",
				   DFile.SectionCount(), DFile.KeyCount());
	SectionItor i = DFile.GetFirstSectionIter();
	//i++;  //���������� ������ ������
	t_Section* section;
	string fileName, method, opt,v_control;
	string sectionName;
	i++;
	while (i!=/* NULL*/DFile.m_Sections.end()/**/) {
		section = (t_Section*)&(*i);
		fileName = localPath+"\\";
		fileName+= section->Keys[0].szValue;//fileName = DFile.GetValue(section->szName, "file");
		t =  Task::loadTask(fileName);
		cout << "Task "<< t->description<< "  loaded" << endl;
		taskList.push_back(t);
		fileName = localPath+"\\";
		fileName += section->Keys[4].szValue;//DFile.GetValue(section->szName, "result");
		out_f.open(fileName.c_str());
		//�������� ����������� ���������� ������
		method =  section->Keys[1].szValue;//DFile.GetValue(section->szName, "method");
		opt =  section->Keys[2].szValue;//DFile.GetValue(section->szName, "search_optimisation");
		v_control =  section->Keys[3].szValue;
		//��������� ���������� ������ ���������� �� �����
		if (opt == "none")
			taskList[j]->perfomance = optNone;
		if (opt == "annealing")
			taskList[j]->perfomance = optAnnealing;
		if (opt == "gradient")  //������ ��� �������� �������� - ������������ ���������
			taskList[j]->perfomance = optGradient;
		taskList[j]->cP->perfomance  = taskList[j]->perfomance;
		taskList[j]->cQ->perfomance  = taskList[j]->perfomance;
		taskList[j]->cM->perfomance  = taskList[j]->perfomance;
		before = clock(); //����� ������� ������ ���������� ��������
		if (method == "pontryagin")
			taskList[j]->TimeCalc_Pontryagin(0);
		if (method == "alt_int")
			taskList[j]->TimeCalc_AltInt(0);
		if (method == "gr1")
			taskList[j]->TimeCalc_AltInt(0);
		if (method == "gr2")
			taskList[j]->TimeCalc_AltInt(0);
		cout << "Time evaluated.." << endl;
		out_f << "Time evaluated.." << endl;
		cout << "time is: " << taskList[j]->tr_s[0].T << endl;
		out_f << "time is: " << taskList[j]->tr_s[0].T << endl;
		if (method == "pontryagin")
			taskList[j]->Control_Pontryagin(0);
		if (method == "alt_int")
			taskList[j]->Control_AltInt(0);
		if (method == "gr1")
			taskList[j]->Control_R1(0);
		if (method == "gr2")
			taskList[j]->Control_R2(0);
		elapsed = clock()-before;  //����� ������� ������ ���������� ��������
		cout << "Control evaluated.." << endl;
		cout << "Traectory: " << endl;
		cout<< *taskList[j]->tr_s[0].x_i;
		cout<< "Evaluation time: "<<elapsed<<endl;   //����� ������� �� �����
		out_f << "Control evaluated.." << endl;
		out_f << "Traectory: " <<  endl;
		out_f << *taskList[j]->tr_s[0].x_i;
		out_f << "Evaluation time: "<<elapsed<<endl; //����� ������� � ����
		out_f.flush();
		out_f.close();
		j++;
		i++;
		//i = DFile.GetNextSectionIter(i);
	}


	for_each(taskList.begin(), taskList.end(), DeleteObj());
}


int _tmain(int argc, _TCHAR* argv[]) {

  #pragma omp parallel
  {
    // Code inside this region runs in parallel.
    printf("Hello!\n");
  }
///
  //	  setlocale( LC_ALL,"Russian.1572" );
  //setlocale( LC_ALL,"RUS" );
  //		SetConsoleCP(1251);
  //	SetConsoleOutputCP(1251);

	string localPath(argv[0]),pathHelper="\\Data";
	cout<<localPath<<endl;
	localPath = getenv("DC_PATH");  //���������� ���������� ��� ��������� ������ � ����. ����� ... � ����� ���������� ��������

	localPath.append(pathHelper);
	cout<<localPath<<endl;
	randomize();
	cout << "Begin.." << endl;
	loadGlobalParameters(localPath);
	cout << "Parameters are loaded" << endl;

	load_and_calc_tasks(localPath);


	return 0;
}
