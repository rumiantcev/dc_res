//---------------------------------------------------------------------------
//#include <stdio.h>
//#include <iostream>
//#include <fstream>
//#include <string>

#include "task.h"
#include "pr_task.h"
#include "pursuerType.h"

#include "CDataFile.h"
#include "general.h"

#ifndef taskloaderH
#define taskloaderH
//---------------------------------------------------------------------------
typedef vector<Task*>VecOfTask;

//extern LDouble _extr_e_param;
//extern LDouble _extr_t0_param;
//extern LDouble _extr_tmin_param;
//extern LDouble _extr_e_val;
//extern LDouble ldZeroDf;

class TaskLoader{
	public:
		string localPath;
		VecOfTask taskList;

		explicit TaskLoader(string);
        virtual ~TaskLoader();

		void loadGlobalParameters();
		void load_and_calc_tasks();
		Task* loadTask(string szFileName);
};
#endif
