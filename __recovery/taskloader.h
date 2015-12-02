//---------------------------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

#include "task.h"
#include "pr_task.h"

#include "CDataFile.h"
#include "general.h"

#ifndef taskloaderH
#define taskloaderH
//---------------------------------------------------------------------------
typedef vector<Task*>VecOfTask;

class TaskLoader{
	public:
		string localPath;
		VecOfTask taskList;

		TaskLoader(string);
        virtual ~TaskLoader();

		void loadGlobalParameters();
		void load_and_calc_tasks();
		Task* loadTask(string szFileName);
};
#endif
