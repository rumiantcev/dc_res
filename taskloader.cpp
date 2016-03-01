//---------------------------------------------------------------------------

#pragma hdrstop

#include "taskloader.h"
#include "environment.h"

//#include <boost/program_options/detail/config_file.hpp>
//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/ini_parser.hpp>
//---------------------------------------------------------------------------
#pragma package(smart_init)
// ------------------------------ constructor --------------------------------//
TaskLoader::TaskLoader(string path):localPath(path) {
	localPath.append("\\Data");
}

// ------------------------------ desttructor --------------------------------//
TaskLoader::~TaskLoader() {
		for_each(taskList.begin(), taskList.end(), DeleteObj());
}

// ------------------------------ загрузка общих параметров ------------------//
 void TaskLoader::loadGlobalParameters()
{
	string fileName = localPath+"\\dconsole.ini";
	cout<<fileName<<endl;
	CDataFile DFile(fileName);
	Report(E_INFO, "[doSomething] The file <dconsole.ini> contains %d sections, & %d keys.",
				   DFile.SectionCount(), DFile.KeyCount());

	Environment::instance()._extr_e_param = DFile.GetFloat("extr_e_param","Параметры метода отжига");
	Environment::instance()._extr_t0_param  = DFile.GetFloat("extr_t0_param", "Параметры метода отжига");

	Environment::instance()._extr_tmin_param  = DFile.GetFloat("extr_tmin_param", "Параметры метода отжига");

	Environment::instance()._extr_e_val = exp(Environment::instance()._extr_e_param);
	/**/
  //  DFile.Clear();
  /*//-------------
	using boost::property_tree::ptree;

	ptree root;

	ptree wave_packet;
	wave_packet.put( "width", "1" );
	wave_packet.put( "position", "2.0" );

	ptree calculation_parameters;
	calculation_parameters.put( "levels", "15" );

	root.push_front(ptree::value_type( "calculation parameters", calculation_parameters ));
	root.push_front(ptree::value_type( "wave packet", wave_packet ));

	write_ini( std::cout, root );
	*/
}
// --------------------- загрузка и запуск задач для расчётов ----------------//

void TaskLoader::load_and_calc_tasks(){

	int j = 0;
	ofstream out_f;
	Task* t = NULL;
	clock_t before;
	double elapsed;
	PR_Task *prTask;

	CDataFile DFile(std::string(localPath+"\\control.ini").c_str());
	Report(E_INFO, "[doSomething] The file <control.ini> contains %d sections, & %d keys.",
				   DFile.SectionCount(), DFile.KeyCount());
	SectionItor i = DFile.GetFirstSectionIter();
	//i++;  //пропускаем пустую секцию
	t_Section* section;
	string fileName, method, opt/*,v_control*/, type;
   //	string sectionName;
	i++;
	while (i!=DFile.m_Sections.end()){
		section = (t_Section*)&(*i);
		fileName = localPath+"\\";
		fileName+= section->Keys[0].szValue;//fileName = DFile.GetValue(section->szName, "file");
		t =  loadTask(fileName);
		cout << "Task "<< t->description<< "  loaded" << endl;
		taskList.push_back(t);
		fileName = localPath+"\\";
		fileName += section->Keys[5].szValue;//DFile.GetValue(section->szName, "result");
		out_f.open(fileName.c_str());
		//загрузка настроечных параметров задачи
		type =  section->Keys[1].szValue;
		method =  section->Keys[2].szValue;//DFile.GetValue(section->szName, "method");
		opt =  section->Keys[3].szValue;//DFile.GetValue(section->szName, "search_optimisation");
	   //	v_control =  section->Keys[4].szValue;
		//Настройка параметров поиска экстремума на сетке
		if (opt == "none")
			taskList[j]->perfomance = optNone;
		if (opt == "annealing")
			taskList[j]->perfomance = optAnnealing;
		if (opt == "gradient")  //только для выпуклых множеств - использовать осторожно
			taskList[j]->perfomance = optGradient;
		if (opt == "TimS")  //с использованием нового метода овыпукления
			taskList[j]->perfomance = optTimS;
		taskList[j]->cP->perfomance  = taskList[j]->perfomance;
		taskList[j]->cQ->perfomance  = taskList[j]->perfomance;
		taskList[j]->cM->perfomance  = taskList[j]->perfomance;
		before = clock(); //замер времени начала выполнения расчётов
		if (type == "pursue"){
			if (method == "pontryagin")
				taskList[j]->TimeCalc_Pontryagin(0);
			if (method == "alt_int")
				taskList[j]->TimeCalc_AltInt(0);
			if (method == "gr1")
				taskList[j]->TimeCalc_AltInt(0);
			if (method == "gr2")
				taskList[j]->TimeCalc_AltInt(0);
		}
		if (type == "pursue_run"){
			prTask = static_cast<PR_Task*>(taskList[j]);
		   	prTask->calcPursuerSets(0);
		   // prTask->Find_Ns(0);
			prTask->TimeCalc_PR(0);
		}
		cout << "Time evaluated.." << endl;
		cout << "time is: " << taskList[j]->tr_s[0].T << endl;
		if (type == "pursue"){
			if (method == "pontryagin")
				taskList[j]->Control_Pontryagin(0);
			if (method == "alt_int")
				taskList[j]->Control_AltInt(0);
			if (method == "gr1")
				taskList[j]->Control_R1(0);
			if (method == "gr2")
				taskList[j]->Control_R2(0);
		}

		if (type == "pursue_run"){
			prTask = static_cast<PR_Task*>(taskList[j]);
			//prTask->Control_PR(0);
			prTask->Control_PR_fullSets(0);
		}
		elapsed = clock()-before;  //замер времени начала выполнения расчётов

		taskList[j]->plot(0);

		out_f << "Time evaluated.." << endl;
		out_f << "time is: " << taskList[j]->tr_s[0].T << endl;
		cout << "Control evaluated.." << endl;
		//cout << "Traectory: " << endl;
		//cout<< *taskList[j]->tr_s[0].x_i;
		cout<< "Evaluation time: "<<elapsed<<endl;   //вывод времени на экран
		out_f << "Control evaluated.." << endl;
		out_f << "Traectory: " <<  endl;
		out_f << *taskList[j]->tr_s[0].x_i;
		out_f << "Evaluation time: "<<elapsed<<endl; //вывод времени в файл
		out_f.flush();
		out_f.close();
		++j;
		++i;
		//i = DFile.GetNextSectionIter(i);
	}
}

// ------------------------------ load Pursue-Running Task -------------------//
Task* TaskLoader::loadTask(string szFileName) {
	char c, buf[20], long_buf[2048];
	int Lev, len, vecs, m, n, dx, du, dv, dm, stp, i, perf, meth, prior;
	LDouble ts, prec, tt, maxT, eps;
	string *tmpStr/*, descr*/;
	vector<Traectory>trs;

	fstream in_f(szFileName.c_str(), ios::in|ios::nocreate);

	if ( in_f.is_open() )
	{
		// Код метода (0-преследование, 1 - преследование-убегание)
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		meth = atoi(buf);
		while (c != ';')
			in_f.get(c);

		// Наименование примера
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);

		// Приоритет задачи (будет работать в многопоточном случае)
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		prior = atoi(buf);
		while (c != ';')
			in_f.get(c);

		// Уровень задачи - состояние  расчётов 0- только загружено, 1 - время рассчитано, 2 - управление построено
		// нужно также для многопоточного исполнения задач
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		Lev = atoi(buf);
		while (c != ';')
			in_f.get(c);

		// Начальное время
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		ts = _atold(buf);
		while (c != ';')
			in_f.get(c);

		// Начальная точка - загружается в начальную точку класса Траектория - Traectory по кол-ву точек
		while (c != '(')
			in_f.get(c);
		while (c != ',')
			in_f.get(c);
		in_f.get(buf, 20, ',');
		len = atoi(buf);
		in_f.get(c);
		in_f.get(buf, 20, ')');
		vecs = atoi(buf);
		for (i = 0; i < vecs; i++) {
			if (i != 0)
				while (c != ')')
					in_f.get(c);
			while (c != '\n')
				in_f.get(c);
			Vector xs(len);
			in_f.get(c);
			in_f >> xs;
			Traectory tr;
			tr.x0 = xs;
			trs.push_back(tr);
		}

		// Матрица A
		while (c != ',')
			in_f.get(c);
		in_f >> m;
		in_f.get(c);
		while (c != ',')
			in_f.get(c);
		in_f >> n;
		Matrix AA(m, n);
		while (c != '\n')
			in_f.get(c);
		in_f >> AA;

		// Матрица B
		while (c != ',')
			in_f.get(c);
		in_f >> m;
		in_f.get(c);
		while (c != ',')
			in_f.get(c);
		in_f >> n;
		Matrix BB(m, n);
		while (c != '\n')
			in_f.get(c);
		in_f >> BB;

		// Матрица C
		while (c != ',')
			in_f.get(c);
		in_f >> m;
		in_f.get(c);
		while (c != ',')
			in_f.get(c);
		in_f >> n;
		Matrix CC(m, n);
		while (c != '\n')
			in_f.get(c);
		in_f >> CC;

		// Матрица Pi
		while (c != ',')
			in_f.get(c);
		in_f >> m;
		in_f.get(c);
		while (c != ',')
			in_f.get(c);
		in_f >> n;
		Matrix PPP(m, n);
		while (c != '\n')
			in_f.get(c);
		in_f >> PPP;

		// Точность epsilon по траектории
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		eps = _atold(buf);
		while (c != ';')
			in_f.get(c);

		// Точность precision расчёта экспоненциала
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		prec = _atold(buf);
		while (c != ';')
			in_f.get(c);

		// Точность tau по  времени
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		tt = _atold(buf);
		while (c != ';')
			in_f.get(c);

		// Максимальное время завершения прелседования
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		maxT = _atold(buf);
		while (c != ';')
			in_f.get(c);
		// Размерность системы
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		dx = atoi(buf);
		while (c != ';')
			in_f.get(c);
		// Размерность управления преследователя
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		du = atoi(buf);
		while (c != ';')
			in_f.get(c);
		// Размерность управления убегающего
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		dv = atoi(buf);
		while (c != ';')
			in_f.get(c);

		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		while (c != ';')
			in_f.get(c);
		// размерность терминального множества
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		dm = atoi(buf);
		while (c != ';')
			in_f.get(c);
		// Шагов сетки на ребре гиперкуба
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		stp = atoi(buf);
		while (c != ';')
			in_f.get(c);
		// Методы расчёта, быстро или медленно
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		perf = atoi(buf);
		while (c != ';')
			in_f.get(c);
		/* */

		Task* res;
		switch (meth) {
			case 1:
				res = new PR_Task(dx, du, dv, dm, ts, prec, eps, tt, maxT, stp, perf,
									Lev, 1);
				break;
		default:
			res = new Task(dx, du, dv, dm, ts, prec, eps, tt, maxT, stp, perf,
									Lev, 1);
		};

		res->A = AA;
		res->B = BB;
		res->C = CC;
		res->PP = PPP;
		res->tr_s = trs;
		res->priority = prior;
		res->method = meth;
		res->description = *tmpStr;
		delete tmpStr;

		// Значение опорной функции преследователя

		while (c != '[')
			in_f.get(c);
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);
		res->setFuncToNetF(*(res->cP), *tmpStr);
		delete tmpStr;

		// Значение опорной функции убегающего
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);
		res->setFuncToNetF(*(res->cQ), *tmpStr);
		delete tmpStr;

		// Значение опорной функции терм. множества
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);
		res->setFuncToNetF(*(res->cM), *tmpStr);
		delete tmpStr;
		/* */

		// Значение опорных функций  множеств преследователей
		if (meth==1) {
			in_f.get(c);
			/*
			while (c != '[')
					in_f.get(c);
			while (!in_f.eof()){
				while ((c != '[')&&(!in_f.eof()))
					in_f.get(c);
				if (in_f.eof())  break;
				in_f.get(long_buf, 256, ']');
				tmpStr = new string(long_buf);
				tmpPRNet  =  new TNetF(*res->cM);
				res->setFuncToNetF(*(tmpPRNet), *tmpStr);
				((PR_Task *)res)->PursuerList.push_back(tmpPRNet);
				delete tmpStr;
				in_f.get(c);
			} */
		  pursuerType *tmpPType;
		  while (!in_f.eof()){
			while ((c != '[')&&(!in_f.eof()))
				in_f.get(c);
			if (in_f.eof())
				break;
			in_f.putback(c);
			tmpPType =  new pursuerType(res->dim_m, res->perfomance, res->steps);
			in_f >> *tmpPType;
			(static_cast<PR_Task *>(res))->pTypes.push_back(tmpPType);
		  }
		}


	return res;
	}
	else
	{
		cout<< " Unable to open file. Does it exist?"<<endl;;
		return false;
	}
}
