//---------------------------------------------------------------------------

#pragma hdrstop

#include "taskloader.h"
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

// ------------------------------ �������� ����� ���������� ------------------//
 void TaskLoader::loadGlobalParameters()
{
	string fileName = localPath+"\\dconsole.ini";
	cout<<fileName<<endl;
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
// --------------------- �������� � ������ ����� ��� �������� ----------------//

void TaskLoader::load_and_calc_tasks(){



















































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





// ------------------------------ load Pursue-Running Task -------------------//
Task* TaskLoader::loadTask(string szFileName) {
	char c, buf[20], long_buf[2048];
	int Lev, len, vecs, m, n, dx, du, dv, dm, stp, i, perf, meth, prior;
	LDouble ts, prec, tt, maxT, eps;
	// TNet *tmpNet;
	// Vector *tmpVec;
	string *tmpStr, descr;
	TNetF* tmpPRNet;
	vector<Traectory>trs;

	fstream in_f(szFileName.c_str(), ios::in|ios::nocreate);

	if ( in_f.is_open() )
	{
		// ��� ������ (0-�������������, 1 - �������������-��������)
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		meth = atoi(buf);
		while (c != ';')
			in_f.get(c);

		// ������������ �������
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);

		// ��������� ������ (����� �������� � ������������� ������)
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		prior = atoi(buf);
		while (c != ';')
			in_f.get(c);

		// ������� ������ - ���������  �������� 0- ������ ���������, 1 - ����� ����������, 2 - ���������� ���������
		// ����� ����� ��� �������������� ���������� �����
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		Lev = atoi(buf);
		while (c != ';')
			in_f.get(c);

		// ��������� �����
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		ts = _atold(buf);
		while (c != ';')
			in_f.get(c);

		// ��������� ����� - ����������� � ��������� ����� ������ ���������� - Traectory �� ���-�� �����
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

		// ������� A
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

		// ������� B
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

		// ������� C
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

		// ������� Pi
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

		// �������� epsilon �� ����������
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		eps = _atold(buf);
		while (c != ';')
			in_f.get(c);

		// �������� precision ������� �������������
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		prec = _atold(buf);
		while (c != ';')
			in_f.get(c);

		// �������� tau ��  �������
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		tt = _atold(buf);
		while (c != ';')
			in_f.get(c);

		// ������������ ����� ���������� �������������
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		maxT = _atold(buf);
		while (c != ';')
			in_f.get(c);
		// ����������� �������
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		dx = atoi(buf);
		while (c != ';')
			in_f.get(c);
		// ����������� ���������� ��������������
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		du = atoi(buf);
		while (c != ';')
			in_f.get(c);
		// ����������� ���������� ����������
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
		// ����������� ������������� ���������
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		dm = atoi(buf);
		while (c != ';')
			in_f.get(c);
		// ����� ����� �� ����� ���������
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		stp = atoi(buf);
		while (c != ';')
			in_f.get(c);
		// ������ �������, ������ ��� ��������
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
			res = new PR_Task(dx, du, dv, dm, ts, prec, eps, tt, maxT, stp, perf,
									Lev, 1);
		}

		res->A = AA;
		res->B = BB;
		res->C = CC;
		res->PP = PPP;
		res->tr_s = trs;
		res->priority = prior;
		res->method = meth;
		res->description = *tmpStr;
		delete tmpStr;

		// �������� ������� ������� ��������������

		while (c != '[')
			in_f.get(c);
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);
		res->setFuncToNetF(*(res->cP), *tmpStr);
		delete tmpStr;

		// �������� ������� ������� ����������
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);
		res->setFuncToNetF(*(res->cQ), *tmpStr);
		delete tmpStr;

		// �������� ������� ������� ����. ���������
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);
		res->setFuncToNetF(*(res->cM), *tmpStr);
		delete tmpStr;
		/* */

		// �������� ������� �������  �������� ���������������
		if (meth==1) {
			in_f.get(c);
			while (c != EOF){
				while (c != '[')
					in_f.get(c);
				in_f.get(long_buf, 256, ']');
				tmpStr = new string(long_buf);
				tmpPRNet  =  new TNetF(*res->cM);
				res->setFuncToNetF(*(tmpPRNet), *tmpStr);
				((PR_Task *)res)->PursuerList.push_back(tmpPRNet);
				delete tmpStr;
				in_f.get(c);
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