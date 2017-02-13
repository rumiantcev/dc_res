//---------------------------------------------------------------------------

#pragma hdrstop

#include "pr_task.h"

#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>
//---------------------------------------------------------------------------
#pragma package(smart_init)
// ------------------------------ desttructor --------------------------------//

PR_Task::~PR_Task() {
	for_each(PursuerList.begin(), PursuerList.end(), DeleteObj());
	for_each(NList.begin(), NList.end(), DeleteObj());
	for_each(pTypes.begin(), pTypes.end(), DeleteObj());
   	for_each(Pursuers.begin(), Pursuers.end(), DeleteObj());
}

// ------------------------------ constructor --------------------------------//
PR_Task::PR_Task(long dimX,  long dimU,  long dimV,  long dimM, LDouble ts, LDouble prec,
	LDouble eps, LDouble delta, LDouble tmax,  long st, short perf, long stat,
	int tr_count):Task(dimX,  dimU,  dimV,  dimM,  ts,  prec,
	 eps,  delta,  tmax,  st,  perf,  stat,  tr_count) {
	}
// ------------------------------ copy constructor --------------------------------//
PR_Task::PR_Task(const PR_Task& ts):Task(ts) {
	}
	//------------------------------расчёт множеств альтернированных интегралов преследователей
void PR_Task::Find_Ns(int trNum) {
	vector<TNetF*> tmplist;     //'input'
	tmplist = PursuerList;

   //	vector<TNetF*> NList;       //'output'// перенесено отсюда в класс, чтобы результат не разрушался при выходе из функции
	//bool *L;
	LDouble t = t0;
	Matrix EtA(A.m()), PiEtA(A.m(), A.n()), PiEtAkB(B.m(), B.n()), PiEtAkC(C.m(), C.n()), intEtA(A.m());
	Matrix PiEtAk(A.m(), A.n());//было забыто
	//Vector PiEtAx0(A.m()), psi(A.v->m), Fs(A.v->m), vmin(A.v->m);

	TNetF /* *c,*/ tmpPNet(*cP), tmpQNet(*cQ), tmpNet(PP.m(), perfomance, steps);

	//long i, j;

	intEtA = (0.5 * tau) * (Exponential(A, 0, precision) + Exponential(A, -tau, precision));

	while (!tmplist.empty()) {

		TNetF NiNet(*(tmplist.back()));
		NiNet.oporn(t0, 1);

	   //	cout<<NiNet<<endl;

		t = t0;   // время  инициализируется для каждого мн-ва
	   //	EtA = Exponential(A, t, precision);  // экспоненциал соотв. образом пересчитывается  поэтому это перенесено в цикл *1
	   //	PiEtA = PP * EtA * intEtA;

		while (1) {
			t += tau;

			EtA = Exponential(A, t, precision);  // *1 вот сюда
			PiEtA = PP * EtA * intEtA;
			PiEtAk = PiEtA * intEtA;

			PiEtAkB = PiEtAk * B;
			tmpPNet = *cP;
			tmpPNet.oporn(t, -1);
			tmpPNet *= PiEtAkB;

			PiEtAkC = PiEtAk * C;
			//c = tmplist.back();
			tmpQNet = *cQ;//*(tmplist.back());
			tmpQNet.oporn(t, -1);
			tmpQNet *= PiEtAkC;

		   //	tmpNet = tmpQNet - tmpPNet;
			NiNet = NiNet + tmpQNet;
			NiNet -= tmpPNet; // геометрическая разность
			cout<<NiNet<<endl;
		  /* */
			cout << NiNet.is_empty << endl;
			if (tmpNet.is_empty) {
				break;
			}

		   //	NiNet += (tau) * tmpNet;
		  //	NiNet.update();
		}
		// L = new bool[NiNet.Count];    //в рамках расчёта геом. разности уже встроено овыпукление
	   //	NiNet.Conv(L);
	   //	delete[]L;
	  	NList.push_back(&NiNet);
	  //	if (c!=NULL) delete c;// c тем, чтобы корректно разрушить то, что иы удаляешь из tmpList
		tmplist.pop_back();

		cout << t << endl;
	}
	//return *this;        //?
		/**/
		//
}

// ------ расчёт множеств альтернированных интегралов преследователей - Костыль --------//
void PR_Task::calcPursuerSets(int trNum){
//	TNetF m2Net(*PursuerList[0]); // сеточная опорная функция терминального множества
// для первого преследователя
  //	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()), PiEtAkC(C.m(), C.n()),
  //		PiEtAkB(A.m(), A.n()), EtA(A.m()), intEtA(A.m());
	// EtA(A.m()) - единичная при t=0
  //	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);
   //	LDouble t = t0, min;
 //	TNetF c(PP.m(), perfomance, steps), tmpPNet(*cP), tmpQNet(*cQ), Net(c),
 //		x0Net(c), tmpNet(c);
    unsigned long /*i,*/j/*, Ind = 0*/;
	pursuerType *pt;

   /*
	//костыль
	//просто рассчитываем опорную функцию, пользуясь тем, что управление константа не зависящая от времени
	for(j=0;j<PursuerList.size();j++){
	   PursuerList[j]->oporn(0.0,1);
	 //  cout<<*PursuerList[j]<<endl;
	}
	*/

   for(j=0;j<pTypes.size();j++){
	   pt = pTypes[j];
	   buildFunnels(pt);
	 //  cout<<*PursuerList[j]<<endl;
	}

}
//--------------------расчёт шага альтернированного интеграла-----------------//
void PR_Task::calcNextAltInt(int trNum, LDouble t){
	TNetF c(PP.m(), perfomance, steps);
	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()),  EtA(A.m()), intEtA(A.m()),
	//PiEtAkC(C.m(), C.n()),
	PiEtAkB(A.m(), A.n());

	c=*tr_s[trNum].NetList[tr_s[trNum].NetList.size()-1];

	intEtA = (0.5 * tau) * (Exponential(A, 0, precision) + Exponential(A, -tau,
		precision));
	intEtA.update();
	EtA = Exponential(A, t, precision);
	PiEtA = PP * EtA;
	PiEtAk = PiEtA * intEtA;

	calcNextAltInt(trNum, t, /*PiEtA,*/ PiEtAk, PiEtAkB/*,PiEtAkC*/,c);

}

void PR_Task::calcNextAltInt(int trNum, LDouble t, /* Matrix& PiEtA,*/ Matrix& PiEtAk, Matrix& PiEtAkB, /* Matrix& PiEtAkC,*/ TNetF& c){

		TNetF  tmpPNet(*cP)/*, tmpQNet(*cQ),*/ ;

		PiEtAkB = PiEtAk * B;
		//PiEtAkC = PiEtAk * C;

		tmpPNet = *cP;
		//cP->update();
		//tmpQNet = *cQ;
		// cQ->update();
		tmpPNet.oporn(t, 1);
		//tmpQNet.oporn(t, 1);
		tmpPNet *= PiEtAkB;
		//tmpQNet *= PiEtAkC;
		tmpPNet.update();
		//tmpQNet.update();

		c = c + tmpPNet;
		//c -= tmpQNet; // геометрическая разность
		c.update();
		tr_s[trNum].NetList.push_back(new TNetF(c));
		// Сохраняем для последующего поиска Psi
}

  // ------------------------------- Time calc--------------------------------//
LDouble PR_Task::TimeCalc_PR(int trNum) {
	TNetF m2Net(*cM); // сеточная опорная функция терминального множества
	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()), //PiEtAkC(C.m(), C.n()),
		PiEtAkB(A.m(), A.n()), EtA(A.m()), intEtA(A.m());
	// EtA(A.m()) - единичная при t=0
	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);
	LDouble t = t0, min;
	TNetF c(PP.m(), perfomance, steps), tmpPNet(*cP), /*tmpQNet(*cQ),*/ Net(c),
		x0Net(c), tmpNet(c);
	unsigned  long i;
	long Ind = 0;

	// интегрирование методом трапеций по двум точкам
	intEtA = (0.5 * tau) * (Exponential(A, 0, precision) + Exponential(A, -tau, precision));
	intEtA.update();

	m2Net.oporn(t0, 1); // опорное значение сетки по M2
	c = m2Net;
	tr_s[trNum].NetList.push_back(new TNetF(c));


	EtA = Exponential(A, t, precision);
	PiEtA = PP * EtA;
	cout << tr_s[0].x0<< endl;

	PiEtAx0 = PiEtA * tr_s[trNum].x0;
	//PiEtAx0.update();

	// опорная функция PiEtA*x0 - строим пользуясь тем что это просто скалярное произведение
	for (i = 0; i < x0Net.Count; i++)
		if (x0Net.f->v->v[i] != pinMark){
			x0Net.f->v->v[i] = scm(i, PiEtAx0, &x0Net,NULL);
		}

	Net = c + (-1) * x0Net;

	Ind = Net.GetExtrGlobal(opMin,  Ind, min);
	cout << min << endl;
	while (min < 0) {
		t += tau;
		min = -1;
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAk = PiEtA * intEtA;
		calcNextAltInt( trNum, t, /*PiEtA,*/ PiEtAk,PiEtAkB,/*PiEtAkC,*/ c);
		//-----------------
		/*
		PiEtAkB = PiEtAk * B;
		//PiEtAkC = PiEtAk * C;
		tmpPNet = *cP;
		//cP->update();
		//tmpQNet = *cQ;
		// cQ->update();
		tmpPNet.oporn(t, 1);
		//tmpQNet.oporn(t, 1);
		tmpPNet *= PiEtAkB;
		//tmpQNet *= PiEtAkC;
		tmpPNet.update();
		//tmpQNet.update();

		c = c + tmpPNet;
		//c -= tmpQNet; // геометрическая разность
		c.update();
		*/
		PiEtAx0 = PiEtA * tr_s[trNum].x0;
		// опорная функция PiEtA*x0 - строим пользуясь тем что это просто скалярное произведение
		for (i = 0; i < x0Net.Count; i++)
			if (x0Net.f->v->v[i] != pinMark){
				x0Net.f->v->v[i] = scm(i, PiEtAx0, &x0Net,NULL);
			}

	   //	tr_s[trNum].NetList.push_back(new TNetF(c));
		// Сохраняем для последующего поиска Psi
		Net = c + (-1) * x0Net;
		//Net.update();
		Ind = Net.GetExtrGlobal(opMin,   Ind, min);
		tr_s[trNum].psiExtr.push_back(Ind);
		cout << min << " : " << Ind << " : " << t << endl;
	}
	tr_s[trNum].T = t;
	/* */
	return t;
}

// -------------------------------------- control finding----------------------//
void PR_Task::Control_PR(int trNum) {
	TNetF c(PP.m(), perfomance, steps), x_Net(c), p_Net(c), px_Net(c);
	TNetF tmpPNet(*cP), tmpQNet(*cQ);
	long  timeSign=-1,  prevSign;
	unsigned long i, j, k, m, l;
	 long Ind = 0,prevInd, k_up;
	Matrix PiEtA(PP.m(), A.n()), PiEtAC(PP.m(), C.n()), PiEtAB(PP.m(), B.n()),
		EtA(A.m()); // EtA(A.m()) - единичная при t=0
	LDouble t, min, absmin=-1, min_x, extr;
	Vector PiEtAx(PP.m()), psi(PP.v->m), extrVec(PP.v->m);
	Vector u_i(cP->Dim), v_i(cQ->Dim), x_i(dim_x);
	bool isCollisionPossible=false, isInit;

	VecOfVec vx_i, vu_i, vv_i;
	Vector *r_i;

	j = 0;
	k = tr_s[trNum].NetList.size()-1;
	k_up=k;


	x_i = tr_s[trNum].x0; // заполняем x_i  начальным значением
	r_i =  new Vector(x_i);
	r_i->detach();
	vx_i.push_back(r_i);//... и сохраняем его для траектории.
	t = tr_s[trNum].T; // значение t = конечному времени;

	c.perfomance = 0;

	while (k >= 1) {

		if (!isCollisionPossible) {     //если шёл ёжижик, шёл и никого не встретил
			k--;
			t -= tau;
			EtA = Exponential(A, t, precision);  //то идём себе дальше
			PiEtA = PP * EtA;
			PiEtAx = PiEtA * x_i;

			for (i = 0; i < x_Net.Count; i++) {// считаем опорную функцию точки PiEtAx
				x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
				c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
			}

			//выбираем psi
			Ind = c.GetExtrGlobal(opMin, 0, min_x);
		} else{
			//если обходили множества преследователей смотрим, а не проскочили ли  чутка
			//и оказались внутри очередного множества управляемости - корректируем время
//
			isInit = true;
			do{
				EtA = Exponential(A, t, precision);
				PiEtA = PP * EtA;
				PiEtAx = PiEtA * x_i;
				for (i = 0; i < x_Net.Count; i++){
					if (x_Net.f->v->v[i] != pinMark){
						x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
						c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
					}
				}

				Ind = c.GetExtrGlobal(opMin, 0, min_x);

				if (min_x>0) //если попали внуть множетсва - выбираем поменьше
					timeSign=-1;
				else       //если оказались снаружи - побольше
					timeSign=1;
				if(isInit){
					prevInd = Ind;
					prevSign = timeSign;
					isInit = false;
				} /**/
				if (prevSign!=timeSign) {  //проверить на возмножность выноса в условие цикла
                    Ind = prevInd;
					break;
				}
				if((static_cast<long>(k)+timeSign)==k_up){
					//пересчитать Альтернированный интеграл и подвинуть вверх k_up
					calcNextAltInt(trNum, t+timeSign*tau);
					k_up = tr_s[trNum].NetList.size()-1;
				}
				k += timeSign;
				t += timeSign*tau;

			} while (true);/**/
		}
		//фиксируем psi  для множества управляемости, ближайшего к таектории
		for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = c.getIJ(Ind, m);

		//выбираем u_i - как если бы не было противников
		PiEtAB = PiEtA * B;
		PiEtAC = PiEtA * C;
		cP->oporn(t, 1);
		cQ->oporn(t, 1);
		extrVec = Transpose(PiEtAB) * psi;
		cP->perfomance = 0;
		Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cP);
		for (m = 0; m < dim_u; m++)
			u_i.v->v[m] = cP->getIJ(Ind, m);
		u_i = cP->getBorderPoint(Ind, u_i);
		for (m = 0; m < dim_v; m++) //и зануляем v если не приблизились к множеcтву откуда возможно преследование
				v_i.v->v[m] =   0.0;  //набо будет подумать v1 - стандартная помеха + v2 - из соотв.
				//множеств  преследователей

		//ввиду того, что кадо каждый раз перепроверять есть или нет пересечение с множествами откуда
		//возможно завершение преследования  то v выбираем по принципу: если траектория приблизилась
		//к множеству преследователей, то выбираем v  соотвествующую ближайшему множеству преследователя
		tmpPNet = *cP;
		tmpPNet *= PiEtAB;
		absmin = 0.0;
		isCollisionPossible = false;
	   for (l=0; l < PursuerList/**//*NList*/.size(); l++) {
		   px_Net = *PursuerList/**//* *NList*/[l]+ (-1)*x_Net;
		   p_Net = px_Net + tmpPNet;   //проверка принадлежности точки траетории сумме
		   //множества достижимости из точки траекторрии и множества соотв. преследователя
		   p_Net.perfomance = 0;
		   Ind = p_Net.GetExtrGlobal(opMin, Ind, min_x);
		   if(min_x>=0){  //если на следующем шаге есть пересечение с множествами откуда
			//возможно завершение преследования то u выбираем так, чтобы оттолкнуться от множества
			if(!isCollisionPossible){
				isCollisionPossible = true;
				absmin=min_x;
			}
			if(min_x<=absmin){  // отталкиваемся от ближайшего из которых возможно преследование
				absmin = min_x;
				px_Net.perfomance = 0;
				Ind = px_Net.GetExtrGlobal(opMin,  0, min);
				for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = px_Net.getIJ(Ind, m);
				//выбираем u_i - надо доработать и сделать проверку на случай если u_i+1 = -u_i c тем, чтобы всегда было отклонение в сторону главной цели
				extrVec = Transpose(PiEtAB) * psi;
				cP->perfomance = 0;
				Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cP);
				for (m = 0; m < dim_u; m++)
					u_i.v->v[m] = cP->getIJ(Ind, m);
				u_i = cP->getBorderPoint(Ind, u_i);
				//выбираем v_i
//				extrVec = Transpose(PiEtAC) * psi;
//				cQ->perfomance = 0;
//				Ind = cQ->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cQ);
//				for (m = 0; m < dim_v; m++)
//					v_i.v->v[m] =   cQ->getIJ(Ind, m);
//				v_i=cQ->getBorderPoint(Ind,v_i);
			}
		   }

		}   /**/

		cout << (j) * tau << " : "<< t << " : " << x_i<< " : " << u_i;
		//находим следующий  x_i  методом  Рунге-Кутты
		x_i = rungeCutt(x_i, u_i, v_i);

//		if(!isCollisionPossible){  //в противном случае надо проверить  не попали ли внуть множетва управляемости
//			t -= tau;
//			k--;
//		}

		//сохраняем результаты расчётов
		r_i =  new Vector(x_i);
		r_i->detach();
		vx_i.push_back(r_i);
		//cout<<x_i;
		r_i =  new Vector(u_i);
		r_i->detach();
		vu_i.push_back(r_i);
		r_i =  new Vector(v_i);
		r_i->detach();
		vv_i.push_back(r_i);

		j++;
	}

	 /* TODO -orum : Сделать отдеьный конструктор из вектора векторов в матрицу */
	k= vx_i.size()-1;
	tr_s[trNum].x_i = new Matrix(k+1, dim_x);
	// массив векторов со значениями  x_i
	tr_s[trNum].u_i = new Matrix(k, dim_u);
	// массив векторов со значениями  u_i
	tr_s[trNum].v_i = new Matrix(k, dim_v);
	// массив векторов со значениями  v_i

	for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[0][m] = vx_i[0]->v->v[m];
	for (j=0; j < k; j++) {
		//сохраняем результаты расчётов
		for (m = 0; m < dim_u; m++)
			tr_s[trNum].u_i->v->v[j][m] = vu_i[j]->v->v[m];
		for (m = 0; m < dim_v; m++)
			tr_s[trNum].v_i->v->v[j][m] = vv_i[j]->v->v[m];
		for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[j + 1][m] = vx_i[j+1]->v->v[m];
	  // cout<<*vx_i[j+1];
	}

	for_each(vx_i.begin(), vx_i.end(), DeleteObj());
	for_each(vu_i.begin(), vu_i.end(), DeleteObj());
	for_each(vv_i.begin(), vv_i.end(), DeleteObj());
	//--------------------------------------------------------------------
	cout<<k;
}

//-------------------постройка множеств откуда возможно преследование----------
void PR_Task::buildFunnels(pursuerType* pT){

	LDouble t,maxRad;//, min;
	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()), PiEtAkC(C.m(), C.n()),
		PiEtAkB(A.m(), A.n()), EtA(A.m()), intEtA(A.m());
	// EtA(A.m()) - единичная при t=0
	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);

	TNetF c(pT->dim, pT->perfomance, pT->steps), tmpPNet(*cP);
	TNetF tmpQNet(c);//, Net(c),  x0Net(c), tmpNet(c);
	unsigned long i,j;//, Ind = 0;

	Pursuer *tmpPr;
	TNetF *qNet;//=pT->funnel[0];
	string func_m, func_v;
	//TNetF m2Net(c);

	for (j = 0; j < pT->centres.size(); j++) { //строим для всех центров, откуда возможно выполнять преследование
		t = t0;
		maxRad = c.maxRad;
		tmpPr=new Pursuer();
		tmpPr->center = new Vector(*pT->centres[j]);
		func_m = pT->func_m;
		for (i = 0; i < pT->dim; i++)  //опорная функция привязанная к конкретной точке начала преследования
			func_m.append("+p"+ intToStr(i) + "*(" + ldToStr(pT->centres[j]->v->v[i])+")");
		cout << func_m<<endl;
		qNet  =  new TNetF(c);
		qNet->SetFunc(func_m);
		qNet->oporn(t0, 1);

		c = *qNet; // сеточная опорная функция терминального множества преследователя
		//cout<<c;
		tmpPr->funnel.push_back(qNet);

		func_v = pT->func_v;
		qNet  =  new TNetF(c);
		qNet->SetFunc(func_v);
		qNet->oporn(t0, 1);

		// интегрирование методом трапеций по двум точкам
		intEtA = (0.5 * tau) * (Exponential(A, 0, precision) + Exponential(A, -tau,	precision));
		intEtA.update();

		//c = m2Net;
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;

		c.is_empty = false;
		c.perfomance = optNone;
		while (!c.is_empty) {
			t += tau;
			//min = -1;
			EtA = Exponential(A, t, precision);
			PiEtA = PP * EtA;
			PiEtAk = PiEtA * intEtA;
			PiEtAkB = PiEtAk * B;
			PiEtAkC = PiEtAk * C;
			tmpPNet = *cP;
			//cP->update();
			tmpQNet = *qNet;//*cQ;
			// cQ->update();
			tmpPNet.oporn(t, 1);
			tmpQNet.oporn(t, 1);
			//cout<<tmpPNet;
			tmpPNet *= PiEtAkB;
			tmpQNet *= PiEtAkC;
			//	tmpPNet.update();
			// 	tmpQNet.update();

			c = c + tmpQNet;
			c -= tmpPNet; // геометрическая разность
			c.update();
			//cout<<c<<endl;
		  /* */
			if(c.maxRad>maxRad)
				maxRad=c.maxRad;
			cout << c.is_empty<<" : "<< t << endl;
			tmpPr->funnel.push_back(new TNetF(c));
		}
		tmpPr->funnelDepth = t;
		tmpPr->radarVisibility = sqrt(t*t+c.Dim*maxRad*maxRad); //проверить расчёт видимости
		Pursuers.push_back(tmpPr);

		delete qNet;
	}
}


// -------------------------------------- control finding----------------------//
void PR_Task::Control_PR_fullSets(int trNum) {
	TNetF c(PP.m(), perfomance, steps), x_Net(c), p_Net(c), px_Net(c);
	TNetF tmpPNet(*cP), tmpQNet(*cQ);
	long  timeSign=-1, prevSign;
	unsigned long i, j, m, n, l, k;
	long    fi, Ind = 0,prevInd, k_up;
	Matrix PiEtA(PP.m(), A.n()), PiEtAC(PP.m(), C.n()), PiEtAB(PP.m(), B.n()),
		EtA(A.m()); // EtA(A.m()) - единичная при t=0
	LDouble t, min, absmin=-1, min_x, extr;
	Vector PiEtAx(PP.m()), psi(PP.v->m), extrVec(PP.v->m);
	Vector u_i(cP->Dim), v_i(cQ->Dim), x_i(dim_x);
	bool isCollisionPossible=false, isInit;
	Pursuer* ps;

	VecOfVec vx_i, vu_i, vv_i;
	Vector *r_i;

	j = 0;
	k = tr_s[trNum].NetList.size()-1;
	k_up=k;


	x_i = tr_s[trNum].x0; // заполняем x_i  начальным значением
	r_i =  new Vector(x_i);
	r_i->detach();
	vx_i.push_back(r_i);//... и сохраняем его для траектории.
	t = tr_s[trNum].T; // значение t = конечному времени;

	c.perfomance = 0;

	while (k >= 1) {

		if (!isCollisionPossible) {     //если шёл ёжик, шёл и никого не встретил
			k--;
			t -= tau;
			EtA = Exponential(A, t, precision);  //то идём себе дальше
			PiEtA = PP * EtA;
			PiEtAx = PiEtA * x_i;

			for (i = 0; i < x_Net.Count; i++) {// считаем опорную функцию точки PiEtAx
				if (x_Net.f->v->v[i] != pinMark){
					x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
					c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
				}
			}

			//выбираем psi
			Ind = c.GetExtrGlobal(opMin, 0, min_x);
		} else{
			//если обходили множества преследователей смотрим, а не проскочили ли  чутка
			//и оказались внутри очередного множества управляемости - корректируем время
//
			isInit = true;
			do{
				EtA = Exponential(A, t, precision);
				PiEtA = PP * EtA;
				PiEtAx = PiEtA * x_i;
				for (i = 0; i < x_Net.Count; i++){
					if (x_Net.f->v->v[i] != pinMark){
						x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
						c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
					}
				}

				Ind = c.GetExtrGlobal(opMin, 0, min_x);

				if (min_x>0) //если попали внуть множества - выбираем  время поменьше
					timeSign=-1;
				else       //если оказались снаружи - побольше
					timeSign=1;
				if(isInit){
					prevInd = Ind;
					prevSign = timeSign;
					isInit = false;
				} /**/
				if(prevSign!=timeSign){      //проверить на возмножность выноса в условие цикла
					Ind = prevInd;
					break;
				}
				if((static_cast<long>(k)+timeSign)==k_up){
					//пересчитать Альтернированный интеграл и подвинуть вверх k_up
					calcNextAltInt(trNum, t+timeSign*tau);
					k_up = tr_s[trNum].NetList.size()-1;
				}
				k += timeSign;
				t += timeSign*tau;

			} while (true);/**/

		}
		//фиксируем psi  для множества управляемости, ближайшего к таектории
		for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = c.getIJ(Ind, m);

		//выбираем u_i - как если бы не было противников
		PiEtAB = PiEtA * B;
		PiEtAC = PiEtA * C;
		cP->oporn(t, 1);
		cQ->oporn(t, 1);
		extrVec = Transpose(PiEtAB) * psi;
		cP->perfomance = 0;
		Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cP);
		for (m = 0; m < dim_u; m++)
			u_i.v->v[m] = cP->getIJ(Ind, m);
		u_i = cP->getBorderPoint(Ind, u_i);
		for (m = 0; m < dim_v; m++) //и зануляем v если не приблизились к множеcтву откуда возможно преследование
				v_i.v->v[m] =   0.0;  //набо будет подумать v1 - стандартная помеха + v2 - из соотв.
				//множеств  преследователей

		//ввиду того, что надо каждый раз перепроверять есть или нет пересечение с множествами откуда
		//возможно завершение преследования  то v выбираем по принципу: если траектория приблизилась
		//к множеству преследователей, то выбираем v  соотвествующую ближайшему множеству преследователя
		tmpPNet = *cP;
		tmpPNet *= PiEtAB;
		absmin = 0.0;
		isCollisionPossible = false;
		for (l=0; l < Pursuers/* *PursuerList*//*NList*/.size(); l++) {
		   ps = Pursuers[l];
		   //проверить радиус видимости
		   cout << eu_dist(*ps->center, x_i)<<endl;
		   if(eu_dist(*ps->center, x_i) <  ps->radarVisibility){
			fi = ps->funnel.size() - 1;
			if(ps->currInd >0 ){   //т.е. если при t-tau уже натыкались на множество откуда возможно преследование и уклонялись, то
				n= ps->currInd;    //при  надо проверить, накнёмся ли при t
				ps->currInd = 0;
			}else                 //иначе считаем, что все проверки с начала "воронки"
				n=0;
			for (/*n = 0*/; n < ps->funnel.size(); n++) {
				px_Net = *ps->funnel[fi-n] + (-1)*x_Net;
				p_Net = px_Net + tmpPNet;   //проверка принадлежности точки траектории сумме
				//множества достижимости из точки траекторрии и множества соотв. преследователя
				p_Net.perfomance = 0;
				Ind = p_Net.GetExtrGlobal(opMin, Ind, min_x);
				if(min_x>=0){  //если на следующем шаге есть пересечение с множествами откуда
					//возможно завершение преследования то u выбираем так, чтобы оттолкнуться от множества
					if(!isCollisionPossible){
						isCollisionPossible = true;
						absmin=min_x;
					}
					if(min_x<=absmin){  // отталкиваемся от ближайшего из которых возможно преследование
						absmin = min_x;
						px_Net.perfomance = 0;
						Ind = px_Net.GetExtrGlobal(opMin,  0, min);
						for (m = 0; m < c.Dim; m++)
						psi.v->v[m] = px_Net.getIJ(Ind, m);
						//выбираем u_i - надо доработать и сделать проверку на случай если u_i+1 = -u_i c тем, чтобы всегда было отклонение в сторону главной цели
						extrVec = Transpose(PiEtAB) * psi;
						cP->perfomance = 0;
						Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cP);
						for (m = 0; m < dim_u; m++)
							u_i.v->v[m] = cP->getIJ(Ind, m);
						u_i = cP->getBorderPoint(Ind, u_i);
						//выбираем v_i
						//extrVec = Transpose(PiEtAC) * psi;
						//cQ->perfomance = 0;
						//Ind = cQ->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cQ);
						//for (m = 0; m < dim_v; m++)
						//v_i.v->v[m] =   cQ->getIJ(Ind, m);
						//v_i=cQ->getBorderPoint(Ind,v_i);
						ps->currInd = n;  //запоминаем время в которое оттолкнулись от множества из котороно возможно преследование
						break;
					}
				}
			}
		   }
		   //подобрать подходящее множество из воронки
		 //


		}   /**/

		cout << (j) * tau << " : "<< t << " : " << x_i<< " : " << u_i;
		//находим следующий  x_i  методом  Рунге-Кутты
		x_i = rungeCutt(x_i, u_i, v_i);


		//сохраняем результаты расчётов
		r_i =  new Vector(x_i);
		r_i->detach();
		vx_i.push_back(r_i);
		//cout<<x_i;
		r_i =  new Vector(u_i);
		r_i->detach();
		vu_i.push_back(r_i);
		r_i =  new Vector(v_i);
		r_i->detach();
		vv_i.push_back(r_i);

		j++;
	}

	tr_s[trNum].T= (j) * tau;
	 /* TODO -orum : Сделать отдеьный конструктор из вектора векторов в матрицу */
	k= vx_i.size()-1;
	tr_s[trNum].x_i = new Matrix(k+1, dim_x);
	// массив векторов со значениями  x_i
	tr_s[trNum].u_i = new Matrix(k, dim_u);
	// массив векторов со значениями  u_i
	tr_s[trNum].v_i = new Matrix(k, dim_v);
	// массив векторов со значениями  v_i

	for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[0][m] = vx_i[0]->v->v[m];
	for (j=0; j < k; j++) {
		//сохраняем результаты расчётов
		for (m = 0; m < dim_u; m++)
			tr_s[trNum].u_i->v->v[j][m] = vu_i[j]->v->v[m];
		for (m = 0; m < dim_v; m++)
			tr_s[trNum].v_i->v->v[j][m] = vv_i[j]->v->v[m];
		for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[j + 1][m] = vx_i[j+1]->v->v[m];
	  // cout<<*vx_i[j+1];
	}

	for_each(vx_i.begin(), vx_i.end(), DeleteObj());
	for_each(vu_i.begin(), vu_i.end(), DeleteObj());
	for_each(vv_i.begin(), vv_i.end(), DeleteObj());
	//--------------------------------------------------------------------
	cout<<k;
}

//------------------------------------------------------------------------
void PR_Task::plot(int trNum){      //пока тест только для двумерных векторов
//------------------------
Gnuplot gp;

	vector<pair<LDouble, LDouble> > xy_pts_A;
	for(long i=0; i<tr_s[trNum].x_i->m(); ++i) {
	   //	LDouble y = x*x*x;
		xy_pts_A.push_back(make_pair(tr_s[trNum].x_i->v->v[i][0], tr_s[trNum].x_i->v->v[i][1]));
	}

	//vector<pair<LDouble, LDouble> > xy_pts_B;
   //	for(LDouble alpha=0; alpha<1; alpha+=1.0/24.0) {
   //		LDouble theta = alpha*2.0*3.14159;
   //		xy_pts_B.push_back(make_pair(cos(theta), sin(theta)));
   //	}

	//gp << "set xrange [-10:10]\nset yrange [-10:10]\n";
	// Data will be sent via a temporary file.  These are erased when you call
	// gp.clearTmpfiles() or when gp goes out of scope.  If you pass a filename
	// (e.g. "gp.file1d(pts, 'mydata.dat')"), then the named file will be created
	// and won't be deleted (this is useful when creating a script).
	gp << "plot" << gp.file1d(xy_pts_A) << "with lines title 'x(t)'"//<<" ,"
		/*<< gp.file1d(xy_pts_B) << "with points title 'circle'"*/ << endl;

#ifdef _WIN32
	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Press enter to exit." << std::endl;
	std::cin.get();
#endif
//------------------------

}

// -------------------------------------- control finding----------------------//
void Control_PR_fullSets_smooth(int trNum, PR_Task& mt) {//mt - mainTask
	TNetF c(mt.PP.m(), mt.perfomance, mt.steps), x_Net(c), p_Net(c), px_Net(c);
	TNetF tmpPNet(*mt.cP), tmpQNet(*mt.cQ);
	long  timeSign=-1, prevSign;
	unsigned long i, j, m, n, l, k;
	long    fi, Ind = 0,prevInd, k_up;
	Matrix PiEtA(mt.PP.m(), mt.A.n()), PiEtAC(mt.PP.m(), mt.C.n()), PiEtAB(mt.PP.m(), mt.B.n()),
		EtA(mt.A.m()); // EtA(A.m()) - единичная при t=0
	LDouble t, min, absmin=-1, min_x, extr;
	Vector PiEtAx(mt.PP.m()), psi(mt.PP.v->m), extrVec(mt.PP.v->m);
	Vector u_i(mt.cP->Dim), v_i(mt.cQ->Dim), x_i(mt.dim_x);
	bool isCollisionPossible=false, isInit;
	Pursuer* ps;

	VecOfVec vx_i, vu_i, vv_i;
	Vector *r_i;

	VecOfLong vec_j, vec_k; // Значения индексов j и k вдоль траектории
	VecOfBool vec_cp; // Значения  признака collisionPossible вдоль траектории
	alphType vec_t;//значения t вдоль траектории



	j = 0;
	k = mt.tr_s[trNum].NetList.size()-1;
	k_up=k; //"наибольшее" удаление


	x_i = mt.tr_s[trNum].x0; // заполняем x_i  начальным значением
	r_i =  new Vector(x_i);
	r_i->detach();
	vx_i.push_back(r_i);//... и сохраняем его для траектории.
	t = mt.tr_s[trNum].T; // значение t = конечному времени;

	c.perfomance = 0;

	while (k >= 1) {

		if (!isCollisionPossible) {     //если шёл ёжик, шёл и никого не встретил,
			// при отсутствии коллизий повторяем обычный алгоритм без сглаживания
			k--;
			t -= mt.tau;
			EtA = Exponential(mt.A, t, mt.precision);  //то идём себе дальше
			PiEtA = mt.PP * EtA;
			PiEtAx = PiEtA * x_i;

			for (i = 0; i < x_Net.Count; i++) {// считаем опорную функцию точки PiEtAx
				if (x_Net.f->v->v[i] != mt.pinMark){
					x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
					c.f->v->v[i] = mt.tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
				}
			}

			//выбираем psi
			Ind = c.GetExtrGlobal(opMin, 0, min_x);
		} else{
			//если обходили множества преследователей смотрим, а не проскочили ли  чутка
			//и оказались внутри очередного множества управляемости - корректируем время
			// при этом сначала моделируем отскок, а потом прицеливаемся на точку в которую отскочили  из x(t_i-tau)
			isInit = true;
			do{
				EtA = Exponential(mt.A, t, mt.precision);
				PiEtA = mt.PP * EtA;
				PiEtAx = PiEtA * x_i;
				//cout << k << endl;
				for (i = 0; i < x_Net.Count; i++){
					if (x_Net.f->v->v[i] != mt.pinMark){
						x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
						c.f->v->v[i] = mt.tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
					}
				}

				Ind = c.GetExtrGlobal(opMin, 0, min_x);

				if (min_x>0) //если попали внуть множества - выбираем  время поменьше
					timeSign=-1;
				else       //если оказались снаружи - побольше
					timeSign=1;
				if(isInit){
					prevInd = Ind;
					prevSign = timeSign;
					isInit = false;
				} /**/
				if(prevSign!=timeSign){ //смена знака при проверке "снаружи" или "внутри" очередного мн-ва альт. инт. находится точка  и есть граница цикла
					Ind = prevInd;
					break;
				}
				if((static_cast<long>(k)+timeSign)==k_up){  //если дошли до верхней границы предрассчитанных альтернированных интегралов
					//пересчитать Альтернированный интеграл и подвинуть вверх k_up
					mt.calcNextAltInt(trNum, t+timeSign * mt.tau);
					k_up = mt.tr_s[trNum].NetList.size()-1;
				}
				k += timeSign;
				t += timeSign * mt.tau;

			} while (true);/**/

		}
		//фиксируем psi  для множества управляемости, ближайшего к таектории
		for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = c.getIJ(Ind, m);

		//выбираем u_i - как если бы не было противников
		PiEtAB = PiEtA * mt.B;
		PiEtAC = PiEtA * mt.C;
		mt.cP->oporn(t, 1);
		mt.cQ->oporn(t, 1);
		extrVec = Transpose(PiEtAB) * psi;
		mt.cP->perfomance = 0;
		Ind = mt.cP->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *mt.cP);
		for (m = 0; m < mt.dim_u; m++)
			u_i.v->v[m] = mt.cP->getIJ(Ind, m);
		u_i = mt.cP->getBorderPoint(Ind, u_i);
		for (m = 0; m < mt.dim_v; m++) //и зануляем v если не приблизились к множеcтву откуда возможно преследование
				v_i.v->v[m] =   0.0;  //набо будет подумать v1 - стандартная помеха + v2 - из соотв.
				//множеств  преследователей

		//ввиду того, что надо каждый раз перепроверять есть или нет пересечение с множествами откуда
		//возможно завершение преследования  то v выбираем по принципу: если траектория приблизилась
		//к множеству преследователей, то выбираем v  соотвествующую ближайшему множеству преследователя
		tmpPNet = *mt.cP;
		tmpPNet *= PiEtAB;
		absmin = 0.0;
		isCollisionPossible = false;
		for (l=0; l < mt.Pursuers/* *PursuerList*//*NList*/.size(); l++) {
		   ps = mt.Pursuers[l];
		   //проверить радиус видимости
		   cout << eu_dist(*ps->center, x_i)<<endl;
		   if(eu_dist(*ps->center, x_i) <  ps->radarVisibility){
			fi = ps->funnel.size() - 1;
			if(ps->currInd >0 ){   //т.е. если при t-tau уже натыкались на множество откуда возможно преследование и уклонялись, то
				n= ps->currInd;    //при  надо проверить, накнёмся ли при t
				ps->currInd = 0;
			}else                 //иначе считаем, что все проверки с начала "воронки"
				n=0;
			for (/*n = 0*/; n < ps->funnel.size(); n++) {
				px_Net = *ps->funnel[fi-n] + (-1)*x_Net;
				p_Net = px_Net + tmpPNet;   //проверка принадлежности точки траектории сумме
				//множества достижимости из точки траекторрии и множества соотв. преследователя
				p_Net.perfomance = 0;
				Ind = p_Net.GetExtrGlobal(opMin, Ind, min_x);
				if(min_x>=0){  //если на следующем шаге есть пересечение с множествами откуда
					//возможно завершение преследования то u выбираем так, чтобы оттолкнуться от множества
					if(!isCollisionPossible){
						isCollisionPossible = true;
						absmin=min_x;
					}
					if(min_x<=absmin){  // отталкиваемся от ближайшего из которых возможно преследование
						absmin = min_x;
						px_Net.perfomance = 0;
						Ind = px_Net.GetExtrGlobal(opMin,  0, min);
						for (m = 0; m < c.Dim; m++)
						psi.v->v[m] = px_Net.getIJ(Ind, m);
						//выбираем u_i - надо доработать и сделать проверку на случай если u_i+1 = -u_i c тем, чтобы всегда было отклонение в сторону главной цели
						extrVec = Transpose(PiEtAB) * psi;
						mt.cP->perfomance = 0;
						Ind = mt.cP->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *mt.cP);
						for (m = 0; m < mt.dim_u; m++)
							u_i.v->v[m] = mt.cP->getIJ(Ind, m);
						u_i = mt.cP->getBorderPoint(Ind, u_i);
						//выбираем v_i
						//extrVec = Transpose(PiEtAC) * psi;
						//cQ->perfomance = 0;
						//Ind = cQ->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cQ);
						//for (m = 0; m < dim_v; m++)
						//v_i.v->v[m] =   cQ->getIJ(Ind, m);
						//v_i=cQ->getBorderPoint(Ind,v_i);
						ps->currInd = n;  //запоминаем время в которое оттолкнулись от множества из котороно возможно преследование
						break;
					}
				}
			}
		   }
		   //подобрать подходящее множество из воронки
		}   /**/

		cout << (j) * mt.tau << " : "<< t << " : " << x_i<< " : " << u_i;


	  if(isCollisionPossible && /*(vec_k.size()>0)*/(vec_cp[vec_cp.size()-1]!=true)){
			//находим следующий  x_i  методом  Рунге-Кутты т.е. место куда "отскочит" траектория под управлением уклонения
			x_i = mt.rungeCutt(x_i, u_i, v_i);
		   // если препятствие есть и есть куда отклониться, а также отклонение не стачала, то создаём копию объкета, даём начальную точку x_{i-2}, конечную точку считаем x_i
			PR_Task mid_step(mt); //для этого создаём новую задачу
			vector<Traectory>ms_tr_s;
			Traectory tr;
			vector<TNetF*>ms_NetList;
		   //	ms_tr_s.clear();
			mid_step.tr_s = ms_tr_s;
			Vector mid_x_0(*vx_i[j-1]);  //в качестве начальной точки выбираем x(t_{i-1})
			Vector mid_x_t(x_i); //в качестве терминального множества  x(t_{i}) + сферка малого радиуса
			tr =  mt.tr_s[trNum];
			tr.x0 = mid_x_0;
			tr.NetList = ms_NetList;
			mid_step.tr_s.push_back(tr);
		   //	cout <<"Начало: "<< mid_x_0 ;
		 //	cout <<"Конец: "<< mid_x_t << endl;
			int ii;
			string str = "";
			for (ii = 0; ii < mid_step.dim_m-1; ii++) {
				str = str + "p"+ldToStr(ii)+"*"+ldToStr(mid_x_t[ii])+"+";
			}
			str = str + "p"+ldToStr(ii)+"*"+ldToStr(mid_x_t[ii]);
			mid_step.setFuncToNetF(*(mid_step.cM), str);
			mid_step.calcPursuerSets(trNum);
			mid_step.TimeCalc_PR(trNum);
			mid_step.Control_PR_fullSets(trNum);
			u_i = mid_step.tr_s[trNum].u_i->GetRow(0);
			x_i = mid_x_0;
			k  = vec_k[vec_k.size()-1];
			vec_k.pop_back();
			j  = vec_j[vec_j.size()-1];
			vec_j.pop_back();
			t  = vec_t[vec_t.size()-1];

			vec_t.pop_back();
			vx_i.pop_back();
			vu_i.pop_back();
			vv_i.pop_back();
			vec_cp.pop_back();

			//cout<<u_i<<x_i<<endl;
			/**/
		}  /**/

		//находим следующий  x_i  методом  Рунге-Кутты
		x_i = mt.rungeCutt(x_i, u_i, v_i);
		//сохраняем результаты расчётов
		r_i =  new Vector(x_i);
		r_i->detach();
		vx_i.push_back(r_i);
		//cout<<x_i;
		r_i =  new Vector(u_i);
		r_i->detach();
		vu_i.push_back(r_i);
		r_i =  new Vector(v_i);
		r_i->detach();
		vv_i.push_back(r_i);

		vec_j.push_back(j);
		vec_k.push_back(k);
		vec_t.push_back(t);
		vec_cp.push_back(isCollisionPossible);

		j++;
	}

	mt.tr_s[trNum].T= (j) * mt.tau;
	 /* TODO -orum : Сделать отдеьный конструктор из вектора векторов в матрицу */
	k= vx_i.size()-1;
	mt.tr_s[trNum].x_i = new Matrix(k+1, mt.dim_x);    	// массив векторов со значениями  x_i
	mt.tr_s[trNum].u_i = new Matrix(k, mt.dim_u);  	// массив векторов со значениями  u_i
	mt.tr_s[trNum].v_i = new Matrix(k, mt.dim_v);  	// массив векторов со значениями  v_i

	for (m = 0; m < mt.dim_x; m++)
			mt.tr_s[trNum].x_i->v->v[0][m] = vx_i[0]->v->v[m];
	for (j=0; j < k; j++) {
		//сохраняем результаты расчётов
		for (m = 0; m < mt.dim_u; m++)
			mt.tr_s[trNum].u_i->v->v[j][m] = vu_i[j]->v->v[m];
		for (m = 0; m < mt.dim_v; m++)
			mt.tr_s[trNum].v_i->v->v[j][m] = vv_i[j]->v->v[m];
		for (m = 0; m < mt.dim_x; m++)
			mt.tr_s[trNum].x_i->v->v[j + 1][m] = vx_i[j+1]->v->v[m];
	  // cout<<*vx_i[j+1];
	}

	for_each(vx_i.begin(), vx_i.end(), DeleteObj());
	for_each(vu_i.begin(), vu_i.end(), DeleteObj());
	for_each(vv_i.begin(), vv_i.end(), DeleteObj());
	//--------------------------------------------------------------------
	cout<<k;
}

