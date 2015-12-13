//---------------------------------------------------------------------------

#pragma hdrstop

#include "pr_task.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
// ------------------------------ desttructor --------------------------------//

PR_Task::~PR_Task() {
	for_each(PursuerList.begin(), PursuerList.end(), DeleteObj());
}

// ------------------------------ constructor --------------------------------//
PR_Task::PR_Task(long dimX, long dimU, long dimV, long dimM, LDouble ts, double prec,
	double eps, LDouble delta, LDouble tmax, long st, short perf, long stat,
	int tr_count):Task(dimX,  dimU,  dimV,  dimM,  ts,  prec,
	 eps,  delta,  tmax,  st,  perf,  stat,  tr_count) {
	}

// ------ расчёт множеств альтернированных интегралов преследователей --------//
LDouble PR_Task::calcPursuerSets(int trNum){
TNetF m2Net(*PursuerList[0]); // сеточная опорная функция терминального множества
// для первого преследователя
	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()), PiEtAkC(C.m(), C.n()),
		PiEtAkB(A.m(), A.n()), EtA(A.m()), intEtA(A.m());
	// EtA(A.m()) - единичная при t=0
	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);
	LDouble t = t0, min;
	TNetF c(PP.m(), perfomance, steps), tmpPNet(*cP), tmpQNet(*cQ), Net(c),
		x0Net(c), tmpNet(c);
	long i,j, Ind = 0;

    //костыль
	//просто рассчитываем опорную функцию, пользуясь тем, что управление константа не зависящая от времени
	for(j=0;j<PursuerList.size();j++){
	   PursuerList[j]->oporn(0.0,1);
	 //  cout<<*PursuerList[j]<<endl;
	}


	/*
    // интегрирование методом трапеций по двум точкам
	intEtA = (0.5 * tau) * (Exponential(A, 0, precision) + Exponential(A, -tau,
		precision));
	intEtA.update();

	m2Net.oporn(t0, 1); // опорное значение сетки по M2
	c = m2Net;
	tr_s[trNum].NetList.push_back(new TNetF(c));


	EtA = Exponential(A, t, precision);
	PiEtA = PP * EtA;
	PiEtAx0 = PiEtA * tr_s[trNum].x0;


	Net = c;

	Ind = Net.GetExtrGlobal(opMax,  Ind, min);
	cout << min << endl;
	while (min < 0) {
		t += tau;
		min = -1;
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAk = PiEtA * intEtA;
		PiEtAkB = PiEtAk * B;
		PiEtAkC = PiEtAk * C;
		PiEtAx0 = PiEtA * tr_s[trNum].x0;
		tmpPNet = *cP;
		// cP->update();
		tmpQNet = *cQ;
		// cQ->update();
		tmpPNet.oporn(t, 1);
		tmpQNet.oporn(t, 1);
		tmpPNet *= PiEtAkB;
		tmpQNet *= PiEtAkC;
		tmpPNet.update();
		tmpQNet.update();

		c = c + tmpPNet;
		c -= tmpQNet; // геометрическая разность
		c.update();

		// опорная функция PiEtA*x0 - строим пользуясь тем что это просто скалярное произведение
		for (i = 0; i < x0Net.Count; i++)
			x0Net.f->v->v[i] = scm(i, PiEtAx0, &x0Net,NULL);

		tr_s[trNum].NetList.push_back(new TNetF(c));
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

  // ------------------------------- Time calc-----------------------------------//
LDouble PR_Task::TimeCalc_PR(int trNum) {
	TNetF m2Net(*cM); // сеточная опорная функция терминального множества
	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()), //PiEtAkC(C.m(), C.n()),
		PiEtAkB(A.m(), A.n()), EtA(A.m()), intEtA(A.m());
	// EtA(A.m()) - единичная при t=0
	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);
	LDouble t = t0, min;
	TNetF c(PP.m(), perfomance, steps), tmpPNet(*cP), /*tmpQNet(*cQ),*/ Net(c),
		x0Net(c), tmpNet(c);
	long i, Ind = 0;

    // интегрирование методом трапеций по двум точкам
	intEtA = (0.5 * tau) * (Exponential(A, 0, precision) + Exponential(A, -tau,
		precision));
	intEtA.update();

	m2Net.oporn(t0, 1); // опорное значение сетки по M2
	c = m2Net;
	tr_s[trNum].NetList.push_back(new TNetF(c));


	EtA = Exponential(A, t, precision);
	PiEtA = PP * EtA;
	PiEtAx0 = PiEtA * tr_s[trNum].x0;
	//PiEtAx0.update();

	// опорная функция PiEtA*x0 - строим пользуясь тем что это просто скалярное произведение
	for (i = 0; i < x0Net.Count; i++)
		x0Net.f->v->v[i] = scm(i, PiEtAx0, &x0Net,NULL);

	Net = c + (-1) * x0Net;

	Ind = Net.GetExtrGlobal(opMin,  Ind, min);
	cout << min << endl;
	while (min < 0) {
		t += tau;
		min = -1;
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAk = PiEtA * intEtA;
		PiEtAkB = PiEtAk * B;
		//PiEtAkC = PiEtAk * C;
		PiEtAx0 = PiEtA * tr_s[trNum].x0;
		tmpPNet = *cP;
		// cP->update();
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
		//c.update();

		// опорная функция PiEtA*x0 - строим пользуясь тем что это просто скалярное произведение
		for (i = 0; i < x0Net.Count; i++)
			x0Net.f->v->v[i] = scm(i, PiEtAx0, &x0Net,NULL);

		tr_s[trNum].NetList.push_back(new TNetF(c));
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
	long i, j, k, m, l, Ind = 0;
	Matrix PiEtA(PP.m(), A.n()), PiEtAC(PP.m(), C.n()), PiEtAB(PP.m(), B.n()),
		EtA(A.m()); // EtA(A.m()) - единичная при t=0
	LDouble t, min, absmin,  extr;
	Vector PiEtAx(PP.m()), psi(PP.v->m), extrVec(PP.v->m);
	Vector u_i(cP->Dim), v_i(cQ->Dim), x_i(dim_x);
	//bool extr_exist;

	VecOfVec vx_i, vu_i, vv_i;
	Vector *r_i;

	j = 0;
	k = tr_s[trNum].NetList.size()-1;


	x_i = tr_s[trNum].x0; // заполняем x_i  начальным значением
	r_i =  new Vector(x_i);
	r_i->detach();
	vx_i.push_back(r_i);//... и сохраняем его для траектории.
	t = tr_s[trNum].T; // значение t = конечному времени;

	while (t >= precision) {
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAx = PiEtA * x_i;

		for (i = 0; i < x_Net.Count; i++) {// считаем опорную функцию точки PiEtAx
				x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
				c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
		}

		//выбираем psi
		c.perfomance = 0;
		Ind = c.GetExtrGlobal(opMin,  0, min);
	   //	for (m = 0; m < c.Dim; m++)
	   //		psi.v->v[m] = c.getIJ(Ind, m);

		while (min>=0){  //если оказались внутри очередного множества управляемости - корректируем время
			k--;
			t -= tau;
			EtA = Exponential(A, t, precision);
			PiEtA = PP * EtA;
			PiEtAx = PiEtA * x_i;
			for (i = 0; i < x_Net.Count; i++){ //смотрим, а не проскочили ли ещё чутка
				x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
				c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
            }
			Ind = c.GetExtrGlobal(opMin,  0, min);
		}

		//фиксируем psi
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
		for (m = 0; m < dim_v; m++) //и зануляем v если не приблизились к множетву откуда возможно преследование
				v_i.v->v[m] =   0.0;

		//ввиду того, что кадо каждый раз перепроверять есть или нет пересечение с множествами откуда
		//возможно завершение преследования  то v выбираем по принципу: если траектория приблизилась
		//к множеству преследователей, то выбираем v  соотвествующую ближайшему множеству преследователя
		tmpPNet = *cP;
		tmpPNet *= PiEtAB;
		absmin = -1.0;
		for (l=0; l < PursuerList.size(); l++) {
		   px_Net = *PursuerList[l];
		   p_Net = px_Net + tmpPNet + (-1)*x_Net;   //проверка принадлежности точки траетории сумме
		   //множества достижимости из точки траекторрии и множества соотв. преследователя
		   Ind = p_Net.GetExtrGlobal(opMin,   Ind, min);
		   if((min>=0)&&(absmin<min)){  //если на следующем шаге есть пересечение с множествами откуда
			//возможно завершение преследования то u выбираем так, чтобы оттолкнуться от множества
			absmin=min; // отталкиваемся от ближайшего из которых возможно преследование
			px_Net.perfomance = 0;
			Ind = px_Net.GetExtrGlobal(opMin,  0, min);
			for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = px_Net.getIJ(Ind, m);
			//выбираем u_i - надо доработать и сделать проверку на случай если u_i+1 = -u_i c тем, чтобы всегда было отклонение в сторону главной цели
			extrVec = Transpose(PiEtAB) * psi;
			cP->perfomance = 0;
			Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cP);
			for (m = 0; m < dim_u; m++)
				u_i.v->v[m] = cP->getIJ(Ind, m);
			u_i = cP->getBorderPoint(Ind, u_i);
			//выбираем v_i
			extrVec = Transpose(PiEtAC) * psi;
			cQ->perfomance = 0;
			Ind = cQ->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cQ);
			for (m = 0; m < dim_v; m++)
				v_i.v->v[m] =   cQ->getIJ(Ind, m);
			v_i=cQ->getBorderPoint(Ind,v_i);
		   }

		}

		cout << (j) * tau << " : "<< t << " : " << x_i;
		//находим следующий  x_i  методом  Рунге-Кутты
		x_i = rungeCutt(x_i, u_i, v_i);

		 if(absmin <0){  //исправлять  т.к. данное условие не надёжно т.е. если препятствия нет то норм.,
		 //а вот если препятствие  есть нужны доп. проверки.
			t -= tau;
			k--;
		}

		// в противном случаенадо проверять в каком "времени" находимся !!!
		// для тестов пока не сделано - потом доделать

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
		//for (m = 0; m < dim_u; m++)
		//	tr_s[trNum].u_i->v->v[j][m] = u_i[m];
		//for (m = 0; m < dim_v; m++)
		//	tr_s[trNum].v_i->v->v[j][m] = v_i[m];
		//for (m = 0; m < dim_x; m++)
		//	tr_s[trNum].x_i->v->v[j + 1][m] = x_i[m];

		j++;
	}

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
