// ---------------------------------------------------------------------------
#pragma hdrstop

#include "task.h"
#include "Netfunc.h"
// ---------------------------------------------------------------------------
#pragma package(smart_init)

// ---------------------------- constructor -----------------------------------//
Task::Task() : A(1, 1), B(1, 1), C(1, 1), PP(1, 1), /* x0(1), */ Level(0) {
	Vector x_0(1);
	Traectory tr;
	tr.x0 = x_0;
	tr_s.push_back(tr);
	tau = 0.1;
	perfomance = 0;
	precision = 0.001;
	epsilon = 0.01;
	t0 = 0;
	maxTime = 10;
	steps = 10;
	v_control = 0;
	// traectoryCount=1;
}
// ---------------------------- constructor -----------------------------------//

Task::Task(long dimX, long dimU, long dimV, long dimM, LDouble ts, double prec,
	double eps, LDouble delta, LDouble tmax, long st, short perf, long stat,
	int tr_count) : A(dimX, dimX), B(dimX, dimU), C(dimX, dimV),
	PP(dimM, dimX) {
	Traectory tr;
	Vector x_0(dimX);
	tr.x0 = x_0;
	for (int i = 0; i < tr_count; i++)
		tr_s.push_back(tr);

	t0 = ts;
	perfomance = perf;
	state = (TaskState)stat;
	precision = prec;
	epsilon = eps;
	tau = delta;
	maxTime = tmax;
	currMu = 0;
	currT = ts;
	dim_x = dimX;
	dim_u = dimU;
	dim_v = dimV;
	dim_m1 = dimM;
	dim_m = dimM;
	steps = st;
	cP = new TNetF(dim_u, perf, st, "0");
	cQ = new TNetF(dim_v, perf, st, "0");
	cM = new TNetF(dim_m, perf, st, "0");
	Level = 2;
	IndI = 1;
	v_control = 0;
	
}

// -----------------------------destructor----------------------------------//
Task::~Task() {
	delete cP;
	cP = NULL;
	delete cQ;
	cQ = NULL;
	delete cM;
	cM = NULL;
}

// ---------------------------- copy constructor ------------------------------//
Task::Task(const Task& ts) : A(ts.A.m(), ts.A.n()), B(ts.B.m(), ts.B.n()),
	C(ts.C.m(), ts.C.n()), PP(ts.PP.m(), ts.PP.n()) {
	long i, size = tr_s.size();

	tr_s.clear();
	for (i = 0; i < (long)ts.tr_s.size(); i++)
		tr_s[i] = (Traectory)ts.tr_s[i];

	perfomance = ts.perfomance;
	state = ts.state;
	A = ts.A;
	B = ts.B;
	C = ts.C;
	PP = ts.PP;
	currT = ts.currT;
	currMu = ts.currMu;
	t0 = ts.t0;
	precision = ts.precision;
	epsilon = ts.epsilon;
	tau = ts.tau;
	maxTime = ts.maxTime;
	dim_x = ts.dim_x;
	dim_u = ts.dim_u;
	dim_v = ts.dim_v;
	dim_m1 = ts.dim_m1;
	dim_m = ts.dim_m;
	steps = ts.steps;
	Level = ts.Level;
	v_control = ts.v_control;
}

// ---------------------------- destructor ------------------------------------//
// ------------------------------- Time calc-----------------------------------//
/*LDouble Task::TimeCalc() {
}
  /**/
// ------------------------------- Time calc-----------------------------------//
LDouble Task::TimeCalc_Pontryagin(int trNum) {
	TNetF m2Net(*cM); // сеточная опорная функция терминального множества
	Matrix PiEtA(A.m(), A.n()), PiEtAC(C.m(), C.n()), PiEtAB(A.m(), A.n()),
		EtA(A.m()); // EtA(A.m()) - единичная при t=0
	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);
	LDouble t = t0, min;
	TNetF c(PP.m(), perfomance, steps), tmpPNet(*cP), tmpQNet(*cQ), Net(c),
		x0Net(c), tmpNet(c);
	long i, Ind = 0;
   //	pathType path;

	m2Net.oporn(t0, 1); // опорное значение сетки по M2
	tr_s[trNum].NetList.push_back(new TNetF(m2Net));

	EtA = Exponential(A, t, precision);
	PiEtA = PP * EtA;
	PiEtAx0 = PiEtA * tr_s[trNum].x0;
	PiEtAx0.update();

	// опорная функция PiEtA*x0 - строим пользуясь тем что это просто скалярное произведение
	for (i = 0; i < x0Net.Count; i++)
		x0Net.f->v->v[i] = scm(i, PiEtAx0, &x0Net,NULL);

	c = m2Net + (-1) * x0Net;

	Ind = c.GetExtrGlobal(opMin, Ind, min);
	cout << min << endl;
	while (min < 0) {
		t += tau;
		min = -1;
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAB = PiEtA * B;
		PiEtAC = PiEtA * C;
		PiEtAx0 = PiEtA * tr_s[trNum].x0;
		tmpPNet = *cP;
		cP->update();
		tmpQNet = *cQ;
		cQ->update();
		tmpPNet.oporn(t, -1);
		tmpQNet.oporn(t, -1);
		tmpPNet *= PiEtAB;
		tmpQNet *= PiEtAC;
		Net = tmpPNet;
		Net -= tmpQNet;
		// геометрическая разность  - пока не оптимально считается
		// - можно писать   tmpPNet - tmpQNet, но так меньше операций копирвоания

		tmpNet += (tau) * Net;
		// накопление интеграла сейчас считается по сути методом Эйлера, можно попробовать методом Рунге-Кутты
	   //	tmpNet *= (-1);
		tmpNet.update();
		tr_s[trNum].NetList.push_back(new TNetF(tmpNet));
		// Сохраняем P-Q для последующего поиска Psi

		m2Net.oporn(t, 1);
		// опорная функция PiEtA*x0 - строим пользуясь тем что это просто скалярное произведение
		for (i = 0; i < x0Net.Count; i++)
			x0Net.f->v->v[i] = scm(i, PiEtAx0, &x0Net,NULL);

		c = m2Net + tmpNet + (-1) * x0Net;
		// Вот на этой сетке и ищем мин, для того, чтобы найти T
		c.update();
		Ind = c.GetExtrGlobal(opMin,  Ind, min);
		tr_s[trNum].psiExtr.push_back(Ind);
		cout << min << " : " << Ind << " : " << t << endl;
	}
	tr_s[trNum].T = t; // запоминаем T

	return t;
}

// ------------------------------- Time calc-----------------------------------//
LDouble Task::TimeCalc_AltInt(int trNum) {
	TNetF m2Net(*cM); // сеточная опорная функция терминального множества
	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()), PiEtAkC(C.m(), C.n()),
		PiEtAkB(A.m(), A.n()), EtA(A.m()), intEtA(A.m());
	// EtA(A.m()) - единичная при t=0
	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);
	LDouble t = t0, min;
	TNetF c(PP.m(), perfomance, steps), tmpPNet(*cP), tmpQNet(*cQ), Net(c),
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
		c -= tmpQNet; // геометрическая разность  - пока не оптимально считается
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

// ------------------------------------- rungeCutt ----------------------------//
//Интегрирование для стандартного шага по времени
Vector Task::rungeCutt(const Vector& xn, const Vector& un, const Vector& vn) {
	Vector K1(A.m()), K2(A.m()), K3(A.m()), K4(A.m()), result(A.m());

	K1 = (A * xn + (-1)* B * un + C * vn);
	K2 = (A * (xn + 0.5 * tau * K1) + (-1)* B * un + C * vn);
	K3 = (A * (xn + 0.5 * tau * K2) + (-1)* B * un + C * vn);
	K4 = (A * (xn + tau * K3) + (-1)* B * un + C * vn);
	result = xn + tau / 6 * (K1 + 2 * (K2 + K3) + K4);
	result.update();
	return result;
	/* */
	// return tau*(A*xn+B*un+C*vn)+xn; //Euler method
}
// ------------------------------------- rungeCutt ----------------------------//
// Интегрирование для НЕстандартного шага по времени
Vector Task::rungeCutt(const Vector& xn, const Vector& un, const Vector& vn, LDouble tau) {
	Vector K1(A.m()), K2(A.m()), K3(A.m()), K4(A.m()), result(A.m());

	K1 = (A * xn + (-1)* B * un + C * vn);
	K2 = (A * (xn + 0.5 * tau * K1) + (-1)* B * un + C * vn);
	K3 = (A * (xn + 0.5 * tau * K2) + (-1)* B * un + C * vn);
	K4 = (A * (xn + tau * K3) + (-1)* B * un + C * vn);
	result = xn + tau / 6 * (K1 + 2 * (K2 + K3) + K4);
	result.update();
	return result;
	/* */
	// return tau*(A*xn+B*un+C*vn)+xn; //Euler method
}

// -------------------------------------- control finding----------------------//
void Task::Control_AltInt(int trNum) {
	TNetF c(PP.m(), perfomance, steps), x_Net(c);
	TNetF tmpPNet(*cP), tmpQNet(*cQ);
	long i, j, k, m, Ind = 0;
	Matrix PiEtA(PP.m(), A.n()), PiEtAC(PP.m(), C.n()), PiEtAB(PP.m(), B.n()),
		EtA(A.m()); // EtA(A.m()) - единичная при t=0
	LDouble t, min,  extr;
	Vector PiEtAx(PP.m()), psi(PP.v->m), extrVec(PP.v->m);
	Vector u_i(cP->Dim), v_i(cQ->Dim), x_i(dim_x);
	//bool extr_exist;

	j = 0;
	k = tr_s[trNum].NetList.size();
	tr_s[trNum].x_i = new Matrix(k, dim_x);
	// массив векторов со значениями  x_i
	tr_s[trNum].u_i = new Matrix(k - 1, dim_u);
	// массив векторов со значениями  u_i
	tr_s[trNum].v_i = new Matrix(k - 1, dim_v);
	// массив векторов со значениями  v_i
	k--;
	for (m = 0; m < dim_x; m++)
		tr_s[trNum].x_i->v->v[0][m] = tr_s[trNum].x0[m];
	x_i = tr_s[trNum].x0; // заполняем x_i  начальным значением
	t = tr_s[trNum].T; // значение t = конечному времени;

	while (t >= precision) {
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAB = PiEtA * B;
		PiEtAC = PiEtA * C;
		PiEtAx = PiEtA * x_i;

		cP->oporn(t, 1);
		cQ->oporn(t, 1);
		for (i = 0; i < x_Net.Count; i++) {// считаем опорную функцию точки PiEtAx
				x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
				c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
		}
		k--;
		//выбираем psi
		c.perfomance = 0;
		Ind = c.GetExtrGlobal(opMin,  0, min);
		for (m = 0; m < c.Dim; m++)
			psi.v->v[m] = c.getIJ(Ind, m);

		//выбираем u_i
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
		cout << (j) * tau << " : " << x_i;

		//находим следующий  x_i  методом  Рунге-Кутты
		x_i = rungeCutt(x_i, u_i, v_i);

		//сохраняем результаты расчётов
		for (m = 0; m < dim_u; m++)
			tr_s[trNum].u_i->v->v[j][m] = u_i[m];
		for (m = 0; m < dim_v; m++)
			tr_s[trNum].v_i->v->v[j][m] = v_i[m];
		for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[j + 1][m] = x_i[m];

		t -= tau;
		j++;
	}

}

// -------------------------------------- control finding----------------------//
void Task::Control_R1(int trNum) {
	TNetF c(PP.m(), perfomance, steps), x_Net(c);
	TNetF tmpPNet(PP.m(), perfomance, steps), tmpQNet(PP.m(), perfomance, steps);
	long i,j,  k, m, Ind = 0;
	Matrix PiEtA(PP.m(), A.n()), PiEtAC(PP.m(), C.n()), PiEtAB(PP.m(), B.n()),
		EtA(A.m()); // EtA(A.m()) - единичная при t=0
	LDouble t, min,  extr;
	Vector PiEtAx(PP.m()), psi(PP.v->m), extrVec(PP.v->m);
	Vector u_i(cP->Dim), v_i(cQ->Dim), x_i(dim_x);
	LDouble  tau_s, tau_delta, xNorm, xExtNorm;
	VecOfVec vx_i, vu_i, vv_i;
	Vector *r_i;

	LDouble tmin = _extr_tmin_param, tmax = _extr_t0_param, tcurr = tmax, p, a =
		LDouble(_lrand() % cP->Count) / cP->Count, ccurr = _extr_e_val;

	long jk=-1,jj, prevInd,indExtr;

	j = 0;
	k = tr_s[trNum].NetList.size()-1;

	x_i = tr_s[trNum].x0; // заполняем x_i  начальным значением
	r_i =  new Vector(x_i);
	r_i->detach();
	vx_i.push_back(r_i);//... и сохраняем его для траектории.

	t = tr_s[trNum].T; // значение t = конечному времени;
	tau_s = tau/(tmpPNet.Count*tmpPNet.Count);//пока так - потом посчитаем на сколько надо делить


	while (t >= precision) {
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAB = PiEtA * B;
		PiEtAC = PiEtA * C;
		PiEtAx = PiEtA * x_i;

		cP->oporn(t, 1);
		cQ->oporn(t, 1);
		for (i = 0; i < x_Net.Count; i++) {// считаем опорную функцию точки PiEtAx
				x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
				c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
		}
		k--;
		//выбираем psi
		c.perfomance = 0;
		Ind = c.GetExtrGlobal(opMin,  0, min);
		for (m = 0; m < c.Dim; m++)
			psi.v->v[m] = c.getIJ(Ind, m);
		//prevInd = Ind;
		//выбираем u_i
		extrVec = Transpose(PiEtAB) * psi;
		cP->perfomance = 0;
		Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cP);
		for (m = 0; m < dim_u; m++)
			u_i.v->v[m] = cP->getIJ(Ind, m);
		u_i = cP->getBorderPoint(Ind, u_i);
		prevInd = Ind;

		cout << Ind << " : " << u_i;

		//выбираем v_i
		extrVec = Transpose(PiEtAC) * psi;
		cQ->perfomance = 0;
		Ind = cQ->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cQ);
		for (m = 0; m < dim_v; m++)
			v_i.v->v[m] =   cQ->getIJ(Ind, m);
		v_i=cQ->getBorderPoint(Ind,v_i);

		if (jk<0) {
            jk = _lrand() % (c.Count);
		}
		xExtNorm = c.f->v->v[jk];
		xNorm =  xExtNorm;
		indExtr = jk;

		 //--------ищем u_i при условии отсутствия информации о системе и наличия данных только о расстоянии до терм. множества
		tau_delta=tau;
		tmin = _extr_tmin_param; tcurr = tmax;
		while (tcurr>tmin) {

			xNorm = c.f->v->v[jk];
			//шаг метода annealing
			if(xNorm < xExtNorm){
				xExtNorm = xNorm;
				indExtr = jk;
				tcurr = tcurr * ccurr;
			}else{
				p = 1.0 / (1.0 + exp(-fabs(xNorm - xExtNorm) /tcurr));
				a = double(_lrand() % c.Count) / (double)c.Count;
				if (a > p) {
					xExtNorm = xNorm;
					indExtr = jk;
					tcurr = tcurr * ccurr;
				}
			}
			a = double(_lrand() % c.Count) / double(c.Count);
			jj = signof(a - 0.5) * tcurr * (pow((1.0 + 1.0 / tcurr), fabs(2.0 * a - 1.0)) -	1.0) * double(c.Count);
			jk = abs(jj + jk) % c.Count;
			for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = c.getIJ(indExtr, m);

			cP->oporn(t, 1);
			extrVec = Transpose(PiEtAB) * psi;
			cP->perfomance = 0;
			Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cP);
			for (m = 0; m < dim_u; m++)
				u_i.v->v[m] = cP->getIJ(Ind, m);
			u_i = cP->getBorderPoint(Ind, u_i);

			x_i = rungeCutt(x_i, u_i, v_i, tau_s);
			t -=tau_s;
			tau_delta -= tau_s;

		   r_i =  new Vector(x_i);
		   r_i->detach();
		   vx_i.push_back(r_i);
		   r_i =  new Vector(u_i);
		   r_i->detach();
		   vu_i.push_back(r_i);
		   r_i =  new Vector(v_i);
		   r_i->detach();
		   vv_i.push_back(r_i);
		}

		for (m = 0; m < dim_v; m++)
			u_i.v->v[m] = cP->getIJ(Ind, m);

		u_i = cP->getBorderPoint(Ind, u_i);
		cout << Ind << " : " << u_i;

		x_i = rungeCutt(x_i, u_i, v_i, tau_delta);
		cout<<xExtNorm<<" : "<< x_i <<endl;
		t -=tau_delta;

	   //----------------------------------------------------------
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
	k= vx_i.size();
	tr_s[trNum].x_i = new Matrix(k+1, dim_x);
	// массив векторов со значениями  x_i
	tr_s[trNum].u_i = new Matrix(k, dim_u);
	// массив векторов со значениями  u_i
	tr_s[trNum].v_i = new Matrix(k, dim_v);
	// массив векторов со значениями  v_i

	for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[0][m] = vx_i[0]->v->v[m];
	for (j=0; j < k-1; j++) {
		//сохраняем результаты расчётов
		for (m = 0; m < dim_u; m++)
			tr_s[trNum].u_i->v->v[j][m] = vu_i[j]->v->v[m];
		for (m = 0; m < dim_v; m++)
			tr_s[trNum].v_i->v->v[j][m] = vv_i[j]->v->v[m];
		for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[j + 1][m] = vx_i[j+1]->v->v[m];
	}
	for_each(vx_i.begin(), vx_i.end(), DeleteObj());
	for_each(vu_i.begin(), vu_i.end(), DeleteObj());
	for_each(vv_i.begin(), vv_i.end(), DeleteObj());
	//--------------------------------------------------------------------
	cout<<k;
}

// -------------------------------------- control finding----------------------//
void Task::Control_R2(int trNum) {
	TNetF c(PP.m(), perfomance, steps), x_Net(c);
	TNetF tmpPNet(PP.m(), perfomance, steps), tmpQNet(PP.m(), perfomance, steps);
	long i,j,  k, m, Ind = 0;
	Matrix PiEtA(PP.m(), A.n()), PiEtAC(PP.m(), C.n()), PiEtAB(PP.m(), B.n()),
		EtA(A.m()); // EtA(A.m()) - единичная при t=0
	LDouble t, min,  extr;
	Vector PiEtAx(PP.m()), psi(PP.v->m), extrVec(PP.v->m);
	Vector u_i(cP->Dim), v_i(cQ->Dim), x_i(dim_x);
	bool extr_exist;

	bool  isNotExtrFound, borderChanged;
	LDouble  tau_s, tau_delta, xNorm, xExtNorm;
	VecOfVec vx_i, vu_i, vv_i;
	Vector *r_i;

	long jk=-1, prevInd,indExtr;
	int currPlateNorm, currDirection=-1, prevDirection = -1, moveSign = 1, exitSign = 0, exitLim = (c.Dim-1)*2;

	j = 0;
	k = tr_s[trNum].NetList.size()-1;

	x_i = tr_s[trNum].x0; // заполняем x_i  начальным значением
	r_i =  new Vector(x_i);
	r_i->detach();
	vx_i.push_back(r_i);//... и сохраняем его для траектории.

	t = tr_s[trNum].T; // значение t = конечному времени;
	tau_s = tau/(tmpPNet.Count*tmpPNet.Count);//пока так - потом посчитаем на сколько надо делить


	while (t >= precision) {
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAB = PiEtA * B;
		PiEtAC = PiEtA * C;
		PiEtAx = PiEtA * x_i;

		cP->oporn(t, 1);
		cQ->oporn(t, 1);
		for (i = 0; i < x_Net.Count; i++) {// считаем опорную функцию точки PiEtAx
				x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
				c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
		}
		k--;
		//выбираем psi
		c.perfomance = 0;
		Ind = c.GetExtrGlobal(opMin,  0, min);
		for (m = 0; m < c.Dim; m++)
			psi.v->v[m] = c.getIJ(Ind, m);
		//выбираем u_i  чисто для контроля работы метода - в принципе данный кусок можно комменитровать
		extrVec = Transpose(PiEtAB) * psi;
		cP->perfomance = 0;
		Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cP);
		for (m = 0; m < dim_u; m++)
			u_i.v->v[m] = cP->getIJ(Ind, m);
		u_i = cP->getBorderPoint(Ind, u_i);
		prevInd = Ind;

		cout << Ind << " : " << u_i;

		//выбираем v_i
		extrVec = Transpose(PiEtAC) * psi;
		cQ->perfomance = 0;
		Ind = cQ->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cQ);
		for (m = 0; m < dim_v; m++)
			v_i.v->v[m] =   cQ->getIJ(Ind, m);
		v_i=cQ->getBorderPoint(Ind,v_i);

		if (jk<0){      //первый узел по любому случаен
			jk = _lrand() % (c.Count);
			indExtr = jk;
			currDirection= rand()%c.Dim;
		}

		xExtNorm = c.f->v->v[indExtr];
		xNorm = xExtNorm;
		jk = indExtr; //по ходу игры начальный узел поиска со следующего шага берём лучший с предыдущего шага
		currPlateNorm = jk / (c.NumOfPoints * 2);
        if (currDirection==currPlateNorm)
			currDirection= (currDirection+1)%c.Dim;
		borderChanged = false;

		 //--------ищем u_i при условии отсутствия информации о системе и наличия данных только о расстоянии до терм. множества
		extr_exist = false;
		tau_delta=tau;

		isNotExtrFound = true;
        exitSign  = 0;
		while (isNotExtrFound) {

			//шаг градиентного метода
			if (exitSign<exitLim) {
               c.wasteCache();
			   prevInd = jk;
				//jk= 67; currDirection = 0;
				jk = c.shift(jk, currDirection, moveSign, borderChanged);  //проверить jk= 66, currDirection = 1;
			   if(borderChanged){ //при переходе на другую грань меняем соотв. образом нормаль
				currPlateNorm = jk / (c.NumOfPoints * 2);
				if (currDirection==currPlateNorm) //Если же следующее направление совпадает с нормалью просто выбираем слеюующее напраавление
						currDirection= (currPlateNorm+1)%(c.Dim);
				exitSign  = 0;
                borderChanged = false;
			   }
			   if (prevInd == jk) {   //попытка сдвига не туда - МЕНЯЕМ НАПРАВЛЕНИЕ
				currDirection= (currDirection+1)%(c.Dim-1);
				if ((prevDirection == currDirection)&&(c.Dim>2))//если в данном направлении уже ходили и есть куда поворачивать (Dim>2) - переходим на следующее
						currDirection = (currDirection +1)%c.Dim;
					else
						moveSign *= -1; //если не ходили или размерность 2, то просто меняем знак  движения
			   }
			   else{
				xNorm = c.f->v->v[jk];

				if (xNorm <= xExtNorm) {
					if (xNorm < xExtNorm){
						xExtNorm = xNorm;
						indExtr = jk;
					}else
						exitSign++;
				  //	exitSign = 0;//если нашли куда куда ходить, то  счётчик сбрасываем
				}
				else{//меняем направление
					if ((prevDirection == currDirection)&&(c.Dim>2))//если в данном направлении уже ходили и есть куда поворачивать (Dim>2) - переходим на следующее
						currDirection = currDirection +1%c.Dim;
					else
						moveSign *= -1; //если не ходили или размерность 2, то просто меняем знак  движения
					if (currDirection==currPlateNorm) //Если же следующее направление совпадает с нормалью просто выбираем слеюующее напраавление
						currDirection= (currPlateNorm+1)%c.Dim;
					exitSign++;
				}
			   }
			  prevDirection = currDirection; //запоминаем предыдущее направление
			}
			else{
				if(exitSign>=exitLim){ //всюду потыкались и не нашли экстремума - выходим
					isNotExtrFound = false;
				   //	seekPath.clear();
				}
			}

			for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = c.getIJ(indExtr, m);

			cP->oporn(t, 1);
			extrVec = Transpose(PiEtAB) * psi;
			cP->perfomance = 0;
			Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMax, nZAware, Ind, NULL,  extr,  *cP);
            /**/
			for (m = 0; m < dim_u; m++)
				u_i.v->v[m] = cP->getIJ(Ind, m);
			u_i = cP->getBorderPoint(Ind, u_i);

			x_i = rungeCutt(x_i, u_i, v_i, tau_s);
		  //	cout<< xExtNorm <<" : "<< x_i <<endl;
			t -=tau_s;
			tau_delta -= tau_s;

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
		}

		for (m = 0; m < dim_v; m++)
			u_i.v->v[m] = cP->getIJ(Ind, m);

		u_i = cP->getBorderPoint(Ind, u_i);
		cout << Ind << " : " << u_i;

		x_i = rungeCutt(x_i, u_i, v_i, tau_delta);
		cout<<xExtNorm<<" : "<< x_i <<endl;
		t -=tau_delta;

	   //----------------------------------------------------------
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

// -------------------------------------- control finding----------------------//
void Task::Control_Pontryagin(int trNum) {
	TNetF c(PP.m(), perfomance, steps), x_Net(c);
	TNetF tmpPNet(*cP), tmpQNet(*cQ);
	long i, j, k, m, Ind = 0;
	Matrix PiEtA(PP.m(), A.n()), PiEtAC(PP.m(), C.n()), PiEtAB(PP.m(), B.n()),
		EtA(A.m()); // EtA(A.m()) - единичная при t=0
	LDouble t, min, val, extr;
	Vector PiEtAx(PP.m()), psi(PP.v->m), extrVec(PP.v->m);
	Vector u_i(cP->Dim), v_i(PP.m()), x_i(dim_x);
	bool extr_exist;

	j = 0;
	k = tr_s[trNum].NetList.size();
	tr_s[trNum].x_i = new Matrix(k, dim_x);
	// массив векторов со значениями  x_i
	tr_s[trNum].u_i = new Matrix(k - 1, dim_u);
	// массив векторов со значениями  u_i
	tr_s[trNum].v_i = new Matrix(k - 1, dim_v);
	// массив векторов со значениями  v_i
	k--;
	for (m = 0; m < dim_x; m++)
		tr_s[trNum].x_i->v->v[0][m] = tr_s[trNum].x0[m];
	x_i = tr_s[trNum].x0; // заполняем x_i  начальным значением
	t = tr_s[trNum].T; // значение t = конечному времени;

	while (t >= precision) {
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAB = PiEtA * B;
		PiEtAC = PiEtA * C;
		PiEtAx = PiEtA * x_i;

		for (i = 0; i < x_Net.Count; i++) // считаем опорную функцию точки PiEtAx
				x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
		c = *tr_s[trNum].NetList[k] + (-1) * x_Net;
		k--;
		 c.perfomance = 0;
		 Ind = c.GetExtrGlobal(opMin,  Ind, min);

		for (m = 0; m < c.Dim; m++)
			psi.v->v[m] = c.getIJ(Ind, m);
		tmpQNet = *cQ;
		tmpQNet.oporn(t, 1);
		tmpQNet *= PiEtAC;
		tmpQNet.update();
		extr_exist = false;
		extrVec = /* Transpose(PiEtAC)* */ psi;
		for (i = 0; i < tmpQNet.Count; i++) {
			val = scm(i, extrVec, &tmpQNet,NULL);
			if (!extr_exist) {
				extr_exist = true;
				extr = val;
				Ind = i;
			}
			else {
				if (val > extr) {
					extr = val;
					Ind = i;
				}
			}
		} /* */
		for (m = 0; m < dim_v; m++)
			v_i.v->v[m] = tmpQNet.getIJ(Ind, m);
		v_i = tmpQNet.getBorderPoint(Ind, v_i);
		cout << Ind << " : " << v_i;
		v_i += psi;

		u_i = Solve(PiEtAB,v_i,epsilon);
		cout << "U_I : " << u_i;

		cout << (j) * tau << " : " << x_i;
		x_i = rungeCutt(x_i, u_i, v_i);


		for (m = 0; m < dim_u; m++)
			tr_s[trNum].u_i->v->v[j][m] = u_i[m];
		for (m = 0; m < dim_v; m++)
			tr_s[trNum].v_i->v->v[j][m] = v_i[m];
		for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[j + 1][m] = x_i[m];

			t -= tau;
		j++;
	}

}

// -------------------------------------- control finding  ----------------------//
  /*
void Task::Control() {

	Level = 2;
}
  /**/
// ------------------------------ saveTask ------------------------------------//
void Task::saveTask(char *fileName) {
	long len, i, size = tr_s.size();
	// TNet Net(dim_u,perfomance);
	ofstream out_f(fileName);

	out_f << "(Method);" << '\n';
	out_f << '[' << priority << ']' << ';' << '\n' << '\n';

	len = description.length();
	out_f << "(Description," << len << ')' << ';' << '\n' << '[';
	out_f << description;
	out_f << ']' << ';' << '\n' << '\n';

	out_f << "(Priority);" << '\n';
	out_f << '[' << priority << ']' << ';' << '\n' << '\n';

	out_f << "(Level);" << '\n';
	out_f << '[' << Level << ']' << ';' << '\n' << '\n';

	out_f << "(t0);" << '\n';
	out_f << '[' << t0 << ']' << ';' << '\n' << '\n';
	for (i = 0; i < size; i++) {
		out_f << "(x0," << tr_s.size() << ',' << tr_s[i].x0.size()
			<< ')' << ';' << '\n';
		out_f << tr_s[i].x0 << '\n';
	}

	out_f << "(A" << ',' << A.m() << ',' << A.n() << ')' << ';' << '\n';
	out_f << A << '\n';

	out_f << "(B" << ',' << B.m() << ',' << B.n() << ')' << ';' << '\n';
	out_f << B << '\n';

	out_f << "(C" << ',' << C.m() << ',' << C.n() << ')' << ';' << '\n';
	out_f << C << '\n';

	out_f << "(PP" << ',' << PP.m() << ',' << PP.n() << ')' << ';' << '\n';
	out_f << PP << '\n';

	out_f << "(precision);" << '\n';
	out_f << '[' << precision << ']' << ';' << '\n' << '\n';

	out_f << "(tau);" << '\n';
	out_f << '[' << tau << ']' << ';' << '\n' << '\n';

	out_f << "(maxTime);" << '\n';
	out_f << '[' << maxTime << ']' << ';' << '\n' << '\n';

	out_f << "(dim_x);" << '\n';
	out_f << '[' << dim_x << ']' << ';' << '\n' << '\n';

	out_f << "(dim_u);" << '\n';
	out_f << '[' << dim_u << ']' << ';' << '\n' << '\n';

	out_f << "(dim_v);" << '\n';
	out_f << '[' << dim_v << ']' << ';' << '\n' << '\n';

	out_f << "(dim_m1);" << '\n';
	out_f << '[' << dim_m1 << ']' << ';' << '\n' << '\n';

	out_f << "(dim_m);" << '\n';
	out_f << '[' << dim_m << ']' << ';' << '\n' << '\n';

	out_f << "(steps);" << '\n';
	out_f << '[' << steps << ']' << ';' << '\n' << '\n';

	// out_f<<"(ResNum);"<<'\n';
	// out_f<<'['<<ResNum<<']'<<';'<<'\n'<<'\n';

	len = cP->fStr.length();
	out_f << "(fP," << len << ')' << ';' << '\n' << '[';
	out_f << cP->fStr;
	out_f << ']' << ';' << '\n' << '\n';

	len = cQ->fStr.length();
	out_f << "(fQ," << len << ')' << ';' << '\n' << '[';
	out_f << cQ->fStr;
	out_f << ']' << ';' << '\n' << '\n';

	len = cM->fStr.length();
	out_f << "(fM," << len << ')' << ';' << '\n' << '[';
	out_f << cM->fStr;
	out_f << ']' << ';' << '\n' << '\n';

	out_f.flush();
	// if(Level>0)
	// {
	/* if (Level>1)
	 {
	 tmpNet=(TNet*)NetList->Items[0];
	 out_f<<"(W,"<<NetList->Count<<','<<tmpNet->BaseDim<<','<<tmpNet->Count<<')'<<'\n';
	 len=NetList->Count;
	 for(i=0;i<len;i++)
	 {
	 tmpNet=(TNet*)NetList->Items[i];
	 Net=(*tmpNet);
	 out_f<<Net;
	 }
	 out_f<<'\n';
	 } */

	/* Vector *vec;
	 vec=(Vector*)selList->Items[0];
	 out_f<<"(Selector,"<<selList->Count<<','<<vec->size<<')'<<'\n';
	 len=selList->Count;
	 for(i=0;i<len;i++)
	 {
	 vec=(Vector*)selList->Items[i];
	 out_f<<(*vec);
	 }
	 out_f<<'\n';

	 vec=(Vector*)contList->Items[0];
	 out_f<<"(u,"<<contList->Count<<','<<vec->size<<')'<<'\n';
	 len=contList->Count;
	 for(i=0;i<len;i++)
	 {
	 vec=(Vector*)contList->Items[i];
	 out_f<<(*vec);
	 }
	 out_f<<'\n';

	 vec=(Vector*)vList->Items[0];
	 out_f<<"(v,"<<vList->Count<<','<<vec->size<<')'<<'\n';
	 len=vList->Count;
	 for(i=0;i<len;i++)
	 {
	 vec=(Vector*)vList->Items[i];
	 out_f<<(*vec);
	 }
	 out_f<<'\n';

	 }
	/* */
	// out_f.flush();
	out_f.close();
}

// ------------------------------ loadTask ------------------------------------//
Task* Task::loadTask(/*istream& in_f*/string szFileName) {
	char c, buf[20], long_buf[2048];
	int Lev, len, vecs, m, n, dx, du, dv, dm, stp, i, perf, meth, prior;
	LDouble ts, prec, tt, maxT, eps;
	// TNet *tmpNet;
	// Vector *tmpVec;
	string *tmpStr, descr;
	vector<Traectory>trs;

	fstream in_f(szFileName.c_str(), ios::in|ios::nocreate);

	if ( in_f.is_open() )
	{
		// Код метода (пока не реализовано, на будущее)
		in_f.get(c);
		while (c != '[')
			in_f.get(c);
		in_f.get(buf, 20, ']');
		meth = atoi(buf);
		while (c != ';')
			in_f.get(c);

		// Наименование примера
		in_f.get(c);
		// while (c != ',')
		// in_f.get(c);
		// in_f >> len;
		while (c != '[')
			in_f.get(c);
		// for (i = 0; i < len; i++)
		// in_f >> long_buf[i];
		// long_buf[i] = '\0';
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

		Task* res = new Task(dx, du, dv, dm, ts, prec, eps, tt, maxT, stp, perf,
			Lev, 1);
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
		// while (c != ',')
		// in_f.get(c);
		// in_f >> len;
		while (c != '[')
			in_f.get(c);
		// for (i = 0; i < len; i++)
		// in_f >> long_buf[i];
		// long_buf[i] = '\0';
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);
		res->setFuncToNetF(*(res->cP), *tmpStr);
		delete tmpStr;

		// Значение опорной функции убегающего
		// while (c != ',')
		in_f.get(c);
		// in_f >> len;
		while (c != '[')
			in_f.get(c);
		// for (i = 0; i < len; i++)
		// in_f >> long_buf[i];
		// long_buf[i] = '\0';
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);
		res->setFuncToNetF(*(res->cQ), *tmpStr);
		delete tmpStr;

		// Значение опорной функции терм. множества
		// while (c != ',')
		in_f.get(c);
		// in_f >> len;
		while (c != '[')
			in_f.get(c);
		// for (i = 0; i < len; i++)
		// in_f >> long_buf[i];
		// long_buf[i] = '\0';
		in_f.get(long_buf, 256, ']');
		tmpStr = new string(long_buf);
		res->setFuncToNetF(*(res->cM), *tmpStr);
		delete tmpStr;
		/* */

	return res;
	}
	else
	{
		cout<< " Unable to open file. Does it exist?"<<endl;;
		return false;
	}
}

/* istream& operator >>(istream& in_f, Task* res)
 {
 return in_f;
 }
/* */
// ---------------------- changes default function zero value--------------------
void __fastcall Task::setFuncToNetF(TNetF& net, string func) {
	net.SetFunc(func);
}

long Task::ImageUp() {
	// TODO: Add your source code here

	return 1;
}
