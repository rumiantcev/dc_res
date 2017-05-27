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
	//------------------------------������ �������� ���������������� ���������� ���������������
void PR_Task::Find_Ns(int trNum) {
	vector<TNetF*> tmplist;     //'input'
	tmplist = PursuerList;

   //	vector<TNetF*> NList;       //'output'// ���������� ������ � �����, ����� ��������� �� ���������� ��� ������ �� �������
	//bool *L;
	LDouble t = t0;
	Matrix EtA(A.m()), PiEtA(A.m(), A.n()), PiEtAkB(B.m(), B.n()), PiEtAkC(C.m(), C.n()), intEtA(A.m());
	Matrix PiEtAk(A.m(), A.n());//���� ������
	//Vector PiEtAx0(A.m()), psi(A.v->m), Fs(A.v->m), vmin(A.v->m);

	TNetF /* *c,*/ tmpPNet(*cP), tmpQNet(*cQ), tmpNet(PP.m(), perfomance, steps);

	//long i, j;

	intEtA = (0.5 * tau) * (Exponential(A, 0, precision) + Exponential(A, -tau, precision));

	while (!tmplist.empty()) {

		TNetF NiNet(*(tmplist.back()));
		NiNet.oporn(t0, 1);

	   //	cout<<NiNet<<endl;

		t = t0;   // �����  ���������������� ��� ������� ��-��
	   //	EtA = Exponential(A, t, precision);  // ������������ �����. ������� ���������������  ������� ��� ���������� � ���� *1
	   //	PiEtA = PP * EtA * intEtA;

		while (1) {
			t += tau;

			EtA = Exponential(A, t, precision);  // *1 ��� ����
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
			NiNet -= tmpPNet; // �������������� ��������
		 //	cout<<NiNet<<endl;
		  /* */
		 //	cout << NiNet.is_empty << endl;
			if (tmpNet.is_empty) {
				break;
			}

		   //	NiNet += (tau) * tmpNet;
		  //	NiNet.update();
		}
		// L = new bool[NiNet.Count];    //� ������ ������� ����. �������� ��� �������� �����������
	   //	NiNet.Conv(L);
	   //	delete[]L;
	  	NList.push_back(&NiNet);
	  //	if (c!=NULL) delete c;// c ���, ����� ��������� ��������� ��, ��� �� �������� �� tmpList
		tmplist.pop_back();

	   //	cout << t << endl;
	}

}

// ------ ������ �������� ���������������� ���������� ��������������� - ������� --------//
void PR_Task::calcPursuerSets(int trNum){
//	TNetF m2Net(*PursuerList[0]); // �������� ������� ������� ������������� ���������
// ��� ������� ��������������
  //	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()), PiEtAkC(C.m(), C.n()),
  //		PiEtAkB(A.m(), A.n()), EtA(A.m()), intEtA(A.m());
	// EtA(A.m()) - ��������� ��� t=0
  //	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);
   //	LDouble t = t0, min;
 //	TNetF c(PP.m(), perfomance, steps), tmpPNet(*cP), tmpQNet(*cQ), Net(c),
 //		x0Net(c), tmpNet(c);
    unsigned long /*i,*/j/*, Ind = 0*/;
	pursuerType *pt;

   /*
	//�������
	//������ ������������ ������� �������, ��������� ���, ��� ���������� ��������� �� ��������� �� �������
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
//--------------------������ ���� ����������������� ���������-----------------//
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
		//c -= tmpQNet; // �������������� ��������
		c.update();
		tr_s[trNum].NetList.push_back(new TNetF(c));
		// ��������� ��� ������������ ������ Psi
}

  // ------------------------------- Time calc--------------------------------//
LDouble PR_Task::TimeCalc_PR(int trNum) {
	TNetF m2Net(*cM); // �������� ������� ������� ������������� ���������
	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()), //PiEtAkC(C.m(), C.n()),
		PiEtAkB(A.m(), A.n()), EtA(A.m()), intEtA(A.m());
	// EtA(A.m()) - ��������� ��� t=0
	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);
	LDouble t = t0, min;
	TNetF c(PP.m(), perfomance, steps), tmpPNet(*cP), /*tmpQNet(*cQ),*/ Net(c),
		x0Net(c), tmpNet(c);
	unsigned  long i;
	long Ind = 0;

	if (tr_s[trNum].NetList.size()>0){// ���� ������� �������� �������� - �������� ���������� ����������
		tr_s[trNum].NetList.clear();
		tr_s[trNum].psiExtr.clear();
	}

	// �������������� ������� �������� �� ���� ������
	intEtA = (0.5 * tau) * (Exponential(A, 0, precision) + Exponential(A, -tau, precision));
	intEtA.update();

	m2Net.oporn(t0, 1); // ������� �������� ����� �� M2
	c = m2Net;
	tr_s[trNum].NetList.push_back(new TNetF(c));


	EtA = Exponential(A, t, precision);
	PiEtA = PP * EtA;
   //	cout << tr_s[0].x0<< endl;

	PiEtAx0 = PiEtA * tr_s[trNum].x0;
	//PiEtAx0.update();

	// ������� ������� PiEtA*x0 - ������ ��������� ��� ��� ��� ������ ��������� ������������
	x0Net = point_oporn(PiEtAx0, x0Net, i);
	Net = c + (-1) * x0Net;

	Ind = Net.GetExtrGlobal(opMin,  Ind, min);
   //	cout << min << endl;
	while (min < 0) {
		t += tau;
		min = -1;
		EtA = Exponential(A, t, precision);
		PiEtA = PP * EtA;
		PiEtAk = PiEtA * intEtA;
		calcNextAltInt( trNum, t, /*PiEtA,*/ PiEtAk,PiEtAkB,/*PiEtAkC,*/ c);
	
		PiEtAx0 = PiEtA * tr_s[trNum].x0;
		// ������� ������� PiEtA*x0 - ������ ��������� ��� ��� ��� ������ ��������� ������������
		x0Net = point_oporn(PiEtAx0, x0Net, i);
		Net = c + (-1) * x0Net;
		//Net.update();
		Ind = Net.GetExtrGlobal(opMin,   Ind, min);
		tr_s[trNum].psiExtr.push_back(Ind);
	 //	cout << min << " : " << Ind << " : " << t << endl;
	}
	if (t <= tau)
		return -1;

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
		EtA(A.m()); // EtA(A.m()) - ��������� ��� t=0
	LDouble t, min, absmin=-1, min_x, extr;
	Vector PiEtAx(PP.m()), psi(PP.v->m), extrVec(PP.v->m);
	Vector u_i(cP->Dim), v_i(cQ->Dim), x_i(dim_x);
	bool isCollisionPossible=false, isInit;

	VecOfVec vx_i, vu_i, vv_i;
	Vector *r_i;

	j = 0;
	k = tr_s[trNum].NetList.size()-1;
	k_up=k;


	x_i = tr_s[trNum].x0; // ��������� x_i  ��������� ���������
	r_i =  new Vector(x_i);
	r_i->detach();
	vx_i.push_back(r_i);//... � ��������� ��� ��� ����������.
	t = tr_s[trNum].T; // �������� t = ��������� �������;

	c.perfomance = 0;

	while (k >= 1) {

		if (!isCollisionPossible) {     //���� ��� ������, ��� � ������ �� ��������
			k--;
			t -= tau;
			EtA = Exponential(A, t, precision);  //�� ��� ���� ������
			PiEtA = PP * EtA;
			PiEtAx = PiEtA * x_i;

			pointGeomDiff(trNum, k, PiEtAx, c, x_Net);// ������� ������� ������� ����� PiEtAx
			//�������� psi
			Ind = c.GetExtrGlobal(opMin, 0, min_x);
		} else{
			//���� �������� ��������� ��������������� �������, � �� ���������� ��  �����
			//� ��������� ������ ���������� ��������� ������������� - ������������ �����
//
			isInit = true;
			do{
				EtA = Exponential(A, t, precision);
				PiEtA = PP * EtA;
				PiEtAx = PiEtA * x_i;
				pointGeomDiff(trNum, k, PiEtAx, c, x_Net);// ������� ������� ������� ����� PiEtAx

				Ind = c.GetExtrGlobal(opMin, 0, min_x);

				if (min_x>0) //���� ������ ����� ��������� - �������� ��������
					timeSign=-1;
				else       //���� ��������� ������� - ��������
					timeSign=1;
				if(isInit){
					prevInd = Ind;
					prevSign = timeSign;
					isInit = false;
				} /**/
				if (prevSign!=timeSign) {  //��������� �� ������������ ������ � ������� �����
                    Ind = prevInd;
					break;
				}
				if((static_cast<long>(k)+timeSign)==k_up){
					//����������� ���������������� �������� � ��������� ����� k_up
					calcNextAltInt(trNum, t+timeSign*tau);
					k_up = tr_s[trNum].NetList.size()-1;
				}
				k += timeSign;
				t += timeSign*tau;

			} while (true);/**/
		}
		//��������� psi  ��� ��������� �������������, ���������� � ���������
		for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = c.getIJ(Ind, m);

		//�������� u_i - ��� ���� �� �� ���� �����������
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
		for (m = 0; m < dim_v; m++) //� �������� v ���� �� ������������ � �����c��� ������ �������� �������������
				v_i.v->v[m] =   0.0;  //���� ����� �������� v1 - ����������� ������ + v2 - �� �����.
				//��������  ���������������

		//����� ����, ��� ���� ������ ��� ������������� ���� ��� ��� ����������� � ����������� ������
		//�������� ���������� �������������  �� v �������� �� ��������: ���� ���������� ������������
		//� ��������� ���������������, �� �������� v  �������������� ���������� ��������� ��������������
		tmpPNet = *cP;
		tmpPNet *= PiEtAB;
		absmin = 0.0;
		isCollisionPossible = false;
	   for (l=0; l < PursuerList/**//*NList*/.size(); l++) {
		   px_Net = *PursuerList/**//* *NList*/[l]+ (-1)*x_Net;
		   p_Net = px_Net + tmpPNet;   //�������� �������������� ����� ��������� �����
		   //��������� ������������ �� ����� ����������� � ��������� �����. ��������������
		   p_Net.perfomance = 0;
		   Ind = p_Net.GetExtrGlobal(opMin, Ind, min_x);
		   if(min_x>=0){  //���� �� ��������� ���� ���� ����������� � ����������� ������
			//�������� ���������� ������������� �� u �������� ���, ����� ������������ �� ���������
			if(!isCollisionPossible){
				isCollisionPossible = true;
				absmin=min_x;
			}
			if(min_x<=absmin){  // ������������� �� ���������� �� ������� �������� �������������
				absmin = min_x;
				px_Net.perfomance = 0;
				Ind = px_Net.GetExtrGlobal(opMin,  0, min);
				for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = px_Net.getIJ(Ind, m);
				//�������� u_i - ���� ���������� � ������� �������� �� ������ ���� u_i+1 = -u_i c ���, ����� ������ ���� ���������� � ������� ������� ����
				extrVec = Transpose(PiEtAB) * psi;
				cP->perfomance = 0;
				Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cP);
				for (m = 0; m < dim_u; m++)
					u_i.v->v[m] = cP->getIJ(Ind, m);
				u_i = cP->getBorderPoint(Ind, u_i);
				//�������� v_i
//				extrVec = Transpose(PiEtAC) * psi;
//				cQ->perfomance = 0;
//				Ind = cQ->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cQ);
//				for (m = 0; m < dim_v; m++)
//					v_i.v->v[m] =   cQ->getIJ(Ind, m);
//				v_i=cQ->getBorderPoint(Ind,v_i);
			}
		   }

		}   /**/

		//cout << (j) * tau << " : "<< t << " : " << x_i<< " : " << u_i;
		//������� ���������  x_i  �������  �����-�����
		x_i = rungeCutt(x_i, u_i, v_i);

		storeLocResults(u_i, v_i, x_i, vx_i, vu_i, vv_i, r_i);//��������� ������������� ���������� ��������

		j++;
	}

	storeResults(trNum, vx_i, vu_i, vv_i); //���������  ���������� ��������
	//--------------------------------------------------------------------
	//cout<<k;
}

//-------------------��������� �������� ������ �������� �������������----------
void PR_Task::buildFunnels(pursuerType* pT){

	LDouble t,maxRad;//, min;
	Matrix PiEtA(A.m(), A.n()), PiEtAk(A.m(), A.n()), PiEtAkC(C.m(), C.n()),
		PiEtAkB(A.m(), A.n()), EtA(A.m()), intEtA(A.m());
	// EtA(A.m()) - ��������� ��� t=0
	Vector PiEtAx0(A.m()), psi(A.v->m), extrVec(PP.v->m);

	TNetF c(pT->dim, pT->perfomance, pT->steps), tmpPNet(*cP);
	TNetF tmpQNet(c);//, Net(c),  x0Net(c), tmpNet(c);
	unsigned long i,j;//, Ind = 0;

	Pursuer *tmpPr;
	TNetF *qNet;//=pT->funnel[0];
	string func_m, func_v;
	//TNetF m2Net(c);

	for (j = 0; j < pT->centres.size(); j++) { //������ ��� ���� �������, ������ �������� ��������� �������������
		t = t0;
		maxRad = c.maxRad;
		tmpPr=new Pursuer();
		tmpPr->center = new Vector(*pT->centres[j]);
		func_m = pT->func_m;
		for (i = 0; i < pT->dim; i++)  //������� ������� ����������� � ���������� ����� ������ �������������
			func_m.append("+p"+ intToStr(i) + "*(" + ldToStr(pT->centres[j]->v->v[i])+")");
	   //	cout << func_m<<endl;
		qNet  =  new TNetF(c);
		qNet->SetFunc(func_m);
		qNet->oporn(t0, 1);

		c = *qNet; // �������� ������� ������� ������������� ��������� ��������������
		//cout<<c;
		tmpPr->funnel.push_back(qNet);

		func_v = pT->func_v;
		qNet  =  new TNetF(c);
		qNet->SetFunc(func_v);
		qNet->oporn(t0, 1);

		// �������������� ������� �������� �� ���� ������
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
			c -= tmpPNet; // �������������� ��������
			c.update();
			//cout<<c<<endl;
		  /* */
			if(c.maxRad>maxRad)
				maxRad=c.maxRad;
		//	cout << c.is_empty<<" : "<< t << endl;
			tmpPr->funnel.push_back(new TNetF(c));
		}
		tmpPr->funnelDepth = t;
		tmpPr->radarVisibility = sqrt(t*t+c.Dim*maxRad*maxRad); //��������� ������ ���������
		Pursuers.push_back(tmpPr);

		delete qNet;
	}
}


// -------------------------------------- control finding----------------------//
int PR_Task::Control_PR_fullSets(int trNum) {
	TNetF c(PP.m(), perfomance, steps), x_Net(c), p_Net(c), px_Net(c);
	TNetF tmpPNet(*cP), tmpQNet(*cQ);
	long  timeSign=-1, prevSign;
	unsigned long i, j, m, n, l, k;
	long    fi, Ind = 0,prevInd, k_up;
	Matrix PiEtA(PP.m(), A.n()), PiEtAC(PP.m(), C.n()), PiEtAB(PP.m(), B.n()),
		EtA(A.m()); // EtA(A.m()) - ��������� ��� t=0
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

	x_i = tr_s[trNum].x0; // ��������� x_i  ��������� ���������
	r_i =  new Vector(x_i);
	r_i->detach();
	vx_i.push_back(r_i);//... � ��������� ��� ��� ����������.
	t = tr_s[trNum].T; // �������� t = ��������� �������;

	c.perfomance = 0;

	while (k >= 1) {

		if (!isCollisionPossible) {     //���� ��� ����, ��� � ������ �� ��������
			k--;
			t -= tau;
			EtA = Exponential(A, t, precision);  //�� ��� ���� ������
			PiEtA = PP * EtA;
			PiEtAx = PiEtA * x_i;

			pointGeomDiff(trNum, k, PiEtAx, c, x_Net);// ������� ������� ������� ����� PiEtAx

			//�������� psi
			Ind = c.GetExtrGlobal(opMin, 0, min_x);
		} else{
			//���� �������� ��������� ��������������� �������, � �� ���������� ��  �����
			//� ��������� ������ ���������� ��������� ������������� - ������������ �����
//
			isInit = true;
			do{
				EtA = Exponential(A, t, precision);
				PiEtA = PP * EtA;
				PiEtAx = PiEtA * x_i;

				pointGeomDiff(trNum, k, PiEtAx, c, x_Net);// ������� ������� ������� ����� PiEtAx

				Ind = c.GetExtrGlobal(opMin, 0, min_x);

				if (min_x>0) //���� ������ ����� ��������� - ��������  ����� ��������
					timeSign=-1;
				else       //���� ��������� ������� - ��������
					timeSign=1;
				if(isInit){
					prevInd = Ind;
					prevSign = timeSign;
					isInit = false;
				} /**/
				if(prevSign!=timeSign){      //��������� �� ������������ ������ � ������� �����
					Ind = prevInd;
					break;
				}
				if((static_cast<long>(k)+timeSign)==k_up){
					//����������� ���������������� �������� � ��������� ����� k_up
					calcNextAltInt(trNum, t+timeSign*tau);
					k_up = tr_s[trNum].NetList.size()-1;
				}
				k += timeSign;
				t += timeSign*tau;

			} while (true);/**/

		}
		//��������� psi  ��� ��������� �������������, ���������� � ���������
		for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = c.getIJ(Ind, m);

		//�������� u_i - ��� ���� �� �� ���� �����������
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
		for (m = 0; m < dim_v; m++) //� �������� v ���� �� ������������ � �����c��� ������ �������� �������������
				v_i.v->v[m] =   0.0;  //���� ����� �������� v1 - ����������� ������ + v2 - �� �����.
				//��������  ���������������

		//����� ����, ��� ���� ������ ��� ������������� ���� ��� ��� ����������� � ����������� ������
		//�������� ���������� �������������  �� v �������� �� ��������: ���� ���������� ������������
		//� ��������� ���������������, �� �������� v  �������������� ���������� ��������� ��������������
		tmpPNet = *cP;
		tmpPNet *= PiEtAB;
		absmin = 0.0;
		isCollisionPossible = false;
		for (l=0; l < Pursuers/* *PursuerList*//*NList*/.size(); l++) {
		   ps = Pursuers[l];
		   //��������� ������ ���������
		  // cout << eu_dist(*ps->center, x_i)<<endl;
		   if(eu_dist(*ps->center, x_i) <  ps->radarVisibility){
			fi = ps->funnel.size() - 1;
			if(ps->currInd >0 ){   //�.�. ���� ��� t-tau ��� ���������� �� ��������� ������ �������� ������������� � ����������, ��
				n= ps->currInd;    //���  ���� ���������, ������� �� ��� t
				ps->currInd = 0;
			}else                 //����� �������, ��� ��� �������� � ������ "�������"
				n=0;
			for (/*n = 0*/; n < ps->funnel.size(); n++) {
				px_Net = *ps->funnel[fi-n] + (-1)*x_Net;
				p_Net = px_Net + tmpPNet;   //�������� �������������� ����� ���������� �����
				//��������� ������������ �� ����� ����������� � ��������� �����. ��������������
				p_Net.perfomance = 0;
				Ind = p_Net.GetExtrGlobal(opMin, Ind, min_x);
				if(min_x>=0){  //���� �� ��������� ���� ���� ����������� � ����������� ������
					//�������� ���������� ������������� �� u �������� ���, ����� ������������ �� ���������
					if(!isCollisionPossible){
						isCollisionPossible = true;
						absmin=min_x;
					}
					if(min_x<=absmin){  // ������������� �� ���������� �� ������� �������� �������������
						absmin = min_x;
						px_Net.perfomance = 0;
						Ind = px_Net.GetExtrGlobal(opMin,  0, min);
						for (m = 0; m < c.Dim; m++)
						psi.v->v[m] = px_Net.getIJ(Ind, m);
						//�������� u_i - ���� ���������� � ������� �������� �� ������ ���� u_i+1 = -u_i c ���, ����� ������ ���� ���������� � ������� ������� ����
						extrVec = Transpose(PiEtAB) * psi;
						cP->perfomance = 0;
						Ind = cP->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cP);
						for (m = 0; m < dim_u; m++)
							u_i.v->v[m] = cP->getIJ(Ind, m);
						u_i = cP->getBorderPoint(Ind, u_i);
						//�������� v_i
						//extrVec = Transpose(PiEtAC) * psi;
						//cQ->perfomance = 0;
						//Ind = cQ->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cQ);
						//for (m = 0; m < dim_v; m++)
						//v_i.v->v[m] =   cQ->getIJ(Ind, m);
						//v_i=cQ->getBorderPoint(Ind,v_i);
						ps->currInd = n;  //���������� ����� � ������� ������������ �� ��������� �� �������� �������� �������������
						break;
					}
				}
			}
		   }
		   //��������� ���������� ��������� �� �������
		}   /**/

	   //	cout << (j) * tau << " : "<< t << " : " << x_i<< " : " << u_i;
		//������� ���������  x_i  �������  �����-�����
		x_i = rungeCutt(x_i, u_i, v_i);

		storeLocResults(u_i, v_i, x_i, vx_i, vu_i, vv_i, r_i);//��������� ������������� ���������� ��������

		j++;
	}
	if ((vx_i.size()-1)>1) {
			storeResults(trNum, vx_i, vu_i, vv_i);  //��������� ���������� ��������
		}else{
            for_each(vx_i.begin(), vx_i.end(), DeleteObj());
			for_each(vu_i.begin(), vu_i.end(), DeleteObj());
			for_each(vv_i.begin(), vv_i.end(), DeleteObj());
			return -1;
        }

	//--------------------------------------------------------------------
   //	cout<<k;
	return 0;
}

//------------------------------------------------------------------------
void PR_Task::plot(int trNum){      //���� ���� ������ ��� ��������� ��������
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
		EtA(mt.A.m()); // EtA(A.m()) - ��������� ��� t=0
	LDouble t, min, absmin=-1, min_x, extr;
	Vector PiEtAx(mt.PP.m()), psi(mt.PP.v->m), extrVec(mt.PP.v->m);
	Vector u_i(mt.cP->Dim), v_i(mt.cQ->Dim), x_i(mt.dim_x);
	bool isCollisionPossible=false,  //������� ����������� ��������� �� ��������� ������ �������� ���������� �������������
		  isInit, //������� ������, ����� ���������� ��� ���������� ����. ��������� ����������
		  isControlChange = false, //������� ������������ ������ ����������  �������������<->��������
		  isCtrlCanChange = false;//������� "�� ������ ������" ����� ����������
	long controlType = -1;//������� ���� ���������� -1- �������������, 0 - ��������, 1 - ��������
	long prevCt,prevCt2;
	//bool isCtrlChangeNextStep=false;
	Pursuer* ps;

	VecOfVec vx_i, vu_i, vv_i;
	Vector *r_i, *prev_u_i;

	VecOfLong vec_j, vec_k; // �������� �������� j � k ����� ����������
	VecOfLong vec_ct; // ��������  �������� controlType  ����� ����������
    vec_ct.push_back(controlType);
	alphType vec_t;//�������� t ����� ����������

	j = 0;
	k = mt.tr_s[trNum].NetList.size()-1;
	k_up=k; //"����������" ��������

	x_i = mt.tr_s[trNum].x0; // ��������� x_i  ��������� ���������
	r_i =  new Vector(x_i);
	r_i->detach();
	vx_i.push_back(r_i);//... � ��������� ��� ��� ����������.
	t = mt.tr_s[trNum].T; // �������� t = ��������� �������;

	c.perfomance = 0;

	while (k >= 1) {

		if (!isCollisionPossible) {     //���� ��� ����, ��� � ������ �� ��������,
			// ��� ���������� �������� ��������� ������� �������� ��� �����������
			k--;
			t -= mt.tau;
			EtA = Exponential(mt.A, t, mt.precision);  //�� ��� ���� ������
			PiEtA = mt.PP * EtA;
			PiEtAx = PiEtA * x_i;

			mt.pointGeomDiff(trNum, k, PiEtAx, c, x_Net);// ������� ������� ������� ����� PiEtAx
			//�������� psi
			Ind = c.GetExtrGlobal(opMin, 0, min_x);
		} else{
			//���� �������� ��������� ��������������� �������, � �� ���������� ��  �����
			//� ��������� ������ ���������� ��������� ������������� - ������������ �����
			// ��� ���� ������� ���������� ������, � ����� ������������� �� ����� � ������� ���������  �� x(t_i-tau)
			isInit = true;
			do{
				EtA = Exponential(mt.A, t, mt.precision);
				PiEtA = mt.PP * EtA;
				PiEtAx = PiEtA * x_i;
				//cout << k << endl;
				mt.pointGeomDiff(trNum, k, PiEtAx, c, x_Net);// ������� ������� ������� ����� PiEtAx

				Ind = c.GetExtrGlobal(opMin, 0, min_x);

				if (min_x>0) //���� ������ ����� ��������� - ��������  ����� ��������
					timeSign=-1;
				else       //���� ��������� ������� - ��������
					timeSign=1;
				if(isInit){
					prevInd = Ind;
					prevSign = timeSign;
					isInit = false;
				} /**/
				if(prevSign!=timeSign){ //����� ����� ��� �������� "�������" ��� "������" ���������� ��-�� ����. ���. ��������� �����  � ���� ������� �����
					Ind = prevInd;
					break;
				}
				if((static_cast<long>(k)+timeSign)==k_up){  //���� ����� �� ������� ������� ���������������� ���������������� ����������
					//����������� ���������������� �������� � ��������� ����� k_up
					mt.calcNextAltInt(trNum, t+timeSign * mt.tau);
					k_up = mt.tr_s[trNum].NetList.size()-1;
				}
				k += timeSign;
				t += timeSign * mt.tau;

			} while (true);/**/

		}
		//��������� psi  ��� ��������� �������������, ���������� � ���������
		for (m = 0; m < c.Dim; m++)
				psi.v->v[m] = c.getIJ(Ind, m);

		//�������� u_i - ��� ���� �� �� ���� �����������
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
		for (m = 0; m < mt.dim_v; m++) //� �������� v ���� �� ������������ � �����c��� ������ �������� �������������
				v_i.v->v[m] =   0.0;  //���� ����� �������� v1 - ����������� ������ + v2 - �� �����.
				//��������  ���������������

		//����� ����, ��� ���� ������ ��� ������������� ���� ��� ��� ����������� � ����������� ������
		//�������� ���������� �������������  �� v �������� �� ��������: ���� ���������� ������������
		//� ��������� ���������������, �� �������� v  �������������� ���������� ��������� ��������������
		tmpPNet = *mt.cP;
		tmpPNet *= PiEtAB;
		absmin = 0.0;
		isCollisionPossible = false;
		controlType  = -1;
		//vec_ct.push_back(controlType);

		for (l=0; l < mt.Pursuers/* *PursuerList*//*NList*/.size(); l++) {
		   ps = mt.Pursuers[l];
		   //��������� ������ ���������
		   //cout << eu_dist(*ps->center, x_i)<<endl;
		   if(eu_dist(*ps->center, x_i) <  ps->radarVisibility){
			fi = ps->funnel.size() - 1;
			if(ps->currInd >0 ){   //�.�. ���� ��� t-tau ��� ���������� �� ��������� ������ �������� ������������� � ����������, ��
				n = ps->currInd;    //���  ���� ���������, ������� �� ��� t
				ps->currInd = 0;
			}else                 //����� �������, ��� ��� �������� � ������ "�������"
				n=0;
			for (/*n = 0*/; n < ps->funnel.size(); n++) {
				px_Net = *ps->funnel[fi-n] + (-1)*x_Net;
				p_Net = px_Net + tmpPNet;   //�������� �������������� ����� ���������� �����
				//��������� ������������ �� ����� ����������� � ��������� �����. ��������������
				p_Net.perfomance = 0;
				Ind = p_Net.GetExtrGlobal(opMin, Ind, min_x);
				if(min_x>=0){  //���� �� ��������� ���� ���� ����������� � ����������� ������
					//�������� ���������� ������������� �� u �������� ���, ����� ������������ �� ���������
					if(!isCollisionPossible){
						isCollisionPossible = true;
						absmin=min_x;
					}
					if(min_x<=absmin){  // ������������� �� ���������� �� ������� �������� �������������
						absmin = min_x;
						px_Net.perfomance = 0;
						Ind = px_Net.GetExtrGlobal(opMin,  0, min);
						for (m = 0; m < c.Dim; m++)
						psi.v->v[m] = px_Net.getIJ(Ind, m);
						//�������� u_i - ���� ���������� � ������� �������� �� ������ ���� u_i+1 = -u_i c ���, ����� ������ ���� ���������� � ������� ������� ����
						extrVec = Transpose(PiEtAB) * psi;
						mt.cP->perfomance = 0;
						Ind = mt.cP->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *mt.cP);
						for (m = 0; m < mt.dim_u; m++)
							u_i.v->v[m] = mt.cP->getIJ(Ind, m);
						u_i = mt.cP->getBorderPoint(Ind, u_i);
						//�������� v_i
						//extrVec = Transpose(PiEtAC) * psi;
						//cQ->perfomance = 0;
						//Ind = cQ->GetExtrDirection(extrVec, scm, convCriteria, opMin, nZAware, Ind, NULL,  extr,  *cQ);
						//for (m = 0; m < dim_v; m++)
						//v_i.v->v[m] =   cQ->getIJ(Ind, m);
						//v_i=cQ->getBorderPoint(Ind,v_i);
						ps->currInd = n;  //���������� ����� � ������� ������������ �� ��������� �� �������� �������� �������������

						controlType  = 1;

						break;
					}
				}
			}
		   }
		   //��������� ���������� ��������� �� �������
		}   /**/


		prevCt = vec_ct[vec_ct.size()-1];

		if (
				(((prevCt == 0)||(prevCt == 1))&&(controlType  == -1))
				||
				(((prevCt == 0)||(prevCt == -1))&&(controlType  == 1))
		   )
			isControlChange = true;
		else
			isControlChange = false;

	   //isControlChange = (prevCt == controlType)?false:true;

		cout<<"From: " << (j) * mt.tau << " : "<< t << " : " << x_i<< " : " << u_i << " : " << vx_i.size()<< endl;

		 if(isControlChange && isCtrlCanChange){ //��� ������������ ������ ���������� (��� �����������)
			recusionStep(trNum, mt, v_i, vx_i, vu_i, vv_i, vec_j, vec_k, vec_t, j, k, t, u_i, x_i);
			cout <<"Back: "<< (j) * mt.tau << " : "<< t << " : " << x_i<< " : " << u_i << " : " << vx_i.size()<< endl;
			controlType = 0;
			//isCtrlChangeNextStep= false;
			isCtrlCanChange= false;
		  //  isControlChange = false;
		}  else
			isCtrlCanChange= true;
		 /**/

		//������� ���������  x_i  �������  �����-�����
		x_i = mt.rungeCutt(x_i, u_i, v_i);
		mt.storeLocResults(u_i, v_i, x_i, vx_i, vu_i, vv_i, r_i);//��������� ������������� ���������� ��������
		cout <<"after: "<< (j) * mt.tau << " : "<< t << " : " << x_i<< " : " << u_i << " : " << vec_t.size()<< endl;

		vec_ct.push_back(controlType);
		vec_j.push_back(j);
		vec_k.push_back(k);
		vec_t.push_back(t);
		j++;
	}

	mt.storeResults(trNum, vx_i, vu_i, vv_i); //��������� ���������� ��������

	//--------------------------------------------------------------------
	cout<<k;
}

// ������� ������� PiEtA*x0 - ������ ��������� ��� ��� ��� ������ ��������� ������������
TNetF &PR_Task::point_oporn(const Vector &PiEtAx0, TNetF &x0Net, unsigned long i) const {
    for (i = 0; i < x0Net.Count; i++)
          if (x0Net.f->v->v[i] != pinMark){
              x0Net.f->v->v[i] = scm(i, PiEtAx0, &x0Net,NULL);
          }
	return x0Net;
}

void PR_Task::pointGeomDiff(int trNum,  unsigned long k, const Vector &PiEtAx, TNetF &c, TNetF &x_Net) const {
    unsigned long i;
    for (i = 0; i < x_Net.Count; i++){
		if (x_Net.f->v->v[i] != pinMark){
			x_Net.f->v->v[i] = scm(i, PiEtAx, &x_Net,NULL);
			c.f->v->v[i] = tr_s[trNum].NetList[k]->f->v->v[i] - x_Net.f->v->v[i];
		}
	}
}

void PR_Task::storeResults(int trNum, VecOfVec &vx_i, VecOfVec &vu_i,  VecOfVec &vv_i) {

	unsigned long j,m, k;

	tr_s[trNum].T= (j) * tau;
    /* TODO -orum : ������� �������� ����������� �� ������� �������� � ������� */
    k= vx_i.size()-1;
    tr_s[trNum].x_i = new Matrix(k + 1, dim_x);
    // ������ �������� �� ����������  x_i
	tr_s[trNum].u_i = new Matrix(k, dim_u);
    // ������ �������� �� ����������  u_i
    tr_s[trNum].v_i = new Matrix(k, dim_v);
    // ������ �������� �� ����������  v_i

    for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[0][m] = vx_i[0]->v->v[m];
    for (j=0; j < k; j++) {
		//��������� ���������� ��������
		for (m = 0; m < dim_u; m++)
			tr_s[trNum].u_i->v->v[j][m] = vu_i[j]->v->v[m];
		for (m = 0; m < dim_v; m++)
			tr_s[trNum].v_i->v->v[j][m] = vv_i[j]->v->v[m];
		for (m = 0; m < dim_x; m++)
			tr_s[trNum].x_i->v->v[j + 1][m] = vx_i[j + 1]->v->v[m];
	  // cout<<*vx_i[j+1];
	}

	for_each(vx_i.begin(), vx_i.end(), DeleteObj());
	for_each(vu_i.begin(), vu_i.end(), DeleteObj());
	for_each(vv_i.begin(), vv_i.end(), DeleteObj());
   // return k;
}

void PR_Task::storeLocResults(const Vector &u_i, const Vector &v_i, const Vector &x_i, VecOfVec &vx_i, VecOfVec &vu_i,
							  VecOfVec &vv_i, Vector *r_i) const {//��������� ���������� ��������
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

void recusionStep(int trNum, PR_Task &mt, const Vector &v_i, VecOfVec &vx_i, VecOfVec &vu_i, VecOfVec &vv_i,
				  VecOfLong &vec_j, VecOfLong &vec_k,   alphType &vec_t, unsigned long &j,
				  unsigned long &k, LDouble &t, Vector &u_i, Vector &x_i) {
	//������� ���������  x_i  �������  �����-����� �.�. ����� ���� "��������"
	//���������� ��� ����������� ���������  ��� �������, ���� �� �����  �� ����������� ���������� ���������

	// ���� ����������� ���� � ���� ���� �����������, � ����� ���������� �� �������, �� ������ ����� �������, ��� ��������� ����� x_{i-2}, �������� ����� ������� x_i
	PR_Task mid_step(mt); //��� ����� ������ ����� ������

	vector<TNetF*>ms_NetList;
	mid_step.tr_s.clear();
	mid_step.PursuerList.clear(); //�.�. � ���� ������� �� �� ���� ���������� �� ���� - �������, ��� ��������������� ���,
	//� ����������� ���� ����� ������� ������ �������������� ������ ��������� �������� �� ������������

	Vector mid_x_0(*vx_i[j-1]);  //� �������� ��������� ����� �������� x(t_{i-1})  vx_i.size()-1

	x_i = mt.rungeCutt(x_i, u_i, v_i);
	Vector mid_x_t(x_i); //� �������� ������������� ���������  x(t_{i})
	LDouble t_0;

	Traectory tr;
	tr =  mt.tr_s[trNum];
	tr.x0 = mid_x_0;
	tr.NetList = ms_NetList;
	mid_step.tr_s.push_back(tr);

	int ii;
	string str = "";
	for (ii = 0; ii < mid_step.dim_m-1; ii++) {
				str = str + "p"+ldToStr(ii)+"*("+ldToStr(mid_x_t[ii])+")+";
			}
	str = str + "p"+ldToStr(ii)+"*("+ldToStr(mid_x_t[ii])+")";
	mid_step.setFuncToNetF(*(mid_step.cM), str);
	mid_step.cM->oporn(mid_step.t0, 1);

   //		cout <<"������: "<< mid_x_0 ;
   //		cout <<"�����: "<< mid_x_t ;

	int kk;  //� ������, ���� ����������� ������� ������ - ���������� ��� �� ������� � �����. ������� ����������� ��������
	do{
		mid_step.calcPursuerSets(trNum);
		t_0 = mid_step.TimeCalc_PR(trNum);
		if (t_0> 0) //mid_step.tau
			kk = mid_step.Control_PR_fullSets(trNum);
		if ((kk<0)||(t_0<0)) {
			mid_step.tau /=2;
			mid_step.precision /=2;
			mid_step.epsilon /=2;
		}
	}while ((kk < 0)||(t_0<0));

   //	cout << *mid_step.tr_s[trNum].u_i;
	u_i = mid_step.tr_s[trNum].u_i->GetRow(0);  //mid_step.tr_s[trNum].u_i->m()-1
  //	cout <<"u_i: "<<  u_i;

	x_i = mid_x_0;
	k  = vec_k[vec_k.size()-1];
	j  = vec_j[vec_j.size()-1];
	t  = vec_t[vec_t.size()-1];

	vec_k.pop_back();
	vec_j.pop_back();
	vec_t.pop_back();

	vx_i.pop_back();
	vu_i.pop_back();
	vv_i.pop_back();

}
