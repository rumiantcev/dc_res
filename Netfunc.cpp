// ---------------------------------------------------------------------------
#pragma hdrstop

#include "Netfunc.h"
// ---------------------------------------------------------------------------
#pragma package(smart_init)

// -----------------------------------------------------------------------------//
LDouble convCriteria(long num, LDouble sk, const TNetF *v) {
	return sk;
};

// -----------------------------------------------------------------------------//
LDouble convCriteria1(long num, LDouble sk, const TNetF *v) {
	return (v->f->v->v[num]) / sk;
};

// -----------------------------------------------------------------------------//
LDouble convCriteria2(long num, LDouble sk, const TNetF *v) {
	return v->f->v->v[num] * sk;
};

// ------------------------------------Create----------------------------------//
/* !inline */
void __fastcall TNetF::create(long Dim /* , long perf, long res */) {
	long i, j;
	LDouble acc;
	NumOfSur = Dim * 2;
	Vector norm(Dim);
	f = new Vector(Count);
	for (i = 0; i < Count; i++) {
		f->v->v[i] = 0;
		acc = 0;
		if (!isVirtual) {
			for (j = 0; j < Dim; j++)
				acc += v->v->v[i][j] * v->v->v[i][j];
			acc = sqrt(acc);
			for (j = 0; j < Dim; j++)
				v->v->v[i][j] /= acc;
		}
	}
	zeroPrecision = 1 / (Count * 100000.0);
	vars = NULL; // new double[Dim];
}

// ------------------------------constructor ----------------------------------//
__fastcall TNetF::TNetF(int dim, long perf, long res)
	: TNet(dim, perf, res, true) {
	initNetFDefault();
	create(dim /* , perf, res */);
	SetFunc("0");
};

// ------------------------------constructor ----------------------------------//
__fastcall TNetF::TNetF(int Dim, long perf, long res, string fstr)
	: TNet(Dim, perf, res, true) {
	/* TODO -orum -caddon : ������� ����������� � ������������ ������������ Virt/Non Virt net */
	initNetFDefault();
	create(Dim /* , perf, res */);
	SetFunc(fstr);
}

// ------------------------------constructor ----------------------------------//
__fastcall TNetF::TNetF(long mm, int nn, long perf, long res)
	: TNet(nn, perf, res, true) {
	initNetFDefault();
	f = new Vector(mm);
	create(nn /* , perf, res */);
	SetFunc("0");
};

// ------------------------copy constructor-----------------------------------//
__fastcall TNetF::TNetF(TNet& Net) : TNet(Net) {
	initNetFDefault();
	create(Dim /* , perfomance, Res */);
	SetFunc("0");

};

/* */// ------------------------copy constructor-----------------------------------//
__fastcall TNetF::TNetF(TNet& Net, const string& fstr) : TNet(Net) {
	initNetFDefault();
	create(Dim /* , perfomance, Res */);
	SetFunc(fstr);
};

// ------------------------copy constructor-----------------------------------//
__fastcall TNetF::TNetF(const TNetF& NetF) : TNet(NetF) {
	initNetFDefault();
	copyNetFFrom(NetF, false);
	SetFunc(NetF.fStr);
};

/* */// ------------------------copy constructor-----------------------------------//
__fastcall TNetF::TNetF(TNetF& Net, const string& fstr) : TNet(Net) {
	initNetFDefault();
	copyNetFFrom(Net, false);
	SetFunc(fstr);
};

// ------------------------------- destructor ---------------------------------//
void TNetF:: operator delete(void *p) {
	TNetF *ptr = (TNetF*)p;

	if (ptr->f == NULL)
		delete(void*) p;
	else
		p = NULL;
};

// ---------------------------- destructor-------------------------------------//
__fastcall TNetF::~TNetF() {

	// if (--f->v->linkCount==0)
	{
		delete f;
		f = NULL;
	}
#ifdef _WIN64
	TDllStdProcV1<TSIC_Data64*>sic_done(*dll, "sic_done");
#else
	TDllStdProcV1<TSIC_Data32*>sic_done(*dll, "sic_done");
#endif
	sic_done(&sic);

	TDllProcV0 sic_fretab(*dll, "sic_fretab");
	sic_fretab();

	if (dll != NULL) {
		delete dll;
		dll = NULL;
	}

	if (vars != NULL) {
		delete[]vars;
		vars = NULL;
	}
	// destroy();
};

// -------------------------------- << ----------------------------------------//
ostream& __fastcall operator << (ostream& out_data, TNetF& C) {
	long SurfaceDim = C.Dim - 1;
	long NumOfSur = C.Dim * 2;
	long NumOfPoints = C.Res;
	long i, j;

	if (!C.updated)
		C.update();
	for (i = 2; i <= SurfaceDim; i++)
		NumOfPoints *= C.Res;

	for (i = 0; i < C.Count; i++) {
		out_data << "[";
		for (j = 0; j < C.Dim; j++)
			out_data << C.getIJ(i, j) << ",";
		out_data << C.f->v->v[i];
		out_data << "]";
		if (i == NumOfPoints*NumOfSur - 1)
			out_data << ";";
		out_data << endl;
	}
	/* */
	return out_data;
};

// ----------------------------------- >> -------------------------------------//
istream& __fastcall operator >> (istream& in_data, TNetF& C) {
	// TNet v=C;
	char c = 0;
	double tempVal;
	long i, j;

	if ((!C.isVirtual) && (C.v == NULL))
		C.v = new Matrix(C.Count, C.Res);
	for (i = 0; i < C.Count;)
		for (j = 0; j < C.Dim + 1;) {
			in_data.get(c);
			switch (c) {
			case '[':
				break;
			case ']':
				break;
			case ',':
				break;
			case '\n':
				break;
			default: {
					in_data.putback(c);
					in_data >> tempVal;
					if (j == C.Dim)
						C.f->v->v[i] = tempVal;

					else if (!C.isVirtual)
						C.v->v->v[i][j] = tempVal;
					j++;
					if (j > C.Dim)
						i++;
				};
			};
		};
	in_data.get(c);
	while (c != ';' && c != EOF)
		in_data.get(c);
	in_data.get(c);
	return in_data;
};
// ----------------------------------- = --------------------------------------//
TNetF& __fastcall TNetF:: operator = (const TNetF & NetF) {
	if (this == &NetF)
		return*this;

	if (!isVirtual)
		if (--v->v->linkCount == 0) {
			delete v;
			v = NULL;
		}
	copyNetFrom(NetF);
	copyNetFFrom(NetF, true);
	return *this;
};

// ----------------------------------- + --------------------------------------//
TNetF& __fastcall TNetF:: operator += (const TNetF& B) {
	long i;
	TNetF* pB = (TNetF*)&B;
	// double coeff = B.upd/upd;
	// if ((v == B.v) && (!isVirtual)) {
	// upd += B.upd;
	// updated = false;
	// }
	// else {
	if (!updated)
		update();
	if (!B.updated)
		pB->update();
	if (virtDim != B.virtDim) {
		Count = B.Count;
		virtDim = B.virtDim;
		if (f != NULL)
			delete f;
		f = new Vector(Count);
	}
	for (i = 0; i < Count; i++)
		f->v->v[i] += pB->f->v->v[i];
	// }

	return *this;
}

// ----------------------------------- + --------------------------------------//
const TNetF __fastcall operator +(const TNetF& A, const TNetF& B) {
	return TNetF(A) += B;
};

// ----------------------------------- �������������� - ----------------------//
TNetF& __fastcall TNetF:: operator -= (const TNetF& B) {
	long i, j, m;
	Vector vec(Dim), st(Dim);
	TNetF st0Net(Dim, perfomance, Res);
	LDouble coeff = (LDouble)Dim / Count, val;
	bool extr_exist, *L;
	FindPath ynPath = perfomance > 1 ? yPath : nPath;

	long k = 1, jj;
	LDouble tmin = _extr_tmin_param, tmax = _extr_t0_param, t = tmax, p, a =
		double(_lrand() % Count) / Count, c = _extr_e_val;

	// ���� ������������- ������� �� ���������

	TNetF* pB = (TNetF*)&B;
	if (!updated)
		update();
	if (!B.updated)
		pB->update();
	if (virtDim != B.virtDim) {
		Count = B.Count;
		virtDim = B.virtDim;
		if (f != NULL)
			delete f;
		f = new Vector(Count);
	}
	for (i = 0; i < Count; i++)
		f->v->v[i] -= pB->f->v->v[i];
	// cout<<*this;
	// ������ ������� ������� ������ ��������
	st = 0;
	for (i = 0; i < Count; /* i=i+2 */ i++) {
		for (m = 0; m < Dim; m++)
			vec.v->v[m] = st0Net.getIJ(i, m);
		st += (f->v->v[i]*vec);
	}
	st *= coeff;

	// �������� ��������� "�� ����� ��������" ��� ����� 0 ��� � ������
	for (i = 0; i < st0Net.Count; i++) {
		st0Net.f->v->v[i] = scm(i, st, &st0Net,0);
		f->v->v[i] -= st0Net.f->v->v[i];
	}
	// ������������ lambda_i
	L = new bool[Count];
	// alpha.clear();
	alpha.reserve(Count);

	makeAlpha(alpha, L, st0Net);
	// ���������� ����������� � ����� ��������� �������
	for (i = 0; i < st0Net.Count; i++) {
		if (L[i]) {
			extr_exist = false;
			for (m = 0; m < Dim; m++)
				vec.v->v[m] = st0Net.getIJ(i, m);

			 j = GetExtrDirection(vec, scm1, convCriteria, opMax, nZAware, i, &alpha,  f->v->v[i],  st0Net);
		  /*	if (perfomance == optNone) {
				// ---------------������� �������
				for (j = 0; j < st0Net.Count; j++) {
					val = alpha[j] * scm(j, vec, &st0Net,NULL);
					if (!extr_exist) {
						extr_exist = true;
						f->v->v[i] = val;
					}
					else if (f->v->v[i] < val)
						f->v->v[i] = val;
				}
			}
			if (perfomance == optAnnealing) {
				// ����� ������� �������� ������
				j = _lrand() % (st0Net.Count);
				while (t > tmin) {
					val = alpha[j] * scm(j, vec, &st0Net,NULL);
					if (!extr_exist) {
						extr_exist = true;
						f->v->v[i] = val;
					}
					else if (f->v->v[i] < val) {
						f->v->v[i] = val;
						t = t * c;
						// t = tmax / k; // ���������������� ���� ����� �� ����
						// k++; // ���������������� ���� ����� �� ����
					}
					else {
						p = 1 / (1 + exp(-abs(alpha[i] - val) / t));
						// p= t/(M_PI*(pow(abs(val-alpha[i])/st0Net.Count,2)+pow(t,2)));  // ���������������� ���� ����� �� ����
						a = double(_lrand() % st0Net.Count) / st0Net.Count;
						if (a > p) {
							f->v->v[i] = val;
							t = t * c;
							// t = tmax / k; // ���������������� ���� ����� �� ����
							// k++; // ���������������� ���� ����� �� ����
						}
					}
					a = double(_lrand() % st0Net.Count) / st0Net.Count;
					jj = signof(a - 0.5) * t *
						(pow((1 + 1 / t), abs(2 * a - 1)) - 1) * st0Net.Count;
					// jj = t * tan(M_PI*(a-0.5))*st0Nett.Count;  // ���������������� ���� ����� �� ����
					j = abs(jj + j) % st0Net.Count;
				}
			}/**/
		}
		f->v->v[i] += st0Net.f->v->v[i];
	}
	/* */
	alpha.clear();
	delete[]L;
	return *this;
}

// ----------------------------------- �������������� - ----------------------//
const TNetF __fastcall operator -(const TNetF& A, const TNetF& B) {
	return TNetF(A) -= B;
};

// ------------------------------------Update---------------------------------//
void __fastcall TNetF::update() {
	TNetF *f1 = NULL;
	Vector *res;
	detach();
	if (!isVirtual) {
		long i, j;
		if (!updated) { // ��������� �� ���� ������� ��� ���������
			updated = true;
			if (!umx) {
				// Vector res(u_mx->v->m);
				if (u_mx->v->m != Dim)
				{ // � ������ ���� ���������� ����������� - ������ ����� �����
					res = new Vector(u_mx->v->m);
					f1 = new TNetF(u_mx->v->m, perfomance, Res);
				}
				else { // ����� ��������� � ������������
					res = new Vector(Dim);
					f1 = this;
				}
			}
			for (i = 0; i < Count; i++) {
				if (!umx) {
					for (j = 0; j < Dim; j++)
						res->v->v[j] = getIJ(i, j);
					(*res) = (*u_mx) * (*res);
					f1->f->v->v[i] = oporn(*res, t, upd);
					for (j = 0; j < Dim; j++)
						f1->v->v->v[i][j] = res->v->v[j];
				}
				else if (upd != 1.0)
					f->v->v[i] *= upd;

			}

			if ((!umx) && (u_mx->v->m != Dim)) {
				*this = *f1;
				delete f1;
			}
			delete res;
		}

	}
	else
		dynUpdate();

	umx = true;
	upd = 1.0;
}

// ------------------------------------Detach----------------------------------//
void __fastcall TNetF::detach() {
	TNet::detach();
	if (f->v->linkCount > 1) {
		long i;
		Vector *vv = NULL;
		vv = Vector::copy(f, vv);
		delete f;
		f = new Vector(vv->v->size);
		for (i = 0; i < vv->v->size; i++)
			f->v->v[i] = vv->v->v[i];
		delete vv;
	}; /* */
}

// ----------------------------------- * --------------------------------------//
TNetF& __fastcall TNetF:: operator *= (const Matrix& A) {
	bool isDimChanged = (Dim != A.n()) ? true : false;

	// if (!uvec)
	// update();
	if (u_mx != NULL)
		* u_mx *= Transpose(A);
	else
		u_mx = new Matrix(Transpose(A));
	Dim = A.n();
	if (isDimChanged) {
		SetFunc(fStr);
		buildPowerVectors(Dim);
	}

	updated = false;
	umx = false;
	return *this;
}

// ----------------------------------- * --------------------------------------//
const TNetF __fastcall operator*(const Matrix &A, const TNetF& B) {
	return TNetF(B) *= A;
};

// ----------------------------------- * --------------------------------------//
TNetF& __fastcall TNetF:: operator *= (const LDouble &a) {
	upd *= a;
	updated = false;
	return *this;
}

// ----------------------------------- + --------------------------------------//
/* TNetF& __fastcall TNetF:: operator += (const double &a) {
 updSum += a;
 updated = false;
 return *this;
 }
/* */
// ----------------------------------- * --------------------------------------//
const TNetF __fastcall operator *(const double a, const TNetF& B) {
	return TNetF(B) *= a;
};

// ----------------------------------- + --------------------------------------//
/* const TNetF __fastcall operator +(const double a, const TNetF& B) {
 return TNetF(B) += a;
 };
/* */
// ------------------------------------AddVariables---------------------------//
void __fastcall TNetF::AddVariables() {
	int i;
	string var;
	ostringstream ss;

	// func->ExtVariables->FreeItems();
	if (vars != NULL) {
		delete[]vars;
		vars = NULL;
	}
	if (vars == NULL)
		vars = new double[Dim];
	/* if(cache!=NULL)
	 {
	 delete cache;
	 cache = new Vector(Dim);
	 }
	/* */
#ifdef _WIN64
	TDllStdProcV3<TSIC_Data64*, char*, double*>sic_avarf(*dll, "sic_avarf");
#else
	TDllStdProcV3<TSIC_Data32*, char*, double*>sic_avarf(*dll, "sic_avarf");
#endif
	for (i = 0; i < Dim; i++) {
		ss << i;
		var = "p";
		var.append(ss.str());
		ss.str("");
		ss.clear();
		sic_avarf(&sic, (char*)var.c_str(), &vars[i]);
	}
	sic_avarf(&sic, "t", &t);
#ifdef _WIN64
	TDllStdProcV1<TSIC_Data64*>sic_patab(*dll, "sic_patab");
#else
	TDllStdProcV1<TSIC_Data32*>sic_patab(*dll, "sic_patab");
#endif
	sic_patab(&sic);
}

// ------------------------------------SetVariables---------------------------//
/* void __fastcall TNetF::SetVariables(const DynamicArray<double> &vv) {
 int i;
 for (i = 0; i < vv.Length - 1; i++)
 vars[i] = vv[i];
 t = vv[vv.Length - 1];
 }
/* */
// ---------------------------------------------------------------------------
LDouble __fastcall TNetF::oporn(const Vector &x, LDouble t, int sign) {
	LDouble result;
	DWORD err;

	long i;
	for (i = 0; i < Dim; i++)
		vars[i] = sign * x.v->v[i];

	this->t = t;
	// result = func->Evaluate();//   ->Value;
#ifdef _WIN64
	TDllStdProc2<double, TSIC_Data64*, DWORD*>sic_exec(*dll, "sic_exec");
#else
	TDllStdProc2<double, TSIC_Data32*, DWORD*>sic_exec(*dll, "sic_exec");
#endif
	result = sic_exec(&sic, &err);
	return result;
};

// ---------------------------------------------------------------------------
LDouble __fastcall TNetF::oporn(const Vector &x, LDouble t,
	const string &funcStr, int sign) {
	SetFunc(funcStr);
	return oporn(x, t, sign);
};

// ---------------------------------------------------------------------------
/* LDouble __fastcall TNetF::oporn(const Vector &x,  LDouble t,
 const TRapidEvaluator &func,int sign)
 {
 SetFunc(string(func.Formula.c_str()));
 return oporn(x,t,sign);
 };
/* */// ---------------------------------------------------------------------------//
void __fastcall TNetF::oporn(LDouble t, int sign) {
	long i, j;
	Vector vv(Dim);

	if ((f != NULL) && (f->v->size != Count)) {
		delete f;
		f = NULL;
	}

	if (f == NULL)
		f = new Vector(Count);

	for (i = 0; i < Count; i++) {
		for (j = 0; j < Dim; j++)
			vv[j] = getIJ(i, j);
		// cout<<vv[0]<<endl;
		f->v->v[i] = oporn(vv, t, sign);
	};
};

// ---------------------------------------------------------------------------//
TNetF __fastcall TNetF::opornx(LDouble t, int sign) {
	oporn(t, sign);
	return *this;
};

// ---------------
LDouble __fastcall TNetF::oporn(const Vector &x, LDouble t, int sign,
	const Matrix &A) {
	umx = true;
	Vector xx = A * x;
	return oporn(xx, t, sign);
};

// ---------------------------------------------------------------------------//
void __fastcall TNetF::oporn(LDouble t, int sign, const Matrix &A) {
	long i, j;
	Vector vv(Dim);

	for (i = 0; i < Count; i++) {
		for (j = 0; j < Dim; j++)
			vv.v[j] = v->v->v[i][j];
		f->v->v[i] = oporn(vv, t, sign, A);
	};
};
// ---------------------------------------------------------------------------//

void /* !inline */ __fastcall TNetF::SetFunc(const string& fstr) {
	detach();
	if (fstr.length() > 0) {
#ifdef _WIN64
		TDllStdProcV1<TSIC_Data64*>sic_done(*dll, "sic_done");
#else
		TDllStdProcV1<TSIC_Data32*>sic_done(*dll, "sic_done");
#endif
		sic_done(&sic);

		TDllProcV0 sic_fretab(*dll, "sic_fretab");
		sic_fretab();

		TDllProcV0 sic_cretab(*dll, "sic_cretab");
		sic_cretab();
#ifdef _WIN64
		TDllStdProcV1<TSIC_Data64*>sic_init(*dll, "sic_init");
#else
		TDllStdProcV1<TSIC_Data32*>sic_init(*dll, "sic_init");
#endif
		sic_init(&sic);

		// if(func!=NULL)
		// {
		// delete func;
		// func=NULL;
		// }

		// if(func==NULL)
		// {
		// func = new TRapidEvaluator(NULL);

		AddVariables();
		// }
		fStr = fstr;
#ifdef _WIN64
		DWORD sop = SIC_OPT_DEFAULT_X64;
#else
		DWORD sop = SIC_OPT_DEFAULT_X32;
#endif
		/* */
#ifdef _WIN64
		TDllStdProcV3<TSIC_Data64*, const char*, DWORD>sic_compile(*dll,
			"sic_compile");
#else
		TDllStdProcV3<TSIC_Data32*, const char*, DWORD>sic_compile(*dll,
			"sic_compile");
#endif
		sic_compile(&sic, fStr.c_str(), sop);
	}
}



// -------------------------- findExtrDirectionSlow ----------------------------------//
long __fastcall TNetF::findExtrSlowDirection(const Vector& vec, scM scmul,
	cCrit crit, ZeroAware isZeroAware, OpType extrOper, long index, alphType* coeff, LDouble &extr,  TNetF& net) {
	// ����� ��������� ���������� ������������ scmul �� ����� ������� ��������� � ����������� ��������� �������

	long i, j = 0;
	LDouble val, sc;
	bool isExtrExist = false, isGrZero, isExtr;

	for (i = 0; i < Count; i++) {
		sc = scmul(i, vec, &net,coeff);
		// j = selectExtrX(vec, scmul, crit, i, j, extr, isMax, isZeroAware,	isExtrExist, index, net);   // �� �������
		// ������� ��������� ������������ �������� �������  ����� (�urrent)  � ��������� (vec)
		isGrZero = sc > zeroPrecision;
		if ((isZeroAware == nZAware) || ((isZeroAware == ZAware) && isGrZero)){
			val = crit(i, sc, this);
			if (!isExtrExist){
				// ���� ����������  ��� �� ���� - ���� ������ �������� � ���������� ��� ���������
				isExtrExist = true;
				extr = val;
				j = i;
			}
			else{
				extrOper == opMax? isExtr = (val > extr) :  isExtr = (val < extr);
				if (isExtr)	{
					extr = val;
                   j=i;
				}
			}
		}

	}
	return j;
} /* */



// -------------------------- findMaxDirectionAnnealing ----------------------------------//
long  __fastcall TNetF::findExtrAnnealingDirection(const Vector& vec, scM scmul,
	cCrit crit, ZeroAware isZeroAware,OpType extrOper, long index, alphType* coeff, LDouble &extr, TNetF& net) {

	long  j, jj;
	LDouble val,sc;
	LDouble tmin = _extr_tmin_param, tmax = _extr_t0_param, t = tmax, p, a, c =
		_extr_e_val;

	bool isGrZero, isExtrExist = false, isExtr;

	//extr = f->v->v[j];
	//val = extr;
	clock_t b0;
	double e0, e1, e2, e3, e4, e5;



	j = _lrand() % (net.Count);
	while (t > tmin) {
	   //	b0 = clock();
		sc = scmul(j, vec, &net, coeff);
	 //	e0 = clock()-b0;
		isGrZero = sc > zeroPrecision;
		if ((isZeroAware == nZAware) || ((isZeroAware == ZAware) && isGrZero)) {
		 //	b0 = clock();
			val = crit(j, sc, this);
		 //	e1 = clock()-b0;
			if (!isExtrExist) {
				// ���� ����������  ��� �� ���� - ���� ������ �������� � ���������� ��� ���������
				isExtrExist = true;
				extr = val;
			}
			else{
				extrOper == opMax? isExtr = (val > extr) :  isExtr = (val < extr);
				if (isExtr) {
					extr = val;
					t = t * c;
				}
				else {
					p = 1 / (1 + exp(-fabs(extr-val) / t));
					a = double(_lrand() % net.Count) / net.Count;
					if (a > p) {
						extr = val;
						t = t * c;
					}
				}

			}
		}
		a = double(_lrand() % net.Count) / net.Count;
		jj = signof(a - 0.5) * t * (pow((1 + 1 / t), fabs(2 * a - 1)) - 1)
				* net.Count;
		j = abs(jj + j) % net.Count;
		//if ((e0>1) || e1>1)
        //    cout << "in call: e0- " << e0<<" e1-"<<e1<<endl;


	}



	return j;
	/**/

  /*	j = _lrand() % (net.Count);
	while (t > tmin) {
		sc = scm(j, vec, &net,NULL);
		if (sc > zeroPrecision) {
			val = f->v->v[j] / sc;//convCriteria1(j,sk,&net);
			if (!isExtrExist) {
				isExtrExist = true;
				extr = val;
			}
			else if (extr >= val) {
				extr = val;
				t = t * c;
			}
			else {
				p = 1 / (1 + exp(-abs(val - extr) / t));
				a = double(_lrand() % net.Count) / net.Count;
				if (a > p) {
					extr = val;
					t = t * c;
				}
			}
		}
		a = double(_lrand() % net.Count) / net.Count;
		jj = signof(a - 0.5) * t * (pow((1 + 1 / t), abs(2 * a - 1)) -1)
			* net.Count;
		j = abs(jj + j) % net.Count;
	}
	/**/
	return j;
} /* */

// -------------------------- findExtrSlowGlobal ----------------------------------//
long __fastcall TNetF::findExtrSlowGlobal(OpType extrOper, LDouble &extr) {
	// ����� �������� ������� ������� �� ����� ������� ���������
	long i, j = 0;
	LDouble val;
	bool isExtr;

	extr = f->v->v[0];
	val = extr;
	for (i = 1; i < Count; i++) {
		val = f->v->v[i];
		extrOper == opMax? isExtr = (val > extr) :  isExtr = (val < extr);
		if (isExtr) {
			extr = val;
			j = i;
		}
	}
	return j;
} /* */

// -------------------------- findExtrAnnealingGlobal ----------------------------------//
long  __fastcall TNetF::findExtrAnnealingGlobal(OpType extrOper, LDouble &extr) {
	// ����� ��������� ������� ������� �� ����� ������� ������
	long  j=_lrand() % (Count), jj;
	LDouble val;
	LDouble tmin = _extr_tmin_param, tmax = _extr_t0_param, t = tmax, p, a, c =
		_extr_e_val;
	bool isExtr;

	extr = f->v->v[j];
	val = extr;

	j = _lrand() % (Count);
	while (t > tmin) {
		val = f->v->v[j];
		extrOper == opMax? isExtr = (val > extr) :  isExtr = (val < extr);
		if (isExtr) {
			extr = val;
			t = t * c;
		}
		else {
			p = 1 / (1 + exp(-fabs(val - extr) / t));
			a = double(_lrand() % Count) / Count;
			if (a > p) {
				extr = val;
				t = t * c;
			}
		}
		a = double(_lrand() % Count) / Count;
		jj = signof(a - 0.5) * t * (pow((1 + 1 / t), fabs(2 * a - 1)) - 1)
			* Count;
		j = abs(jj + j) % Count;
	}

	return j;
} /* */

// -------------------------- selectExtrX ------------------------------------//
long /* !inline */ TNetF::selectExtrX(const Vector& vec, scM scmul, cCrit crit,
	long current, long& result, LDouble &extr, OpType isMax,
	ZeroAware isZeroAware, bool &isExtrExist, long index, TNetF& net) {
	LDouble sc, val = extr;
	bool isProxy, isGrZero;

	sc = (*scmul)(current, vec, &net,0);
	// ������� ��������� ������������ �������� �������  ����� (�urrent)  � ��������� (vec)
	isGrZero = sc > zeroPrecision;
	if ((isZeroAware == nZAware) || ((isZeroAware == ZAware) && isGrZero)) {
		val = (*crit)(current, sc, this);
		if (!isExtrExist)
			// ���� ����������  ��� �� ���� - ���� ������ �������� � ���������� ��� ���������
		{
			isExtrExist = true;
			extr = val;
			result = current;
		}
		else {
			isMax == opMax ? isProxy = val >= extr : isProxy = val <= extr;
			// ���� ���� ��������, �� ��������� �� �� ��� ����� �������� >= ������������, ���� ������� -<= ������������
			if (isProxy) {
				extr = val;
				result = current;
			}
		}
		return result;
	}
	// extr = val;
	result = current;
	return current;
	// return  result;
}

// -------------------------- findExtrFastX ----------------------------------//
long __fastcall TNetF::findExtrFastXDirection(const Vector& vec, scM scmul,
	cCrit crit, OpType isMax, ZeroAware isZeroAware, long index, alphType* coeff, LDouble &extr,
	TNetF& net){
	// ����� ���������� ���������� ������������ scmul �� ����� � ����������� ��������� ������� ��������� ���, ��� ������� ��������� �������� ����� ��������

	long i, j = index, je, k1, k2;
	// LDouble       extr;
	bool isExtrExist = false, isNotExtrFound = true, borderChanged;
	// , maxOper;
	seekType seekPath;

	seekPath.clear();

	while (isNotExtrFound) {
		/* DONE -orum -crepair : ���������� �������� ������ �� ����� */
		if (!seekPath[j]) {
			// prevInd=j;
			je = j;
			seekPath[j] = true;
			// }
			for (i = 0; i < Dim; i++) {
				try {
					k1 = shift(je, i, 1, borderChanged);
				}
				catch (exInvalidMoveDirection) { // do nothing;
				}
				if (k1 != je)
					j = selectExtrX(vec, scmul, crit, k1, j, extr, isMax,
					isZeroAware, isExtrExist, index, net);
				try {
					k2 = shift(je, i, -1, borderChanged);
				}
				catch (exInvalidMoveDirection) { // do nothing;
				}
				if (k2 != je)
					j = selectExtrX(vec, scmul, crit, k2, j, extr, isMax,
					isZeroAware, isExtrExist, index, net);
			}
			if (j < 0)
				j = shift(je, 0, 1, borderChanged);

		}
		else {
			// if(prevInd==j)
			isNotExtrFound = false;
			seekPath.clear();
		}

	}
	return j; /* */
}

// -------------------------- findExtrFastXGlobal ----------------------------------//
long __fastcall TNetF::findExtrFastXGlobal(OpType isMax, long index,
	LDouble &extr)
	// ����� ���������� ������� ������� ��������� ���, ��� ������� ��������� �������� ����� ��������
{
	long i, j = index, je, k1, k2;
	LDouble val;
	bool isExtrExist = false, isNotExtrFound = true, borderChanged, isProxy;
	seekType seekPath;

	seekPath.clear();

	while (isNotExtrFound) {
		/* DONE -orum -crepair : ���������� �������� ������ �� ����� */
		if (!seekPath[j]) {
			// prevInd=j;
			je = j;
			seekPath[j] = true;
			// }
			for (i = 0; i < Dim; i++) {
				try {
					k1 = shift(je, i, 1, borderChanged);
				}
				catch (exInvalidMoveDirection) { // do nothing;
				}
				if (k1 != je) {
					val = f->v->v[k1];
					if (!isExtrExist) {
						isExtrExist = true;
						extr = val;
					}
					else {
						isMax == opMax ? isProxy = val >= extr : isProxy =
							val <= extr;
						if (isProxy) {
							extr = val;
							je = k1;
							j = je;
						}
					}
				}
				try {
					k2 = shift(je, i, -1, borderChanged);
				}
				catch (exInvalidMoveDirection) { // do nothing;
				}
				if (k2 != je)
					if (k2 != je) {
						val = f->v->v[k2];
						if (!isExtrExist) {
							isExtrExist = true;
							extr = val;
						}
						else {
							isMax == opMax ? isProxy = val >= extr : isProxy =
								val <= extr;
							if (isProxy) {
								extr = val;
								je = k2;
								j = je;
							}
						}
					}
			}
			if (j < 0)
				j = shift(je, 0, 1, borderChanged);

		}
		else {
			// if(prevInd==j)
			isNotExtrFound = false;
			seekPath.clear();
		}

	}
	return j; /* */
}

// ----------------------------- GetExtr � ����������� ��������� �������------//
long __fastcall TNetF::GetExtrDirection(const Vector& vec, scM scmul,
	cCrit crit, OpType extrOper, ZeroAware isZeroAware, long index, alphType* coeff, LDouble& extr,
	TNetF& net) {

	if (perfomance == optNone)
		return findExtrSlowDirection(vec, scmul, crit, isZeroAware, extrOper, index, coeff,
			extr, net);

	if (perfomance == optAnnealing)
			return findExtrAnnealingDirection(vec, scmul, crit, isZeroAware,  extrOper,
			index, coeff, extr, net);

	if (perfomance == optGradient)
		return findExtrFastXDirection(vec, scmul, crit, extrOper, isZeroAware,
			index, coeff, extr, net);
		/* */ // Fast Evaluator

	return -1;
};

// ----------------------------- GetExtrGlobal---------------------------------------//
long __fastcall TNetF::GetExtrGlobal(OpType extrOper, long index, LDouble& extr) {
	if (perfomance == optNone)
			return findExtrSlowGlobal(extrOper,extr);

	if (perfomance == optAnnealing)
			return findExtrAnnealingGlobal(extrOper,extr);

	if (perfomance == optGradient) {
		// if (isMax == opMax)
		// return findMaxSlowGlobal(extr);
		// else
		return findExtrFastXGlobal(extrOper, index, extr); /* */ // Fast Evaluator
	}
	// else {

	// };
	return -1;
};

// -----------------------------------------------------------------------------//
LDouble  __fastcall scm(long num, const Vector &vec, TNetF *v, alphType* coeff) {
	LDouble sc = 0;
	int m;
	for (m = 0; m < vec.v->size; m++)
		sc += v->getIJ(num, m) * vec[m];

	return sc;
}

// -----------------------------------------------------------------------------//
LDouble __fastcall scm1(long num, const Vector &vec, TNetF *v, alphType* coeff) {
	LDouble sc = 0;
	int m;
	for (m = 0; m < vec.v->size; m++)
		sc += v->getIJ(num, m) * vec[m];
	return  (*coeff)[num]*sc;
}
// -----------------------------------------------------------------------------//
LDouble __fastcall scm2(long num, const Vector &vec, TNetF *v, alphType* coeff) {
	LDouble sc = 0;
	int m;
	for (m = 0; m < vec.v->size; m++)
		sc += (-v-> /* TNet:: */ getIJ(num, m)) * vec[m];
	return v->f->v->v[num] + sc;
}

// -----------------------------------------------------------------------------//
/* !inline */
/* DONE : ������������ ���������������� ������ */
void __fastcall TNetF::makeAlpha(alphType& alpha, bool* L, TNetF &net) {
	long i, j, m;
	LDouble extr, sk; // , *chk_alpha, val;
	bool extr_exist;
	Vector vec(Dim);

	long k = 1, jj;
	LDouble tmin = _extr_tmin_param, tmax = _extr_t0_param, t = tmax, p, a =
		LDouble(_lrand() % net.Count) / net.Count, c = _extr_e_val;

	// chk_alpha = new LDouble[Count];
	// cout<<alpha[i]<<" : "<<chk_alpha[i]<<endl;
	alpha.reserve(Count);
	for (i = 0; i < Count; i++) {
		if (L != NULL)
			L[i] = false;
		for (m = 0; m < Dim; m++)
			vec[m] = net.getIJ(i, m);
		// ------------GetMin by criteria--------------------------------------

		alpha[i] = 0.0;
	  /*	j = GetExtrDirection(vec, scm, convCriteria1, opMin, ZAware, i, NULL,  extr,  net);
		 alpha[i] = extr;
		 if (L != NULL)
			L[i] = true;
		/**/

		extr_exist = false;
		if (perfomance == optNone) {
			// ����� ���������������� ���������
			for (j = 0; j < net.Count; j++) {
				sk = scm(j, vec, &net,NULL);
				if (sk > 0.0) {
					extr = f->v->v[j] / sk;
					if (!extr_exist) {
						extr_exist = true;
						alpha[i] = extr;
						if (L != NULL)
							L[i] = true;
					}
					else if (alpha[i] > extr)
						alpha[i] = extr;
				}
			}
		}
		if (perfomance == optAnnealing) {
			// ����� ������� �������� ������
  //			clock_t before;
  //			double elapsed;
  //			before = clock();

			j = _lrand() % (net.Count);
			while (t > tmin) {
				sk = scm(j, vec, &net,NULL);
				if (sk > zeroPrecision) {
					extr = f->v->v[j] / sk;//convCriteria1(j,sk,&net);
					if (!extr_exist) {
						extr_exist = true;
						alpha[i] = extr;
						if (L != NULL)
							L[i] = true;
					}
					else if (alpha[i] >= extr) {
						alpha[i] = extr;
						t = t * c;
						// t = tmax / k; // ���������������� ���� ����� �� ����
						//k++;  // ���������������� ���� ����� �� ����
					}
					else {
						p = 1 / (1 + exp(-fabs(alpha[i] - extr) / t));
						// p= t/(M_PI*(pow(fabs(extr-alpha[i])/net.Count,2)+pow(t,2)));  // ���������������� ���� ����� �� ����
						a = double(_lrand() % net.Count) / net.Count;
						if (a > p) {
							alpha[i] = extr;
							t = t * c;
							// t = tmax / k; // ���������������� ���� ����� �� ����
							//k++;    // ���������������� ���� ����� �� ����
						}
					}
				}
				a = double(_lrand() % net.Count) / net.Count;
				jj = signof(a - 0.5) * t * (pow((1 + 1 / t), abs(2 * a - 1)) -
					1) * net.Count;
				// jj = t * tan(M_PI*(a-0.5))*net.Count;  // ���������������� ���� ����� �� ����
				j = abs(jj + j) % net.Count;
			}

//			elapsed = clock()-before;
//			cout << "in code:" << elapsed<<endl;
		}
		/**/

	}
}

// -----------------------------------------------------------------------------//

void __fastcall TNetF::Conv(bool *L) {
	long i, j, m, k;
	LDouble coeff = (LDouble)Dim / Count;
	TNetF st0Net(Dim, perfomance, NumOfPoints);
	Vector vec(Dim), st(Dim);
	// char dc=FormatSettings.DecimalSeparator; //��������� �� ����, ��� ��� �������������� ����� � ������ ������ '.' ����� ',' ��� �����������
	// FormatSettings.DecimalSeparator = '.';
	bool extr_exist;
	LDouble extr;
	pathType path;
	pathType::iterator pitr;
	alphType::iterator where;
	FindPath ynPath = perfomance > 1 ? yPath : nPath;

	long jj;
	LDouble tmin = _extr_tmin_param, tmax = _extr_t0_param, t = tmax, p, a =
		double(_lrand() % st0Net.Count) / st0Net.Count, c = _extr_e_val;
	k = 1;
	// ������ ������� ������� ������ ��������
	st = 0;
	for (i = 0; i < Count; i++) {
		for (m = 0; m < Dim; m++)
			vec.v->v[m] = st0Net.getIJ(i, m);
		st += (f->v->v[i]*vec);
	}
	st *= coeff;

	// �������� ��������� "�� ����� ��������" ��� ����� 0 ��� � ������
	for (i = 0; i < st0Net.Count; i++) {
		st0Net.f->v->v[i] = scm(i, st, &st0Net,0);
		f->v->v[i] -= st0Net.f->v->v[i];
	}

	// ������������ lambda_i
	makeAlpha(alpha, L, st0Net);


	// ���������� ���������� � �������� �����
	for (i = 0; i < st0Net.Count; i++) {
		if (L[i]) {
			extr_exist = false;
			for (m = 0; m < Dim; m++)
				vec.v->v[m] = st0Net.getIJ(i, m);

			 j = GetExtrDirection(vec, scm1, convCriteria, opMax, nZAware, i, &alpha,  f->v->v[i],  st0Net);
			/*
			if (perfomance == optNone) {
				// ������� �������
				for (j = 0; j < st0Net.Count; j++) {
					extr = alpha[j] * scm(j, vec, &st0Net,0);
					if (!extr_exist) {
						extr_exist = true;
						f->v->v[i] = extr;
					}
					else if (f->v->v[i] < extr)
						f->v->v[i] = extr;
				}
			}

			if (perfomance == optAnnealing) {
				// ����� ������� �������� ������
				j = _lrand() % (st0Net.Count);
				while (t > tmin) {
					extr = alpha[j] * scm(j, vec, &st0Net,0);
					if (!extr_exist) {
						extr_exist = true;
						f->v->v[i] = extr;
					}
					else if (f->v->v[i] <= extr) {
						f->v->v[i] = extr;
						t = t * c;
						// t = tmax / k; // ���������������� ���� ����� �� ����
						// k++; // ���������������� ���� ����� �� ����
					}
					else {
						p = 1 / (1 + exp(-abs(alpha[i] - extr) / t));
						// p= t/(M_PI*(pow(abs(extr-alpha[i])/st0Net.Count,2)+pow(t,2)));  // ���������������� ���� ����� �� ����
						a = double(_lrand() % st0Net.Count) / st0Net.Count;
						if (a > p) {
							f->v->v[i] = extr;
							t = t * c;
							// t = tmax / k; // ���������������� ���� ����� �� ����
							// k++; // ���������������� ���� ����� �� ����
						}
					}
					a = double(_lrand() % st0Net.Count) / st0Net.Count;
					jj = signof(a - 0.5) * t *
						(pow((1 + 1 / t), abs(2 * a - 1)) - 1) * st0Net.Count;
					// jj = t * tan(M_PI*(a-0.5))*st0Nett.Count;  // ���������������� ���� ����� �� ����
					j = abs(jj + j) % st0Net.Count;
				}
			}
			/**/
		}
		f->v->v[i] += st0Net.f->v->v[i];
		// ����� ��������� �������� �������� �� �������� �����
	}

	// FormatSettings.DecimalSeparator = dc;
	alpha.clear();
}

// -----------------------------------------------------------------------------//
TNet __fastcall TNetF::Points(bool compactPoints) {
	long i, j;
	pathType path;
	// if(!alpha.size())
	makeAlpha(alpha, 0, *this);

	TNet result(alpha.size(), Dim, perfomance, Res);
	result.copyNetFrom(*this);
	result.isVirtual = false;
	// alphType::iterator where;
	i = 0;
	for (i = 0; i < Count; i++) {
		for (j = 0; j < Dim; j++)
			result.v->v->v[i][j] = alpha[i] * getIJ(i, j);
		i++;
	};

	return result;
};

// -----------------------------------------------------------------------------//
Vector __fastcall TNetF::getBorderPoint(long index, const Vector& psi) {
	LDouble extr, val, sk;
	long i, Ind;
  	bool extr_exist = false;
  //	pathType path;
	Vector result = *getVecAt(index);
	result.update();
   	//perfomance = 0;
	 GetExtrDirection(psi, scm, convCriteria1, opMin, ZAware, index, NULL, extr, *this);
	// cout<<result;

  /*	for (i = 0; i < Count; i++) {
		sk = scm(i, psi, this,0);
		if (sk > 0.0) {
			val = f->v->v[i] / sk;
			if (!extr_exist) {
				extr_exist = true;
				extr = val;
				// L[i] = true;
			}
			else if (val < extr)
				extr = val;
		}
	}

	/**/
	result *= extr;
   	result.update();
	return result;

}

// -----------------------------------------------------------------------------//
void __fastcall TNetF::saveAsVrml(string) {
	//
}

// ------------------------------ Clear ---------------------------------------//
// ��������� �������������� ����
void __fastcall TNetF::Clear() {
	long i; // ,j;
	TNet::Clear();
	if (f != NULL)
		for (i = 0; i < Count; i++)
			f->v->v[i] = 0;
}

// ------------------------------copy variables -------------------------//
void /* !inline */ __fastcall TNetF::copyNetFFrom(const TNetF& NetF,
	bool isBuildFunc) {
	alpha = NetF.alpha;
	f = Vector::copy(NetF.f, f);
	fStr = NetF.fStr;
	// vars=NULL;
	if (isBuildFunc)
		SetFunc(fStr);
	/* if(f!=NULL)
	 f->v->linkCount++;
	/* */
	zeroPrecision = NetF.zeroPrecision;

}

// ------------------------------inintialize pointers -------------------------//
void /* !inline */ __fastcall TNetF::initNetFDefault() {
	// func=NULL;
	vars = NULL;
	// cache = NULL;
	f = NULL;
	// updSum = 0;
	alphaMode = false;
	t = 0;
#ifdef _WIN64
	lib_name = "SICx64.dll";
#else
	lib_name = "SICx32.dll";
#endif
	dll = new TDll(lib_name.c_str());

	TDllProcV0 sic_cretab(*dll, "sic_cretab");
	sic_cretab();
#ifdef _WIN64
	TDllStdProcV1<TSIC_Data64*>sic_init(*dll, "sic_init");
#else
	TDllStdProcV1<TSIC_Data32*>sic_init(*dll, "sic_init");
#endif
	sic_init(&sic);

}

// ------------------------------------Dynamically Update----------------------------------//
void __fastcall TNetF::dynUpdate() {
	int i, j;
	TNetF *f1 = NULL;
	Vector *res = NULL;
	if (cache == NULL) {
		cache = new Vector(Dim);
		cacheCurrent = -1;
	}
	if (!updated) {
		updated = true;
		if (!umx) {
			// Vector res(u_mx->v->m);
			if (u_mx->v->m != Dim)
			{ // � ������ ���� ���������� ����������� - ������ ����� �����
				res = new Vector(u_mx->v->m);
				f1 = new TNetF(u_mx->v->m, perfomance, Res);
			}
			else { // ����� ��������� � ������������
				res = new Vector(Dim);
				f1 = this;
			}
		}
		for (i = 0; i < Count; i++) {
			if (!umx) {
				for (j = 0; j < Dim; j++)
					res->v->v[j] = getIJ(i, j);
				(*res) = (*u_mx) * (*res);
				f1->f->v->v[i] = oporn(*res, t, upd);
			}
			else if (upd != 1.0)
				f->v->v[i] *= upd;

		}

		if ((!umx) && (u_mx->v->m != Dim)) {
			*this = *f1;
			delete f1;
		}
		delete res;
	}
	cacheCurrent = -1;
}

// ------------------------------------Get IJ----------------------------------//
inline double __fastcall TNetF::getIJ(long current, int coordNumber) {
	double res;
	if (isVirtual) {
		checkCacheVector();
		if (current != cacheCurrent) {
			cache = parseCoordinate(current);
			cacheCurrent = current;
			cache->norm(halfRes);
		}

		if (!updated)
			dynUpdate();
		res = cache->v->v[coordNumber];
	}
	else {
		if (!updated)
			update();
		res = v->v->v[current][coordNumber];
	}

	return ((alphaMode) && ((long)alpha.size() == Count)) ?
		res * alpha[current] : res;
}

// ------------------------------------Smoothing-------------------------------//
void __fastcall TNetF::smoothFunction(double epsilon) {
	int i;
	update();
	for (i = 0; i < Count; i++) {
		f->v->v[i] = (f->v->v[i] + epsilon) / 2 +
			sqrt(pow(epsilon, 2) + pow(((f->v->v[i] - epsilon) / 2), 2));
	}
}

// ------------------------------- getVecAt -----------------------------------//

Vector* __fastcall TNetF::getVecAt(long i) {
	int j;
	for (j = 0; j < Dim; j++)
		cache->v->v[j] = getIJ(i, j);
	return cache;
}
