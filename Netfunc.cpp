// ---------------------------------------------------------------------------
#pragma hdrstop

#include "Netfunc.h"
// ---------------------------------------------------------------------------
#pragma package(smart_init)

// ----------------------------------------------------------------------------//
LDouble convCriteria( long num, LDouble sk, const TNetF *v) {
	return sk;
}

// ----------------------------------------------------------------------------//
LDouble convCriteria1(long num, LDouble sk, const TNetF *v) {
	return (v->f->v->v[num]) / sk;
}

// ----------------------------------------------------------------------------//
LDouble convCriteria2(long num, LDouble sk, const TNetF *v) {
	return v->f->v->v[num] * sk;
}

// ------------------------------------Create----------------------------------//
/* !inline */
void __fastcall TNetF::create( unsigned long Dim /* , long perf, long res */) {
	unsigned long i, j;
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
	vars = NULL;
}

// ------------------------------constructor ----------------------------------//
__fastcall TNetF::TNetF(unsigned long dim, long perf, unsigned long res)
	: TNet(dim, perf, res, true) {
	initNetFDefault();
	create(dim );
	SetFunc("0");
}

// ------------------------------constructor ----------------------------------//
__fastcall TNetF::TNetF(unsigned  long Dim, long perf, unsigned  long res, string fstr)
	: TNet(Dim, perf, res, true) {
	/* TODO -orum -caddon : Сделать конструктор с возможностью переключения Virt/Non Virt net */
	initNetFDefault();
	create(Dim );
	SetFunc(fstr);
}

// ------------------------------constructor ----------------------------------//
__fastcall TNetF::TNetF(unsigned long mm, unsigned long nn, long perf, unsigned long res)
	: TNet(nn, perf, res, true) {
	initNetFDefault();
	f = new Vector(mm);
	create(nn );
	SetFunc("0");
}

// ------------------------copy constructor-----------------------------|------//
__fastcall TNetF::TNetF(TNet& Net) : TNet(Net) {
	initNetFDefault();
	create(Dim );
	SetFunc("0");

}

/* */// ------------------------copy constructor-------------------------------//
__fastcall TNetF::TNetF(TNet& Net, const string& fstr) : TNet(Net) {
	initNetFDefault();
	create(Dim);
	SetFunc(fstr);
}

// ------------------------copy constructor----------------------------|-------//
__fastcall TNetF::TNetF(const TNetF& NetF) : TNet(NetF) {
	initNetFDefault();
	copyNetFFrom(NetF, false);
	SetFunc(NetF.fStr);
	is_empty = NetF.is_empty;
	maxRad = NetF.maxRad;
}

/* */// ------------------------copy constructor-------------------------------//
__fastcall TNetF::TNetF(TNetF& Net, const string& fstr) : TNet(Net) {
	initNetFDefault();
	copyNetFFrom(Net, false);
	SetFunc(fstr);
	is_empty = Net.is_empty;
	maxRad = Net.maxRad;
}

// ------------------------------- destructor ---------------------------------//
void TNetF:: operator delete(void *p) {
	if (p == 0) return;

	TNetF *ptr = static_cast<TNetF*>(p);

	if (ptr->f == NULL)
		delete(void*) p;
	else
		p = NULL;
}

// ---------------------------- destructor-------------------------------------//
__fastcall TNetF::~TNetF() {
		delete f;
		f = NULL;
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
}

// -------------------------------- << ----------------------------------------//
ostream& __fastcall operator << (ostream& out_data, TNetF& C) {
	unsigned long SurfaceDim = C.Dim - 1;
	unsigned long NumOfSur = C.Dim * 2;
	unsigned long NumOfPoints = C.Res;
	unsigned long i, j;

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
}

// ----------------------------------- >> -------------------------------------//
istream& __fastcall operator >> (istream& in_data, TNetF& C) {
	char c = 0;
	LDouble tempVal;
	unsigned long i, j;

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
				}
			}
		}
	in_data.get(c);
	while (c != ';' && c != EOF)
		in_data.get(c);
	in_data.get(c);
	return in_data;
}
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
}

// ----------------------------------- + --------------------------------------//
TNetF& __fastcall TNetF:: operator += (const TNetF& B) {
	unsigned long i;
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
		f->v->v[i] += pB->f->v->v[i];
	return *this;
}

// ----------------------------------- + --------------------------------------//
const TNetF __fastcall operator +(const TNetF& A, const TNetF& B) {
	return TNetF(A) += B;
}

// ----------------------------------- Геометрическая - -----------------------//
TNetF& __fastcall TNetF:: operator -= (const TNetF& B) {
	unsigned long i, /*j, */ m;
	Vector vec(Dim), st(Dim);
	TNetF st0Net(Dim, perfomance, Res);
	LDouble coeff = (LDouble)Dim / Count;
	bool /*extr_exist,*/ *L;
   //	FindPath ynPath = perfomance > 1 ? yPath : nPath;

   //	long k = 1;
   //	LDouble tmin = _extr_tmin_param, tmax = _extr_t0_param, t = tmax, a =
   //		LDouble(_lrand() % Count) / Count, c = _extr_e_val;

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
	// расчёт опорной функции центра Штейнера
	st = 0;
	for (i = 0; i < Count; /* i=i+2 */ i++) {
		for (m = 0; m < Dim; m++)
			vec.v->v[m] = st0Net.getIJ(i, m);
		st += (f->v->v[i]*vec);
	}
	st *= coeff;

	// сдвигаем множество "на центр Штейнера" так чтобы 0 был в центре
	for (i = 0; i < st0Net.Count; i++) {
		st0Net.f->v->v[i] = scm(i, st, &st0Net,0);
		f->v->v[i] -= st0Net.f->v->v[i];
	}
	// рассчитываем lambda_i
	L = new bool[Count];
	// alpha.clear();
	alpha.reserve(Count);

	//makeAlpha(alpha, L, st0Net); //тут по ходу была ошибка
	makeAlpha(alpha, L, *this);

	// Собственно овыпукление и сдвиг множества обратно
	for (i = 0; i < st0Net.Count; i++) {
		if (L[i]) {
			//extr_exist = false;
			for (m = 0; m < Dim; m++)
				vec.v->v[m] = st0Net.getIJ(i, m);

			/*j =*/ GetExtrDirection(vec, scm1, convCriteria, opMax, nZAware, i, &alpha,  f->v->v[i],  st0Net);
		}
		f->v->v[i] += st0Net.f->v->v[i];
	}
	alpha.clear();
	delete[]L;
	return *this;
}

// ----------------------------------- Геометрическая - -----------------------//
const TNetF __fastcall operator -(const TNetF& A, const TNetF& B) {
	return TNetF(A) -= B;
}

// ------------------------------------Update----------------------------------//
void __fastcall TNetF::update() {
	detach();
	if (!isVirtual) {
		if (!updated) { // обновляем за счет матрицы или множителя
			unsigned long i, j;
			TNetF *f1 = NULL;
			Vector *res = NULL;
			updated = true;
			if (!umx) {
				if (u_mx->v->m != Dim)
				{ // в случае если поменялась размерность - создаём новую сетку
					res = new Vector(u_mx->v->m);
					f1 = new TNetF(u_mx->v->m, perfomance, Res);
				}
				else { // иначе оперируем с существующей
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
			if (res != NULL)
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
		unsigned long i;
		Vector *vv = NULL;
		vv = Vector::copy(f, vv);
		delete f;
		f = new Vector(vv->v->size);
		for (i = 0; i < vv->v->size; i++)
			f->v->v[i] = vv->v->v[i];
		delete vv;
	}
}

// ----------------------------------- * --------------------------------------//
TNetF& __fastcall TNetF:: operator *= (const Matrix& A) {
	bool isDimChanged = (Dim != A.n()) ? true : false;

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
}

// ----------------------------------- * --------------------------------------//
TNetF& __fastcall TNetF:: operator *= (const LDouble &a) {
	upd *= a;
	updated = false;
	return *this;
}

// ----------------------------------- + --------------------------------------//
/* TNetF& __fastcall TNetF:: operator += (const LDouble &a) {
 updSum += a;
 updated = false;
 return *this;
 }
/* */
// ----------------------------------- * --------------------------------------//
const TNetF __fastcall operator *(const LDouble a, const TNetF& B) {
	return TNetF(B) *= a;
}

// ----------------------------------- + --------------------------------------//
/* const TNetF __fastcall operator +(const LDouble a, const TNetF& B) {
 return TNetF(B) += a;
 }
/* */
// ------------------------------------AddVariables----------------------------//
void __fastcall TNetF::AddVariables() {
	unsigned long i;
	string var;
	ostringstream ss;

	if (vars != NULL) {
		delete[]vars;
		vars = NULL;
	}
	if (vars == NULL)
		vars = new LDouble[Dim];

#ifdef _WIN64
	TDllStdProcV3<TSIC_Data64*, char*, LDouble*>sic_avarf(*dll, "sic_avarf");
#else
	TDllStdProcV3<TSIC_Data32*, char*, LDouble*>sic_avarf(*dll, "sic_avarf");
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


// ----------------------------------------------------------------------------//
LDouble __fastcall TNetF::oporn(const Vector &x, LDouble t, int sign) {
	LDouble result;
	DWORD err;

	unsigned long i;
	for (i = 0; i < Dim; i++)
		vars[i] = sign * x.v->v[i];

	this->t = t;
#ifdef _WIN64
	TDllStdProc2<LDouble, TSIC_Data64*, DWORD*>sic_exec(*dll, "sic_exec");
#else
	TDllStdProc2<LDouble, TSIC_Data32*, DWORD*>sic_exec(*dll, "sic_exec");
#endif
	result = sic_exec(&sic, &err);
	return result;
}

// ----------------------------------------------------------------------------//
LDouble __fastcall TNetF::oporn(const Vector &x, LDouble t,
	const string &funcStr, int sign) {
	SetFunc(funcStr);
	return oporn(x, t, sign);
}


// ----------------------------------------------------------------------------//
void __fastcall TNetF::oporn(LDouble t, int sign) {
	unsigned long i, j;
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
		f->v->v[i] = oporn(vv, t, sign);
	}
}

// ----------------------------------------------------------------------------//
TNetF __fastcall TNetF::opornx(LDouble t, int sign) {
	oporn(t, sign);
	return *this;
}

// ----------------------------------------------------------------------------//
LDouble __fastcall TNetF::oporn(const Vector &x, LDouble t, int sign,
	const Matrix &A) {
	umx = true;
	Vector xx = A * x;
	return oporn(xx, t, sign);
}

// ----------------------------------------------------------------------------//
void __fastcall TNetF::oporn(LDouble t, int sign, const Matrix &A) {
	unsigned long i, j;
	Vector vv(Dim);

	for (i = 0; i < Count; i++) {
		for (j = 0; j < Dim; j++)
			vv.v->v[j] = v->v->v[i][j];
		f->v->v[i] = oporn(vv, t, sign, A);
	}
}
// ----------------------------------------------------------------------------//

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
		AddVariables();
		fStr = fstr;
#ifdef _WIN64
		DWORD sop = SIC_OPT_DEFAULT_X64;
#else
		DWORD sop = SIC_OPT_DEFAULT_X32;
#endif
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



// -------------------------- findExtrDirectionSlow ---------------------------//
unsigned long __fastcall TNetF::findExtrSlowDirection(const Vector& vec, scM scmul,
	cCrit crit, ZeroAware isZeroAware, OpType extrOper, long index, alphType* coeff, LDouble &extr,  TNetF& net) {
	// Поиск максимума скалярного произведения scmul на сетке простым перебором в направлении заданного вектора

	unsigned long i, j = 0;
	LDouble val, sc;
	bool isExtrExist = false, isGrZero, isExtr;

	for (i = 0; i < Count; i++) {
		sc = scmul(i, vec, &net,coeff);
		// на будущее
		// считаем скалярное произведение текущего вектора  сетки (сurrent)  и заданного (vec)
		isGrZero = sc > zeroPrecision;
		if ((isZeroAware == nZAware) || ((isZeroAware == ZAware) && isGrZero)){
			val = crit(i, sc, this);
			if (!isExtrExist){
				// Если экстремума  ещё не было - берём первое значение и запоминаем как экстремум
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



// -------------------------- findMaxDirectionAnnealing -----------------------//
unsigned long  __fastcall TNetF::findExtrAnnealingDirection(const Vector& vec, scM scmul,
	cCrit crit, ZeroAware isZeroAware,OpType extrOper,  long index, alphType* coeff, LDouble &extr, TNetF& net) {

	long j, jj;
	LDouble val,sc;
	LDouble tmin = Environment::instance()._extr_tmin_param, tmax = Environment::instance()._extr_t0_param, t = tmax, p, a, c =
		Environment::instance()._extr_e_val;

	bool isGrZero, isExtrExist = false, isExtr;

	j = _lrand() % (net.Count);
	while (t > tmin) {
		sc = scmul(j, vec, &net, coeff);
		isGrZero = sc > zeroPrecision;
		if ((isZeroAware == nZAware) || ((isZeroAware == ZAware) && isGrZero)) {
			val = crit(j, sc, this);
			if (!isExtrExist) {
				// Если экстремума  ещё не было - берём первое значение и запоминаем как экстремум
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
					a = LDouble(_lrand() % net.Count) / net.Count;
					if (a > p) {
						extr = val;
						t = t * c;
					}
				}

			}
		}
		a = LDouble(_lrand() % net.Count) / net.Count;
		jj = signof(a - 0.5) * t * (pow((1 + 1 / t), fabs(2 * a - 1)) - 1)
				* net.Count;
		j = abs(jj + j) % net.Count;
	}
	return j;
}

// -------------------------- findExtrSlowGlobal ------------------------------//
long __fastcall TNetF::findExtrSlowGlobal(OpType extrOper, LDouble &extr) {
	// Поиск минимума опорной функции на сетке простым перебором
	unsigned long i, j = 0;
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

// -------------------------- findExtrAnnealingGlobal -------------------------//
long  __fastcall TNetF::findExtrAnnealingGlobal(OpType extrOper, LDouble &extr) {
	// Поиск максимума опорной функции на сетке методом отжига
	long  j=_lrand() % (Count), jj;
	LDouble val;
	LDouble tmin = Environment::instance()._extr_tmin_param, tmax = Environment::instance()._extr_t0_param, t = tmax, p, a, c =
		Environment::instance()._extr_e_val;
	bool isExtr;

	extr = f->v->v[j];
	val = extr;

   //	j = _lrand() % (Count);
	while (t > tmin) {
		val = f->v->v[j];
		extrOper == opMax? isExtr = (val > extr) :  isExtr = (val < extr);
		if (isExtr) {
			extr = val;
			t = t * c;
		}
		else {
			p = 1 / (1 + exp(-fabs(val - extr) / t));
			a = LDouble(_lrand() % Count) / Count;
			if (a > p) {
				extr = val;
				t = t * c;
			}
		}
		a = LDouble(_lrand() % Count) / Count;
		jj = signof(a - 0.5) * t * (pow((1 + 1 / t), fabs(2 * a - 1)) - 1)
			* Count;
		j = abs(jj + j) % Count;
	}

	return j;
}

// -------------------------- selectExtrX -------------------------------------//
unsigned long /* !inline */ TNetF::selectExtrX(const Vector& vec, scM scmul, cCrit crit,
	 long current,  long& result, LDouble &extr, OpType isMax,
	ZeroAware isZeroAware, bool &isExtrExist, TNetF& net) {
	LDouble sc, val = extr;
	bool  isGrZero;

	sc = (*scmul)(current, vec, &net,0);
	// считаем скалярное произведение текущего вектора  сетки (сurrent)  и заданного (vec)
	isGrZero = sc > zeroPrecision;
	if ((isZeroAware == nZAware) || ((isZeroAware == ZAware) && isGrZero)) {
		val = (*crit)(current, sc, this);
		if (!isExtrExist){
			// Если экстремума  ещё не было - берём первое значение и запоминаем как экстремум
			isExtrExist = true;
			extr = val;
			result = current;
		}
		else {
			bool isProxy;
			isMax == opMax ? isProxy = val >= extr : isProxy = val <= extr;
			// если ищем максимум, то проверяем на то что новое значение >= запомненному, если минимум -<= запомненному
			if (isProxy) {
				extr = val;
				result = current;
			}
		}
		return result;
	}
	result = current;
	return current;
}

// -------------------------- findExtrFastX -----------------------------------//
unsigned long __fastcall TNetF::findExtrFastXDirection(const Vector& vec, scM scmul,
	cCrit crit, OpType isMax, ZeroAware isZeroAware, long index, alphType* coeff, LDouble &extr,
	TNetF& net){
	// Поиск экстремума скалярного произведения scmul на сетке в направлении заданного вектора пользуясь тем, что функция выпуклого коспакта также выпуклая

	unsigned	long i,  je, k1, k2;
	long  j = index;
	bool isExtrExist = false, isNotExtrFound = true, borderChanged;
	seekType seekPath;

	seekPath.clear();

	while (isNotExtrFound) {
		/* DONE -orum -crepair : Переделать критерий выхода из цикла */
		if (!seekPath[j]) {
			je = j;
			seekPath[j] = true;
			for (i = 0; i < Dim; i++) {
				try {
					k1 = shift(je, i, 1, borderChanged);
				}
				catch (exInvalidMoveDirection) { // do nothing;
				}
				if (k1 != je)
					j = selectExtrX(vec, scmul, crit, k1, j, extr, isMax,
					isZeroAware, isExtrExist,  net);
				try {
					k2 = shift(je, i, -1, borderChanged);
				}
				catch (exInvalidMoveDirection) { // do nothing;
				}
				if (k2 != je)
					j = selectExtrX(vec, scmul, crit, k2, j, extr, isMax,
					isZeroAware, isExtrExist,  net);
			}
			if (j < 0)
				j = shift(je, 0, 1, borderChanged);

		}
		else {
			isNotExtrFound = false;
			seekPath.clear();
		}

	}
	return j;
}

// -------------------------- findExtrFastXGlobal -----------------------------//
unsigned long __fastcall TNetF::findExtrFastXGlobal(OpType isMax, unsigned long index,
	LDouble &extr)
	// Поиск экстремума опорной функции пользуясь тем, что функция выпуклого коспакта также выпуклая
{
	unsigned long i,  je, k1, k2;
	long j = index;
	LDouble val;
	bool isExtrExist = false, isNotExtrFound = true, borderChanged, isProxy;
	seekType seekPath;

	seekPath.clear();

	while (isNotExtrFound) {
		/* DONE -orum -crepair : Переделать критерий выхода из цикла */
		if (!seekPath[j]) {
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
		} else {
			isNotExtrFound = false;
			seekPath.clear();
		}

	}
	return j;
}

// ----------------------------- GetExtr в направлении заданного вектора-------//
long __fastcall TNetF::GetExtrDirection(const Vector& vec, scM scmul,
	cCrit crit, OpType extrOper, ZeroAware isZeroAware, long index, alphType* coeff, LDouble& extr,
	TNetF& net) {

	if (perfomance == optNone)
		return findExtrSlowDirection(vec, scmul, crit, isZeroAware, extrOper, index, coeff,
			extr, net);

	if (perfomance == optAnnealing)
			return findExtrAnnealingDirection(vec, scmul, crit, isZeroAware,  extrOper,
			index, coeff, extr, net);

	if ((perfomance == optGradient)||(perfomance == optTimS))
		return findExtrFastXDirection(vec, scmul, crit, extrOper, isZeroAware,
			index, coeff, extr, net);
		 // Fast Evaluator

	return -1;
}

// ----------------------------- GetExtrGlobal---------------------------------//
long __fastcall TNetF::GetExtrGlobal(OpType extrOper, long index, LDouble& extr) {
	if (perfomance == optNone)
			return findExtrSlowGlobal(extrOper,extr);

	if (perfomance == optAnnealing)
			return findExtrAnnealingGlobal(extrOper,extr);

	if ((perfomance == optGradient)||(perfomance == optTimS)) {
		return findExtrFastXGlobal(extrOper, index, extr); /* */ // Fast Evaluator
	}
	return -1;
}

// ----------------------------------------------------------------------------//
LDouble  __fastcall scm(long num, const Vector &vec, TNetF *v, alphType* coeff) {
	LDouble sc = 0;
	unsigned long m;
	for (m = 0; m < vec.v->size; m++)
		sc += v->getIJ(num, m) * vec[m];
	return sc;
}

// ----------------------------------------------------------------------------//
LDouble __fastcall scm1(long num, const Vector &vec, TNetF *v, alphType* coeff) {
	LDouble sc = 0;
	unsigned long m;
	for (m = 0; m < vec.v->size; m++)
		sc += v->getIJ(num, m) * vec[m];
	return  (*coeff)[num]*sc;
}
// ----------------------------------------------------------------------------//
LDouble __fastcall scm2(long num, const Vector &vec, TNetF *v, alphType* coeff) {
	LDouble sc = 0;
	unsigned long m;
	for (m = 0; m < vec.v->size; m++)
		sc += (-v-> /* TNet:: */ getIJ(num, m)) * vec[m];
	return v->f->v->v[num] + sc;
}

// ----------------------------------------------------------------------------//
void __fastcall TNetF::makeAlpha(alphType& alpha, bool* L, TNetF &net) {
	unsigned long i, k, m;
	LDouble extr, sk;
	bool extr_exist;
	Vector vec(Dim);

	long /*k = 1,*/j, jj, ind=0;
	LDouble tmin = Environment::instance()._extr_tmin_param, tmax = Environment::instance()._extr_t0_param, t = tmax, p, a =
		LDouble(_lrand() % net.Count) / net.Count, c = Environment::instance()._extr_e_val;
	LDouble ldZeroDf = Environment::instance().ldZeroDf;
	alpha.reserve(Count);
	net.is_empty = true;
	net.maxRad = 0;
	for (i = 0; i < Count; i++) {
		if (L != NULL)
			L[i] = false;
		for (m = 0; m < Dim; m++)
			vec[m] = net.getIJ(i, m);
		// ------------GetMin by criteria---------------------------------------

		alpha[i] = 0.0;

		extr_exist = false;
		if (perfomance == optNone) {
			// поиск последовательным перебором
			for (k = 0; k < (net.Count); k++) {
				sk = scm(k, vec, &net,NULL);
				if (sk > 0.0) {
					extr = f->v->v[k] / sk;
					if (!extr_exist) {
						extr_exist = true;
						alpha[i] = extr;
						//if (L != NULL)
						//	L[i] = true;
						ind = k;
					}
					else if (alpha[i] > extr){
						alpha[i] = extr;
						ind = k;
					}
				}
			}
			if (L != NULL)
				L[ind] = true;
		}

		if (perfomance == optAnnealing) {
			// поиск методом эмуляции отжига
			j = _lrand() % (net.Count);
			while (t > tmin) {
				sk = scm(j, vec, &net,NULL);
				if (sk > zeroPrecision) {
					extr = f->v->v[j] / sk;//convCriteria1(j,sk,&net);
					if (!extr_exist) {
						extr_exist = true;
						alpha[i] = extr;
					   //	if (L != NULL)
					   //		L[i] = true;
						ind = j;
					}
					else if (alpha[i] >= extr) {
						alpha[i] = extr;
						t = t * c;
						ind = j;
						// t = tmax / k; // раскомменировать Если отжиг по Коши
						//k++;  // раскомменировать Если отжиг по Коши
					}
					else {
						p = 1 / (1 + exp(-fabs(alpha[i] - extr) / t));
						// p= t/(M_PI*(pow(fabs(extr-alpha[i])/net.Count,2)+pow(t,2)));  // раскомменировать Если отжиг по Коши
						a = LDouble(_lrand() % net.Count) / net.Count;
						if (a > p) {
							alpha[i] = extr;
							t = t * c;
							ind = j;
							// t = tmax / k; // раскомменировать Если отжиг по Коши
							//k++;    // раскомменировать Если отжиг по Коши
						}
					}
				}
				a = LDouble(_lrand() % net.Count) / net.Count;
				jj = signof(a - 0.5) * t * (pow((1 + 1 / t), abs(2 * a - 1)) - 1) * net.Count;
				// jj = t * tan(M_PI*(a-0.5))*net.Count;  // раскомменировать Если отжиг по Коши
				j = abs(jj + j) % net.Count;
			}
			if (L != NULL)
				L[ind] = true;
		}

		if (alpha[i]> ldZeroDf)
			net.is_empty = false;
		if (alpha[i]>net.maxRad)
			net.maxRad = alpha[i];
	}
}

// ----------------------------------------------------------------------------//

void __fastcall TNetF::Conv(bool *L) {
	unsigned long i, /*j,*/ m/*, k*/;
	LDouble coeff = (LDouble)Dim / Count;
	TNetF st0Net(Dim, perfomance, NumOfPoints);
	Vector vec(Dim), st(Dim);

	// расчёт опорной функции центра Штейнера
	st = 0;
	for (i = 0; i < Count; i++) {
		for (m = 0; m < Dim; m++)
			vec.v->v[m] = st0Net.getIJ(i, m);
		st += (f->v->v[i]*vec);
	}
	st *= coeff;

	// сдвигаем множество "на центр Штайнера" так чтобы 0 был в центре
	for (i = 0; i < st0Net.Count; i++) {
		st0Net.f->v->v[i] = scm(i, st, &st0Net,0);
		f->v->v[i] -= st0Net.f->v->v[i];
	}

	// рассчитываем lambda_i
   //	makeAlpha(alpha, L, st0Net);   //тут по ходу была ошибка
	makeAlpha(alpha, L, *this);

	// собственно овыпукляем и сдвигаем назад
	for (i = 0; i < Count; i++) {
		if (!L[i]) {  //была ошибка
			//extr_exist = false;
			for (m = 0; m < Dim; m++)
				vec.v->v[m] = /*st0Net.*/getIJ(i, m);   // была ошибка
			/* j = */GetExtrDirection(vec, scm1, convCriteria, opMax, nZAware, i, &alpha,  f->v->v[i],  *this);    // была ошибка
		}
		f->v->v[i] += st0Net.f->v->v[i];
		// сдвиг множества выпуклой оболочки на исходное место
	}
	alpha.clear();
}

// ----------------------------------------------------------------------------//

void __fastcall TNetF::ConvTimS(bool *L) {
	unsigned long i,j, k, m, lInd, rInd, lIndEx, rIndEx, pr_rInd, pr_lInd, indNormal, ind;
	LDouble coeff = (LDouble)Dim / Count;
	LDouble sk,extr, pr_extr;
	LDouble ldZeroDf = Environment::instance().ldZeroDf;
	TNetF st0Net(Dim, perfomance, NumOfPoints);
	Vector vec(Dim), st(Dim);
	//вектора "левых" и "правых" соседей для каждой точки
	//long lrV[Count][Dim*2];
	t_Mx<long>  lrV(Count,Dim*2);

	int  normSign;
	bool isBC, isExtr, isLExtr, isRExtr, extr_exist;

	// расчёт опорной функции центра Штейнера
	st = 0;
	for (i = 0; i < Count; i++) {
		for (m = 0; m < Dim; m++)
			vec.v->v[m] = st0Net.getIJ(i, m);
		st += (f->v->v[i]*vec);
	}
	st *= coeff;

	// сдвигаем множество "на центр Штайнера" так чтобы 0 был в центре
	for (i = 0; i < st0Net.Count; i++) {
		st0Net.f->v->v[i] = scm(i, st, &st0Net,0);
		f->v->v[i] -= st0Net.f->v->v[i];
	}

	for (i = 0; i < Count; i++) {
		//ищем нормаль и её знак - т.к. по координате нормали нет необходимости искать левого и правого соседей
		indNormal = static_cast<long>(i/NumOfPoints); //замена trunc через приведение типа
		normSign = (!!(indNormal & 1)) ? 1:-1;
		indNormal = static_cast<long>(indNormal/2);  //замена trunc через приведение типа
		//и проверяем на локальную выпуклость
		for (m = 0; m < Dim; m++)
			if (m==indNormal)
				lrV[i][m]=lrV[i][m+Dim] = -1;
			else {
				lInd = shift(i,m,-1,isBC);
				rInd = shift(i,m, 1,isBC);
				//!!!уточнить неравенство Йенсена верно всюду или только для положительно определённых функкций
				//f[i]; - как rvalue не работает, только (*f)[i] или f->v->v[i]
				if (f->v->v[i] >= ((f->v->v[lInd]+f->v->v[rInd])/2))//есть локальная выпуклость
					lrV[i][m]=lrV[i][m+Dim] = 0;
				else {//нашли локальную вогнутость
					isExtr = isLExtr = isRExtr = false;
					pr_lInd = pr_rInd = i;
					lIndEx = rIndEx = 0;
					while (isExtr){
						//шаг влево
						pr_lInd = lInd;
						lInd = shift(pr_lInd,m,-1,isBC);
						if(lrV[lInd][m] > 0){//нашли существующий путь
							lIndEx = lrV[lInd][m];
							rIndEx = lrV[lInd][m+Dim];
							isExtr = true;
							continue;
						}
						if (f->v->v[pr_lInd] >= ((f->v->v[lInd]+f->v->v[rInd])/2))//шаг влево. есть локальная выпуклость слева
							isLExtr = true;
						else
							lrV[lInd][m]= pr_lInd;
						if (f->v->v[pr_rInd] >= ((f->v->v[lInd]+f->v->v[rInd])/2))//шаг влево. есть локальная выпуклость справа
							isRExtr = true;
						else
							lrV[rInd][m+Dim]= pr_rInd;

						if(isLExtr && isRExtr){//если выпукло с всех сторон
							isExtr = true;
							continue;
						}

						//шаг вправо
						pr_rInd = rInd;
						rInd = shift(pr_lInd,m, 1,isBC);
						if(lrV[rInd][m] > 0){//нашли существующий путь
							lIndEx = lrV[rInd][m];
							rIndEx = lrV[rInd][m+Dim];
							isExtr = true;
							continue;
						}
                    	if (f->v->v[pr_lInd] >= ((f->v->v[lInd]+f->v->v[rInd])/2))//шаг влево. есть локальная выпуклость слева
							isLExtr = true;
						else
							lrV[lInd][m]= pr_lInd;

						if (f->v->v[pr_rInd] >= ((f->v->v[lInd]+f->v->v[rInd])/2))//шаг вправо. есть локальная выпуклость справа
							isRExtr = true;
						else
							lrV[rInd][m+Dim]= pr_rInd;

						if(isLExtr && isRExtr)//если выпукло с всех сторон
							isExtr = true;
					}

					if (lIndEx == 0) { //на имеющийся путь не натыкались
						lIndEx = lInd;
						rIndEx = rInd;
					}
					//обратный ход слева  с прометкой левых и правых соседей
					j= lIndEx;
					lrV[j][m]= lIndEx;
					lrV[j][m+Dim]= rIndEx;
					j = pr_lInd;
					while (j != i){ //идём слева к i
						pr_lInd = lrV[j][m];
						lrV[j][m]= lIndEx;
						lrV[j][m+Dim]= rIndEx;
						j = pr_lInd;
					}

					//обратный ход справа с прометкой левых и правых соседей
					j= rIndEx;
					lrV[j][m]= lIndEx;
					lrV[j][m+Dim]= rIndEx;
					j = pr_rInd;
					while (j != i){ //идём справа к i
						pr_rInd = lrV[j][m+Dim];
						lrV[j][m]= lIndEx;
						lrV[j][m+Dim]= rIndEx;
						j = pr_rInd;
					}

				}


			}

	}
	//на выходе имеем массив, в которой для каждой точки указаны номера 2*(Dim-1) соседей в которых удовлетворяется неравенство Йенсена
	// рассчитываем lambda_i
	is_empty = true;
	maxRad = 0;

	for (i = 0; i < Count; i++) {
		if (L != NULL)
			L[i] = false;
		for (m = 0; m < Dim; m++)
			vec[m] = getIJ(i, m);
		alpha[i] = 0.0;

		extr_exist = isExtr = false;
		 // поиск extr градиентным методом (покоординатного спуска)
		k = i;
		j = -1;
		ind = -1;
		while (!isExtr){
			for (m = 0; m < Dim;) {
				if(lrV[k][m]>0)   //на невыпуклом участке берём соотв левого и правого выпуклого соседа
				  j == -1? k = lrV[k][m]: k = lrV[k][m+Dim];
				else
					if(lrV[k][m] == 0) //на выпуклом участке берём левого и правого соседа
						k = shift(k,m, j,isBC);
					else
						continue; //вне поверхности не ищем ничего, просто проскакиваем

				if (j == 1 ) {//цикл проходим дважды, по левым и правым соседям
				  j=-1;
				  m++;
				} else
					j = 1;

				sk = scm(k, vec, this,NULL);
				if (sk > 0.0) {
					sk = f->v->v[k] / sk;
					if (!extr_exist) {
						extr_exist = true;
						alpha[i] = pr_extr  = sk;
						ind = k;
					}else
						if (sk < alpha[i]){
							pr_extr  = alpha[i];
							alpha[i] = sk;
							ind = k;
						}
				}
			}

			if (alpha[i] < 0){ //по идее ситуация невозможная, т.к. обходим только выпуклые области, поэтому вставлена диагностика
			   isExtr = true;
			   cout << "Error: i="<< i <<", k="<<k<<endl;
			} else
				if (pr_extr - alpha[i] < ldZeroDf)  //по построению pr_extr - alpha[i] д.б. > 0
				   isExtr = true;
		}
		if (L != NULL)
			L[ind] = true;

		if (alpha[i]>= ldZeroDf)
			is_empty = false;
		if (alpha[i]> maxRad)
			maxRad = alpha[i];
    }

	// собственно овыпукляем и сдвигаем назад
	for (i = 0; i < Count; i++) {
		if (!L[i]) {
			extr_exist = false;
			for (m = 0; m < Dim; m++)
				vec.v->v[m] = /*st0Net.*/getIJ(i, m);
			extr_exist = isExtr = false;
			k = i;
			j = -1;
			ind = -1;
			while (!isExtr){
				if(lrV[k][m]>0)   //на невыпуклом участке берём соотв левого и правого выпуклого соседа
				  j == -1? k = lrV[k][m]: k = lrV[k][m+Dim];  //по идее опять же всегда должны быть левые и правые соседи т.к. смотрим не точки где зафиксирована локальная невыпуклость
				else
					if(lrV[k][m] == 0) //Если так, то ошибка
						cout << "Error: i="<< i <<", k="<<k<<endl;
					else
						continue; //вне поверхности не ищем ничего, просто проскакиваем

				if (j == 1 ) {//цикл проходим дважды, по левым и правым соседям
				  j=-1;
				  m++;
				} else
					j = 1;

				sk = alpha[k]*scm(k, vec, this, NULL);   //расчёт co_\Psi(c_\psi)
				if (!extr_exist) {
					extr_exist = true;
					pr_extr = extr = sk;
					//ind = k;
				}else
					if (sk > extr){
						pr_extr  = extr;
						extr = sk;
						//ind = k;
					}
				if (extr - pr_extr < ldZeroDf)  //по построению extr - pr_extr д.б. > 0
				   isExtr = true;
			}
			 //j = GetExtrDirection(vec, scm1, convCriteria, opMax, nZAware, i, &alpha,  f->v->v[i], *this);
		}
		f->v->v[i] += st0Net.f->v->v[i];  // сдвиг множества выпуклой оболочки на исходное место
	}
}
// ----------------------------------------------------------------------------//
TNet __fastcall TNetF::Points(/*bool compactPoints*/) {
	unsigned long i, j;
  //	pathType path;
	makeAlpha(alpha, 0, *this);

	TNet result(alpha.size(), Dim, perfomance, Res);
	result.copyNetFrom(*this);
	result.isVirtual = false;
	//i = 0;
	for (i = 0; i < Count; i++) {
		for (j = 0; j < Dim; j++)
			result.v->v->v[i][j] = alpha[i] * getIJ(i, j);
		i++;
	}
	return result;
}

// ----------------------------------------------------------------------------//
Vector __fastcall TNetF::getBorderPoint(unsigned long index, const Vector& psi) {
	LDouble extr;
	//bool extr_exist = false;
	Vector result = *getVecAt(index);
	result.update();
	 GetExtrDirection(psi, scm, convCriteria1, opMin, ZAware, index, NULL, extr, *this);
	result *= extr;
	result.update();
	return result;

}

// ----------------------------------------------------------------------------//
void __fastcall TNetF::saveAsVrml(string) {
	//
}

// ------------------------------ Clear ---------------------------------------//
// Зануление функциональной сети
void __fastcall TNetF::Clear() {
	TNet::Clear();
	if (f != NULL){
		unsigned long i; // ,j;
		for (i = 0; i < Count; i++)
			f->v->v[i] = 0;
	}
}

// ------------------------------copy variables -------------------------------//
void /* !inline */ __fastcall TNetF::copyNetFFrom(const TNetF& NetF,
	bool isBuildFunc) {
	alpha = NetF.alpha;
	f = Vector::copy(NetF.f, f);
	fStr = NetF.fStr;
	if (isBuildFunc)
		SetFunc(fStr);
	zeroPrecision = NetF.zeroPrecision;

}

// ------------------------------inintialize pointers -------------------------//
void /* !inline */ __fastcall TNetF::initNetFDefault() {
	vars = NULL;
	f = NULL;
	alphaMode = false;
	is_empty = false;
	maxRad = 0.0;
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

// ------------------------------------Dynamically Update----------------------//
void __fastcall TNetF::dynUpdate() {

	if (cache == NULL) {
		cache = new Vector(Dim);
		cacheCurrent = -1;
	}
	if (!updated) {
		unsigned long i, j;
		TNetF *f1 = NULL;
		Vector *res = NULL;
		updated = true;
		if (!umx) {
			if (u_mx->v->m != Dim)
			{ // в случае если поменялась размерность - создаём новую сетку
				res = new Vector(u_mx->v->m);
				f1 = new TNetF(u_mx->v->m, perfomance, Res);
			}
			else { // иначе оперируем с существующей
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
inline LDouble __fastcall TNetF::getIJ(unsigned long current, unsigned long coordNumber) {
	LDouble res;
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
	return ((alphaMode) && ((unsigned long)alpha.size() == Count)) ?
		res * alpha[current] : res;
}

// ------------------------------------Smoothing-------------------------------//
void __fastcall TNetF::smoothFunction(LDouble epsilon) {
	unsigned long i;
	update();
	for (i = 0; i < Count; i++) {
		f->v->v[i] = (f->v->v[i] + epsilon) / 2 +
			sqrt(pow(epsilon, 2) + pow(((f->v->v[i] - epsilon) / 2), 2));
	}
}

// ------------------------------- getVecAt -----------------------------------//

Vector* __fastcall TNetF::getVecAt(unsigned long i) {
	unsigned long j;
	for (j = 0; j < Dim; j++)
		cache->v->v[j] = getIJ(i, j);
	return cache;
}
