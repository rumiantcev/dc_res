// ---------------------------------------------------------------------------
#ifndef traectoryH
#define traectoryH

#include <algorithm>
#include<vector>
#include"Vector.h"
#include "Matrix.h"
#include"Netfunc.h"
#include"general.h"
using namespace std;


typedef vector<Vector*>VecOfVec;


// ---------------------------------------------------------------------------
class Traectory {
public:
	Vector x0; // ,x;
	Vector psi;
	LDouble mu, T;
	long ind;
	vector<TNetF*>NetList;
	VecOfLong psiExtr; // �������� psi �� ������� ����������� extr
  //	VecOfVec uu_i, vv_i, xx_i; //��������� � ���������� ��� ���������������� �������
	Matrix *u_i, *v_i, *x_i; //��������� � ���������� ��� �������������� �������
	VecOfVec vecBrige;
	VecOfVec w;

	~Traectory();
	Traectory();

};
// ---------------------------------------------------------------------------
#endif
