//---------------------------------------------------------------------------

#pragma hdrstop

#include "pursuerType.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)


//----------------------------------------------------------------------------
pursuerType::pursuerType( unsigned long d, short p,unsigned long s):
	dim(d), perfomance(p), steps(s){

}
//----------------------------------------------------------------------------
pursuerType::~pursuerType() {
	for_each(centres.begin(), centres.end(), DeleteObj());
	//for_each(funnel.begin(), funnel.end(), DeleteObj());
}

// ----------------------------------- >> -------------------------------------//
istream& __fastcall operator >> (istream& in_data, pursuerType& C) {
	char c;
	char long_buf[2048];
	Vector *xs;

		in_data.get(c);
		while ((c != '['))  //ищем начало строчки с догоняющим
			in_data.get(c);
		while ((c != '='))  //ищем начало функции терм. мн-ва помехи
			in_data.get(c);
		in_data.get(long_buf, 256, ',');
		C.func_m.assign(long_buf);
		cout<<C.func_m<<endl;
		in_data.get(c);

		while ((c != '='))  //ищем начало функции мн-ва управления помехи
			in_data.get(c);
		in_data.get(long_buf, 256, ',');
		C.func_v.assign(long_buf);
		cout<<C.func_v<<endl;
		in_data.get(c);

		do{
			while ((c != '[') &&(c != ';'))
				in_data.get(c);
			if (c == ';')
				break;
			in_data.putback(c);
			xs=new Vector(C.dim);
			in_data.get(c);
			in_data >> *xs;
			C.centres.push_back(xs);
			in_data.get(c);
			cout << *xs;
		}while ((c != ';'));
	do{
		in_data.get(c);
	}while ((c != '\n')&&(!in_data.eof()));

	return in_data;
}


