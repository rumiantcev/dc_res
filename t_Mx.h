/*
* Copyright (c) 2001-2017, Alexey Rumyantsev.
* e-mail: rumiantcev@yandex.ru
* All rights reserved.
*
*/
//---------------------------------------------------------------------------

#ifndef t_MxH
#define t_MxH
//---------------------------------------------------------------------------
#include<assert.h>

template <typename T> class t_Mx {
protected:
	long * v;
	long m,n;

	class pr_Mx {
	public:
		pr_Mx(T * p, long rows):p_n(rows),p_v(p){
		};
		~pr_Mx(){
		};

		T operator [] (long row) const{
			assert ((long)row < (long)p_n);
			return p_v[row];
		};

		T& operator [] (long row){
			assert ((long)row < (long)p_n);
			return p_v[row];
		};
	private:
		T * p_v;
		long p_n;
	};

public:
	t_Mx(long rows, long cols):n(cols),m(rows){
		v = new T [n * m];
	};
	~t_Mx(){
		delete [] v;
	};

	pr_Mx operator [] (long col){
		assert ((long)col < (long)m);
		return pr_Mx(v + col * n, m);
	};
};

#endif
