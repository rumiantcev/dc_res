// ---------------------------------------------------------------------------
#ifndef generalH
#define generalH

#include <string>
#include <sstream>
#include <map>
#include <list>
#include <vector>
#include <math.h>
//#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <locale.h>
#include <windows.h>
#include<time.h>
#include <io.h>
#include<iterator>
#include <stdlib>
#include <sstream>
#include <algorithm>
//#include <system.hpp>




typedef double LDouble;
const std::string states[4] = {
	"Выбор метода", "Ввод параметров", "Решение", "Анализ результатов"};

typedef enum OpType {
	opMin = 0, opMax = 1, opZero = 2
} OpType;

typedef enum ZeroAware {
	nZAware = false, ZAware = true
} ZeroAware;

typedef enum FindBase {
	nFind = false, yFind = true
} FindBase;

typedef enum FindPath {
	nPath = false, yPath = true
} FindPath;

typedef enum optSearch {
	optNone = 0, optAnnealing = 1, optGradient = 2
} optSearch;

typedef std::vector<LDouble>alphType;
typedef alphType::value_type aplhVal;
typedef std::map<long, bool>seekType;
typedef std::list<long>pathType;

inline int signof(LDouble a) { return (a == 0.0) ? 0 : (a<0 ? -1 : 1); }

//глобальные переменные используемые в методе отжига
static LDouble _extr_e_param  = -0.8;
static LDouble _extr_t0_param  = 1.0;
static LDouble _extr_tmin_param  = 0.0001;
static LDouble _extr_e_val = exp(_extr_e_param);

const LDouble ldZeroDf = 0.000000001f; //допустимая погрешность при сравнении с нулём для чисел с плавающей точкой

//очистка ссылок в векторах объектов
struct DeleteObj {
	template<typename T>
	void operator()(const T* ptr) const {
		if (ptr != NULL)
			delete ptr;
	}
};

//---- преобразования в std::string
inline std::string intToStr(int i){
	std::string str;
	std::stringstream oss;
	oss << i;
	oss >> str;
	return str;
}

inline std::string ldToStr(LDouble ld){
	std::string str;
	std::stringstream oss;
	oss << ld;
	oss >> str;
	return str;
}


#endif
