#pragma once

#include <vector>
#include <iostream>
#include "Gauss.h"

typedef double (*fn)(double, double);

class Element_Uni
{
public:
	std::vector<std::vector<double>> mat_dksi;
	std::vector<std::vector<double>> mat_ddeta;
	std::vector<std::vector<std::vector<double>>> mat_bc;

	int N, P;

	double(*dd[4])(double, double);
	double(*dk[4])(double, double);

	Element_Uni(int Ns, int Points,
		double(*tdd[4])(double, double),
		double(*tdk[4])(double, double));
	
	void setup_matrix();
	void print();
};