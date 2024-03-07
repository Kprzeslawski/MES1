#pragma once
#include "ElementBaz.h"
//#define COMPILE
class Element
{
public:
	int ID[4];
	double* detJ;
	Element(int i1, int i2, int i3, int i4,int PC=2);
	Element() {};

	double** CalculateH(std::vector<std::vector<double>> mat_dksi,
		std::vector<std::vector<double>> mat_deta,
		double** locmat,
		int* BC,
		int PCAL,
		double COND,
		double ALF);
	double** CalculateC(
		int PC,
		double DENS,
		double SHEAT);
	double* CalculateVectP(
		double** locmat,
		int* BC,
		int PC,
		double AMB_TEMP,
		double ALF);

	int* getNodeIds();
};

