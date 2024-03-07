#pragma once
#include "Element.h"
#include "Node.h"
#include <vector>
//#define COMPILE 
class Grid
{
public:
	int nn, ne;
	std::vector<Node> ND;
	std::vector<Element> ELEM;
	Element_Uni EU;
	double** globalH;
	double** globalC;
	double* globalP;

	//double** EQ1; // (H+C/dt)^-1

	Grid(int nn, int ne, Element_Uni eu);
	void setupMat(double DENS, double S_HEAT, double AMB_TEMP, double ALF,double COND);
	void divH(double dT);
	double* tempVectCalcByStepUsingGauss(double* t0);
	
	/*
	void setupEqSolver(
		double dt
	);
	double* tempVectCalcByStep(double* t0);
	double calculateDET(double** mat,int dim);
	double** calculateDopMat(double** mat, int dim);*/
};

