#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#pragma warning(disable : 4996)

#include "GlobalData.h"
#include "Grid.h"
#include "Gauss.h"
#include "ElementBaz.h"

//#define COMPILE 

#define PunktyCa 3
using namespace std;


double f1dd(double ksi, double deta) {
	return -0.25 * (1. - ksi);
}
double f2dd(double ksi, double deta) {
	return 0.25 * (1. - ksi);
}
double f3dd(double ksi, double deta) {
	return 0.25 * (1. + ksi);
}
double f4dd(double ksi, double deta) {
	return -0.25 * (1. + ksi);
}

double f1dk(double ksi, double deta) {
	return -0.25 * (1. - deta);
}
double f2dk(double ksi, double deta) {
	return -0.25 * (1. + deta);
}
double f3dk(double ksi, double deta) {
	return 0.25 * (1. + deta);
}
double f4dk(double ksi, double deta) {
	return 0.25 * (1. - deta);
}

int main() {
	
	double (*dd[4])(double, double) = { f1dd,f2dd,f3dd,f4dd };
	double (*dk[4])(double, double) = { f1dk,f2dk,f3dk,f4dk };

	Element_Uni e = Element_Uni(4, PunktyCa, dd, dk);
	e.setup_matrix();

	fstream file;
	//file.open("Test1_4_4.txt");
	file.open("Test2_4_4_mg.txt");
	//file.open("Test3_31_31.txt");
	string line;
	int gda[10];
	
	for (int i = 0; i < 10; i++)
	{
		getline(file, line);
		gda[i] = stoi(line.substr(line.find_last_of(' '), line.back()- line.find_last_of(' ')));
	}

	GlobalData GD(gda[0], gda[1], gda[2], gda[3], gda[4], gda[5], gda[6], gda[7], gda[8], gda[9]);
	Grid GR(GD.Nodes_number, GD.Elements_number, e);

	getline(file, line); // *Nody
	float x, y;
	for (int i = 0; i < GR.nn; i++)
	{
		file >> line;
		file >> line; //x
		x = stod(line); 
		file >> line; //y
		y = stod(line);
		//cout << x << " " << y << endl;
		GR.ND[i] = Node(x, y, GD.InitialTemp);
	}

	getline(file, line);
	getline(file, line); // *Element

	int t1, t2, t3, t4;
	
	for (int i = 0; i < GR.ne; i++)
	{
		file >> line;
		file >> line; //1
		t1 = stoi(line);
		file >> line; //2
		t2 = stoi(line);
		file >> line; //3
		t3 = stoi(line);
		file >> line; //4
		t4 = stoi(line);

		GR.ELEM[i] = Element(t1,t2,t3,t4,PunktyCa);
	}

	getline(file, line);
	getline(file, line); // *BC

	while (file >> line) {
		GR.ND[stoi(line)-1].set_bc();
	}

	GR.setupMat(GD.Density,GD.SpecificHeat,GD.Tot,GD.Alfa,GD.Conductivity);

	//GR.divH(1);
	GR.divH(GD.SimulationStepTime);

	double* vecT = new double[GR.nn];
	for (int i = 0; i < GR.nn; i++) {
		vecT[i] = GD.InitialTemp;
	}
	for (int i = 0; i < GD.SimulationTime / GD.SimulationStepTime; i++) {
		vecT = GR.tempVectCalcByStepUsingGauss(vecT);
		double max = vecT[0], min = vecT[0];
		for (int i2 = 1; i2 < GR.nn;i2++) {
			if (vecT[i2] > max)max = vecT[i2];
			if (vecT[i2] < min)min = vecT[i2];
		}
		std::cout << "Time " << (i + 1) * GD.SimulationStepTime 
			<< "s: MinT: " << min << " MaxT: " << max << "\n";
	}
	
	/*
	GR.setupEqSolver(50);

	double* vecT = new double[GR.nn] {100.};

	for (int i = 0; i < 10; i++) {
	
		vecT = GR.tempVectCalcByStep(vecT);
		for (int i2 = 0; i2 < GR.nn; i2++) {
			std::cout << vecT[i] << " ";
		}
		std::cout << std::endl;
	
	}
	*/
	//std::cout << G.calculate(f1, G.base_nodes[0], nullptr, G.base_weights[0], 2, 1) << endl;
	//std::cout << G.calculate(f1, G.base_nodes[1], nullptr, G.base_weights[1], 3, 1) << endl;

	//std::cout << G.calculate(f2, G.base_nodes[0], G.base_nodes[0], G.base_weights[0], 2, 2) << endl;
	//std::cout << G.calculate(f2, G.base_nodes[1], G.base_nodes[1], G.base_weights[1], 3, 2) << endl;
	
	//printf("%d", GD.SimulationTime);
	/*
	Element el = Element(0, 1, 2, 3);

	double** locmat = new double* [2];
	locmat[0] = new double[4] {0.};
	locmat[1] = new double[4] {0.};

	locmat[0][1] = 0.025;
	locmat[0][2] = 0.025;
	locmat[1][2] = 0.025;
	locmat[1][3] = 0.025;
	*/
	//double** Hlok = el.Calculate(e.mat_dksi, e.mat_ddeta, locmat);
	
}
	

