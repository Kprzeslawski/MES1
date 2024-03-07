#include "Grid.h"

Grid::Grid(int nn, int ne,Element_Uni eu) : nn(nn), ne(ne),EU(eu) {

	ND.resize(nn);
	ELEM.resize(ne);
	globalH = new double* [nn];
	globalC = new double* [nn];
	globalP = new double [nn] {0.};
	//EQ1 = new double* [nn];

	for (int i = 0; i < nn; i++) {
		globalH[i] = new double [nn] {0.};
		globalC[i] = new double [nn] {0.};
		//EQ1[i] = new double [nn] {0.};
	}

}

void Grid::setupMat(double DENS,double S_HEAT,double AMB_TEMP,double ALF,double COND) {
	for (int i = 0; i < ne; i++) {
#ifdef COMPILE
		std::cout << i << std::endl;
#endif
		double** locmat = new double* [2];
		locmat[0] = new double[4] {0.};
		locmat[1] = new double[4] {0.};
		for (int i2 = 0; i2 < 4; i2++) {
			locmat[0][i2] = ND[ELEM[i].ID[i2]].get_x();
			locmat[1][i2] = ND[ELEM[i].ID[i2]].get_y();
#ifdef COMPILE
			std::cout << locmat[0][i2] << " " << locmat[1][i2] <<std::endl;
#endif
		}
		int* BC = new int[5] {
			ND[ELEM[i].ID[0]].bc(),
			ND[ELEM[i].ID[1]].bc(),
			ND[ELEM[i].ID[2]].bc(),
			ND[ELEM[i].ID[3]].bc(),	

			ND[ELEM[i].ID[0]].bc()
		};
		double** Hlok = ELEM[i].CalculateH(EU.mat_dksi,EU.mat_ddeta,locmat,BC,EU.P,COND,ALF);
		double** Clok = ELEM[i].CalculateC(EU.P, DENS, S_HEAT);
		double* Plok = ELEM[i].CalculateVectP(locmat, BC,EU.P,AMB_TEMP,ALF); 

		for (int i2 = 0; i2 < 4; i2++)for (int i3 = 0; i3 < 4; i3++) {
			globalH[ELEM[i].ID[i2]][ELEM[i].ID[i3]] += Hlok[i2][i3];
			globalC[ELEM[i].ID[i2]][ELEM[i].ID[i3]] += Clok[i2][i3];
		}
		for (int i2 = 0; i2 < 4; i2++)
			globalP[ELEM[i].ID[i2]] += Plok[i2];
	}
#ifdef COMPILE
	for (int i = 0; i < nn; i++) {
		for (int i2 = 0; i2 < nn; i2++)
			std::cout << globalH[i][i2] << " ";
		std::cout << std::endl;
	}std::cout << std::endl;
	for (int i = 0; i < nn; i++) {
		for (int i2 = 0; i2 < nn; i2++)
			std::cout << globalC[i][i2] << " ";
		std::cout << std::endl;
	}std::cout << std::endl;
	for (int i = 0; i < nn; i++) {
			std::cout << globalP[i] << " ";
	}std::cout << std::endl;
	std::cout << std::endl;
#endif
}
/*
void Grid::setupEqSolver(
	double dt
) {

	double dt_mn = 1. / dt;
	double** HCdt_min_1 = new double* [nn];

	for (int i = 0; i < nn; i++) {
		HCdt_min_1[i] = new double[nn];
		for (int i2 = 0; i2 < nn; i2++) {
			globalC[i][i2] *= dt_mn;
			HCdt_min_1[i][i2] = globalH[i][i2] + globalC[i][i2];
		}
	}
	//liczymy HC ^ -1
	double detHC = calculateDET(HCdt_min_1, nn);

	std::cout << "\nPOLICZONE\n";

	HCdt_min_1 = calculateDopMat(HCdt_min_1, nn);
	
	double detHC_mn = 1. / detHC;
	for (int i = 0; i < nn; i++)
		for (int i2 = 0; i2 < nn; i2++) {
			HCdt_min_1[i][i2] *= detHC_mn;
			EQ1[i2][i] = HCdt_min_1[i][i2];
		}
}
double Grid::calculateDET(double** mat, int dim) {
	if (dim == 1) return mat[0][0];
	if (dim == 2) return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
	if (dim == 3) return mat[0][0] * mat[1][1] * mat[2][2]
		+ mat[1][0] * mat[2][1] * mat[0][2]
		+ mat[2][0] * mat[0][1] * mat[1][2]
		- mat[2][0] * mat[1][1] * mat[0][2]
		- mat[1][0] * mat[0][1] * mat[2][2]
		- mat[0][0] * mat[2][1] * mat[1][2];

	int res = 0;
	double** matT = new double* [dim - 1];
	for (int i = 0; i < dim-1; i++)matT[i] = new double[dim - 1];

	for (int i = 0; i < dim; i++) {
		if (mat[0][i] == 0)continue;
		
		for (int i2 = 0; i2 < dim - 1; i2++)for (int i3 = 0; i3 < dim - 1; i3++) {
			if (i > i2)matT[i3][i2] = mat[i3+1][i2];
			else matT[i3][i2] = mat[i3 + 1][i2+1];
		}
		res += mat[0][i] * pow(-1,i) * calculateDET(matT, dim - 1);
	}
	delete[] matT;
	return res;
}

double** Grid::calculateDopMat(double** mat, int dim) {
	double** dopmat = new double* [dim];
	for (int i = 0; i < dim; i++)dopmat[i] = new double[dim];

	double** matT = new double* [dim - 1];
	for (int i = 0; i < dim - 1; i++)matT[i] = new double[dim - 1];

	for (int i = 0; i < dim; i++) {
		for (int i2 = 0; i2 < dim; i2++) {
			for (int i3 = 0; i3 < dim - 1; i3++)
				for (int i4 = 0; i4 < dim - 1; i4++) {
					int at=i3, at2=i4;
					if (i < i3)at++;
					if (i2 < i4)at2++;
					matT[i3][i4] = mat[at][at2];

				}
			dopmat[i][i2] = pow(-1, i2+i) * calculateDET(matT, dim - 1);
		}
	}
	return dopmat;
}

double* Grid::tempVectCalcByStep(double* t0) {
	double* t1 = new double[nn] {0.};
	double* matR = new double[nn] {0.};

	for (int i = 0; i < nn; i++) {
		for (int i2 = 0; i2 < nn; i2++)
			matR[i] += globalC[i][i2] * t0[i2];
		matR[i] -= globalP[i];
	}

	for (int i = 0; i < nn; i++)
		for (int i2 = 0; i2 < nn; i2++)
			t1[i] += EQ1[i][i2] * matR[i2];

	return t1;
}
*/
void Grid::divH(double dT) {
	double dt_mn = 1. / dT;

	for (int i = 0; i < nn; i++) {
		for (int i2 = 0; i2 < nn; i2++) {
			globalC[i][i2] *= dt_mn;
#ifdef COMPILE
			std::cout << globalC[i][i2] << " ";
#endif
		}
#ifdef COMPILE
		std::cout << "\n";
#endif
	}

	
}

double* Grid::tempVectCalcByStepUsingGauss(double* t0) {
	double** matToSolve = new double* [nn + 1];
	for (int i = 0; i < nn + 1; i++)
		matToSolve[i] = new double[nn] {0.};
	double* tRES = new double[nn] {0.};

	for (int i = 0; i < nn; i++)for (int i2 = 0; i2 < nn; i2++)
		matToSolve[i][i2] =  globalH[i][i2] + globalC[i][i2];

	for (int i = 0; i < nn; i++) {
		for (int i2 = 0; i2 < nn; i2++)
			matToSolve[nn][i] += globalC[i2][i] * t0[i2];
		matToSolve[nn][i] += globalP[i];
	}
#ifdef COMPILE
	std::cout << "Macierz " << std::endl;

	for (int i2 = 0; i2 < nn; i2++) {
		for (int i3 = 0; i3 < nn + 1; i3++)
			std::cout << matToSolve[i3][i2] << " ";
		std::cout << std::endl;
	}
#endif
	double m;
	for (int i5 = 0; i5 < nn; i5++)
		for (int i4 = i5 + 1; i4 < nn; i4++) {

			if (matToSolve[i5][i4]) {
				m = matToSolve[i5][i4] / matToSolve[i5][i5];

				for (int i6 = i5; i6 < nn+1; i6++)
					matToSolve[i6][i4] -= m * matToSolve[i6][i5];
			}
		}
#ifdef COMPILE

	std::cout << "Macierz po obliczeniu post. prostym" << std::endl;

	for (int i2 = 0; i2 < nn; i2++) {
		for (int i3 = 0; i3 < nn+1; i3++)
			std::cout << matToSolve[i3][i2] << " ";
		std::cout << std::endl;
	}
	
	std::cout << "Rozwiazany uklad" << std::endl;
#endif
	for (int j = nn-1; j >= 0; j--) {
		tRES[j] = matToSolve[nn][j];
		for (int j2 = j + 1; j2 < nn; j2++)
			tRES[j] -= matToSolve[j2][j] * tRES[j2];


		tRES[j] /= matToSolve[j][j];
		//std::cout << tRES[j] << " ";

	}
	//std::cout << "\n";
	return tRES;
	//solve by gauss 

}