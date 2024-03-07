#include "Element.h"
#ifndef PunktyCa
#define PunktyCa 2
#endif // PunktyCa

Element::Element(int i1, int i2, int i3, int i4,int PC) {
	ID[0] = i1-1;
	ID[1] = i2-1;
	ID[2] = i3-1;
	ID[3] = i4-1;
	detJ = new double[PC*PC] {0.};//PC*PC
}

double** Element::CalculateH(
	std::vector<std::vector<double>> mat_dksi,
	std::vector<std::vector<double>> mat_deta,
	double** locmat,
	int* BC,
	int PCAL,
	double COND,
	double ALF) {
	Gauss G = Gauss();
	
	double*** PC = new double** [PCAL*PCAL];
	double*** dNdxy = new double** [PCAL * PCAL];
	double*** Hi = new double** [PCAL * PCAL];
	double** H = new double* [4];

	for (int i = 0; i < 4; i++)H[i] = new double[4] {0.};

	for (int i = 0; i < PCAL*PCAL; i++) {
		PC[i] = new double*[2];
		
		PC[i][0] = new double [2] {
			locmat[1][0] * mat_deta[0][i]+ locmat[1][1] * mat_deta[1][i]+ locmat[1][2] * mat_deta[2][i]+ locmat[1][3] * mat_deta[3][i],
			-(locmat[1][0] * mat_dksi[0][i] + locmat[1][1] * mat_dksi[1][i] + locmat[1][2] * mat_dksi[2][i] + locmat[1][3] * mat_dksi[3][i]) };
		PC[i][1] = new double [2] {
			-(locmat[0][0] * mat_deta[0][i] + locmat[0][1] * mat_deta[1][i] + locmat[0][2] * mat_deta[2][i] + locmat[0][3] * mat_deta[3][i]),
				locmat[0][0] * mat_dksi[0][i] + locmat[0][1] * mat_dksi[1][i] + locmat[0][2] * mat_dksi[2][i] + locmat[0][3] * mat_dksi[3][i] };
			//std::cout << PC[i][0][0] << " " << PC[i][0][1] << std::endl;
			//std::cout << PC[i][1][0] << " " << PC[i][1][1] << std::endl;
		
		double det =(PC[i][0][0] * PC[i][1][1] - PC[i][1][0] * PC[i][0][1]);
		detJ[i] = det;
		det = 1 / det;
		//std::cout << det << std::endl;
		for (int i2 = 0; i2 < 2; i2++)
			for (int i3 = 0; i3 < 2; i3++)
			PC[i][i2][i3] *= det;

		//std::cout << PC[i][0][0] << " " << PC[i][0][1] << std::endl;
		//std::cout << PC[i][1][0] << " " << PC[i][1][1] << std::endl;
		
		dNdxy[i] = new double* [4];
		for (int i2 = 0; i2 < 4; i2++) {
			dNdxy[i][i2] = new double [2] {
			PC[i][0][0]*mat_dksi[i2][i] + PC[i][0][1] * mat_deta[i2][i],
				PC[i][1][0] * mat_dksi[i2][i] + PC[i][1][1] * mat_deta[i2][i]};

			//std::cout << dNdxy[i][i2][1] << " ";
		}
		//std::cout << std::endl;
		
		det = 1. / det;
		//std::cout << det << std::endl;
		Hi[i] = new double* [4];
		//tu git
		for (int i2 = 0; i2 < 4; i2++) {
			Hi[i][i2] = new double[4];
			for (int i3 = 0; i3 < 4; i3++) {
				Hi[i][i2][i3] = COND * det* (dNdxy[i][i2][0] * dNdxy[i][i3][0] + dNdxy[i][i2][1] * dNdxy[i][i3][1])
					* G.base_weights[PCAL - 2][i / PCAL] * G.base_weights[PCAL - 2][i % PCAL];
				//std::cout << Hi[i][i2][i3] << " ";
				H[i2][i3] += Hi[i][i2][i3];
			}
			//std::cout << std::endl;
		}
		//std::cout << std::endl;
	}
#ifdef COMPILE
	for (int i = 0; i < 4; i++) {
		for (int i2 = 0; i2 < 4; i2++)std::cout << H[i][i2] << " ";
		std::cout << std::endl;
	}
#endif
	for (int i = 0; i < 4; i++) { // 4 scianki

		if (BC[i] == 0 || BC[i + 1] == 0)continue; // no bc

		double* bc = new double[4] {0.};
		double detj = 0.5*sqrt(pow(locmat[0][i]-locmat[0][(i+1)%4], 2) + pow(locmat[1][i]-locmat[1][(i+1)%4], 2));

		for (int i2 = 0; i2 < PCAL; i2++) {
			for (int i3 = 0; i3 < 4; i3++)bc[i3] = 0.;
			if(i==0){
				bc[0] = 0.5 * (1 - G.base_nodes[PCAL-2][i2]);
				bc[1] = 0.5 * (1 + G.base_nodes[PCAL-2][i2]);
			}									
			else if(i==1){						
				bc[1] = 0.5 * (1 - G.base_nodes[PCAL - 2][i2]);
				bc[2] = 0.5 * (1 + G.base_nodes[PCAL - 2][i2]);
			}									
			else if(i==2){						
				bc[2] = 0.5 * (1 - G.base_nodes[PCAL - 2][i2]);
				bc[3] = 0.5 * (1 + G.base_nodes[PCAL - 2][i2]);
			}									
			else if(i==3){						
				bc[3] = 0.5 * (1 - G.base_nodes[PCAL - 2][i2]);
				bc[0] = 0.5 * (1 + G.base_nodes[PCAL - 2][i2]);
			}

			for (int i3 = 0; i3 < 4; i3++)
				for (int i4 = 0; i4 < 4; i4++) {
					//std::cout << "val: " << bc[i3] * bc[i4] * detj << std::endl;
					H[i3][i4] += bc[i3] * bc[i4] * detj * ALF
						* G.base_weights[PCAL - 2][i2];
				}// *G.base_weights[2 - 2][] * G.base_weights[2 - 2][]; //*waga
		}

	}
#ifdef COMPILE
	for (int i = 0; i < 4; i++) {
		for (int i2 = 0; i2 < 4; i2++)std::cout << H[i][i2] << " ";
		std::cout << std::endl;
	}
#endif
	return H;
}

double** Element::CalculateC(
	int PC,
	double DENS,
	double SHEAT) {

	Gauss G = Gauss();

	double** CLok = new double* [4];
	for (int i = 0; i < 4; i++)
		CLok[i] = new double[4] {0.};
	double* N = new double[4] {0.};
#ifdef COMPILE
	std::cout << "vectN ";
#endif // COMPILE
	for (int i = 0; i < PC * PC; i++) {

		N[0] = 0.25 * (1. - G.base_nodes[PC-2][i%PC]) * (1. - G.base_nodes[PC-2][i/PC]);
		N[1] = 0.25 * (1. + G.base_nodes[PC-2][i%PC]) * (1. - G.base_nodes[PC-2][i/PC]);
		N[2] = 0.25 * (1. + G.base_nodes[PC-2][i%PC]) * (1. + G.base_nodes[PC-2][i/PC]);
		N[3] = 0.25 * (1. - G.base_nodes[PC-2][i%PC]) * (1. + G.base_nodes[PC-2][i/PC]);
#ifdef COMPILE
		for (int i = 0; i < 4; i++)std::cout << N[i] << " ";
		std::cout << std::endl;
#endif // COMPILE

		
		for (int i3 = 0; i3 < 4; i3++)
			for (int i4 = 0; i4 < 4; i4++)
				CLok[i3][i4] += N[i3] * N[i4] 
				* G.base_weights[PC - 2][i / PC] * G.base_weights[PC - 2][i % PC]
				*detJ[i] * DENS * SHEAT; 
	}
#ifdef COMPILE
	std::cout << "\nCLOK";

	for (int i = 0; i < 4; i++) {
		for (int i2 = 0; i2 < 4; i2++)std::cout << CLok[i][i2] << " ";
		std::cout << std::endl;
	}
#endif // COMPILE
	return CLok;
}

double* Element::CalculateVectP(
	double** locmat,
	int* BC,
	int PC,
	double AMB_TEMP,
	double ALF
) {

	Gauss G = Gauss();
	double* vectP = new double[4] {0.};
	for (int i = 0; i < 4; i++) { // 4 scianki

		if (BC[i] == 0 || BC[i + 1] == 0)continue; // no bc

		double* bc = new double[4] {0.};
		double detj = 0.5 * sqrt(pow(locmat[0][i] - locmat[0][(i + 1) % 4], 2) + pow(locmat[1][i] - locmat[1][(i + 1) % 4], 2));

		for (int i2 = 0; i2 < PC; i2++) {
			if (i == 0) {
				bc[0] = 0.5 * (1 - G.base_nodes[PC - 2][i2]);
				bc[1] = 0.5 * (1 + G.base_nodes[PC - 2][i2]);
			}									
			else if (i == 1) {					
				bc[1] = 0.5 * (1 - G.base_nodes[PC - 2][i2]);
				bc[2] = 0.5 * (1 + G.base_nodes[PC - 2][i2]);
			}									
			else if (i == 2) {					
				bc[2] = 0.5 * (1 - G.base_nodes[PC - 2][i2]);
				bc[3] = 0.5 * (1 + G.base_nodes[PC - 2][i2]);
			}									
			else if (i == 3) {					
				bc[3] = 0.5 * (1 - G.base_nodes[PC - 2][i2]);
				bc[0] = 0.5 * (1 + G.base_nodes[PC - 2][i2]);
			}
			/*
			bc[0] = 0.25 * (1. - mat_bc[i][i2][0]) * (1. - mat_bc[i][i2][1]);
			bc[1] = 0.25 * (1. + mat_bc[i][i2][0]) * (1. - mat_bc[i][i2][1]);
			bc[2] = 0.25 * (1. + mat_bc[i][i2][0]) * (1. + mat_bc[i][i2][1]);
			bc[3] = 0.25 * (1. - mat_bc[i][i2][0]) * (1. + mat_bc[i][i2][1]);
			
			//for (int i3 = 0; i3 < 4; i3++)std::cout << bc[i] << " ";
			//std::cout << std::endl;
			*/
			for (int p = 0; p < 4; p++)vectP[p] += bc[p] * AMB_TEMP * ALF * detj * G.base_weights[PC - 2][i2];
		}

	}
#ifdef COMPILE
	for (int i = 0; i < 4; i++)std::cout << vectP[i] << " ";
	std::cout << std::endl;
#endif // COMPILE

	return vectP;
}

int* Element::getNodeIds() {
	return new int[4] {ID[0], ID[1], ID[2], ID[3]};
}