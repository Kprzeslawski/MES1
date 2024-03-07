#include "ElementBaz.h"

Element_Uni::Element_Uni(int Ns, int Points,
	double(*tdd[4])(double, double),
	double(*tdk[4])(double, double)
	) :N(Ns), P(Points) {

	mat_dksi.resize(Ns);
	mat_ddeta.resize(Ns);

	for (int i = 0; i < Ns; i++) {
		mat_dksi[i].resize(Points * Points);
		mat_ddeta[i].resize(Points * Points);
		dd[i] = tdd[i];
		dk[i] = tdk[i];
	}	
}

void Element_Uni::setup_matrix() {
	Gauss G = Gauss();
	for (int i = 0; i < N; i++) {
		for (int i2 = 0; i2 < P * P; i2++) {
			mat_dksi[i][i2] =  dd[i](G.base_nodes[P - 2][i2 / P], 0);
			mat_ddeta[i][i2] =   dk[i](0, G.base_nodes[P - 2][i2 % P]);
		}
	}
}
/*
void Element_Uni::setup_BC_matrix() {
	Gauss G = Gauss();
	for (int i = 0; i < 2; i++) {
		mat_bc[0][i][0] = -1.;
		mat_bc[0][i][1] = G.base_nodes[0][i];

		mat_bc[1][i][0] = G.base_nodes[0][i];
		mat_bc[1][i][1] = 1.;

		mat_bc[2][i][0] = 1.;
		mat_bc[2][i][1] = G.base_nodes[0][i];

		mat_bc[3][i][0] = G.base_nodes[0][i];
		mat_bc[3][i][1] = - 1.;
	}

	for (int i = 0; i < 4; i++) {
		for (int i2 = 0; i2 < 2; i2++) {
			for (int i3 = 0; i3 < 2; i3++)
				std::cout << mat_bc[i][i2][i3] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}*/

void Element_Uni::print() {
	std::cout << std::endl;
	for (int i = 0; i < P*P; i++) {
		for (int i2 = 0; i2 < N; i2++) {
			std::cout << mat_dksi[i2][i] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	for (int i = 0; i < P*P; i++) {
		for (int i2 = 0; i2 < N; i2++) {
			std::cout << mat_ddeta[i2][i] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

}


