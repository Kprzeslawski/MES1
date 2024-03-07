#pragma once
#include <math.h>
class Gauss
{
public:

	double** base_weights = new double* [4];
	double** base_nodes = new double* [4];

	Gauss() {
		base_weights[0] = new double[2] {1.,1.};
		base_weights[1] = new double[3] {5./9.,8./9.,5./9.};
		base_weights[2] = new double[4] {0.5-sqrt(30)/36, 0.5 + sqrt(30) / 36, 0.5 + sqrt(30) / 36, 0.5 - sqrt(30) / 36};
		base_weights[3] = new double[5] {(322.-13*sqrt(70.))/900., (322. + 13 * sqrt(70.)) / 900.,128./225., (322. + 13 * sqrt(70.)) / 900., (322. - 13 * sqrt(70.)) / 900.};

		base_nodes[0] = new double[2] {-1./sqrt(3), 1. / sqrt(3)};
		base_nodes[1] = new double[3] { -sqrt(0.6), 0, sqrt(0.6)};
		base_nodes[2] = new double[4] {-sqrt(3./7. + 2./7.*sqrt(1.2)), -sqrt(3. / 7. - 2. / 7. * sqrt(1.2)), sqrt(3. / 7. - 2. / 7. * sqrt(1.2)), sqrt(3. / 7. + 2. / 7. * sqrt(1.2))};
		base_nodes[3] = new double[5] {-1./3.* sqrt(5.+ 2.* sqrt(10./7.)), -1. / 3. * sqrt(5. - 2. * sqrt(10. / 7.)), 0, 1. / 3. * sqrt(5. - 2. * sqrt(10. / 7.)), 1. / 3. * sqrt(5. + 2. * sqrt(10. / 7.))};
	};
	
	double calculate(double (* f)(double, double), double* ksi, double* fi, double* waga, int il_el, int dim) {
		double res = 0;

		for (int i = 0; i < il_el; i++) {
			if (dim == 1) {
				res += f(ksi[i], 0) * waga[i];
			}
			if (dim == 2)
				for (int i2 = 0; i2 < il_el; i2++) {
					res += f(ksi[i], fi[i2]) * waga[i] * waga[i2];
				}
		}
		return res;
	}
};

