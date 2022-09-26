// wobble.cpp
// Main file for wobble simulator.

#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include "toml.hpp"

typedef Eigen::Matrix<std::complex<double>, 2, 1> Vector2c;
typedef Eigen::Matrix<std::complex<double>, 2, 2> Matrix2c;

const double PI = std::acos(-1);

Vector2c bloch_state(double theta, double phi) {
	// Make a state vector on the Bloch sphere.
	Vector2c psi;
	psi(0) = std::cos(theta / 2.);
	psi(1) = std::polar(std::sin(theta / 2.), phi);
	return psi;
}

std::vector<Matrix2c> jump_sim (
	Vector2c psi0, Matrix2c ham, Matrix2c jump,
	double t_int, int steps, int runs, std::mt19937 gen,
	std::uniform_real_distribution<> dis
){
	std::vector<Matrix2c> rho;
	Matrix2c h_eff;
	Vector2c psi;
	double pjump, rnd, norm_psi;

	// Assign a matrix of zeros to each element of rho.
	for (int i = 0; i <= steps; i++) {
		rho.push_back(Matrix2c::Zero());
	}
	// Calculate effective Hamiltonian.
	h_eff = ham - std::complex(0., 1.) * jump.adjoint() * jump;
	// Simulate each trajectory, adding the state vectors
	// to the elements of rho.
	for (int i = 0; i < runs; i++) {
		// Set initial wave function
		// And starting density matrix.
		psi = psi0;
		rho[0] += psi * psi.adjoint();
		// Step the state forward in time. 
		for (int j = 1; j < steps + 1; j++) {
			// Generate jump probability.
			pjump = t_int / steps * std::abs(psi.conjugate().dot(
					jump.adjoint() * jump * psi));
			// Evolve under either hamiltonian or
			// jump operator. Generate a random number
			// and compare with pjump.
			rnd = dis(gen);
			if (rnd <= pjump) {
				// Collapse the state.
				psi = jump * psi;
			} else {
				// Evolve under h_eff.
				psi -= std::complex(0., t_int / double(steps)) * h_eff * psi;
			}
			// Normalize the state vector.
			psi = psi / std::sqrt(std::abs(psi.conjugate().dot(psi)));
			// Push to the accumulator.
			rho[j] += psi * psi.adjoint();
		}
	}
	// Divide each element of rho by the number of runs, so that 
	// it is an average.
	for (int i = 0; i <= steps; i++) {
		rho[i] = rho[i] / runs;
	}
	return rho;
}

int main() {
	Vector2c ket_g, ket_e;
	Matrix2c ham, jump;
	double omega0, gamma, t_int;
	double prob_g;
	int steps, runs;
	std::vector<Matrix2c> rho;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0., 1.);
	std::ofstream outfile;
	toml::table tbl;

	ket_e = bloch_state(0., 0.);
	ket_g = bloch_state(PI, 0.);
	// Read constants from file.
	tbl = toml::parse_file("input.toml");
	omega0 = tbl["atom"]["omega0"].value_exact <double> ().value();
	gamma = tbl["atom"]["gamma"].value_exact <double> ().value();
	t_int = tbl["integration"]["t_int"].value_exact <double> ().value();
	steps = tbl["integration"]["steps"].value_exact <int64_t> ().value();
	runs = tbl["integration"]["runs"].value_exact <int64_t> ().value();
	// Define operators.
	ham = omega0 / 2. 
		* (ket_e * ket_e.adjoint()
		  - ket_g * ket_g.adjoint());
	jump = sqrt(gamma) * ket_g * ket_e.adjoint();
	// Run the simulation.
	rho = jump_sim(ket_e, ham, jump, t_int, steps, runs, gen, dis);
	// Print data to file.
	outfile.open("out.csv");
	outfile << "i,rho_e,rho_g\n";
	for (int i = 0; i <= steps; i++) {
		outfile << i << ',' << rho[i](0,0).real() << ',' << rho[i](1,1).real() << std::endl;
	}
	outfile.close();
	return 0;
}
