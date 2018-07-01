// Includes

#include "stdafx.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

// Namespace

using namespace std;

// Definitions

using vec = vector<double>;
using matrix = vector<vec>;

// -- == ~~ EXERCISE SHEET 7 ~~ == -------------------------------------------------------------------------------------

void JacobiMethod(matrix &A, int dim) {

	/// Jacobi Method for determining eigenvalues and eigenvectors of a Matrix A of dimension 'dim'

	// Declaration and Initialization of Variables

	double S; // Deviation from diagonal matrix
	const double epsilon = 1.0e-20; // bound for S 

	double Theta;
	double tan_phi;

	double sin_phi;
	double cos_phi;

	double Tau;

	double A_pp;
	double A_qq;

	vec A_rp(dim);
	vec A_rq(dim);

	vec V_rp(dim);
	vec V_rq(dim);

	// == Initialize Matrix V of eigenvectors ==

	matrix V(dim); // matrix of eigenvectors
	for (int i = 0; i < dim; i++) {
		V[i].resize(dim, 0);
	}

	for (int row = 0; row < dim; row++) {

		V[row][row] = 1.;

	}

	// == Jacobi Method ==

	while (1) {

		S = 0.;

		// Calculate Deviation from Diagonal Matrix

		for (int row = 0; row < dim; row++) {

			for (int col = 0; col < row; col++) {

				S += pow(A[row][col], 2.);

			}

		}

		S *= 2;

		if (S <= epsilon) { break; }// Matrix is diagonal enough

		// Sweep once over all Off-Diagonal Elements

		for (int p = 0; p < dim; p++) {

			for (int q = 0; q < dim; q++) {

				if (fabs(A[p][q] < epsilon / dim * dim)) continue; // Avoid devision by almost 0

				// == Determine Theta ==

				Theta = 0.5*(A[q][q] - A[p][p]) / A[p][q];

				// == Determine tan_phi ==

				tan_phi = 1. / (fabs(Theta) + sqrt(pow(Theta, 2.) + 1.));
				if (Theta < 0) tan_phi = -tan_phi;

				// == Determine cos_phi and sin_phi ==

				cos_phi = 1. / sqrt(pow(tan_phi, 2.) + 1.);
				sin_phi = tan_phi * cos_phi;

				// == Determine Tau ==

				Tau = sin_phi / (1. + cos_phi);


				// --- Jacobi Rotation ---

				// == Matrix A ==

				// p q Elements

				A_pp = A[p][p] - tan_phi * A[p][q];
				A_qq = A[q][q] + tan_phi * A[p][q];

				// Non p q Elements

				for (int k = 0; k < dim; k++) {

					if (k != p && k != q) {

						A_rp[k] = A[k][p] - sin_phi * (A[k][q] + Tau * A[k][p]);
						A_rq[k] = A[k][q] + sin_phi * (A[k][p] - Tau * A[k][q]);

					}

				}

				// Update Matrix

				A[p][q] = 0.;
				A[q][p] = 0.;

				A[p][p] = A_pp;
				A[q][q] = A_qq;

				for (int k = 0; k < dim; k++) {

					if (k != p && k != q) {

						A[k][p] = A_rp[k];
						A[p][k] = A_rp[k];

						A[k][q] = A_rq[k];
						A[q][k] = A_rq[k];

					}

				}

				// == Matrix V ==

				for (int k = 0; k < dim; k++) {

					V_rp[k] = V[k][p] - sin_phi * (V[k][q] + Tau * V[k][p]);
					V_rq[k] = V[k][q] + sin_phi * (V[k][p] - Tau * V[k][q]);

				}

				for (int k = 0; k < dim; k++) {

					V[k][p] = V_rp[k];
					V[k][q] = V_rq[k];

				}

			}

		}

	}

	// == Output ==

	for (int row = 0; row < dim; row++) {

		cout << "| ";

		for (int col = 0; col < dim; col++) {

			cout << V[row][col] << " | ";

		}

		cout << endl;

	}

	cout << endl << "| ";

	for (int i = 0; i < dim; i++) {

		cout << A[i][i] << " | ";

	}

	cout << endl;

}

void Exercise_7_1(void) {

	cout << "---- Exercise 7: ---- " << endl << endl;

	int dim = 3;

	matrix A(dim);
	for (int i = 0; i < dim; i++) {
		A[i].resize(dim, 0);
	}

	A[0][0] = 0;
	A[0][1] = 2;
	A[0][2] = -1;
	A[1][0] = 2;
	A[1][1] = -1;
	A[1][2] = 1;
	A[2][0] = 2;
	A[2][1] = -1;
	A[2][2] = -3;


	JacobiMethod(A, dim);

	


}

