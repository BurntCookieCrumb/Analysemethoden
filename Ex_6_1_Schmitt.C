
/* ------------------------------------------------------------------------------------------ */

// Analysemethoden der Hochenergiephysik
// Uebungsblatt 6

// bearbeitet von Kristina Schmitt

/* ------------------------------------------------------------------------------------------ */


// Includes

#include <iostream>
#include <iomanip>
#include <cmath>

#include <stdlib.h>
#include <time.h>

#include "TF1.h"

// Namespace

using namespace std;

/*-------------------------------------------------------------------------------------------------------*/
//		Exercise 1
/*-------------------------------------------------------------------------------------------------------*/

// ***** General Definitions and Functions **********************************

int sign() {

    if (rand()%2 == 0) return +1;
    else return -1;

}


// ****** Classes ***********************************************************

// --- == Nucleon == --------------------------------------------------------

class Nucleon {

public:

	Nucleon(double x_pos, double y_pos);

	double GetX() const;
	double GetY() const;

	int GetNColl() const;

	void ShiftX(double shift);
	void AddCollision();

private:

	double x;
	double y;

	int N_Coll;


};

// == public ==

Nucleon::Nucleon(double x_pos, double y_pos)
	: x(x_pos), y(y_pos), N_Coll(0)
{

}

double Nucleon::GetX() const {

	/// Get x coordinate

	return x;

}

double Nucleon::GetY() const {

	/// Get y coordinate

	return y;

}

int Nucleon::GetNColl() const {

    return N_Coll;

}

void Nucleon::ShiftX(double shift) {

	/// Shift Nucleon in x direction by 'shift'

	x += shift;

}

void Nucleon::AddCollision() {

    N_Coll++;

}


// ---  == Nucleus == ---------------------------------------------------------

class Nucleus {

public:

	Nucleus(int mass_number, int proton_number, double radius, double r_N, double d);
	Nucleus(bool random);

	int GetMassNumber() const;
	int GetProtonNumber() const;

	Nucleon** const * GetContent() const;

	void ShiftInX(double offset);
	void PrintNucleusContent();

	int GetNCollisions(); /// !!!
	int GetNParticipants();

	~Nucleus();


private:

	int mass_number;
	int proton_number;

	double radius;
	double r_N;
	double d;

	Nucleon** content = NULL;

	void FillNucleus(double r, double r_N, double d);



};

// == public ==

Nucleus::Nucleus(int n_mass, int n_proton, double r /*fm*/, double r_N, double d)
	: mass_number(n_mass), proton_number(n_proton), radius(r), r_N(r_N), d(d)
{

	/// Full Constructor with all options

	content = new Nucleon*[n_mass];
	FillNucleus(r, r_N, d);

}

Nucleus::Nucleus(bool Pb)
	: mass_number(208), proton_number(82), radius(7), r_N(6.62), d(0.54)
{

	content = new Nucleon*[208];
	FillNucleus(7, 6.62, 0.54);

}

int Nucleus::GetMassNumber() const {

	return mass_number;

}

int Nucleus::GetProtonNumber() const {

	return proton_number;

}

Nucleon** const* Nucleus::GetContent() const {

	return &content;

}

void Nucleus::ShiftInX(double offset) {

	/// realize shift of Nucleons in x direction by 'offset'

	for (int nucleon = 0; nucleon < mass_number; nucleon++) {

		content[nucleon]->ShiftX(offset);

	}

}

void Nucleus::PrintNucleusContent() {

	/// Print Nucleon content of Nucleus

	for (int nucleon = 0; nucleon < mass_number; nucleon++) {
		cout << left << showpos << "Nukleon " << setw(4) << nucleon + 1
			<< " x: " << setw(8) << setprecision(5) << content[nucleon]->GetX()
			<< " y: " << setw(8) << setprecision(5) << content[nucleon]->GetY()
			<< endl;
	}

	cout << endl;

}

int Nucleus::GetNCollisions() {

    int coll = 0;

    for (int nucleon = 0; nucleon < mass_number; nucleon++) {
		 coll += content[nucleon]->GetNColl();
	}

	return coll;

}

// **
// -- Determination of #Participants ---------------------------------------------

int Nucleus::GetNParticipants(){

    int part = 0;

    for (int nucleon = 0; nucleon < mass_number; nucleon++) {
		 if (content[nucleon]->GetNColl() > 0) part++;
	}

	return part;

}

// ------------------------------------------------------------------------------
// **

Nucleus::~Nucleus() {

	for (int nucleon = 0; nucleon < mass_number; nucleon++) {
		delete content[nucleon];
	}

	delete[] content;
	content = NULL;

}

// == private ==

// **
// -- Fill Nucleus following the Woods-Saxon Potential -------------------------

void Nucleus::FillNucleus(double r, double r_N, double d) {

	/// Fill Nucleus randomly with Nucleons

	double x;
	double y;

	TF1* WoodsSaxon = new TF1("WoodsSaxon", "1. / (1. + exp((x - [0])/[1]))", -r, r);
	WoodsSaxon->SetParameters(0, r_N);
	WoodsSaxon->SetParameters(1, d);

	for (int nucleon = 0; nucleon < mass_number; nucleon++) {

		x = WoodsSaxon->GetRandom(-r, r);
		y = WoodsSaxon->GetRandom(-r, r);

		content[nucleon] = new Nucleon(x, y);

	}

    delete WoodsSaxon;

}

// -----------------------------------------------------------------------------
// **

// ***** Main Programm *********************************************************

// --- == Simulation == --------------------------------------------------------

int GlauberMonteCarlo(int collisions, double impct_prmtr, double sigma_mb = 60,
	int A1 = 208, int A2 = 208, int Z1 = 82, int Z2 = 82, double r1 = 7.0 /*fm*/,
	double r2 = 7.0 /*fm*/, double r_N1 = 6.62 /*fm*/, double r_N2 = 6.62 /*fm*/,
	double d1 = 0.54 /*fm*/, double d2 = 0.54 /*fm*/)
{

	/** Simulates 'coll' Nucleus-Nucleus collisions and returns the mean number of
	Nucleon-Nucleon collisions
	**/

	// Parameters


	int N_part = 0;
	int N_coll = 0;
	double sigma = sigma_mb / 10.; // 1 fm^2 = 10 mb

	double x_diff;
	double y_diff;

	double d_max = sqrt(sigma / M_PI);

	// Collisions

	for (int coll = 1; coll <= collisions; coll++) {

		// Initialize Nuclei

		Nucleus N1(true);
		Nucleus N2(true); //A2, Z2, r2, r_N2, d2

		N2.ShiftInX(impct_prmtr); // Shift corresponding to impact parameter

        // Check for hitting Nucleons

		for (int n1 = 0; n1 < A1; n1++) {
			for (int n2 = 0; n2 < A2; n2++) {

				x_diff = (*N1.GetContent())[n1]->GetX() - (*N2.GetContent())[n2]->GetX();
				y_diff = (*N1.GetContent())[n1]->GetY() - (*N2.GetContent())[n2]->GetY();

				if (sqrt(pow(x_diff, 2) + pow(y_diff, 2)) <= d_max) {

                    // **
				    // -- new determination of N_coll with class functions --

				    (*N1.GetContent())[n1]->AddCollision();
				    (*N2.GetContent())[n2]->AddCollision();

				    // ------------------------------------------------------
				    // **

				    N_coll ++;

				    }

			}
		}

		//N_coll += N1.GetNCollisions(); // adding all Nucleon-Nucleon collisions
		N_part += N1.GetNParticipants() + N2.GetNParticipants(); // adding number participants

	}

    cout << N_part / collisions << endl;
	return N_coll / collisions; // return mean

}

void Exercise_6_1(void) {

	cout << endl << " -- Exercise 1: ---------------------- " << endl << endl;

	cout << "Stossparameter b = 0fm -> <N_coll> = " << GlauberMonteCarlo(1000, 0) << endl;
	cout << "Stossparameter b = 4fm -> <N_coll> = " << GlauberMonteCarlo(1000, 4) << endl;
	cout << "Stossparameter b = 8fm -> <N_coll> = " << GlauberMonteCarlo(1000, 8) << endl;

	cout << endl;

}

void Ex_6_1_Schmitt(){

    srand((int)time(NULL));

    Exercise_6_1();

    //Nucleus N1(true);
    //N1.PrintNucleusContent();

}
