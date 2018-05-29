

/* ------------------------------------------------------------------------------------------ */

// Analysemethoden der Hochenergiephysik
// Uebungsblatt 4

// bearbeitet von Kristina Schmitt

/* ------------------------------------------------------------------------------------------ */


// Includes

#include <iostream>
#include <iomanip>
#include <cmath>

#include <stdlib.h>
#include <time.h>

// Namespace

using namespace std;

/*-------------------------------------------------------------------------------------------------------*/
//		Exercise 1
/*-------------------------------------------------------------------------------------------------------*/

// --- == Nucleon == --------------------------------------------------------

class Nucleon {

    public:

        Nucleon(double x_pos, double y_pos);

        double GetX() const;
        double GetY() const;

        void ShiftX(double shift);

    private:

        double x;
        double y;


};

Nucleon::Nucleon(double x_pos, double y_pos)
: x(x_pos), y(y_pos)
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

void Nucleon::ShiftX(double shift){

    /// Shift Nucleon in x direction by 'shift'

    x += shift;

}


// ---  == Nucleus == ---------------------------------------------------------

class Nucleus{

    public:

        Nucleus(int mass_number, int proton_number, double radius);

        int GetMassNumber() const;
        int GetProtonNumber() const;

        Nucleon** const * GetContent() const;

        void ShiftInX(double offset);
        void PrintNucleusContent();

        ~Nucleus();


    private:

        int mass_number;
        int proton_number;

        double radius;

        Nucleon** content = NULL;

        void FillNucleus();



};

Nucleus::Nucleus(int n_mass, int n_proton, double r /*fm*/)
: mass_number(n_mass), proton_number(n_proton), radius(r)
{

    content = new Nucleon*[n_mass];
    FillNucleus();

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

    for (int nucleon = 0; nucleon < mass_number; nucleon++){

        content[nucleon]->ShiftX(offset);

    }

}

void Nucleus::FillNucleus(){

    /// Fill Nucleus randomly with Nucleons

    double x;
    double y;

    for (int nucleon = 0; nucleon < mass_number; nucleon++){

        x = ((double)rand()/RAND_MAX)*2*radius - radius;
        y = ((double)rand()/RAND_MAX)*2*radius - radius;

        content[nucleon] = new Nucleon(x,y);

    }

}

void Nucleus::PrintNucleusContent(){

    /// Print Nucleon content of Nucleus

    for (int nucleon = 0; nucleon < mass_number; nucleon++){
        cout << left << showpos << "Nukleon " << setw(4) << nucleon + 1
         << " x: " << setw(8) << setprecision(5) << content[nucleon]->GetX()
         << " y: " << setw(8) << setprecision(5) << content[nucleon]->GetY()
         << endl;
    }

    cout << endl;

}

Nucleus::~Nucleus(){

    for (int nucleon = 0; nucleon < mass_number; nucleon++){
        delete content[nucleon];
    }

    delete [] content;
    content = NULL;

}


// --- == Simulation == --------------------------------------------------------

int GlauberMonteCarlo(int collisions, double impct_prmtr, double sigma_mb = 60,
                      int A1 = 208, int A2 = 208, int Z1 = 82, int Z2 = 82,
                      double r1 = 7.0 /*fm*/, double r2 = 7.0 /*fm*/)
{

    /** Simulates 'coll' Nucleus-Nucleus collisions and returns the mean number of
        Nucleon-Nucleon collisions
    **/

    // Parameters

    int N_coll = 0;
    double sigma = sigma_mb / 10.; // 1 fm^2 = 10 mb

    double x_diff;
    double y_diff;

    double d_max = sqrt(sigma / M_PI);

    // Collisions

    for (int coll = 1; coll <= collisions; coll++){

        // Initialize Nuclei

        Nucleus N1(A1, Z1, r1);
        Nucleus N2(A2, Z2, r2);

        N2.ShiftInX(impct_prmtr); // Shift corresponding to impact parameter

        // Check for hitting Nucleons

        for (int n1 = 0; n1 < A1; n1++){
            for (int n2 = 0; n2 < A2; n2++){

                x_diff = (*N1.GetContent())[n1]->GetX() - (*N2.GetContent())[n2]->GetX();
                y_diff = (*N1.GetContent())[n1]->GetY() - (*N2.GetContent())[n2]->GetY();

                if (sqrt(pow(x_diff, 2) + pow(y_diff, 2)) <= d_max) N_coll++;

            }
        }
    }

    return N_coll / collisions; // return mean

}

void Exercise_5_1(void){

    cout << endl << " -- Exercise 1: ---------------------- " << endl << endl;

    cout << "Stossparameter b = 0fm -> N_coll = " << GlauberMonteCarlo(1000, 0) << endl;
    cout << "Stossparameter b = 4fm -> N_coll = " << GlauberMonteCarlo(1000, 4) << endl;
    cout << "Stossparameter b = 8fm -> N_coll = " << GlauberMonteCarlo(1000, 8) << endl;





}


int main(){

    srand((int)time(NULL));

    Exercise_5_1();

    return 0;

}
