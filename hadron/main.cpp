
// g++ main.cpp -I`./lhapdf-config --incdir` -L`./lhapdf-config --libdir` -static -lLHAPDF

#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>
#include <iterator>
#include <random>
#include <cmath>       /* floor */
#include <stdio.h>
#include <signal.h>
#include "LHAPDF/LHAPDF.h"

using namespace std;
using namespace LHAPDF;

const double ALPHA_S 		 = 0.118;              // coupling constant alpha_s
const double TOP_QUARK_MASS  = 172.5;              // in GeV
const double ENERGY 		 = 6500;               // in GeV, energy of single proton
const double S               = 4 * pow(ENERGY, 2); // Mandelstam variable for CoM massless scattering
const double GEV_TO_PB 		 = 2.56819e-9;         // GeV^-2 to 1 pb (picobarn) conversion factor
const int    PRECISION		 = 10;

const double UPPER_BOUND_INC = 0.1;                // default: 0.1
const double UPPER_BOUND_DEC = 0.05;               // default: 0.05

const bool   CONVERT_TO_PB   = true;
const bool   GENERATE_POINT  = true;

const string FILENAME        = "output.dat";
const string PDF_NAME 		 = "NNPDF31_nnlo_as_0118";


// for reactivating cursor
void signal_callback_handler(int signum) {
   system("setterm -cursor on");
   exit(signum);
}

// class to use ofstream with the possibility to retrive filename
class myOfstream : public ofstream {
    public:
        string filename;

        myOfstream(string filename) {
            this->filename = filename;
            open(filename);
        }
};

// #region class IntegrationGeneration

class IntegrationGeneration {
public:
	IntegrationGeneration(double (f(double *)), const int Ndim, const int NumIntervals) {
		this->f = f;
		this->NumIntervals = NumIntervals;
		this->Ndim = Ndim;

		xi = new double * [Ndim];
		for (int k = 0; k < Ndim; k++) {
			xi[k] = new double[NumIntervals + 1];
			for (int j = 0; j <= NumIntervals; j++) xi[k][j] = ((double) j) / NumIntervals;
		}

		Si = new double * [Ndim];
		for (int k = 0; k < Ndim; k++) {
			Si[k] = new double[NumIntervals + 1];
			for (int j = 0; j <= NumIntervals; j++) Si[k][j] = 0;
		}

		Ti = new double * [Ndim];
		for (int k = 0; k < Ndim; k++) {
			Ti[k] = new double[NumIntervals + 1];
			for (int j = 0; j <= NumIntervals; j++) Ti[k][j] = 0;
		}

		Ii = new double * [Ndim];
		Ni = new long * [Ndim];
		for (int k = 0; k < Ndim; k++){
			Ii[k] = new double[NumIntervals + 1];
			Ni[k] = new long[NumIntervals + 1];
			for (int j = 0; j <= NumIntervals; j++) { Ii[k][j] = 0; Ni[k][j] = 0; }
		}
	}

	~IntegrationGeneration() {
		for (int k=0; k<Ndim; k++) {
			delete xi[k];
			delete Si[k];
			delete Ti[k];
			delete Ii[k];
			delete Ni[k];
		}
		delete xi;
		delete Si;
		delete Ti;
		delete Ii;
		delete Ni;
	}

	void refineGrid(long N);
    void multipleRefineGrid(int nTimes, long N);
	void computeIntegralAndError(long, double &, double &);
	void generatePoint(double *, int &, int &);
    void generatePointsToFile(long N, myOfstream &filestream, int &miss, int &ub_violation, string header);

private:
	double (*f)(double *);
	uniform_real_distribution<double> unif;
	default_random_engine re;
	int NumIntervals;
	int Ndim;
	double **xi;
	double **Si;
	double **Ti;
	double **Ii;
	long   **Ni;
	bool integrationDone = false;
	void findIndexAndx(double * y, int * m, double * x) {
		for (int k = 0; k < Ndim; k++) {
			int intpart_y = floor(y[k] * NumIntervals);
			int l = intpart_y + 1;
			m[k] = l;
			x[k] = xi[k][l-1] + (xi[k][l] - xi[k][l-1]) * (y[k] * NumIntervals - intpart_y);
		}
	}
};


// #region Implementation

// refine the grid performing N calls
void IntegrationGeneration::refineGrid(long N) {
	for (int k = 0; k < Ndim; k++) {
		for (int j = 0; j <= NumIntervals; j++) { Ii[k][j] = 0; Ni[k][j] = 0; }
	}

	system("setterm -cursor off");
	cout << fixed;
	cout << setprecision(2);
	for (long j = 0; j < N; j++) { // call the function and fill the Ii and Ni
		cout << "\r" << "Refining grid with N = " << N << "... " << (double)(j+1) * 100 / N << "%        ";

		double y[Ndim];
		for (int k = 0; k < Ndim; k++) y[k] = unif(re);
		int l[Ndim];
		double x[Ndim];
		findIndexAndx(y, l, x);
		double xval = f(x);
		for (int k = 0; k < Ndim; k++) {
			Ni[k][l[k]]++;
			Ii[k][l[k]] += xval * ( xi[k][l[k]] - xi[k][l[k]-1] );
		}
	} // end of loop on function calls
	cout << endl;
	system("setterm -cursor on");
	cout << defaultfloat;
	cout << setprecision(PRECISION);

	for (int k = 0; k < Ndim; k++) {
		for (int l = 1; l <= NumIntervals; l++) { // Compute the cumulated integral
			Ii[k][l] = Ii[k][l-1] + Ii[k][l] / Ni[k][l];
		}
		for (int l = 1; l <= NumIntervals; l++) { // Normalize it to 1
			Ii[k][l] = Ii[k][l] / Ii[k][NumIntervals];
		}
		// compute the new grid
		double xbar[NumIntervals];
		xbar[0] = 0;
		xbar[NumIntervals] = 1;
		int lp = 1;
		for (int l = 1; l < NumIntervals; l++) {
			double xl = ((double) l) / NumIntervals;
			while (Ii[k][lp] < xl) lp++;
			xbar[l] = xi[k][lp-1] + ( xl - Ii[k][lp-1] ) * ( xi[k][lp] - xi[k][lp-1] ) / ( Ii[k][lp] - Ii[k][lp-1] );
		}
		for (int l = 1; l < NumIntervals; l++) xi[k][l] = xbar[l];
	}
}

void IntegrationGeneration::multipleRefineGrid(int nTimes, long N) {
    system("setterm -cursor off");
    cout << fixed;
    cout << setprecision(2);
    for (int n = 0; n < nTimes; n++) {
        for (int k = 0; k < Ndim; k++) {
    		for (int j = 0; j <= NumIntervals; j++) { Ii[k][j] = 0; Ni[k][j] = 0; }
    	}

     	for (long j = 0; j < N; j++) { // call the function and fill the Ii and Ni
            cout << "\r" << "Refining grid " << (n+1) << "/" << nTimes << " with N = " << N << "... "
                 << (double)(j+1) * 100 / N << "% - TOT: " << (double)(j+1+N*n) * 100 / (N*nTimes) << "%       ";

    		double y[Ndim];
    		for (int k = 0; k < Ndim; k++) y[k] = unif(re);
    		int l[Ndim];
    		double x[Ndim];
    		findIndexAndx(y, l, x);
    		double xval = f(x);
    		for (int k = 0; k < Ndim; k++) {
    			Ni[k][l[k]]++;
    			Ii[k][l[k]] += xval * ( xi[k][l[k]] - xi[k][l[k]-1] );
    		}
    	} // end of loop on function calls

    	for (int k = 0; k < Ndim; k++) {
    		for (int l = 1; l <= NumIntervals; l++) { // Compute the cumulated integral
    			Ii[k][l] = Ii[k][l-1] + Ii[k][l] / Ni[k][l];
    		}
    		for (int l = 1; l <= NumIntervals; l++) { // Normalize it to 1
    			Ii[k][l] = Ii[k][l] / Ii[k][NumIntervals];
    		}
    		// compute the new grid
    		double xbar[NumIntervals];
    		xbar[0] = 0;
    		xbar[NumIntervals] = 1;
    		int lp = 1;
    		for (int l = 1; l < NumIntervals; l++) {
    			double xl = ((double) l) / NumIntervals;
    			while (Ii[k][lp] < xl) lp++;
    			xbar[l] = xi[k][lp-1] + ( xl - Ii[k][lp-1] ) * ( xi[k][lp] - xi[k][lp-1] ) / ( Ii[k][lp] - Ii[k][lp-1] );
    		}
    		for (int l = 1; l < NumIntervals; l++) xi[k][l] = xbar[l];
    	}
    }
    system("setterm -cursor on");
    cout << defaultfloat;
    cout << setprecision(PRECISION);
}


void IntegrationGeneration::computeIntegralAndError(long N, double &integ, double &err) {
	integ = 0;
	err = 0;

	system("setterm -cursor off");
	cout << fixed;
	cout << setprecision(2);
	for (long j = 0; j < N; j++) { // call the function and fill the Ii and Ni
		cout << "\r" << "Computing integral with N = " << N << "... " << (double)(j+1) * 100 / N << "%        ";

		double y[Ndim];
		for (int k = 0; k < Ndim; k++) y[k] = unif(re);
		int l[Ndim];
		double x[Ndim];
		findIndexAndx(y, l, x);
		double fval = f(x);
		double dint = fval;

		for (int k = 0; k < Ndim; k++) dint *= ( xi[k][l[k]] - xi[k][l[k]-1] ) * NumIntervals;
		integ += dint;
		err   += dint * dint;

		// Computation of the upper bound. Here we must use fantasy
		double current_ub = 1;
		for (int k = 0; k < Ndim; k++) current_ub *= Si[k][l[k]];
		if(current_ub < fval) {
			if(current_ub == 0) {
				double vvv = pow(fval, 1.0 / Ndim);
				for (int k = 0; k < Ndim; k++) {
					Si[k][l[k]] = vvv;
				}
			}
			else {
				double inc = 1 + UPPER_BOUND_INC / (Ndim); //increment upper bound
				for (int k = 0; k < Ndim; k++) Si[k][l[k]] *= inc;
			}
		}
		else {
			double inc = 1 - UPPER_BOUND_DEC / Ndim;
			for (int k = 0; k < Ndim; k++) Si[k][l[k]] *= inc;
		}
	}
	cout << endl;
	system("setterm -cursor on");
	cout << defaultfloat;
	cout << setprecision(PRECISION);

	integ = integ / N;
	err   = sqrt( (err / N - integ * integ) / N );

	// compute cumulant of Si //
	for (int k = 0; k < Ndim; k++) {
		for (int j = 1; j <= NumIntervals; j++) Ti[k][j] = Ti[k][j-1] + Si[k][j] * ( xi[k][j] - xi[k][j-1] );
		for (int j = 1; j <= NumIntervals; j++) Ti[k][j] = Ti[k][j] / Ti[k][NumIntervals]; // Normalize
	}

	ofstream cumulant;
	cumulant.open("cumulant1.dat");
	for (int k = 0; k < Ndim; k++) {
		for (int j = 0; j <= NumIntervals; j++) cumulant << ((double) j) / NumIntervals << ' ' << xi[k][j] << ' ' << Ti[k][j] << endl;
		cumulant << endl;
	}
	cumulant.close();

	integrationDone = true;
}

void IntegrationGeneration::generatePoint(double * xval, int &miss, int &ub_violation) {
	if (!integrationDone) {
		cout << "You are calling generatePoint before having done the integration!" << endl;
		exit(-1);
	}

 	int l[Ndim];
	double r[Ndim];
	double rr;
	double fval;
	double upper_b;

    for (int k = 0; k < Ndim; k++) r[k] = unif(re);
    for (int k = 0; k < Ndim; k++) {
        for (int j = 1; j <= NumIntervals; j++) {
            if (Ti[k][j] > r[k]) { l[k] = j; break; }
        }
        xval[k] = xi[k][l[k]-1] + ( r[k] - Ti[k][l[k]-1] ) / ( Ti[k][l[k]] - Ti[k][l[k]-1] ) * ( xi[k][l[k]] - xi[k][l[k]-1] );
    }
    fval = f(xval);   // Now the Hit-And-Miss
    rr = unif(re);
    upper_b = 1.0;
    for (int k = 0; k < Ndim; k++) upper_b *= Si[k][l[k]];
    if (fval > upper_b) ub_violation++;

    while (rr > fval / upper_b) {
        miss++;
		for (int k = 0; k < Ndim; k++) r[k] = unif(re);
		for (int k = 0; k < Ndim; k++) {
			for (int j = 1; j <= NumIntervals; j++) {
				if (Ti[k][j] > r[k]) { l[k] = j; break; }
			}
			xval[k] = xi[k][l[k]-1] + ( r[k] - Ti[k][l[k]-1] ) / ( Ti[k][l[k]] - Ti[k][l[k]-1] ) * ( xi[k][l[k]] - xi[k][l[k]-1] );
		}

		fval = f(xval);   // Now the Hit-And-Miss
		rr = unif(re);
		upper_b = 1;
		for (int k = 0; k < Ndim; k++) upper_b *= Si[k][l[k]];
		if (fval > upper_b) ub_violation++;
    }
}

void IntegrationGeneration::generatePointsToFile(long N, myOfstream &filestream, int &miss, int &ub_violation, string header = "") {
    miss = 0;
    ub_violation = 0; // upper bound violation

    if (filestream.is_open()) {
        if (header != "") filestream << header << endl;

        system("setterm -cursor off");
    	cout << fixed;
    	cout << setprecision(2);

        for (int j = 0; j < N; j++) {
            cout << "\r" << "Generating point to " << filestream.filename << " with N = " << N << "... " << (double)(j+1) * 100 / N << "%        ";
            double xval[Ndim];
            generatePoint(xval, miss, ub_violation);
            for (int i = 0; i < Ndim - 1; i++) {
                filestream << xval[i] << ",";
            }
            filestream << xval[Ndim - 1] << endl;
        }

        filestream.close();

        system("setterm -cursor on");
    	cout << defaultfloat;
    	cout << setprecision(PRECISION);
    }
    else {
        cout << "Error printing points to file" << endl;
        exit(1);
    }
}

// #endregion Implementation

// #endregion class IntegrationGeneration

// print info to the first two lines of file
void printInfoToFile(myOfstream &filestream) {
    if (filestream.is_open()) {
        filestream << "alpha_s,top_quark_mass,energy,upper_bound_inc,upper_bound_dec" << endl;
        filestream << ALPHA_S << "," << TOP_QUARK_MASS << "," << ENERGY << "," << UPPER_BOUND_INC << "," << UPPER_BOUND_DEC << endl;
        filestream << endl; // leave empty line
    }
    else {
        cout << "Error printing info to file" << endl;
        exit(1);
    }
}

// mQ: mass of produced heavy particle
// cross section formula for quark-quark scattering
double qqCrossSection(double s, double t, double u, double mQ, double alpha = ALPHA_S) {
	double m2 = pow(mQ, 2);
	if (s < 4*m2) return 0.0;
	else {
		double factor = pow(alpha, 2) / (9 * pow(s, 3)) * sqrt(1 - (4 * m2) / s);
		return factor * ( pow(m2 - t, 2) + pow(m2 - u, 2) + 2*m2*s );
	}
}

// mQ: mass of produced heavy particle
// cross section formula for gluon-gluon scattering
double ggCrossSection(double s, double t, double u, double mQ, double alpha = ALPHA_S) {
	double m2 = pow(mQ, 2);
	if (s < 4*m2) return 0.0;
	else {
		double factor = pow(alpha, 2) / (32 * s) * sqrt(1 - (4 * m2) / s);
		return factor * ( 6 * (m2 - t) * (m2 - u) / pow(s, 2) -
	 					  m2 * (s - 4*m2) / (3 * (m2 - t) * (m2 - u)) +
					      4 * ((m2 - t) * (m2 - u) - 2*m2 * (m2 + t)) / (3 * pow((m2 - t), 2)) +
					  	  4 * ((m2 - t) * (m2 - u) - 2*m2 * (m2 + u)) / (3 * pow((m2 - u), 2)) -
					      3 * ((m2 - t) * (m2 - u) + m2 * (u - t)) / (s * (m2 - t)) -
					      3 * ((m2 - t) * (m2 - u) + m2 * (t - u)) / (s * (m2 - u)) );
	}
}

// function to integrate; args are in range (0, 1)
// arg[0] = theta
// arg[1] = x1
// arg[2] = x2
double TotalCrossSection(double * arg) {
	double mQ = TOP_QUARK_MASS;
    double m2 = pow(mQ, 2);

	double th_a = 0;     // min range of theta
	double th_b = M_PI;	 // max range of theta
	double th   = th_a + (th_b - th_a) * arg[0]; // theta in (th_a, th_b)

	double x1 = arg[1];
	double x2 = arg[2];

    // in the center of mass
    double s = x1 * x2 * S;
    if (s < 4*m2) return 0.0;
    double t = m2 - (s/2) * ( 1 - sqrt(1 - 4*m2/s)*cos(th) );
    double u = 2*m2 - t - s;

	double res = 0.0;
	double temp;

	// loop through all pair of light quarks and gluons
	// fl = -2 <--> ~u  u
	// fl = -2 <--> ~d  d
	// fl =  0 <-->  g  g
	// fl =  1 <-->  d ~d
	// fl =  2 <-->  u ~u
	for (int fl = -2; fl <= 2; fl++) {
		if (fl == 0) temp = ggCrossSection(s, t, u, mQ);
		else temp = qqCrossSection(s, t, u, mQ);

        // xfxM() returns x*f(x, mQ, fl), so it must be divided by x
        temp *= xfxM(1, x1, mQ, fl) * xfxM(1, x2, mQ, -fl) / (x1 * x2);
		res += temp;
	}

	// there's a 2pi factor due to the integration in phi,
	// a (th_b - th_a) factor due to change of integration variable from (0, 1) to (th_a, th_b)
	// and sin(theta) because of d_omega = sin(theta) d_theta d_phi
	res *= 2 * M_PI * (th_b - th_a) * sin(th);
	if (CONVERT_TO_PB) return res / GEV_TO_PB; // converted from GeV^-2 to pb (picobarn)
    else return res;
}

int main() {
    string filename = FILENAME;

    string unit;
    if (CONVERT_TO_PB) unit = "pb";
    else unit = "GeV^-2";

    signal(SIGINT, signal_callback_handler);
	initPDFSetM(1, PDF_NAME, LHGRID);
	initPDFM(1, 0);
	cout << endl;

    int Ngrid = 100000;        // N used for refineGrid()
    int numberOfRefine = 5;    // number of times to call refineGrid()
    int Nint  = 1000000;       // N used for calculating integral
    int Num_Intervals = 1000;  // number of intervals used to calculate integral
    int Ndim = 3;              // number of dimension of function to integrate

    cout << "Ngrid = " << Ngrid << endl;
    cout << "numberOfRefine = " << numberOfRefine << endl;
    cout << "Nint = " << Nint << endl;
    cout << "Num_Intervals = " << Num_Intervals << endl;
    cout << endl;

	IntegrationGeneration intf(TotalCrossSection, Ndim, Num_Intervals);

    double integ, err;

    intf.multipleRefineGrid(numberOfRefine, Ngrid);
    cout << endl;

	intf.computeIntegralAndError(Nint, integ, err);

    if (CONVERT_TO_PB) cout << fixed << setprecision(PRECISION) << endl;
    else cout << scientific;
	cout << "Integral = " << integ << " " << unit << endl;
	cout << "Error    = " << err << " " << unit << endl << endl;

    if (GENERATE_POINT) {
        myOfstream filestream(filename); // also opens the file
        printInfoToFile(filestream);

        int miss, ub_violation;
        int Npoints = 100000;
        string header = "theta,x1,x2";
        intf.generatePointsToFile(Npoints, filestream, miss, ub_violation, header);

        filestream.close();

        cout << fixed << setprecision(PRECISION) << endl << endl;
        cout << "# miss = " << (double)miss / Npoints << endl;
	    cout << "# upper bound violation = " << (double)ub_violation / Npoints << endl;
    }

    cout << endl;
	return 0;
}
