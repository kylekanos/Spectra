/*
 * Read & process Average*.curve files that contain information on spectra
 *   Kyle Kanos, 2015
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;

void Spectra(double *, double *, double *, double* );

const int nvar=8+3;            // total number of variables
const int ncell=500;           // spatial cells
const int np=16;               // particle spectrum
const int nb=122;              // photon spectrum 
const double ylow=-9.5;
const double yhigh=+9.5;
const double dlp0=(yhigh-ylow)/(double)(np-1);
const double pscale=3.164311153797431e-5;
const double c=2.9979e10;
const double mc2=1.67e-24*c*c;

// @brief return a padded-zero number
std::string Padding(int num){
    std::ostringstream ss;
    ss << std::setw(5) << std::setfill('0') << num;
    return ss.str();
};

// @brief the pion object
class PionDef{
  public:
    void FillPion();
    double *g;
    double **M;
    double *beta;
    double *sigma;
};

PionDef Pion;

// @brief the spectrum object
class Spectrum{
  public:
    void InitSpectrum();
    void ComputeSpectrum();
    void ShockSpectrum(int);
    void WriteSpectrum(int);
  private:
    double *x;   // pointer array to positions
    double *p;   // pointer array to momenta
    double **q;  // pointer array to variables
    double **n;  // pointer array to particle spectrum
    
};

int main(){
    Spectrum x;
    Pion.FillPion();
    x.InitSpectrum();

    for(int n = 0; n < 81; n++){
        x.FillSpectrum(n);
        x.ComputeSpectrum();
        x.ShockSpectrum(n);
        x.WriteSpectrum(n);
    }
    return 0;
}


// @brief initialize the pion structure
void PionDef::FillPion(){
  // allocate memory space
    M = new double *[nb];
    for(int ip = 0; ip < nb; ip++){
        M[ip] = new double [np];
    }
    sigma = new double [np];
    beta = new double [np];
    g = new double [nb];
    
  // open pion.matrix.data and read in file
    ifstream fp;
    fp.open("pion.matrix.data",ifstream::in);
  // read line of beta"s
    for(int i = 0; i < np; i++){
        fp >> beta[i];
    }
    for(int i = 0; i < np; i++){
        if(i < 5){
            sigma[i] = 0.0;
        }else{
            fp >> sigma[i];
        }
    }
  // read in matrix
    for(int i = 0; i < nb; i++){
        fp >> g[i];
        for(int j = 0; j < np; j++){
            if(j < 5){
                M[i][j] = 0.0;
            }else{
                fp >> M[i][j];
            }
        }
    }    
}

// @brief Initialize the spectrum
void Spectrum::InitSpectrum(){
    x = new double [ncell];
    p = new double [np];
    q = new double *[ncell];
    n = new double *[nb];
    for(int ic = 0; ic < ncell; ic++){
        q[ic] = new double [nvar];
    }
    for(int ib = 0; ib < nb; ib++){
        n[ib] = new double [ncell]; // n => n[1->nb][1->ncell]
        for(int nc = 0; nc < ncell; nc++){
            n[ib][nc] = 0.0;
        }
    }
    p[0] = exp(ylow);
    for(int ip = 1; ip < np; ip++){
        p[ip] = p[ip-1]*exp(dlp0);
    }
}

// @brief compute the spectrum (all the dirty work here)
void Spectrum::ComputeSpectrum(){
    double *Ne, pp;
    Ne = new double [np]; // temporary helpers
    for(int i=0; i<ncell; i++){
        
        
   // your spectrum from your data goes here!
        Ne = 0.0;
        
        double dE = exp(dlp0);
      // loop over photon spectrum
        for(int ib = 0; ib < nb; ib++){
            double Fnu = 0.0;
            for(int ip = 0; ip < np; ip++){
                Fnu += dE*Pion.M[ib][ip]*Pion.sigma[ip]*Pion.beta[ip]*Ne[ip];
            }
            n[ib][i] += q[i][0]*Fnu;
        }
    }
}

// @brief This fills the x & q variable arrays
void Spectrum::FillSpectrum(int nframe){
    int ic,nv;
    std::string filename, fnum, line, hash;
    hash = '#';
  // write the file
    fnum = Padding(nframe);
    filename = "out/Average" + fnum + ".curve";
    cout << "Reading file " << filename << endl;
    
    std::ifstream fp;
    fp.open(filename.c_str(),ifstream::in);
    
    nv = -1;
    ic = ncell+1;
    while(!fp.eof()){
        getline(fp, line);
      // if we've hit more than enough cells, we need to do some resetting
        if(ic >= ncell){
            ic = 0;
            nv++;
            getline(fp, line);
        }
        std::istringstream iss(line);
        double position, variable;
        sscanf(line.c_str(), "%lf %lf", &position, &variable);
        if(variable != variable){
            ic++;
            continue;
        }
        x[ic] = position;
        q[ic][nv] = variable;
        ic++;
    }
    fp.close();
}

// @brief write the spectrum (just use GNUPLOT compatible for now)
void Spectrum::WriteSpectrum(int nframe){
    std::string filename, fnum;
    fnum = Padding(nframe);
    filename = "out/Spectrum"+fnum+".dat";
    cout << "Writing file " << filename << endl << endl;
    ofstream fp;
    fp.open(filename.c_str(), ofstream::out);
    fp << scientific;
    for(int ib = 0; ib < nb; ib++){
        for(int ic = 0; ic < ncell; ic++){
            fp << Pion.g[ib] << "\t" << x[ic] << "\t" << n[ib][ic] << endl;
        }
        fp << endl;
    }
    fp.close();
}



    
        
    
