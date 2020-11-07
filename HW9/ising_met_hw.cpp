#include <fstream> 
#include <iostream> 
#include <math.h>
#include <string>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include "MersenneTwister.h"

using namespace std;

class PARAMS
{
public:
  int Nlin; // linear size of lattice
  double beta; // 1/T
  int Neql;// eql sweeps
  int Nmcs;// sweeps in a bin
  int Nbin;// number of bin
  int SEED; // for MT
  string latt_;// lattice kind
  PARAMS();//constructor
};

class LATTICE
{
public:
  LATTICE(const PARAMS&);//constructor
  int Nsite;// number of lattice sites
  int Lx,Ly;
  vector<vector<int> > nrnbrs;
  void print ();
  string latt_;
};

class ISING_CONF
{
public:
  ISING_CONF(const PARAMS&, const LATTICE&, MTRand&);//constructor
  void conf_write(const PARAMS&, const LATTICE& latt);
  void sweep(const PARAMS&,const LATTICE&,MTRand&);
  void meas_clear(const PARAMS&,const LATTICE&,MTRand&);
  void meas(const PARAMS&,const LATTICE&,MTRand&);
  void binwrite(const PARAMS&,const LATTICE&,MTRand&);
  ~ISING_CONF();// destructor
private:
  vector <int> spin;
  ofstream dfout;
  vector <vector<vector<double> > > wght_tbl;
  // these are the observables
  double energy;
  double energy_sq;
  double mag;
  double mag_sq;
};

int main(void)
{
  PARAMS p;
  LATTICE latt(p);
  MTRand ran1(p.SEED);
  ISING_CONF ising(p,latt,ran1);

  // latt.print();
  //EQUILIBRATE
   for(int eql=0;eql<p.Neql;eql++)
     ising.sweep(p,latt,ran1);
  
  //PRODUCTION
   for (int bin=0;bin<p.Nbin;bin++)
     {
       ising.meas_clear(p,latt,ran1);
       for(int mcs=0;mcs<p.Nmcs;mcs++)
	 {
	    ising.sweep(p,latt,ran1);
	    ising.meas(p,latt,ran1);
	 }
       ising.binwrite(p,latt,ran1);
     }
}

PARAMS::PARAMS(){
 // ---------COMPLETE CODE HERE--------!
};//constructor


LATTICE::LATTICE (const PARAMS& p)
{
 // ---------COMPLETE CODE HERE--------!
}

void LATTICE::print()
{
 // ---------COMPLETE CODE HERE--------!
}


ISING_CONF::ISING_CONF(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  // INITIALIZE SPIN CONFIG
  // ---------COMPLETE CODE HERE--------!

  //CREATE WEIGHT TABLE
  if(p.latt_=="sqlatt_PBC")
    {
  // ---------COMPLETE CODE HERE--------!
    }
  else
    {cout <<"NEED TO CODE ALL LATTICE OPTIONS"<<endl;}
}

ISING_CONF::~ISING_CONF()
{
  dfout.close();
}

void ISING_CONF::conf_write(const PARAMS& p, const LATTICE& latt)
{
 // ---------COMPLETE CODE HERE--------!
}


void ISING_CONF::sweep(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  /*SWEEP THROUGH THE LATTICE*/
 // ---------COMPLETE CODE HERE--------!
}

void ISING_CONF::meas_clear(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
// ---------COMPLETE CODE HERE--------!
}

void ISING_CONF::meas(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
 // ---------COMPLETE CODE HERE--------!
}

void ISING_CONF::binwrite(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
// ---------COMPLETE CODE HERE--------!
}

