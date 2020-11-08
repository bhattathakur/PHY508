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
  //initializes commonly used parameters from a file
  //int Nlin; // linear size of lattice
  //double beta; // 1/T
  //int Neql;// eql sweeps
  //int Nmcs;// sweeps in a bin
  //int Nbin;// number of bin
  //int SEED; // for MT
  //string latt_;// lattice kind
  ifstream pfin;
  pfin.open("param.dat");  
  if (pfin.is_open()) { 
    pfin >> Nlin;
    pfin >> beta;
    pfin >> Neql;
    pfin >> Nmcs;
    pfin >> Nbin;
    pfin >> SEED;
    pfin >> latt_;
  }
  else
  {cout << "No input file to read ... exiting!"<<endl;exit(1);}
  pfin.close();
  // print out all parameters for record
  cout << "---Input Parameters---"<<endl; 
  cout <<"Nlin = "<<Nlin<<"; beta= "<<beta<<endl;
  cout <<"Neql= "<<Neql<<"\t Nmcs\t"<<Nmcs<<"\t; Number of bins = "<<Nbin<<endl;
  cout <<"RNG will be given SEED of = "<<SEED<<endl;
  cout <<"Percolation problem on lattice --> "<<latt_<<endl;
};//constructor

LATTICE::LATTICE (const PARAMS& p)
{
  if(p.latt_=="sqlatt_PBC")
  {
    Lx=p.Nlin;Ly=p.Nlin;
    Nsite=Lx*Ly;
    for(int site=0;site<Nsite;site++)  //site runs from 0 to Nlin*Nlin	      
    {
      int x=site%Lx;                  //gives x of site
      int y=site/Ly;                  //gives y of site
      vector <int> nrn; 		//vector to hold neighbour elements
      (x+1<Lx)?nrn.push_back(site+1):nrn.push_back(site+1-Lx);  //use of ternary operator
      (y+1<Ly)?nrn.push_back(site+Ly):nrn.push_back(site-y*Ly);
      (x-1>=0)?nrn.push_back(site-1):nrn.push_back(site-1+Lx);
      (y-1>=0)?nrn.push_back(site-Ly):nrn.push_back(Nsite-Ly+site);
      nrnbrs.push_back(nrn); 		//neighbour vector is added to the main vector nrnbrs
    }
  }

  else if(p.latt_=="sqlatt_OBC")
  {
    Lx=p.Nlin;Ly=p.Nlin;
    Nsite=Lx*Ly;
    for(int site=0;site<Nsite;site++) 
    {
      int x=site%Lx;                   
      int y=site/Ly;		        
      vector <int> nrn;                
      if(x+1<Lx)nrn.push_back(site+1);
      if(y+1<Ly)nrn.push_back(site+Ly);
      if(x-1>=0)nrn.push_back(site-1);
      if(y-1>=0)nrn.push_back(site-Ly);
      nrnbrs.push_back(nrn);         
    }
  }
  else
  {cout <<"Dont know your option for lattice in param.dat .. exiting"<<endl;exit(1);}
}
void LATTICE::print()
{
  //THIS FUNCTIONS MAY BE CALLED DURING DEBUGGING TO MAKE SURE LATTICE HAS BEEN DEFINED CORRECTLY
  cout <<"---printing out properties of lattice ---"<<endl;
  cout<<"size is  "<<Lx<<"x"<<Ly<<endl;
  cout <<"neighbors are"<<endl;
  for (int site=0;site<Nsite;site++)
  {
    cout <<site<<" : ";
    for (int nn=0;nn<nrnbrs.at(site).size();nn++)
      cout<<nrnbrs.at(site).at(nn)<<" ";
    cout <<endl;
  }
}

ISING_CONF::ISING_CONF(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  // INITIALIZE SPIN CONFIG
  // ---------COMPLETE CODE HERE--------!

  //size of spin conf vector is Nsite
  spin.resize(latt.Nsite);

  //initilize the spin with 0 or 1
  for (int i=0;i<latt.Nsite;i++)spin[i]=ran1.randInt(1);    //randInt(n)->0,..,n

  //CREATE WEIGHT TABLE
  //possible sum of spin of neighbours of a site is 5 (0,1,2,3,4) 0->all down, 4->all up
  //new and old both spin configuration have 2 possible values 0, 1
  //define 2*2 vector for the new and old spin configuration
  //vector <vector<double> > oldnewspin(2,vector<double>(2,0)); //2*2 vector of doubles initilized with 0

  //resize wght_tbl the 5*oldnewspin vector named wght_tbl

  //wght_tbl.resize(5,oldnewspin);
  if(p.latt_=="sqlatt_PBC")
  {
    vector <vector<double> > oldnewspin(2,vector<double>(2,0)); //2*2 vector of doubles initilized with 0
    wght_tbl.resize(5,oldnewspin); //resize wght_tbl
    for (int i=0;i<5;i++){
      for (int j=0;j<2;j++){
        for (int k=0;k<2;k++){
          wght_tbl[i][j][k]=exp(2.*p.beta*1.0*(j-k)*(2.0*i-4));
        }
      }
    }

  // ---------COMPLETE CODE HERE--------!
  }
  else
  {cout <<"NEED TO CODE ALL LATTICE OPTIONS"<<endl;}
}

ISING_CONF::~ISING_CONF()
{
  dfout.close();
}
//writes the configuration in a file
void ISING_CONF::conf_write(const PARAMS& p, const LATTICE& latt)
{
  ofstream outputfile;
  string filename="isingconf.dat"; //stores configuration
  outputfile.open(filename);
  if(outputfile.is_open()){
    for(int i=0;i<latt.Nsite;i++)outputfile<<spin[i]<<endl;
    outputfile.close();
  }
  else cout<<filename<<"\t doesnot exist"<<endl;

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

