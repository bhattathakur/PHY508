#include <fstream> 
#include <iostream> 
#include <math.h>
#include <string>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include "../HW9/MersenneTwister.h"

using namespace std;

class PARAMS
{
public:
  int Nlin; // linear size of lattice
  double betaJ; // Ising coupling J/T; J>0: ferro
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
  int Nnn;
};

class ISING_CONF
{
public:
  ISING_CONF(const PARAMS&, const LATTICE&, MTRand&);//constructor
  void conf_write(const PARAMS&, const LATTICE& latt);
  void sweep(const PARAMS&,const LATTICE&,MTRand&);
  void wolff_update(const PARAMS&, const LATTICE&,MTRand&);
  void meas_clear(const PARAMS&,const LATTICE&,MTRand&);
  void meas(const PARAMS&,const LATTICE&,MTRand&);
  void binwrite(const PARAMS&,const LATTICE&,MTRand&);
  ~ISING_CONF();// destructor
  vector <int> spin;
private:

  ofstream dfout;
  vector <vector<vector<double> > > wght_tbl;
  // these are the observables
  double energy;
  double energy_sq;
  double mag;
  double mag_sq;
  double pr;
  double wolff_size;
};

int main(void)
{
  PARAMS p;
  LATTICE latt(p);
  MTRand ran1(p.SEED);
  ISING_CONF ising(p,latt,ran1);

  //EQUILIBRATE
   for(int eql=0;eql<p.Neql;eql++)
    ising.wolff_update(p,latt,ran1);
// ising.sweep(p,latt,ran1);
  //PRODUCTION
   for(int bin=0;bin<p.Nbin;bin++)
     {
       ising.meas_clear(p,latt,ran1);
       for(int mcs=0;mcs<p.Nmcs;mcs++)
       {
         ising.wolff_update(p,latt,ran1);
         // ising.sweep(p,latt,ran1);
         ising.meas(p,latt,ran1);
       }
       ising.binwrite(p,latt,ran1);
     }
   // ising.conf_write(p,latt);
}

PARAMS::PARAMS(){
  //initializes commonly used parameters from a file
  ifstream pfin;
  pfin.open("param.dat");  
  if (pfin.is_open()) { 
    pfin >> Nlin;
    pfin >> betaJ;
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
  cout << "--- Parameters at input for single site metropolis on Ising model ---"<<endl; 
  cout <<"Nlin = "<<Nlin<<"; J/T = "<<betaJ<<endl;
  cout <<"Number of eql sweeps = "<<Neql<<"; Number of sweeps in a bin = "<<Nmcs<<"; Number of bins = "<<Nbin<<endl;
  cout <<"RNG will be given SEED of = "<<SEED<<endl;
  cout <<"Ising model on lattice --> "<<latt_<<endl;
};//constructor


LATTICE::LATTICE (const PARAMS& p)
{
  latt_=p.latt_;
  if(p.latt_=="sqlatt_PBC")
    {
      Lx=p.Nlin;Ly=p.Nlin;
      Nsite=Lx*Ly;
      nrnbrs.resize(Nsite,vector<int>(4,0));
      Nnn=4;
      int px=0; int py=1;
      int mx=2; int my=3;
      for (int site=0;site<Nsite;site++)
	{
	  nrnbrs[site][px]=(site/Lx)*Lx + ((site%Lx+1)%Lx);
	  nrnbrs[site][py]=((site/Lx+1)%Lx)*Lx + (site)%Lx;
	  nrnbrs[site][mx]=(site/Lx)*Lx + ((site%Lx-1+Lx)%Lx);
	  nrnbrs[site][my]=((site/Lx-1+Lx)%Lx)*Lx + (site)%Lx;
	}
    }
  else if(p.latt_=="trilatt_PBC")
    {
      Lx=p.Nlin;Ly=p.Nlin;
      Nsite=Lx*Ly;
      nrnbrs.resize(Nsite,vector<int>(6,0));
      Nnn=6;
      int px=0; int py=1;
      int mx=2; int my=3;
      int pxpy=4;int mxmy=5;
      for (int site=0;site<Nsite;site++)
	{
	  nrnbrs[site][px]=(site/Lx)*Lx + ((site%Lx+1)%Lx);
	  nrnbrs[site][py]=((site/Lx+1)%Lx)*Lx + (site)%Lx;
	  nrnbrs[site][mx]=(site/Lx)*Lx + ((site%Lx-1+Lx)%Lx);
	  nrnbrs[site][my]=((site/Lx-1+Lx)%Lx)*Lx + (site)%Lx;
	  nrnbrs[site][pxpy]=((site/Lx+1)%Lx)*Lx +  ((site%Lx+1)%Lx);
	  nrnbrs[site][mxmy]=((site/Lx-1+Lx)%Lx)*Lx +  ((site%Lx-1+Lx)%Lx);
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
      cout <<site<<" : ["<<nrnbrs.at(site).size()<<"]  ";
      for (int nn=0;nn<nrnbrs.at(site).size();nn++)
	cout<<nrnbrs.at(site).at(nn)<<" ";
      cout <<endl;
    }
}


ISING_CONF::ISING_CONF(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  // INITIALIZE SPIN CONFIG 
  spin.resize(latt.Nsite);
  if(p.Neql!=0)
    {
      cout << "Initializing randomly"<<endl;
      for (int site=0;site<latt.Nsite;site++)
	spin.at(site)=ran1.randInt(1);//random 0s and 1s
    }
  else
    {
      cout << "Initializing from conf.in"<<endl;
      ifstream cfin;
      cfin.open("conf.in");
      if (!cfin.is_open())
	{cout <<"CANT FIND conf FILE"<<endl;exit(1);}
      else 
	{ 
	  for (int site=0;site<latt.Nsite;site++)
	    cfin >>spin.at(site);
	}
      cfin.close();
    }

  //CREATE WEIGHT TABLE
  vector<vector<double> > tmp;
  tmp.resize(2,vector<double>(2));
  wght_tbl.resize(latt.Nnn+1,tmp);
  
  for(int mag=0;mag<(latt.Nnn+1);mag++)
    for(int os=0;os<2;os++)
      for(int ns=0;ns<2;ns++)
	wght_tbl.at(mag).at(ns).at(os)=exp(2.*p.betaJ* (double)((ns-os)*(2*mag-latt.Nnn)));
  
}

ISING_CONF::~ISING_CONF()
{

  dfout.close();
}

void ISING_CONF::conf_write(const PARAMS& p, const LATTICE& latt)
{
  ofstream cfout;
  cfout.open("ssfdata.out");
  for (int y=0;y<latt.Ly;y++)
    {
      for (int x=0;x<latt.Lx;x++)
	cfout<<spin.at(x+y*latt.Lx)<<" ";
      cfout<<endl;
    }
  cfout.close();
}


void ISING_CONF::sweep(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  // Nsite SSF updates
  for(int repeat=0;repeat<latt.Nsite;repeat++)
    {
      int site=ran1.randInt(latt.Nsite-1);
      int os=spin[site];
      int ns=1-os;
     
      int mag=0;
      for(int nbr=0;nbr<latt.nrnbrs[site].size();nbr++)
	mag+=spin[latt.nrnbrs[site][nbr]]; 

      if(ran1.randDblExc()<wght_tbl[mag][ns][os])
	spin[site]=ns;
    }// end loop on sites

}
//flips the single cluster based on the Wolff algorithm
void ISING_CONF::wolff_update(const PARAMS& p, const LATTICE& latt,MTRand& ran)
{

  //go to random site of the cluster
  int ransite=ran.randInt(latt.Nsite-1); //randInt(n) -> [0,n] 

  //creat a stack vector to store the sites which are asked and accepted
  vector <int> stk;

  //create a vector which stores the sites which are asked
  //at the begining all the elements are assigned to 0, asked will be updated to 1.

  vector <int> asked(latt.Nsite,0);//size same as the size of the lattice

  //update stack vector and asked vector with the ransite
  stk.push_back(ransite); //random site added to stack
  asked[ransite]=1;       //value of random site in asked vector is updated to 1

  //get the spin of random site
  int ran_spin=spin[ransite];
  
  //check the nearest neighbors of values in the stack
  //if the neighbour have same spin as the site and its probability less than activation probability it is 
  //added in the stack and asked vector is also updated

  //activation probabiliy from defination of Wolff algorithm
  pr=1.0-exp(-2.0*p.betaJ);
  
  for (int site=0;site<stk.size();site++)       //loop in the stack
  {
    int stk_value=stk[site];                    //get the stack value
    int nbrsize=latt.nrnbrs[stk_value].size();  //size of nearest neighbours
    //cout<<"size of the nearest neighbour\t"<<nbrsize<<endl;
    for(int nb=0;nb<nbrsize;nb++)               //loop over neighbors
    {
      int test=latt.nrnbrs[stk_value][nb];      //value of neighbor
      if(asked[test]==0)                        //check if already asked
      {
        if(spin[test]==ran_spin)                //check if spin are same
        {
          if(ran.rand()<pr)                     //double rand(); real number in [0,1]
          {
            //cout<<"inside core test\t"<<endl;
            asked[test]=1;                      //add to the asked list
            stk.push_back(test);                //add to the stack list
          }
        }
      }
    }
    
//get the size of wolff cluster
wolff_size+=stk.size();

//flip the spins of the sites in the stack
for (int i=0;i<stk.size();i++)spin[stk[i]]=1-spin[stk[i]];
  }
}


void ISING_CONF::meas_clear(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  energy=0.;
  energy_sq=0.;
  mag=0.;
  mag_sq=0.;
  wolff_size=0.;
}

void ISING_CONF::meas(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  int ener=0;
  for (int site=0;site<latt.Nsite;site++)
    for (int nbr=0;nbr<latt.nrnbrs[site].size();nbr++)
      ener+=(2*spin[site]-1)*(2*spin[latt.nrnbrs[site][nbr]]-1);
  ener=ener/2;
  
  int mag_tmp=0;
  for (int site=0;site<latt.Nsite;site++)
    mag_tmp+=2*spin[site]-1;

  energy+=((double)ener/(double)latt.Nsite);
  energy_sq+=((double)ener)*((double)ener)/((double)latt.Nsite*latt.Nsite);
  mag+= ((double) mag_tmp)/((double)latt.Nsite);
  mag_sq+=((double) mag_tmp)*((double) mag_tmp)/((double)latt.Nsite*latt.Nsite);

}

void ISING_CONF::binwrite(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  //dfout.open("data.out");
  dfout.open("ssfdata.out");
  //dfout.open("ssfdata.out",ios::app);
  dfout << p.betaJ*energy/((double)(p.Nmcs))<<" "<< p.betaJ*p.betaJ*energy_sq/((double)(p.Nmcs))<<" "<<mag/((double)(p.Nmcs))<<" "<<mag_sq/((double)(p.Nmcs))<<" "<<wolff_size/((double)p.Nmcs*latt.Nsite)<<endl;

  dfout.close();
}

