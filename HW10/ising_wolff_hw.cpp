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
  vector<vector<int> > nrnbrs;//vector of nearest neighbors
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

    //Add a function for Wolff sinlge cluster algorithm
    void wolff_update(const PARAMS&, const LATTICE&,MTRand&);
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
    double wolff_size;
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
    ising.wolff_update(p,latt,ran1);
    //ising.sweep(p,latt,ran1);

  //PRODUCTION
  for (int bin=0;bin<p.Nbin;bin++)
  {
    cout<<"bin\t"<<bin<<endl;
    ising.meas_clear(p,latt,ran1);
    for(int mcs=0;mcs<p.Nmcs;mcs++)
    {
      ising.wolff_update(p,latt,ran1);
      //ising.sweep(p,latt,ran1);
      ising.meas(p,latt,ran1);
    }
    ising.binwrite(p,latt,ran1);
  }
}

PARAMS::PARAMS(){
  ifstream pfin;
  pfin.open("param.dat");  //file stores input parameters
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
      vector <int> nrn; 		//vector to hold nearest neighbour elements
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
  //cout<<"inside ising_conf"<<endl;

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

  }
  else
  {cout <<"NEED TO CODE ALL LATTICE OPTIONS"<<endl;}

  //open the file to write output
  dfout.open("outputdata.out");
}

ISING_CONF::~ISING_CONF()
{
  dfout.close();
}
//writes the configuration in a file
void ISING_CONF::conf_write(const PARAMS& p, const LATTICE& latt)
{
  //cout<<"inside conf_write"<<endl;
  ofstream outputfile;
  string filename="isingconf.dat"; //stores configuration
  outputfile.open(filename);
  if(outputfile.is_open()){
    for(int i=0;i<latt.Nsite;i++)outputfile<<spin[i]<<endl;
    outputfile.close();
  }
  else cout<<filename<<"\t doesnot exist"<<endl;

}

//flips the single spin in the entire lattice based on the probability
void ISING_CONF::sweep(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  for (int site=0;site<latt.Nsite;site++){ //go through all the sites 
    int neighborsum=0;          //sum of neighbor spins
    int spin_new=1-spin[site];  //flip the spin
    for(int nbr=0;nbr<latt.nrnbrs[site].size();nbr++){//spin sum of neighbors
      neighborsum+=spin[latt.nrnbrs[site][nbr]];
    } 
    //randDblExc()-> real number in (0,1)
    //accept if ratio>1

    if(wght_tbl[neighborsum][spin_new][spin[site]]>ran1.randDblExc())spin[site]=spin_new;
  }
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
  double pr=1.0-exp(-2.0*p.beta);
  
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

//resets the values of observable energy, magnetization and their sqare to 0
void ISING_CONF::meas_clear(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  //cout<<"inside meas_clear"<<endl;
  // reset the values to 0
  energy=0.;
  energy_sq=0.;
  mag=0;
  mag_sq=0;
  wolff_size=0.;
}

//calcualte the value of observables
void ISING_CONF::meas(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  //cout<<"inside meas"<<endl;
  //magnetization
  //magnetization=sum(s_i) s_i=(+/-)1
  //our spin is in the form of sigma=0,1
  //s=2*sigma-1
  double magnetization=0.;
  for (int site=0;site<latt.Nsite;site++){
    magnetization+=2.*spin[site]-1;
  }
  //average magnetization
  magnetization=1.0*magnetization/latt.Nsite;

  //energy
  //energy=-sum(s_i*s_j) <ij> are the nearest neighbors
  double ener=0.;
  for(int site=0;site<latt.Nsite;site++)
  {
    for (int i=0;i<latt.nrnbrs[site].size();i++){
      ener+=-(2.0*spin[site]-1)*(2.0*spin[latt.nrnbrs[site][i]]-1);
    }
  }
  ener=ener/2.;//counted twice in loop
  ener=1.0*ener/latt.Nsite; //average energy

  //calulation of energy, magnetization and squares
  mag+=magnetization;
  mag_sq+=magnetization*magnetization;
  energy+=ener;
  energy_sq+=ener*ener;
}

void ISING_CONF::binwrite(const PARAMS& p, const LATTICE& latt, MTRand& ran1)
{
  //write the output in the output file
  //dfout<<1.0*energy/p.Nmcs<<"\t"<<1.0*energy_sq/p.Nmcs<<"\t"<<1.0*mag/p.Nmcs<<"\t"<<1.0*mag_sq/p.Nmcs<<endl;
  dfout.precision(5);
  dfout<<fixed<<1.0*energy/double(p.Nmcs)<<setw(15)<<1.0*energy_sq/double(p.Nmcs)<<setw(15)<<1.0*mag/double(p.Nmcs)<<setw(15)<<1.0*mag_sq/double(p.Nmcs)<<setw(15)<<1.0*wolff_size/double(p.Nmcs)<<endl;
 }

