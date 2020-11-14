#include <fstream> 
#include <iostream> 
#include <math.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include "MersenneTwister.h"

using namespace std;

class PARAMS
{
public:
  int Nlin; // linear size of lattice
  double pr; // probability for a site
  double Nclust; // number of clusters in a bin
  double Nbin; // number of bins of data to output
  int SEED; // seed for mersenne twister
  string latt_; // which lattice 
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

};

class CLUSTER
{
public:
  CLUSTER(const PARAMS&, const LATTICE&);//constructor
  void grow(const LATTICE&, MTRand&);    //grows based on the neighbor, stores the size of the neighbor vector 
  void meas_clear(const LATTICE&);	 //clears the sum of avg_size
  void meas(const LATTICE&);		 //stores the sum of sizes in avg_size
  void binwrite(const PARAMS&, const LATTICE&); //gives the mean of sizes avg_size/Ncluster
  void print(const LATTICE& latt, int index);
  ~CLUSTER();// destructor
private:
  int size; //size of the neighbor vector
  vector <int> conf;
  vector <int> stack;
  double pr;
  int stck_pnt,stck_end;
  double avg_size;//sum of size
  ofstream dfout;
 
};

int main(void)
{
  PARAMS p;
  LATTICE latt(p);
  CLUSTER cluster(p,latt);
  MTRand ran1(p.SEED);
  //MTRand ran2;

  //cout<<"ran1\t"<<ran1<<endl;
 // cout<<"ran1\t"<<ran1.rand()<<endl;
  //latt.print();
  //cluster.print(latt,1);
 //cluster.grow(latt,ran1);
  
  for (int bin=0;bin<p.Nbin;bin++)
    {
      cluster.meas_clear(latt);
      for(int clust=0;clust<p.Nclust;clust++)
	{
	  cluster.grow(latt,ran1);
	  cluster.meas(latt);

	}
      cluster.binwrite(p,latt);
    }
  cluster.print(latt,1);
}

PARAMS::PARAMS(){
  //initializes commonly used parameters from a file
  ifstream pfin;
  pfin.open("param.dat");  
  if (pfin.is_open()) { 
    pfin >> Nlin;
    pfin >> pr;
    pfin >> Nclust;
    pfin >> Nbin;
    pfin >> SEED;
    pfin >> latt_;
  }
  else
    {cout << "No input file to read ... exiting!"<<endl;exit(1);}
  pfin.close();
  // print out all parameters for record
  cout << "--- Parameters at input for percolation problem ---"<<endl; 
  cout <<"Nlin = "<<Nlin<<"; prob. of site = "<<pr<<endl;
  cout <<"Number of clusters in a bin = "<<Nclust<<"; Number of bins = "<<Nbin<<endl;
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


CLUSTER::CLUSTER(const PARAMS& p, const LATTICE& latt)
{
  conf.resize(latt.Nsite);
  stack.resize(latt.Nsite);
  pr=p.pr;// store prob in a private member of cluster
  dfout.open("data.out");
}

CLUSTER::~CLUSTER()
{
  dfout.close();
}

void CLUSTER::grow(const LATTICE& latt, MTRand& ran1)
{
  // COMPLETE HERE !!
  // conf vector of size Nsite each 0, 0->not asked, 1->asked joined, 2->asked not joined
  // stack of size Nsite, if site joins the cluster add 1 
  //stk-pnt index of the members who are asked and joined
  //stck-end->size of number of elements in the stack
  //1conf[site]=0, stck-end=0
  //2go to the random site and 'site', conf[site]=1, add it 
  //3to the stack
  //4stck-pnt=0, stack-end=1
  //5check neighbrs of site at stck-pnt
  //6If not asked before ask them if they want to joint the cluster based on the pr
  //update conf, add the new site to stack if it joined,update stack-end
  //increase the stak-pt
	
	//Initialize the conf and stack
	fill(conf.begin(),conf.end(),0);  //fill each element of conf with 0
	stack={};                         //initilization to get rid of junk addresses
	stck_end=0;                       //0 to stck_end 
	stck_pnt=0;                       //tracks the each value in the stack

	int choice=ran1.randInt(latt.Nsite-1); //random pick 0 to Nsite-1
	//int choice=int(ran1.rand()*latt.Nsite); //random pick 0 to Nsite-1
	//cout<<"random pick\t"<<choice<<endl;
	stack.push_back(choice);           //add a choice to the stack
	conf[choice]=1;                    //Update conf
	stck_end=1;                        //Update stak_end
	
	//for (int i=0;i<conf.size();i++)cout<<"conf "<<i<<"\t"<<conf[i]<<endl;
	//for (int i=0;i<stack.size();i++)cout<<"stack "<<i<<"\t"<<stack[i]<<endl;
	
	while(stck_pnt<stck_end)
	{
	        //cout<<"========================"<<endl;
		int stack_value=stack[stck_pnt];  //accss value at stak_pnt
		vector <int> neighbors=latt.nrnbrs.at(stack_value); //neighbor vector for stack value
		
		//cout<<"neighbors for \t"<<stack_value<<"->\t";
		//for (int i=0;i<neighbors.size();i++)cout<<neighbors[i]<<'\t';
		//cout<<"\n.........................."<<endl;
		//cout<<endl;
		//for(int i=0;i<latt.nrnbrs.at(stck_pnt).size();i++) //running through all the neighbours of choice
		for(int i=0;i<neighbors.size();i++)
		{
			//int check=latt.nrnbrs.at(stck_pnt).at(i);
			//cout<<"This number in neighbour vector\t"<<neighbors[i]<<endl;
			int check=neighbors[i];
			//cout<<"value being checked...\t"<<check<<endl;
			if(conf[check]==0) //not alreay asked
			{
			//	cout<<check<<"\t not already checked\t"<<endl;
				double prob=ran1.rand();//random prob value
			//	cout<<"ran1\t"<<prob<<endl;
				if(prob>pr)
				{
				//	cout<<"refused\t"<<endl;
					conf[check]=2;//refused
				}
				else
				{
				//	cout<<"accepted\t"<<endl;
					conf[check]=1;//accepted
					stack.push_back(check);
					stck_end=stck_end+1; //increase the length of stck
				}	
			}
		}
	stck_pnt=stck_pnt+1;  //update the stck_pnt
	//for(int i=0;i<stack.size();i++){cout<<"stack element inside "<<i<<"\t"<<stack[i]<<endl;}
	}
	//cout<<"stack elements:\t";
	//for(int i=0;i<stack.size();i++){cout<<stack[i]<<"\t";}
	//cout<<endl;
	//store the size of the lattice in the varialbe size
	size=stack.size();
	//cout<<"size\t"<<size<<endl;
	//cout<<"stack_end\t"<<stck_end<<endl;
}

void CLUSTER::print(const LATTICE& latt, int index)
{
  stringstream ss;
  string file_name;
  ss<<index<<".clust";
  file_name=ss.str();

  ofstream clout;
  clout.open(file_name.c_str());
  clout <<"#"<<latt.Lx<<" x "<<latt.Ly<<endl;
 
  for (int y=0;y<latt.Ly;y++)
    {
      for (int x=0;x<latt.Lx;x++)
	clout<<conf[x+y*latt.Lx]<<" ";
      clout<<endl;
    }
      
  clout.close();
}

void CLUSTER::meas(const LATTICE& latt)
{
  avg_size+=(double)size;
}


void CLUSTER::meas_clear(const LATTICE& latt)
{
  avg_size=0.;
}


void CLUSTER::binwrite(const PARAMS& p, const LATTICE& latt)
{
  dfout << avg_size/((double)p.Nclust)<<endl;
}
