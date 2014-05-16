// not to complicate the issue but I need a different variant object to handle populations. 

#ifndef __VAR_H
#define __VAR_H

#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <stdio.h>      
#include <stdlib.h>
#include "split.h"

using namespace std;

class zvar{
public:

  int npop;

  string seqid;
  long int pos;

  double nalt ;
  double nref ;
  double af   ;

  double alpha; 
  double beta ;

  virtual void loadPop(vector< map< string, vector<string> > >& group, string seqid, long int position) = 0;
  virtual void estimatePosterior() = 0 ;

};

class genotype : public zvar {
  
public:

  double nhomr;
  double nhoma;
  double nhet ;
  double ngeno;
  double fis  ;
  
  vector<int> genoIndex;
  vector< vector < double > > genoLikelihoods;
  vector< vector < double > > genoLikelihoodsCDF;

  virtual double unphred(map< string, vector<string> > & geno, int index) = 0; 
  virtual void loadPop(vector< map< string, vector<string> > >& group, string seqid, long int position);
  void estimatePosterior();

};

class pooled : public zvar{
public:


  double ntot  ;
  double afsum ; 
    
  vector<double> nalts;
  vector<double> nrefs;
  vector<double> afs  ; 

  void loadPop(vector< map< string, vector<string> > >& group, string seqid, long int position);
  void estimatePosterior();
  double bound(double v);

  pooled(void);
  
};

class gl : public genotype{
public:
  gl(void);
  double unphred(map< string, vector<string> > & geno, int index);
};

class pl : public genotype{
public:
  pl(void);
  double unphred(map< string, vector<string> > & geno, int index);
}; 

#endif 
