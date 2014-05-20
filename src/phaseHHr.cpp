#include "Variant.h"
#include "split.h"
#include "cdflib.hpp"
#include "pdflib.hpp"

#include <string>
#include <iostream>
#include <math.h>  
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <getopt.h>
#include "var.hpp"

using namespace std;
using namespace vcf;


void printVersion(void){
  cerr << "INFO: version 1.0.0 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu " << endl << endl;
}

void printHelp(void){

  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "     phaseHHr is a phasing stragety that maximizes haplotype homozygosity (EHH) in a 15 SNP window with 5 SNP overlap." << endl ;
  cerr << "     A stochastic search is carried out 1000 iteratoration within each window searching attempting to maximize HH.    " << endl ;
  cerr << "     Only bi-allelic sites are considered in phaseHHr.                                                                " << endl << endl ;
   
  cerr << "Output: A phased vcf file.  Note: phaseHHr does not modify the INFO field, meaning AF AN AC are no longer valid       " << endl << endl ; 
 
  cerr << "INFO usage: phaseHHr --file my.vcf --region chr1 --type PL" << endl << endl;

  cerr << "INFO required: f,file     -- a valid VCF file containing a PL or GL field                  " << endl;
  cerr << "INFO required: t,type     -- the genotype likelihood format [PL,GL]                        " << endl;
  cerr << "INFO optional: r,region   -- a region in tabix compliant format: e.g: chr1 or chr1:1-100000" << endl << endl;
  
  printVersion();

}

double EHH(string haplotypes[][2], int nhaps){

  map<string , int> hapcounts;

  for(int i = 0; i < nhaps; i++){
    hapcounts[ haplotypes[i][0] ]++;
    hapcounts[ haplotypes[i][1] ]++;
  }
  
  double sum = 0;
  double nh  = 0;

  for( map<string, int>::iterator it = hapcounts.begin(); it != hapcounts.end(); it++){
    nh  += it->second; 
    sum += r8_choose(it->second, 2);
  }

  double max = (sum /  r8_choose(nh, 2));
  
  return max;

}

void clearHaplotypes(string haplotypes[][2], int ntarget){
  for(int i= 0; i < ntarget; i++){
    haplotypes[i][0].clear();
    haplotypes[i][1].clear();
  }
}

void loadImprovement(string tmpHaplotypes[][2], string haplotypes[][2], int ntarget){

   for(int i= 0; i < ntarget; i++){
    haplotypes[i][0] = tmpHaplotypes[i][0];
    haplotypes[i][1] = tmpHaplotypes[i][1];
  }

}

void appendHaplotypes(string tmpHaplotypes[][2], string haplotypes[][2], int ntarget){ 

  int totalL = 0;

  if(tmpHaplotypes[0][0].size() > 15){
    totalL = 15;
  }
  else{
    totalL = tmpHaplotypes[0][0].size();
  }

  for(int i= 0; i < ntarget; i++){
    haplotypes[i][0].append(tmpHaplotypes[i][0].substr(5,totalL));
    haplotypes[i][1].append(tmpHaplotypes[i][1].substr(5,totalL));
  }
}

void localPhase(string haplotypes[][2], list<genotype*> & window){
 
  int ntarget =  window.front()->genoLikelihoods.size();
  double ehhmax  =  -1;
  string totalHaplotypes[ntarget][2];
  
  for(int k = 0; k < 1000; k++){   

    string tempHaplotypes[ntarget][2];

    int tlength = haplotypes[0][0].size();

    for(int nt = 0; nt < ntarget ; nt++){
      if(tlength > 0 ){
        tempHaplotypes[nt][0] = haplotypes[nt][0].substr(tlength - 5, 5);
        tempHaplotypes[nt][1] = haplotypes[nt][1].substr(tlength - 5, 5);
      }
      else{
        tempHaplotypes[nt][0] = "00000";
        tempHaplotypes[nt][1] = "00000";
      }
    }
    

    for(list< genotype* >::iterator pos = window.begin(); pos != window.end(); pos++){    
      
      int indIndex = 0;

      for(vector<int>::iterator ind = (*pos)->genoIndex.begin(); ind != (*pos)->genoIndex.end(); ind++){      
	int g = (*pos)->genoIndex[indIndex];
	if(g == -1){
	  g = rand() % 3;
	}

	//double rang  = ((double)rand() / (double)(RAND_MAX));
	//
	////cerr <<  "g:0 " << pos->genoLikelihoodsCDF[indIndex][0] << endl;
	////cerr <<  "g:1 " << pos->genoLikelihoodsCDF[indIndex][1] << endl;
	////cerr <<  "g:2 " << pos->genoLikelihoodsCDF[indIndex][2] << endl;
	//
	//if(rang < (*pos)->genoLikelihoodsCDF[indIndex][0]){
	//  g = 0;
	//}
	//else if(rang < (*pos)->genoLikelihoodsCDF[indIndex][1]){
	//  g = 1;
	//}
	//else{
	//  g = 2;
	//}	
	if(g == 0 ){
          tempHaplotypes[indIndex][0].append("0");
          tempHaplotypes[indIndex][1].append("0");
        }
	if( g == 2 ){	
	  tempHaplotypes[indIndex][0].append("1");
	  tempHaplotypes[indIndex][1].append("1");
	}
	if(g == 1){
	  double ranh  = ((double)rand() / (double)(RAND_MAX)); 
	  if(ranh < 0.5){
	    tempHaplotypes[indIndex][0].append("0");
	    tempHaplotypes[indIndex][1].append("1");
	  }
	  else{
	    tempHaplotypes[indIndex][0].append("1");
	    tempHaplotypes[indIndex][1].append("0");
	  }
	}	
	indIndex += 1;
      }
      
    }

    double ehh = EHH(tempHaplotypes, ntarget);    
    if(ehh > ehhmax){
      ehhmax = ehh;
      loadImprovement(tempHaplotypes, totalHaplotypes, ntarget);
    }
  }
  appendHaplotypes(totalHaplotypes, haplotypes, ntarget);
}

  
void printVCF(string haplotypes[][2], list< Variant > & vars, int nsamples, int len){
  
  int step = 0;
  
  for(list < Variant >::iterator it = vars.begin(); it != vars.end(); it++){
    
    // i is the individuals' index 
    for(int i = 0; i < nsamples; i++){
      
      string genotype = "0|0";
   
      //      cerr << "info: " << haplotypes[i][0].substr(len+step, 1) << haplotypes[i][1].substr(len+step, 1) << endl;
   
      while(1){
	if(haplotypes[i][0].substr(len+step, 1) == "0" && haplotypes[i][1].substr(len+step, 1) == "1"){
	  genotype = "0|1";
	  break;
	}
	if(haplotypes[i][0].substr(len+step, 1) == "1" && haplotypes[i][1].substr(len+step, 1) == "0"){
	  genotype = "1|0";
	  break;
	}
	if(haplotypes[i][0].substr(len+step, 1) == "1" && haplotypes[i][1].substr(len+step, 1) == "1"){
	  genotype = "1|1";
	  break;
	}
	break;
      }
      it->samples[it->sampleNames[i]]["GT"][0] = genotype;
    }
    cout << (*it) << endl;
    step += 1;
  }
}

int main(int argc, char** argv) {

  // set the random seed for MCMC

  srand((unsigned)time(NULL));

  // the filename

  string filename = "NA";

  // set region to scaffold

  string region = "NA"; 

  // using vcflib; thanks to Erik Garrison 

  VariantCallFile variantFile;

  // zero based index for the target and background indivudals 
  
  map<int, int> it, ib;
  
  // deltaaf is the difference of allele frequency we bother to look at 

  // ancestral state is set to zero by default

  string type = "NA";

  int counts = 0;
  
  // phased 

  int phased = 0;

    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"region"    , 1, 0, 'r'},
	{"type"      , 1, 0, 't'},
	{0,0,0,0}
      };

    int findex;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "f:r:t:hv", longopts, &findex);
	
	switch (iarg)
	  {
	  case 'h' :
	    printHelp();
	    return 0;
	  case 'v' :
	    printVersion();
	    return 0;
	  case 'r':
	    region = optarg;
	    cerr << "INFO: region has been set to: " << region << endl; 
	    break;
	  case 'f' :
	    filename = optarg;
	    cerr << "INFO: file set to: " << filename << endl;
	    break;
	  case 't':
	    type = optarg;
	    cerr << "INFO: setting genotype likelihood class to: " << "type" << endl;
	    break;
	  default:
	    break;
	  }
      }

    if(filename == "NA"){
      cerr << "FATAL: did not specify a file" << endl;
      cerr << "INFO: please use phaseHHr --help" << endl;
      return(1);
    }

    variantFile.open(filename);
   
    if(type == "NA"){
      cerr << "FATAL: must specify a genotype likelihood format : PL or GL " << endl;
      return 1;
    }
    map<string , int> okayGenotypeLikelihood;
    okayGenotypeLikelihood["GL"] = 1;
    okayGenotypeLikelihood["PL"] = 1;
    
    if(okayGenotypeLikelihood.find(type) == okayGenotypeLikelihood.end()){
      cerr << "FATAL: specified genotype likelihood is incorrectly formatted. Only use \"PL\" or \"GL\"!" << endl;
      return(1);
    }

    if (!variantFile.is_open()) {
      cerr << "FATAL: could not open VCF for reading" << endl;
      return 1;
    }
    if(region != "NA"){
      variantFile.setRegion(region);
    }
    
    variantFile.addHeaderLine("##INFO<ID=EH,Number=1,Type=Float,Description=\"Extended haplotype homozygosity at beginning of window\">");
    variantFile.addHeaderLine("##source=\"GPAT++ phaseEHHr v.1.0.0\"");
    cout << variantFile.header << endl;
    
    Variant var(variantFile);
    
    list< Variant > windowVcfData;
    list< genotype * > windowGenotypeData; 
    
    int nsamples = variantFile.sampleNames.size();
    
    if(nsamples < 2){
      cerr << "FATAL: too few individual to phase " << endl;
      return 1;
    }
    else{
      cerr << "INFO: phasing " << nsamples << " samples" << endl; 
    }
    
    string haplotypes [nsamples][2];
    
    int nsnps = 0;
    
    string currentSeqid = "NA";
    
    while (variantFile.getNextVariant(var)) {
      
      if(var.alt.size() > 1){
	continue;
      }
      
      if(currentSeqid != var.sequenceName ){
	if( windowGenotypeData.size() > 1){
	  localPhase(haplotypes, windowGenotypeData  );
	  printVCF(haplotypes,   windowVcfData, nsamples, nsnps);
	}
	clearHaplotypes(haplotypes, nsamples);
	windowVcfData.clear();
	windowGenotypeData.clear();
	currentSeqid = var.sequenceName;
	nsnps = 0;
      }
      
      map<string, map<string, vector<string> > >::iterator s     = var.samples.begin(); 
      map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();
      
      vector < map< string, vector<string> > >  total;
      
      for (; s != sEnd; s++) {	  
	
	map<string, vector<string> >& sample = s->second;
	
	total.push_back(sample);	
      }
      
      genotype * population;
      if(type == "PL"){
	population = new pl();
      }
      else{
	population = new gl();
      }
      
     population->loadPop(total, var.sequenceName, var.position);
     
     windowVcfData.push_back(var);
     
     windowGenotypeData.push_back(population);
     
     while(windowGenotypeData.size() >= 15 && !windowGenotypeData.empty()){
       localPhase(haplotypes, windowGenotypeData  );
       printVCF(haplotypes,   windowVcfData, nsamples, nsnps);
       nsnps+=15;
       while(!windowGenotypeData.empty()){
	 windowGenotypeData.pop_front();
	 windowVcfData.pop_front();
       }
     }
    }
   localPhase(haplotypes, windowGenotypeData  );
   printVCF(haplotypes,   windowVcfData, nsamples, nsnps);
   
   return 0;		    
}
