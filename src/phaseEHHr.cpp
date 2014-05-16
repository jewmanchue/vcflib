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

void loadIndices(map<int, int> & index, string set){
  
  vector<string>  indviduals = split(set, ",");

  vector<string>::iterator it = indviduals.begin();
  
  for(; it != indviduals.end(); it++){
    index[ atoi( (*it).c_str() ) ] = 1;
  }
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


void printHaplotypes(string haps[][2], int ntargets){
  for(int snp = 0; snp < haps[0][1].length(); snp++){
    for(int ind = 0; ind < ntargets; ind++){
      cout << haps[ind][0].substr(snp , 1) << "\t";
      cout << haps[ind][1].substr(snp , 1) << "\t";
    }
    cout << endl;
  }
}

void loadImprovement(string tmpHaplotypes[][2], string haplotypes[][2], int ntarget){

  for(int i= 0; i < ntarget; i++){
    haplotypes[i][0] = tmpHaplotypes[i][0];
    haplotypes[i][1] = tmpHaplotypes[i][1];
  }

}

void appendHaplotypes(string tmpHaplotypes[][2], string haplotypes[][2], int ntarget){ 
  for(int i= 0; i < ntarget; i++){
    haplotypes[i][0].append(tmpHaplotypes[i][0].substr(5,15));
    haplotypes[i][1].append(tmpHaplotypes[i][1].substr(5,15));
  }
}

void localPhase(string haplotypes[][2], list<pl> & window){
 
  int ntarget =  window.front().genoLikelihoods.size();
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
    

    for(list<pl>::iterator pos = window.begin(); pos != window.end(); pos++){    
      
      int indIndex = 0;

      for(vector<int>::iterator ind = pos->genoIndex.begin(); ind != pos->genoIndex.end(); ind++){      
	int g = pos->genoIndex[indIndex];
	double rang  = ((double)rand() / (double)(RAND_MAX));

	//cerr <<  "g:0 " << pos->genoLikelihoodsCDF[indIndex][0] << endl;
	//cerr <<  "g:1 " << pos->genoLikelihoodsCDF[indIndex][1] << endl;
	//cerr <<  "g:2 " << pos->genoLikelihoodsCDF[indIndex][2] << endl;

	if(rang < pos->genoLikelihoodsCDF[indIndex][0]){
	  g = 0;
	}
	else if(rang < pos->genoLikelihoodsCDF[indIndex][1]){
	  g = 1;
	}
	else{
	  g = 2;
	}	
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

  string mut = "1";

  int counts = 0;
  
  // phased 

  int phased = 0;

    const struct option longopts[] = 
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"background", 1, 0, 'b'},
	{"deltaaf"   , 1, 0, 'd'},
	{"region"    , 1, 0, 'r'},
	{"mutation"  , 1, 0, 'm'},
	{"phased"    , 1, 0, 'p'},
	{0,0,0,0}
      };

    int findex;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "p:m:r:d:t:b:f:hv", longopts, &findex);
	
	switch (iarg)
	  {
	  case 'h':
	    cerr << endl << endl;
	    cerr << "INFO: help" << endl;
	    cerr << "INFO: description:" << endl;
            cerr << "     pl-XPEHH estimates haplotype decay between the target and background populations.  SNVs are integrated                           " << endl;
	    cerr << "     until EHH in the target is less than 0.05.  The reported score is the itegrated EHH (target) / integrated EHH (background).	   " << endl;
	    cerr << "     pl-XPEHH does NOT integrate over genetic distance, as genetic maps are not availible for most non-model organisms. 		   " << endl;
	    cerr << "     pl-XPEHH phases genotypes, imuputes missing genotypes, and changes poor quality genotypes. Phasing is done in a sliding window   " << endl;
	    cerr << "     with a stochastic search, therefore, every time pl-XPEHH is run it will generate slightly different results.                     " << endl;

	    cerr << "Output : 4 columns :     "    << endl;
	    cerr << "     1. seqid            "    << endl;
	    cerr << "     2. position         "    << endl;
	    cerr << "     3. xp-ehh           "    << endl;
	    cerr << "     4. iHS              "    << endl  << endl;

	    cerr << "INFO: pl-XPEHH  --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --file my.vcf --deltaaf 0.1 --ancestral 0        " << endl;
	    cerr << endl;
	    cerr << "INFO: required: r,region     -- a genomice range to calculate pl-XPEHH on in the format : \"seqid:start-end]\" or \"seqid\" " << endl;
	    cerr << "INFO: required: t,target     -- a zero base comma seperated list of target individuals corrisponding to VCF columns        " << endl;
	    cerr << "INFO: required: b,background -- a zero base comma seperated list of background individuals corrisponding to VCF columns    " << endl;
	    cerr << "INFO: required: f,file a     -- proper formatted VCF.  the FORMAT field MUST contain \"PL\" if option phased == 0           " << endl; 
	    cerr << "INFO: optional: m,mutation   -- which state is derived in vcf [0,1] default is 1                                            " << endl;
	    cerr << "INFO: optional: p,phased     -- phasing flag [0,1] 0 = phase vcf, 1 = vcf is already phased                                 " << endl;
	    cerr << endl; 
	    cerr << "INFO: version 1.0.1 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu " << endl;
	    cerr << endl << endl;
	    return 0;
	  case 'v':
	    cerr << endl << endl;
	    cerr << "INFO: version 1.0.1 ; date: April 2014 ; author: Zev Kronenberg; email : zev.kronenberg@utah.edu "  << endl;
	    return 0;
	  case 'p':
	    phased = atoi(optarg);
	    cerr << "INFO: setting phase to: " << phased << endl;
	    break;
	  case 'm':
	    mut = optarg;
	    cerr << "INFO: derived state set to " << mut << endl;
	    break;
	  case 't':
	    loadIndices(it, optarg);
	    cerr << "INFO: there are " << it.size() << " individuals in the target" << endl;
	    cerr << "INFO: target ids: " << optarg << endl;
	    break;
	  case 'b':
	    loadIndices(ib, optarg);
	    cerr << "INFO: there are " << ib.size() << " individuals in the background" << endl;
	    cerr << "INFO: background ids: " << optarg << endl;
	    break;
	  case 'f':
	    cerr << "INFO: file: " << optarg  <<  endl;
	    filename = optarg;
	    break;
	  case 'r':
            cerr << "INFO: set seqid region to : " << optarg << endl;
	    region = optarg; 
	    break;
	  default:
	    break;
	  }

      }

    if(filename == "NA"){
      cerr << "FATAL: did not specify a file" << endl;
      cerr << "INFO: please use pl-XPEHH --help" << endl;
      return(1);
    }

    variantFile.open(filename);
    
    if(region == "NA"){
      cerr << "FATAL: did not specify a region"  << endl;
      cerr << "INFO: please use pl-XPEHH --help" << endl;
    }

   if(region != "NA"){
     variantFile.setRegion(region); 
   }
   
   if (!variantFile.is_open()) {
     return 1;
   }
   
    variantFile.addHeaderLine("##INFO<ID=EH,Number=1,Type=Float,Description=\"Extended haplotype homozygosity at beginning of window\">");
    variantFile.addHeaderLine("##source=\"GPAT++ phaseEHHr v.1.0.0\"");
    cout << variantFile.header << endl;
    
   Variant var(variantFile);
   
   list< Variant > windowVcfData;
   list< pl > windowGenotypeData; 
   
   int nsamples = variantFile.sampleNames.size();

   string haplotypes [nsamples][2];
   
   int nsnps = 0;
   
   while (variantFile.getNextVariant(var)) {

     if(var.alt.size() > 1){
       continue;
     }

     map<string, map<string, vector<string> > >::iterator s     = var.samples.begin(); 
     map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();
           
     vector < map< string, vector<string> > >  total;
      
     for (; s != sEnd; s++) {	  
       
       map<string, vector<string> >& sample = s->second;
       
       total.push_back(sample);	
     }

     pl population;

     population.loadPop(total, var.sequenceName, var.position);

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
