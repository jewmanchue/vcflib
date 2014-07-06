#include "Variant.h"
#include "split.h"
#include "cdflib.hpp"
#include "pdflib.hpp"
#include "var.hpp"

#include <string>
#include <iostream>
#include <math.h>  
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;
using namespace vcf;

int main(int argc, char** argv) {

  string filename = "/home/zkronenb/jDavis_CQF_8_libraries.Unified_Genotyper.vcf.gz";

  VariantCallFile variantFile;
  
  variantFile.open(filename);
  
  //  variantFile.setRegion("scaffold_39");
  
  Variant var(variantFile);

  vector< pooled *> ncData;
  vector< pooled *> scData;
  vector< pooled *> rsData;
  
  while (variantFile.getNextVariant(var)) {
    
    if(var.alt.size() > 1){
      continue;
    }
    
    map<string, map<string, vector<string> > >::iterator s     = var.samples.begin(); 
    map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();
    
    vector < map< string, vector<string> > > northCarolina, southCarolina, resistant, nonResistant;
    
    northCarolina.push_back(var.samples["jDavis_CQF_library_NC_parental"]);
    southCarolina.push_back(var.samples["jDavis_CQF_library_SC_parental"]);
    
    nonResistant.push_back(var.samples["jDavis_CQF_library_4"]);
    nonResistant.push_back(var.samples["jDavis_CQF_library_5"]);
    nonResistant.push_back(var.samples["jDavis_CQF_library_6"]);

    resistant.push_back(var.samples["jDavis_CQF_library_1"]);
    resistant.push_back(var.samples["jDavis_CQF_library_2"]);
    resistant.push_back(var.samples["jDavis_CQF_library_3"]);

    pooled * nc;
    pooled * sc;
    pooled * rs;
    pooled * ns;

    nc = new pooled;
    sc = new pooled;
    rs = new pooled;
    ns = new pooled;

    nc->loadPop(northCarolina, var.sequenceName, var.position);
    sc->loadPop(southCarolina, var.sequenceName, var.position);
    rs->loadPop(resistant,     var.sequenceName, var.position);
    ns->loadPop(nonResistant,  var.sequenceName, var.position);
 
    // cout << var.sequenceName << "\t" << var.position << "\t" << rs->af << "\t" << ns->af << endl;

    if(nc->af == -1 || sc->af == -1 || rs->af == -1 || ns->af == -1){
      continue ; 
    }
    
    if(rs->npop < 2){
      continue;
    }
    if(ns->npop < 2){
      continue;
    }
    
    rs->estimatePosterior();
    ns->estimatePosterior();
    
    

    // if library variance is above 0.1 allele frequency

   if(rs->pvar > 0.1){
     continue;
   }
   if(ns->pvar > 0.1){
     continue;
   }


    double afRs = rs->af;
    double afNs = ns->af;
    double afSc = sc->af;
    double afNc = nc->af;

    
    // the two libraries need to have an allele frequency difference greater than 0.5
    
    if(abs(afRs - afNs) < 0.4){
      continue;
    }

    // if the resistant library has an allele frequency above 0.5 it must also have an af above 0.98

    if(afRs > 0.5 && afRs < 0.99){
      continue;
    }

    // if the resistand library has a frequency below 0.5 is must also have an af below 0.02

    if(afRs < 0.5 && afRs > 0.01){
      continue;
    }
    
    // must be in the non-resistant library

    if(afNs > 0.99 || afNs < 0.01){
      continue;
    }


    // the frequencies need to be different between the north Carolina and south Carolina

    if(abs(afNc - afSc) < 0.05){
      continue;
    }
    
    // allele must be present in both resistant and NC

    if(afRs > 0.99 && afNc < 0.01){
      continue;
    }
    if(afRs < 0.01 && afNc > 0.99){
      continue;
    }

    

    cout << var.sequenceName << "\t" << var.position << "\t" << afRs << "\t" << afNc << "\t" << afNs << "\t" << afSc << "\t" <<  endl;
    
    delete nc;
    delete sc;
    delete rs;
    delete ns;    

  }
  return 0;		    
}
