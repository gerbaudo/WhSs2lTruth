#include "WhSs2lTruth/TSelector_SusyNtuple_Truth.h"
#include "SusyNtuple/ChainHelper.h"

#include "TChain.h"
#include "Cintex/Cintex.h"

#include <iostream>
#include <string>

using namespace std;

//----------------------------------------------------------
//
// run TSelector_SusyNtuple_Truth
//
//----------------------------------------------------------
void printHelp(const char *exeName)
{
  cout<<"Usage :"<<endl
      <<exeName<<endl
      <<"\t -i inputFile.root"<<endl
      <<"\t -s samplename"<<endl
      <<endl;
}

int main(int argc, char **argv)
{

  ROOT::Cintex::Cintex::Enable();
  string sampleName;
  string inputRootFname;
  bool verbose(false);

  int optind(1);
  while ((optind < argc)) {
    if(argv[optind][0]!='-'){optind++; continue;}
    string sw = argv[optind];
    if     (sw == "-h"){ printHelp(argv[0]); return 0; }
    else if(sw == "-i"){ optind++; inputRootFname = argv[optind]; }
    else if(sw == "-s"){ optind++; sampleName     = argv[optind]; }
    else if(sw == "-v"){ verbose = true; }
    else if(argv[optind][0]=='-') cout<<"Unknown switch "<<sw<<endl;
    optind++;
  } // end if(optind<argc)
  cout<<"Using the following options:"<<endl
      <<"inputRootFname : "<<inputRootFname<<endl
      <<"sample         : "<<sampleName<<endl
      <<endl;

  TSelector_SusyNtuple_Truth analysis;
  if(verbose) analysis.setDebug(1);

  TChain chain("susyNt");
  ChainHelper::addFile(&chain, inputRootFname);
  Long64_t nEntries = chain.GetEntries();

  chain.Process(&analysis, sampleName.c_str()); 
  cout<<"Total entries:   "<<nEntries<<endl;
  return 0;
}
