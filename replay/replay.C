#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"
//#pragma cling load("libsbs.so");
//#include "SBSGEMStand.h"
//#include "SBSBigBite.h"

// Replay script for SLAC beam test, 2018

void replay(int run_num){
  //  Steering script for Hall A analyzer demo

    gSystem->Load("libsbs.so");

    // Set up the equipment to be analyzed.    
    SBSBigBite   *sbs = new SBSBigBite("sbs", "Generic apparatus");
    SBSGEMStand *gems = new SBSGEMStand("gems", "Collection of GEMs in stand");
    THaSBUScint *det = new THaSBUScint("sbuscint", "sbu scint");
    sbs->AddDetector(gems);
    sbs->AddDetector(det);

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Apparatus's and PhysicsModules,
  // and executes the output routines.
  THaAnalyzer* analyzer = new THaAnalyzer;
  gHaApps->Add(sbs);

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;
  
  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  THaRun* run = new THaRun(Form("/home/yetao/workarea/slac_beam/raw_data/run_%d.dat",run_num));
  run->SetLastEvent(-1);
  run->SetDataRequired(0);
  run->SetDate(TDatime());

  analyzer->SetVerbosity(3);
  
  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile(Form("/home/yetao/workarea/slac_beam/rootfiles/run_%d.root",run_num));
  // File to record cuts accounting information
  analyzer->SetSummaryFile("summary_example.log"); // optional
  //analyzer->SetCompressionLevel(0); // turn off compression
  analyzer->Process(run);     // start the actual analysis
}
