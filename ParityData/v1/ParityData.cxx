//*-- Author :    Bob Michaels,  Aug 2009

//////////////////////////////////////////////////////////////////////////
//
// ParityData
//
// HAPPEX data from spectrometer DAQ.
// 1. Cavity XYQ Monitors
// 2. HAPPEX electron detectors, alignment.
// 3. UMass quartz and scintillator
// 4. Kinematics (Q^2, etc).
//
// Minimal alignment procedure.  output.def :
// a) Alignment histograms (L-arm example):
//   Cut onetrkL L.tr.n==1
//   Cut happexL L.tr.n==1&&P.hapadcl1>1800
//   TH2F Lali1 'Left arm alignment' L.tr.x L.tr.y 100 -0.8 0.8 100 -0.2 0.2  onetrkL
//   TH2F Lali2 'Left arm alignment (HAP det)' L.tr.x L.tr.y 100 -0.8 0.8 100 -0.2 0.2  happexL
// b) Find detector using field too low by ~4%. Whole
//    detector illuminated by rad. tail.
//    --> Defines "box" of detector.
// c) Put field at 0%, check if tracks fit in "box".
//
//////////////////////////////////////////////////////////////////////////

//#define WITH_DEBUG 1

//#define DUMP_SCALER

#include "ParityData.h"
#include "THaVarList.h"
#include "THaVar.h"
#include "THaGlobals.h"
#include "THaEvData.h"
#include "THaDetMap.h"
// scalers are not done this way anymore
//#include "THaScalerGroup.h"
//#include "THaScaler.h"
#include "THaAnalyzer.h"
#include "THaString.h"
#include "TMath.h"
#include "TDatime.h"
#include "TH1.h"
#include "VarDef.h"
#include <fstream>
#include <iostream>

#define NCAV 6
#define NSTR 12

using namespace std;
using namespace THaString;

typedef vector<PdataLoc*>::iterator Iter_t;

//_____________________________________________________________________________
ParityData::ParityData( const char* name, const char* descript ) : 
  THaApparatus( name, descript )
{
  lscaler = 0;
  rscaler = 0;
  modern_times = 1;  // 0 = 2008, 1 = now
  trigcnt = new Int_t[12];
  cped = new Double_t[NCAV];
  cfb  = new Double_t[NCAV];
  sped = new Double_t[NSTR];
  sfb  = new Double_t[NSTR];
  for (int i = 0; i < 12; i++) trigcnt[i]=0;
  for (int i = 0; i < NCAV; i++) {
        cped[i]=0;
        cfb[i]=0;
  }
  for (int i = 0; i < NSTR; i++) {
        sped[i]=0;
        sfb[i]=0;
  }
  Clear();
}

//_____________________________________________________________________________
ParityData::~ParityData()
{
  if (trigcnt) delete [] trigcnt;
  if (cped) delete [] cped;
  if (cfb) delete [] cfb; 
  if (sped) delete [] sped;
  if (sfb) delete [] sfb;
  SetupParData( NULL, kDelete ); 
  for( Iter_t p = fWordLoc.begin();  p != fWordLoc.end(); p++ )  delete *p;
  for( Iter_t p = fCrateLoc.begin(); p != fCrateLoc.end(); p++ ) delete *p;
}

//_____________________________________________________________________________
void ParityData::Clear( Option_t* opt ) 
{
  evtypebits = 0;
  evtype     = 0;
  s0La = 0; s0Lb = 0; s0Ra = 0;  s0Rb = 0;
  upQadcL =0; upQadcR = 0; // Ahmed.
  loQadcL =0; loQadcR = 0; // Ahmed.
  upQtdcL =0; upQtdcR = 0; // Ahmed.
  loQtdcL =0; loQtdcR = 0; // Ahmed.
  atQadcL =0; atQadcR = 0; // Ahmed.
  atQtdcL =0; atQtdcR = 0; // Ahmed.
  upSadcL =0; upSadcR = 0; // Ahmed.
  upStdcL =0; upStdcR = 0; // Ahmed.
  loSadcL =0; loSadcR = 0; // Ahmed.
  loStdcL =0; loStdcR = 0; // Ahmed.
  amp1=0; // Ahmed.
  bamp=0; // Ahmed.
  hapadcL = 0; hapadcR = 0;
  haptdcL = 0; haptdcR = 0;
  Aquartz = 0;  Tquartz = 0;
  Apscint = 0;  Tpscint = 0;
  for (int i = 0; i < NCAV; i++) cfb[i]=0;
  for (int i = 0; i < NSTR; i++) sfb[i]=0;
  xcav4a = 0; ycav4a = 0; qcav4a = 0;
  xcav4b = 0; ycav4b = 0; qcav4b = 0;
  xstrip4a = 0; ystrip4a = 0;
  xstrip4b = 0; ystrip4b = 0;
  for( Iter_t p = fWordLoc.begin();  p != fWordLoc.end(); p++)  (*p)->Clear();
  for( Iter_t p = fCrateLoc.begin(); p != fCrateLoc.end(); p++) (*p)->Clear();
}

//_____________________________________________________________________________
Int_t ParityData::SetupParData( const TDatime* run_time, EMode mode )
{

  Int_t retval = 0;

  RVarDef vars[] = {
    { "upQadcL",   "PREX Upper Quartz ADC LHRS",  "upQadcL" }, // Ahmed.
    { "upQadcR",   "PREX Upper Quartz ADC RHRS",  "upQadcR" }, // Ahmed.
    { "loQadcL",   "PREX Lower Quartz ADC LHRS",  "loQadcL" }, // Ahmed.
    { "loQadcR",   "PREX Lower Quartz ADC RHRS",  "loQadcR" }, // Ahmed.
    { "upQtdcL",   "PREX Upper Quartz TDC LHRS",  "upQtdcL" }, // Ahmed.
    { "upQtdcR",   "PREX Upper Quartz TDC RHRS",  "upQtdcR" }, // Ahmed.
    { "loQtdcL",   "PREX Lower Quartz TDC LHRS",  "loQtdcL" }, // Ahmed.
    { "loQtdcR",   "PREX Lower Quartz TDC RHRS",  "loQtdcR" }, // Ahmed.
    { "atQadcL",   "PREX A_T Quartz ADC LHRS",  "atQadcL" },   // Ahmed.
    { "atQadcR",   "PREX A_T Quartz ADC RHRS",  "atQadcR" },   // Ahmed.
    { "atQtdcL",   "PREX A_T Quartz TDC LHRS",  "atQtdcL" },   // Ahmed.
    { "atQtdcR",   "PREX A_T Quartz TDC RHRS",  "atQtdcR" },   // Ahmed.
    { "upSadcL",   "PREX Upper Scint ADC LHRS",   "upSadcL" }, // Ahmed.
    { "upSadcR",   "PREX Upper Scint ADC RHRS",   "upSadcR" }, // Ahmed.
    { "loSadcL",   "PREX Lower Scint ADC LHRS",   "loSadcL" }, // Ahmed.
    { "loSadcR",   "PREX Lower Scint ADC RHRS",   "loSadcR" }, // Ahmed.
    { "upStdcL",   "PREX Upper Scint TDC LHRS",   "upStdcL" }, // Ahmed.
    { "upStdcR",   "PREX Upper Scint TDC RHRS",   "upStdcR" }, // Ahmed.
    { "loStdcL",   "PREX Lower Scint TDC LHRS",   "loStdcL" }, // Ahmed.
    { "loStdcR",   "PREX Lower Scint TDC RHRS",   "loStdcR" }, // Ahmed.
    { "amp1",   "PREX Upper Quartz ADC RHRS",  "amp1" },       // Ahmed.
    { "bamp",   "PREX Upper Quartz ADC RHRS",  "bamp" },       // Ahmed
    { "evtypebits", "event type bit pattern",   "evtypebits" },  
    { "evtype",     "event type from bit pattern", "evtype" },  
    { "hapadcL",   "HAPPEX ADC Left HRS",     "hapadcL" },  
    { "hapadcR",   "HAPPEX ADC Right HRS",    "hapadcR" },  
    { "haptdcL",   "HAPPEX TDC Left HRS",     "haptdcL" },  
    { "haptdcR",   "HAPPEX TDC Right HRS",    "haptdcR" },  
    { "s0La",      "S0 Left HRS PMT A",       "s0La" },
    { "s0Lb",      "S0 Left HRS PMT B",       "s0Lb" },
    { "s0Ra",      "S0 Right HRS PMT A",      "s0Ra" },
    { "s0Rb",      "S0 Right HRS PMT B",      "s0Rb" },
    { "Aquartz",   "PREX Quartz ADC",  "Aquartz" },
    { "Tquartz",   "PREX Quartz TDC",  "Tquartz" },
    { "Apscint",   "PREX Scint ADC",   "Apscint" },
    { "Tpscint",   "PREX Scint TDC",   "Tpscint" },
    { "xcav4a",     "X Cavity 4A",             "xcav4a" },  
    { "ycav4a",     "Y Cavity 4A",             "ycav4a" },  
    { "qcav4a",     "Q Cavity 4A",             "qcav4a" },  
    { "xcav4b",     "X Cavity 4B",             "xcav4b" },  
    { "ycav4b",     "Y Cavity 4B",             "ycav4b" },  
    { "qcav4b",     "Q Cavity 4B",             "qcav4b" },  
    { "xstrip4a",   "X Stripline 4A",          "xstrip4a" },  
    { "ystrip4a",   "Y Stripline 4A",          "ystrip4a" },  
    { "xstrip4b",   "X Stripline 4B",          "xstrip4b" },  
    { "ystrip4b",   "Y Stripline 4B",          "ystrip4b" },  
    { "q2L",   "Qsq on Left arm",             "q2L" },  
    { "mmL",   "Missing mass on Left arm",    "mmL" },  
    { "q2R",   "Qsq on Right arm",            "q2R" },  
    { "mmR",   "Missing mass on Right arm",   "mmR" },  
    { 0 }
  };

  Int_t nvar = sizeof(vars)/sizeof(RVarDef);

  if( mode != kDefine || !fIsSetup )
    retval = DefineVarsFromList( vars, mode );

  fIsSetup = ( mode == kDefine );

  if( mode != kDefine ) {   // cleanup the dynamically built list
    for (unsigned int i=0; i<fCrateLoc.size(); i++)
      DefineChannel(fCrateLoc[i],mode);

    for (unsigned int i=0; i<fWordLoc.size(); i++)
      DefineChannel(fWordLoc[i],mode);

    return retval;
  }
  

  fCrateLoc.clear();   
  fWordLoc.clear();   

  ifstream pardatafile;

  const char* const here = "SetupParData()";
  const char* name = GetDBFileName();

  TDatime date;
  if( run_time ) date = *run_time;
  vector<string> fnames = GetDBFileList( name, date, Here(here));
  // always look for 'pardata.map' in the current directory first.
  fnames.insert(fnames.begin(),string("pardata.map"));
  if( !fnames.empty() ) {
    vector<string>::iterator it = fnames.begin();
    do {
#ifdef WITH_DEBUG
      if( fDebug>0 ) {
	cout << "<" << Here(here) << ">: Opening database file " << *it;
      }
#endif
      pardatafile.clear();  // Forget previous failures before attempt
      pardatafile.open((*it).c_str());

#ifdef WITH_DEBUG
      if (pardatafile) {
	cout << "Opening file "<<*it<<endl;
      } else {
        cout << "cannot open file "<<*it<<endl;
      }
      if( fDebug>0 ) 
	if( !pardatafile ) cout << " ... failed" << endl;
	else               cout << " ... ok" << endl;
#endif
    } while ( !pardatafile && ++it != fnames.end() );
  }
  if( fnames.empty() || !pardatafile )
    if (PARDATA_VERBOSE) {
      cout << "WARNING:: ParityData: File db_"<<name<<".dat not found."<<endl;
      cout << "An example of this file should be in the examples directory."<<endl;
      cout << "Will proceed with default mapping for ParityData."<<endl;
      Int_t statm = DefaultMap();
      PrintMap();
      return statm;
    }

  string sinput;
  const string comment = "#";
  while (getline(pardatafile, sinput)) {
#ifdef WITH_DEBUG
    cout << "sinput "<<sinput<<endl;
#endif
    vector<string> strvect( vsplit(sinput) );
    if (strvect.size() < 5 || strvect[0] == comment) continue;
    Bool_t found = kFALSE;
    for (int i = 0; i < nvar; i++) {
      if (vars[i].name && strvect[0] == vars[i].name) found = kTRUE;
    }
// !found may be ok, but might be a typo error too, so I print to warn you.
    if ( !found && PARDATA_VERBOSE ) 
      cout << "ParityData: new variable "<<strvect[0]<<" will become global"<<endl;
    Int_t crate = (Int_t)atoi(strvect[2].c_str());  // crate #
    PdataLoc *b = 0;
    if (strvect[1] == "crate") {  // Crate data ?
      Int_t slot = (Int_t)atoi(strvect[3].c_str());
      Int_t chan = (Int_t)atoi(strvect[4].c_str());
      b = new PdataLoc(strvect[0].c_str(), crate, slot, chan); 
      fCrateLoc.push_back(b);
    } else {         // Data is relative to a header
      UInt_t header = header_str_to_base16(strvect[3].c_str());
      Int_t skip = (Int_t)atoi(strvect[4].c_str());
      b = new PdataLoc(strvect[0].c_str(), crate, header, skip);
      fWordLoc.push_back(b);
    }

    if (!found) {
      // a new variable to add to our dynamic list
      DefineChannel(b,mode);
    }
  }
  PrintMap(1);
  return retval;
}

//_____________________________________________________________________________
void ParityData::PrintMap(Int_t flag) {
  cout << "Map for Parity Data "<<endl;
  if (flag == 1) cout << "Map read from file "<<endl;
  for( Iter_t p = fCrateLoc.begin(); p != fCrateLoc.end(); p++) {
    PdataLoc *dataloc = *p;
    dataloc->Print();
  }
  for (Iter_t p = fWordLoc.begin(); p != fWordLoc.end(); p++) {
    PdataLoc *dataloc = *p;
    dataloc->Print();
  }  
}

//_____________________________________________________________________________
PdataLoc* ParityData::DefineChannel(PdataLoc *b, EMode mode, const char* desc) {
  string nm(fPrefix + b->name);
  if (mode==kDefine) {
    if (gHaVars) gHaVars->Define(nm.c_str(),desc,b->rdata[0],&(b->ndata));
  } else {
    if (gHaVars) gHaVars->RemoveName(nm.c_str());
  }
  return b;
}
  
//_____________________________________________________________________________
Int_t ParityData::End( THaRunBase* run ) 
{
  WriteHist();
  return 0;
}

//_____________________________________________________________________________
  void ParityData::WriteHist()
{
  for (UInt_t i = 0; i < hist.size(); i++) hist[i]->Write();
}

//_____________________________________________________________________________
  void ParityData::BookHist()
{
  hist.clear();
  hist.push_back(new TH1F("Cav1Q","Cavity1 Q",1200,-800,17000));
  hist.push_back(new TH1F("Cav1X","Cavity1 X",1200,-800,17000));
  hist.push_back(new TH1F("Cav1Y","Cavity1 Y",1200,-800,17000));
  hist.push_back(new TH1F("Cav2Q","Cavity2 Q",1200,-800,17000));
  hist.push_back(new TH1F("Cav2X","Cavity2 X",1200,-800,17000));
  hist.push_back(new TH1F("Cav2Y","Cavity2 Y",1200,-800,17000));
  hist.push_back(new TH1F("h4axp","4axp raw",1000,-500,4000));
  hist.push_back(new TH1F("h4axm","4axm raw",1000,-500,4000));
  hist.push_back(new TH1F("h4ayp","4ayp raw",1000,-500,4000));
  hist.push_back(new TH1F("h4aym","4ay, raw",1000,-500,4000));
  hist.push_back(new TH1F("h4bxp","4bxp raw",1000,-500,4000));
  hist.push_back(new TH1F("h4bxm","4bxm raw",1000,-500,4000));
  hist.push_back(new TH1F("h4byp","4byp raw",1000,-500,4000));
  hist.push_back(new TH1F("h4bym","4bym raw",1000,-500,4000));  //14 histos


}


//_____________________________________________________________________________
THaAnalysisObject::EStatus ParityData::Init( const TDatime& run_time ) 
{

  fStatus = kNotinit;
  MakePrefix();
  BookHist();
  InitCalib();
  fStatus = static_cast<EStatus>( SetupParData( &run_time ) );

#ifdef OLDSCALER
  // Get scalers.  They must be initialized in the
  // analyzer. 
  THaAnalyzer* theAnalyzer = THaAnalyzer::GetInstance();
  TList* scalerList = theAnalyzer->GetScalers();
  TIter next(scalerList);
  while( THaScalerGroup* tscalgrp = static_cast<THaScalerGroup*>( next() )) {
    THaScaler *scaler = tscalgrp->GetScalerObj();
    string lbank("Left");
    if (CmpNoCase(lbank,scaler->GetName()) == 0) {
         lscaler = scaler; 
    }
    string rbank("Right");
    if (CmpNoCase(rbank,scaler->GetName()) == 0) {
         rscaler = scaler; 
    }
  }
#endif
  return fStatus;
}


//_____________________________________________________________________________
Int_t ParityData::InitCalib() {
// Calibration initialization

// 568  518  427  448  1.0  1.0  18.87  1.11  1.07   -0.47   0.67
// 490  426  553  520  1.0  1.0  18.87  0.95  0.99    0.32   0.34
 
//%%%%%%%%%%%  Calibration for Striplines  %%%%%%%%%%%%
//In each line above, first four entries are pedestal.
//5th, 6th entries are alpha_x, alpha_y gain factors.
//By historical convention they are always 1.0 now.
//7th entry is kappa in units of millimeters.
//8th,9th entries are ''attenuation factors'', x_att, y_att
//Last two entries are x_offset,y_offset.  Again millimeters.
//Note, the offset and attenuation is applied after rotation.
//The formulas:
//  xp = xp_adc - xp_ped (attenna signal above ped); similarly xm
//  x_rot = kappa*(xp - alpha_x*xm)/(xp + alpha_x*xm)
//  y_rot  similar
//  x = -1*x_offset + (x_rot - y_rot)*x_att / sqrt(2)
//  y = y_offset + (x_rot + y_rot)*y_att / sqrt(2)
//Sign convention should agree with EPICS

// -------------------------------------------

  
// cavities I,X,Y for 2 monitors
// (One issue: unplugging is different from no signal 
//  from floor of hall.  presumably a ground diff
//  leading to offset due to differential amplifier)

//    cped[0] = 550.2;
//    cped[1] = 503.6;
//    cped[2] = 514.0;
//    cped[3] = 514.9;
//    cped[4] = 618.6;
//    cped[5] = 587.6;

    aqcav4a = 1;
    axcav4a = 1;
    aycav4a = 1;
    aqcav4b = 1;
    axcav4b = 1;
    aycav4b = 1;
    bqcav4a = 0;
    bxcav4a = 0;
    bycav4a = 0;
    bqcav4b = 0;
    bxcav4b = 0;
    bycav4b = 0;

    // Striplines
    alpha4ax = 1.0;
    alpha4ay = 1.0;
    alpha4bx = 1.0;
    alpha4by = 1.0;
    kappa4a  = 18.87;
    kappa4b  = 18.87;
    att4ax   = 1.11;
    att4ay   = 0.97;
    att4bx   = 0.98;
    att4by   = 0.95;
    off4ax   = -1.358;
    off4ay   = 0.53;
    off4bx   = 0.359;
    off4by   = 1.09;


    //    sped[0] = 652.4;
    //    sped[1] = 441.3;
    //    sped[2] = 292.8;
    //    sped[3] = 225.3;
    //   sped[4] = 364.7;
    //    sped[5] = 236.1;
    //    sped[6] = 269;
    //   sped[7] = 389;
    sped[8] = 0;
    sped[9] = 0;
    sped[10] = 0;
    sped[11] = 0;


    return 0;

}

//_____________________________________________________________________________
Int_t ParityData::DefaultMap() {
// Default setup of mapping of data in 
// this class to locations in the raw data.

// Bit pattern for trigger definition

   for (UInt_t i = 0; i < bits.GetNbits(); i++) {
     fCrateLoc.push_back(new PdataLoc(Form("bit%d",i+1), 3, (Int_t) 11, 48+i));
   }

   fCrateLoc.push_back(new PdataLoc("upQadcL", 3, 23, 16)); 
   fCrateLoc.push_back(new PdataLoc("upQadcR", 1, 25, 40)); 
   fCrateLoc.push_back(new PdataLoc("loQadcL", 3, 23, 20)); 
   fCrateLoc.push_back(new PdataLoc("loQadcR", 1, 25, 36)); 

   fCrateLoc.push_back(new PdataLoc("upQtdcL", 4, 13, 84)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("upQtdcR", 2, 14, 84)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("loQtdcL", 4, 13, 80)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("loQtdcR", 2, 14, 80)); // Ahmed.

   fCrateLoc.push_back(new PdataLoc("atQadcL", 3, 23, 22)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("atQadcR", 1, 25, 39)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("atQtdcL", 4, 13, 86)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("atQtdcR", 2, 14, 85)); // Ahmed.

   fCrateLoc.push_back(new PdataLoc("hapadcL", 3, 25, 38));  // 7th chan.
   fCrateLoc.push_back(new PdataLoc("hapadcR", 1, 25, 32));  // 1st chan.
   fCrateLoc.push_back(new PdataLoc("haptdcL", 4, 13, 86));
   fCrateLoc.push_back(new PdataLoc("haptdcR", 1, 14, 80));

   fCrateLoc.push_back(new PdataLoc("upSadcL", 3, 23, 18)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("upSadcR", 1, 25, 34)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("upStdcL", 4, 13, 82)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("upStdcR", 2, 14, 82)); // Ahmed.

   fCrateLoc.push_back(new PdataLoc("loSadcL", 3, 23, 19)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("loSadcR", 1, 25, 35)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("loStdcL", 4, 13, 83)); // Ahmed.
   fCrateLoc.push_back(new PdataLoc("loStdcR", 2, 14, 83)); // Ahmed.

   fCrateLoc.push_back(new PdataLoc("amp1", 1, 25, 1));
   fCrateLoc.push_back(new PdataLoc("bamp", 1, 25, 15));

//    fCrateLoc.push_back(new PdataLoc("s0La", 3, 23, 16));
//    fCrateLoc.push_back(new PdataLoc("s0Lb", 3, 23, 22));
//    fCrateLoc.push_back(new PdataLoc("s0Ra", 1, 25, 0));
//    fCrateLoc.push_back(new PdataLoc("s0Rb", 1, 25, 6));

   fCrateLoc.push_back(new PdataLoc("s0La", 3, 25, 16));
   fCrateLoc.push_back(new PdataLoc("s0Lb", 3, 25, 22));
   fCrateLoc.push_back(new PdataLoc("s0Ra", 1, 25, 0));
   fCrateLoc.push_back(new PdataLoc("s0Rb", 1, 25, 8));

   // Quartz moved from 2nd to 4th chan, Aug 24
   fCrateLoc.push_back(new PdataLoc("Apscint", 1, 25, 34));
   fCrateLoc.push_back(new PdataLoc("Aquartz", 1, 25, 35));
   fCrateLoc.push_back(new PdataLoc("Tpscint", 1, 14, 82));
   fCrateLoc.push_back(new PdataLoc("Tquartz", 1, 14, 83));


   return 0;
}


//_____________________________________________________________________________
Int_t ParityData::Decode(const THaEvData& evdata)
{

  Int_t i;
  Int_t ldebug=1;

  Clear();


  for( Iter_t p = fCrateLoc.begin(); p != fCrateLoc.end(); p++) {
    PdataLoc *dataloc = *p;
    if ( dataloc->IsSlot() ) {  
      for (i = 0; i < evdata.GetNumHits(dataloc->crate, 
		         dataloc->slot, dataloc->chan); i++) {
	dataloc->Load(evdata.GetData(dataloc->crate, dataloc->slot, 
				     dataloc->chan, i));
      }
    }
  }
  

  for (i = 0; i < evdata.GetEvLength(); i++) {
    for (Iter_t p = fWordLoc.begin(); p != fWordLoc.end(); p++) {
      PdataLoc *dataloc = *p;
      if ( dataloc->DidLoad() || dataloc->IsSlot() ) continue;
      if ( evdata.InCrate(dataloc->crate, i) ) {
        if ((UInt_t)evdata.GetRawData(i) == dataloc->header) {
            dataloc->Load(evdata.GetRawData(i + dataloc->ntoskip));
        }
      }
    }
  }

  evtype = evdata.GetEvType();   // CODA event type 

  for( Iter_t p = fCrateLoc.begin(); p != fCrateLoc.end(); p++) {
    PdataLoc *dataloc = *p;

// bit pattern of triggers
    for (UInt_t i = 0; i < bits.GetNbits(); i++) {
      if ( dataloc->ThisIs(Form("bit%d",i+1)) ) TrigBits(i+1,dataloc);
    }

    if ( dataloc->ThisIs("upQadcL") ) upQadcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("upQadcR") ) upQadcR  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("loQadcL") ) loQadcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("loQadcR") ) loQadcR  = dataloc->Get(); // Ahmed.

    if ( dataloc->ThisIs("upQtdcL") ) upQtdcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("upQtdcR") ) upQtdcR  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("loQtdcL") ) loQtdcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("loQtdcR") ) loQtdcR  = dataloc->Get(); // Ahmed.

    if ( dataloc->ThisIs("atQadcL") ) atQadcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("atQadcR") ) atQadcR  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("atQtdcL") ) atQtdcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("atQtdcR") ) atQtdcR  = dataloc->Get(); // Ahmed.

    if ( dataloc->ThisIs("upSadcL") ) upSadcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("upSadcR") ) upSadcR  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("upStdcL") ) upStdcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("upStdcR") ) upStdcR  = dataloc->Get(); // Ahmed.

    if ( dataloc->ThisIs("loSadcL") ) loSadcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("loSadcR") ) loSadcR  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("loStdcL") ) loStdcL  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("loStdcR") ) loStdcR  = dataloc->Get(); // Ahmed.

    if ( dataloc->ThisIs("amp1") ) amp1  = dataloc->Get(); // Ahmed.
    if ( dataloc->ThisIs("bamp") ) bamp  = dataloc->Get(); // Ahmed.

    if ( dataloc->ThisIs("hapadcL") ) hapadcL  = dataloc->Get();
    if ( dataloc->ThisIs("hapadcR") ) hapadcR  = dataloc->Get();
    if ( dataloc->ThisIs("haptdcL") ) haptdcL  = dataloc->Get();
    if ( dataloc->ThisIs("haptdcR") ) haptdcR  = dataloc->Get();

    if ( dataloc->ThisIs("s0La") ) s0La  = dataloc->Get();
    if ( dataloc->ThisIs("s0Lb") ) s0Lb  = dataloc->Get();
    if ( dataloc->ThisIs("s0Ra") ) s0Ra  = dataloc->Get();
    if ( dataloc->ThisIs("s0Rb") ) s0Rb  = dataloc->Get();

    if ( dataloc->ThisIs("Aquartz") ) Aquartz  = dataloc->Get();
    if ( dataloc->ThisIs("Tquartz") ) Tquartz  = dataloc->Get();
    if ( dataloc->ThisIs("Apscint") )  {
         Apscint  = dataloc->Get();
	 //       cout << "Apscint = "<<Apscint<<endl;
    }
    if ( dataloc->ThisIs("Tpscint") ) {
         Tpscint  = dataloc->Get();
	 //         cout << "Tpscint = "<<Tpscint<<endl;
    }

  }

  // Striplines and Cavities on R-HRS.

  Int_t runnum = evdata.GetRunNum();

  // The following is for 2009 - 2010

  if (modern_times) {
   if (runnum < -24000) {  // obsolete cut
     Int_t str_roc=2;
     Int_t str_slot=23;
     Int_t str_chan0=0;
     for (int i=0; i<NSTR; i++) {
       sfb[i] = evdata.GetData(str_roc,str_slot,
			       str_chan0+i, 0);         
     }

  // Cavities
     Int_t cav_roc=2;
     Int_t cav_slot=6;
     Int_t cav_chan0=32;
     for (int i=0; i<NCAV; i++) {
            cfb[i] = evdata.GetData(cav_roc,cav_slot,
                           cav_chan0+i, 0);          
            if (ldebug) printf("cavities in roc %d   slot %d    chan %d  data = %f\n",
			 cav_roc, cav_slot, cav_chan0+i, cfb[i]);
     }
   } else {
     if (ldebug) cout << endl << endl;
     Int_t str_roc=4;
     Int_t str_slot=18;
     Int_t str_chan0=0;
     for (int i=0; i<NSTR; i++) {
       sfb[i] = evdata.GetData(str_roc,str_slot, str_chan0+i, 0);
       if (ldebug) cout << "Left HRS spot++ raw data "<<i<<"  "<<sfb[i]<<endl;
     }  // cut on run number
   }

  } else { 

  // The following is for 2008 test run

  // Striplines, 2008
     Int_t str_roc=1;
     Int_t str_slot=23;
     Int_t str_chan0=0;
     for (int i=0; i<NSTR; i++) {
        sfb[i] = evdata.GetData(str_roc,str_slot,
                 str_chan0+i, 0);         
     }

  // Cavities, 2008
     Int_t cav_roc=1;
     Int_t cav_slot=23;
     Int_t cav_chan0=16;
     for (int i=0; i<NCAV; i++) {
      cfb[i] = evdata.GetData(cav_roc,cav_slot,
                 cav_chan0+i, 0);          
     if (ldebug) printf("cavities in roc %d   slot %d    chan %d  data = %f\n",
			cav_roc, cav_slot, cav_chan0+i, cfb[i]);

     }

  }

  DoBpm();

  DoKine();

  if (PARDATA_PRINT) Print();

  return 0;
}


//_____________________________________________________________________________
Int_t ParityData::DoBpm( ) {
// The 'stripline' BPM code is redundant to the 
// analyzer.  The 'cavity' BPM is new.

  Double_t x_rot, y_rot;
  static Double_t sqrt2 = 1.41421;
  Int_t loc_debug = 0;

  xp4a = sfb[0]-sped[0];  // "s" = stripline
  xm4a = sfb[1]-sped[1];
  yp4a = sfb[2]-sped[2];
  ym4a = sfb[3]-sped[3];
  xp4b = sfb[4]-sped[4];
  xm4b = sfb[5]-sped[5];
  yp4b = sfb[6]-sped[6];
  ym4b = sfb[7]-sped[7];

  hist[6]->Fill(xp4a);
  hist[7]->Fill(xm4a);
  hist[8]->Fill(yp4a);
  hist[9]->Fill(ym4a);
  hist[10]->Fill(xp4b);
  hist[11]->Fill(xm4b);
  hist[12]->Fill(yp4b);
  hist[13]->Fill(ym4b);

  x_rot = kappa4a*(xp4a-alpha4ax*xm4a)/(xp4a+alpha4ax*xm4a);
  y_rot = kappa4a*(yp4a-alpha4ay*ym4a)/(yp4a+alpha4ay*ym4a);
  xstrip4a = -1*off4ax + (x_rot - y_rot)*att4ax / sqrt2;
  ystrip4a = off4ay + (x_rot + y_rot)*att4ay / sqrt2;

  if (loc_debug) cout << "strip 4a tst "<<x_rot<<"  "<<y_rot<<"  "<<xstrip4a<<"  "<<ystrip4a<<endl;


  x_rot = kappa4b*(xp4b-alpha4bx*xm4b)/(xp4b+alpha4bx*xm4b);
  y_rot = kappa4b*(yp4b-alpha4by*ym4b)/(yp4b+alpha4by*ym4b);
  xstrip4b = -1*off4bx + (x_rot - y_rot)*att4bx / sqrt2;
  ystrip4b = off4by + (x_rot + y_rot)*att4by / sqrt2;

  if (loc_debug) cout << "strip 4b tst "<<x_rot<<"  "<<y_rot<<"  "<<xstrip4b<<"  "<<ystrip4b<<endl;

  // cavity monitor QXY.  Calibration here:
  qcav4a = aqcav4a*(cfb[0]-cped[0])+bqcav4a;
  xcav4a = axcav4a*(cfb[1]-cped[1])+bxcav4a;
  ycav4a = aycav4a*(cfb[2]-cped[2])+bycav4a;
  qcav4b = aqcav4b*(cfb[3]-cped[3])+bqcav4b;
  xcav4b = axcav4b*(cfb[4]-cped[4])+bxcav4b;
  ycav4b = aycav4b*(cfb[5]-cped[5])+bycav4b;

  hist[0]->Fill(qcav4a);
  hist[1]->Fill(xcav4a);
  hist[2]->Fill(ycav4a);
  hist[3]->Fill(qcav4b);
  hist[4]->Fill(xcav4b);
  hist[5]->Fill(ycav4b);

  // Build a line, make residuals

  if (qcav4a<1150&&qcav4b<600)  {
  //Convert all cavity readings to mm
  Float_t xcav4amm = -.009416729*xcav4a+28.52658;
  Float_t ycav4amm = -.008864325*ycav4a+24.09716;
  Float_t xcav4bmm = -.002436467*xcav4b+7.310761;
  Float_t ycav4bmm = -.00191471*ycav4b+5.004565;

  //z values of each monitor
  Float_t Zxcav4a = 6.83;
  Float_t Zycav4a = 6.53;
  Float_t Zxcav4b = 3.78;
  Float_t Zycav4b = 3.49;

  Float_t Zstrip4a = 7.52;
  Float_t Zstrip4b = 2.38;

//4A, X

  //Sums 

  Float_t sumxz4AX = xcav4bmm*Zxcav4b+xstrip4a*Zstrip4a+xstrip4b*Zstrip4b;
  Float_t sumz4AX = Zxcav4b+Zstrip4a+Zstrip4b;
  Float_t sumx4AX = xcav4bmm+xstrip4a+xstrip4b;
  Float_t sumz24AX = Zxcav4b*Zxcav4b+Zstrip4a*Zstrip4a+Zstrip4b*Zstrip4b;
  
  //Linear fit

  Float_t m4AX = (3*sumxz4AX-sumz4AX*sumx4AX)/(3*sumz24AX-sumz4AX*sumz4AX);
  Float_t b4AX = (sumx4AX-m4AX*sumz4AX)/3;

  //check fits
  //cout << xstrip4a << " (" << Zstrip4a << ") " << endl;
  //cout << xstrip4b << " (" << Zstrip4b << ") " << endl;
  //cout << xcav4bmm << " (" << Zxcav4b << ") " << endl;
  //cout << endl;

  //cout << "xcav4amm=" << xcav4amm << endl;
  //cout << "m4AX = " << m4AX << endl;
  //cout << "b4AX = " << b4AX << endl;
  //cout << "*" << endl;

  //Residuals

  Float_t xpred4AX = b4AX + Zxcav4a*m4AX;
  Float_t res4AX = xcav4amm-xpred4AX;

//4B, X

  //Sums 

  Float_t sumxz4BX = xcav4amm*Zxcav4a+xstrip4a*Zstrip4a+xstrip4b*Zstrip4b;
  Float_t sumz4BX = Zxcav4a+Zstrip4a+Zstrip4b;
  Float_t sumx4BX = xcav4amm+xstrip4a+xstrip4b;
  Float_t sumz24BX = Zxcav4a*Zxcav4a+Zstrip4a*Zstrip4a+Zstrip4b*Zstrip4b;
  
  //Linear fit

  Float_t m4BX = (3*sumxz4BX-sumz4BX*sumx4BX)/(3*sumz24BX-sumz4BX*sumz4BX);
  Float_t b4BX = (sumx4BX-m4BX*sumz4BX)/3;

  //Residuals

  Float_t xpred4BX = b4BX + Zxcav4b*m4BX;
  Float_t res4BX = xcav4bmm-xpred4BX;

//4A, Y

  //Sums 

  Float_t sumyz4AY = ycav4bmm*Zycav4b+ystrip4a*Zstrip4a+ystrip4b*Zstrip4b;
  Float_t sumz4AY = Zycav4b+Zstrip4a+Zstrip4b;
  Float_t sumy4AY = ycav4bmm+ystrip4a+ystrip4b;
  Float_t sumz24AY = Zycav4b*Zycav4b+Zstrip4a*Zstrip4a+Zstrip4b*Zstrip4b;
  
  //Linear fit

  Float_t m4AY = (3*sumyz4AY-sumz4AY*sumy4AY)/(3*sumz24AY-sumz4AY*sumz4AY);
  Float_t b4AY = (sumy4AY-m4AY*sumz4AY)/3;

  //Residuals

  Float_t ypred4AY = b4AY + Zycav4a*m4AY;
  Float_t res4AY = ycav4amm-ypred4AY;

//4B, Y

  //Sums 

  Float_t sumyz4BY = ycav4amm*Zycav4a+ystrip4a*Zstrip4a+ystrip4b*Zstrip4b;
  Float_t sumz4BY = Zycav4a+Zstrip4a+Zstrip4b;
  Float_t sumy4BY = ycav4amm+ystrip4a+ystrip4b;
  Float_t sumz24BY = Zycav4a*Zycav4a+Zstrip4a*Zstrip4a+Zstrip4b*Zstrip4b;
  
  //Linear fit

  Float_t m4BY = (3*sumyz4BY-sumz4BY*sumy4BY)/(3*sumz24BY-sumz4BY*sumz4BY);
  Float_t b4BY = (sumy4BY-m4BY*sumz4BY)/3;

  //Residuals

  Float_t ypred4BY = b4BY + Zycav4b*m4BY;
  Float_t res4BY = ycav4bmm-ypred4BY;
  }

  if (loc_debug) {
    cout << "sfb ---> "<<endl;
    for (int i = 0; i < 8; i++) {
      cout << "   sig = "<<sfb[i]<<"   ped = "<<sped[i]<<endl;
    }
    cout << "calib --> "<<endl;
    cout << "alphas "<<alpha4ax<<"  "<<alpha4ay;
    cout << "  "<<alpha4bx<<"  "<<alpha4by<<endl;
    cout << "offsets "<<off4ax<<"  "<<off4ay;
    cout << "  "<<off4bx<<"  "<<off4by<<endl;
    cout << "kappa "<<kappa4a<<"  "<<kappa4b<<endl;
    cout << "xstrip4a, etc "<<xstrip4a<<"  "<<ystrip4a;
    cout << "  "<<xstrip4b<<"  "<<ystrip4b<<endl;
    cout << "\n\ncfb --> "<<endl;
    for (int i = 0; i < NCAV; i++) {
      cout << "  sig = "<<cfb[i]<<"   ped = "<<cped[i]<<endl;
    }
    cout << "qcav4a = "<<qcav4a<<endl;
    cout << "xcav4a = "<<xcav4a<<endl;
    cout << "ycav4a = "<<ycav4a<<endl;
    cout << "qcav4b = "<<qcav4b<<endl;
    cout << "xcav4b = "<<xcav4b<<endl;
    cout << "ycav4b = "<<ycav4b<<endl;
  }

  return 1;
}


//_____________________________________________________________________________
Int_t ParityData::DoKine( ) {
// Calculate Q^2 and missing mass squared.
// I realized later that this routine obtains
// the same result as THaElectronKine, but that's
// a good check.


  Double_t pL,thL,phiL;
  Double_t pR,thR,phiR;
  THaVar *pvar;

  Double_t ebeam = 3.48;     // E_beam - < dE/dx >
  Double_t mprot = 0.93827;

  Double_t theta_L = 14.0;    // central angle, degrees
  Double_t theta_R = 14.0;    // central angle, degrees
  Double_t thoff_L = 0;       // offset to play with
  Double_t thoff_R = 0;       // offset "   "    "

  Double_t pscale  = 1.00;    // momentum scale factor
                              // to play with

  Double_t pi = 3.1415926;
  Double_t thL0,thR0,costhL,costhR,sinthL,sinthR;

  thL0 = (theta_L+thoff_L)*pi/180;
  costhL = TMath::Cos(thL0);
  sinthL = TMath::Sin(thL0);
  thR0 = (theta_R+thoff_R)*pi/180;
  costhR = TMath::Cos(thR0);
  sinthR = TMath::Sin(thR0);

  pL   = -999;
  thL  = -999;
  phiL = -999;
  pR   = -999;
  thR  = -999;
  phiR = -999;
  Int_t okL = 0;  
  Int_t okR = 0;

// Must have the THaGoldenTrack module invoked
  
  pvar = gHaVars->Find("L.gold.ok");
  if (pvar) okL = (Int_t) pvar->GetValue();
  if (okL) {
    pvar = gHaVars->Find("L.gold.p");
    if (pvar) pL = pvar->GetValue();
    pvar = gHaVars->Find("L.gold.th");
    if (pvar) thL = pvar->GetValue();
    pvar = gHaVars->Find("L.gold.ph");
    if (pvar) phiL = pvar->GetValue();
  }

  pvar = gHaVars->Find("R.gold.ok");
  if (pvar) okR = (Int_t) pvar->GetValue();
  if (okR) {
    pvar = gHaVars->Find("R.gold.p");
    if (pvar) pR = pvar->GetValue();
    pvar = gHaVars->Find("R.gold.th");
    if (pvar) thR = pvar->GetValue();
    pvar = gHaVars->Find("R.gold.ph");
    if (pvar) phiR = pvar->GetValue();
  }  

  q2L = -999;
  mmL = -999;

  if (pL > -999 && thL > -999 && phiL > -999 &&
      pL < 1e32 && thL < 1e32 && phiL < 1e32) {
    q2L = 2*ebeam*(pscale*pL)*(1-((costhL-sinthL*phiL)/TMath::Sqrt(1+thL*thL+phiL*phiL)));
    mmL = 2*mprot*(ebeam-(pscale*pL))-q2L;
  }

  q2R = -999;
  mmR = -999;
  if (pR > -999 && thR > -999 && phiR > -999 &&
      pR < 1e32 && thR < 1e32 && phiR < 1e32) {
// note : sign convention requires costhR + sinthR*phiR 
//                                    not -
    q2R = 2*ebeam*(pscale*pR)*(1-((costhR+sinthR*phiR)/TMath::Sqrt(1+thR*thR+phiR*phiR)));
    mmR = 2*mprot*(ebeam-(pscale*pR))-q2R;
  }

  return 1;

}


//_____________________________________________________________________________
void ParityData::Print( Option_t* opt ) const {
// Dump the data for purpose of debugging.
  cout << "Dump of ParityData "<<endl;
  cout << "event pattern bits : ";
   for (UInt_t i = 0; i < bits.GetNbits(); i++) 
    cout << " "<<i<<" = "<< bits.TestBitNumber(i)<<"  | ";
  cout << endl;
  cout << "PREX Upper Quartz ADC L, R "<<upQadcL<<"  "<<upQadcR<<endl; // Ahmed.
  cout << "PREX Lower Quartz ADC L, R "<<loQadcL<<"  "<<loQadcR<<endl; // Ahmed.
  cout << "PREX Upper Quartz TDC L, R "<<upQtdcL<<"  "<<upQtdcR<<endl; // Ahmed.
  cout << "PREX Lower Quartz ADC L, R "<<loQadcL<<"  "<<loQadcR<<endl; // Ahmed.
  cout << "PREX A_T Quartz ADC L, R "<<atQadcL<<"  "<<atQadcR<<endl;   // Ahmed.
  cout << "PREX A_T Quartz TDC L, R "<<atQtdcL<<"  "<<atQtdcR<<endl;   // Ahmed.
  cout << "PREX Upper Scint ADC L, R "<<upSadcL<<"  "<<upSadcR<<endl;  // Ahmed.
  cout << "PREX Upper Scint TDC L, R "<<upStdcL<<"  "<<upStdcR<<endl;  // Ahmed.
  cout << "PREX Lower Scint ADC L, R "<<loSadcL<<"  "<<loSadcR<<endl;  // Ahmed.
  cout << "PREX Lower Scint TDC L, R "<<loStdcL<<"  "<<loStdcR<<endl;  // Ahmed.
  cout << "event types,  CODA = "<<evtype;
  cout << "   bit pattern = "<<evtypebits<<endl;
  cout << "HAPPEX ADC L, R "<<hapadcL<<"  "<<hapadcR<<endl;
  cout << "HAPPEX TDC L, R "<<haptdcL<<"  "<<haptdcR<<endl;
  cout << "S0 detectors Left HRS  "<<s0La<<"  "<<s0Lb<<endl;
  cout << "S0 detectors Right HRS  "<<s0Ra<<"  "<<s0Rb<<endl;
  cout << "Quartz detector ADC "<<Aquartz<<"   TDC "<<Tquartz<<endl;
  cout << "PREX scint. det ADC "<<Apscint<<"   TDC "<<Tpscint<<endl;
  cout << "Cavity at 4A  X/Y/Q "<<
    xcav4a<<"  "<<ycav4a<<"  "<<qcav4a<<endl;
  cout << "Cavity at 4B  X/Y/Q "<<
    xcav4b<<"  "<<ycav4b<<"  "<<qcav4b<<endl;
  cout << "Stripline at 4A  X/Y "<<
    xstrip4a<<"  "<<ystrip4a<<endl;
  cout << "Stripline at 4B  X/Y "<<
    xstrip4a<<"  "<<ystrip4b<<endl;
  cout << "Left arm q2 "<<q2L<<"   miss mass sq "<<mmL<<endl;
  cout << "Right arm q2 "<<q2R<<"   miss mass sq "<<mmR<<endl;
#ifdef CHECK1
  cout << "Accepted Trigger Counts "<<endl;
  for (Int_t i = 0; i < 12; i++) {
    cout << "     Trig "<<i+1<<"   count = "<<trigcnt[i]<<endl;
  }
#endif
}


//_____________________________________________________________________________
vector<string> ParityData::vsplit(const string& s) {
// split a string into whitespace-separated strings
  vector<string> ret;
  typedef string::size_type ssiz_t;
  ssiz_t i = 0;
  while ( i != s.size()) {
    while (i != s.size() && isspace(s[i])) ++i;
      ssiz_t j = i;
      while (j != s.size() && !isspace(s[j])) ++j;
      if (i != j) {
         ret.push_back(s.substr(i, j-i));
         i = j;
      }
  }
  return ret;
}

//_____________________________________________________________________________
UInt_t ParityData::header_str_to_base16(const char* hdr) {
// Utility to convert string header to base 16 integer
  static const char chex[] = "0123456789abcdef";
  if( !hdr ) return 0;
  const char* p = hdr+strlen(hdr);
  UInt_t result = 0;  UInt_t power = 1;
  while( p-- != hdr ) {
    const char* q = strchr(chex,tolower(*p));
    if( q ) {
      result += (q-chex)*power; 
      power *= 16;
    }
  }
  return result;
};

//_____________________________________________________________________________
void ParityData::TrigBits(UInt_t ibit, PdataLoc *dataloc) {
// Figure out which triggers got a hit.  These are multihit TDCs, so we
// need to sort out which hit we want to take by applying cuts.

  if( ibit >= kBitsPerByte*sizeof(UInt_t) ) return; //Limit of evtypebits
  bits.ResetBitNumber(ibit);

  static const UInt_t cutlo = 10;
  static const UInt_t cuthi = 1000;

  //  cout << "Bit TDC num hits "<<dataloc->NumHits()<<endl;
    for (int ihit = 0; ihit < dataloc->NumHits(); ihit++) {
      //            cout << "TDC data " << ibit<<"  "<<dataloc->Get(ihit)<<endl;

  if (dataloc->Get(ihit) > cutlo && dataloc->Get(ihit) < cuthi) {
      bits.SetBitNumber(ibit);
      evtypebits |= BIT(ibit);
    }
  }

}
//_____________________________________________________________________________
ClassImp(ParityData)

