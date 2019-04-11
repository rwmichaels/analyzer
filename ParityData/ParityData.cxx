//*-- Author :    Bob Michaels,  April 2019

//////////////////////////////////////////////////////////////////////////
//
// ParityData
//
// PREX data from spectrometer DAQ.
//   Quartz detectors
//   S0 and S2 scintillators
//   Gather data and tracking results for Q2
//   Looking for GEMs ?  it's not here
//
//////////////////////////////////////////////////////////////////////////

//#define WITH_DEBUG 1

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

#define NSTR 12

using namespace std;
using namespace THaString;

typedef vector<PdataLoc*>::iterator Iter_t;

//_____________________________________________________________________________
ParityData::ParityData( const char* name, const char* descript ) : 
  THaApparatus( name, descript )
{
  Clear();
}

//_____________________________________________________________________________
ParityData::~ParityData()
{
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
  s2La = 0; s2Lb = 0; s2Ra = 0;  s2Rb = 0;
  upQadcL =0; upQadcR = 0;  
  loQadcL =0; loQadcR = 0;  
  upQtdcL =0; upQtdcR = 0;  
  loQtdcL =0; loQtdcR = 0;  
  atQadcL =0; atQadcR = 0;  
  atQtdcL =0; atQtdcR = 0;  
  upSadcL =0; upSadcR = 0;  
  upStdcL =0; upStdcR = 0;  
  loSadcL =0; loSadcR = 0;  
  loStdcL =0; loStdcR = 0;  
  for( Iter_t p = fWordLoc.begin();  p != fWordLoc.end(); p++)  (*p)->Clear();
  for( Iter_t p = fCrateLoc.begin(); p != fCrateLoc.end(); p++) (*p)->Clear();
}

//_____________________________________________________________________________
Int_t ParityData::SetupParData( const TDatime* run_time, EMode mode )
{

  Int_t retval = 0;

  // Creating global variables
  RVarDef vars[] = {
    { "upQadcL",   "PREX Upper Quartz ADC LHRS",  "upQadcL" },  
    { "upQadcR",   "PREX Upper Quartz ADC RHRS",  "upQadcR" },  
    { "loQadcL",   "PREX Lower Quartz ADC LHRS",  "loQadcL" },  
    { "loQadcR",   "PREX Lower Quartz ADC RHRS",  "loQadcR" },  
    { "upQtdcL",   "PREX Upper Quartz TDC LHRS",  "upQtdcL" },  
    { "upQtdcR",   "PREX Upper Quartz TDC RHRS",  "upQtdcR" },  
    { "loQtdcL",   "PREX Lower Quartz TDC LHRS",  "loQtdcL" },  
    { "loQtdcR",   "PREX Lower Quartz TDC RHRS",  "loQtdcR" },  
    { "atQadcL",   "PREX A_T Quartz ADC LHRS",  "atQadcL" },    
    { "atQadcR",   "PREX A_T Quartz ADC RHRS",  "atQadcR" },    
    { "atQtdcL",   "PREX A_T Quartz TDC LHRS",  "atQtdcL" },    
    { "atQtdcR",   "PREX A_T Quartz TDC RHRS",  "atQtdcR" },    
    { "upSadcL",   "PREX Upper Scint ADC LHRS",   "upSadcL" },  
    { "upSadcR",   "PREX Upper Scint ADC RHRS",   "upSadcR" },  
    { "loSadcL",   "PREX Lower Scint ADC LHRS",   "loSadcL" },  
    { "loSadcR",   "PREX Lower Scint ADC RHRS",   "loSadcR" },  
    { "upStdcL",   "PREX Upper Scint TDC LHRS",   "upStdcL" },  
    { "upStdcR",   "PREX Upper Scint TDC RHRS",   "upStdcR" },  
    { "loStdcL",   "PREX Lower Scint TDC LHRS",   "loStdcL" },  
    { "loStdcR",   "PREX Lower Scint TDC RHRS",   "loStdcR" },  
    { "evtypebits", "event type bit pattern",   "evtypebits" },  
    { "evtype",     "event type from bit pattern", "evtype" },  
    { "s0La",      "S0 Left HRS PMT A",       "s0La" },
    { "s0Lb",      "S0 Left HRS PMT B",       "s0Lb" },
    { "s0Ra",      "S0 Right HRS PMT A",      "s0Ra" },
    { "s0Rb",      "S0 Right HRS PMT B",      "s0Rb" },
    { "s2La",      "S2 Left HRS PMT A",       "s2La" },
    { "s2Lb",      "S2 Left HRS PMT B",       "s2Lb" },
    { "s2Ra",      "S2 Right HRS PMT A",      "s2Ra" },
    { "s2Rb",      "S2 Right HRS PMT B",      "s2Rb" },
    { "q2L",   "Qsq on Left arm",             "q2L" },  
    { "q2R",   "Qsq on Right arm",            "q2R" },  
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
      if( fDebug>0 ) {
	if( !pardatafile ) cout << " ... failed" << endl;
	else               cout << " ... ok" << endl;
      }
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
  hist.push_back(new TH1F("test1","test histo 1",1200,-800,17000));
  hist.push_back(new TH1F("test2","test histo 2",1200,-800,17000));

}


//_____________________________________________________________________________
THaAnalysisObject::EStatus ParityData::Init( const TDatime& run_time ) 
{

  fStatus = kNotinit;
  MakePrefix();
  BookHist();
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
Int_t ParityData::DefaultMap() {
// Default setup of mapping of data in 
// this class to locations in the raw data.
// The code comes here if db_P.dat is not found in $DB_DIR.

// Bit pattern for trigger definition

   cout << "ParityData:: Warning:  Using the default map"<<endl;

   for (UInt_t i = 0; i < bits.GetNbits(); i++) {
     fCrateLoc.push_back(new PdataLoc(Form("bit%d",i+1), 3, (Int_t) 11, 48+i));
   }

   fCrateLoc.push_back(new PdataLoc("upQadcL", 3, 23, 16)); 
   fCrateLoc.push_back(new PdataLoc("upQadcR", 1, 25, 40)); 
   fCrateLoc.push_back(new PdataLoc("loQadcL", 3, 23, 20)); 
   fCrateLoc.push_back(new PdataLoc("loQadcR", 1, 25, 36)); 

   fCrateLoc.push_back(new PdataLoc("upQtdcL", 4, 13, 84));  
   fCrateLoc.push_back(new PdataLoc("upQtdcR", 2, 14, 84));  
   fCrateLoc.push_back(new PdataLoc("loQtdcL", 4, 13, 80));  
   fCrateLoc.push_back(new PdataLoc("loQtdcR", 2, 14, 80));  

   fCrateLoc.push_back(new PdataLoc("atQadcL", 3, 23, 22));  
   fCrateLoc.push_back(new PdataLoc("atQadcR", 1, 25, 39));  
   fCrateLoc.push_back(new PdataLoc("atQtdcL", 4, 13, 86));  
   fCrateLoc.push_back(new PdataLoc("atQtdcR", 2, 14, 85));  

   fCrateLoc.push_back(new PdataLoc("upSadcL", 3, 23, 18));  
   fCrateLoc.push_back(new PdataLoc("upSadcR", 1, 25, 34));  
   fCrateLoc.push_back(new PdataLoc("upStdcL", 4, 13, 82));  
   fCrateLoc.push_back(new PdataLoc("upStdcR", 2, 14, 82));  

   fCrateLoc.push_back(new PdataLoc("loSadcL", 3, 23, 19));  
   fCrateLoc.push_back(new PdataLoc("loSadcR", 1, 25, 35));  
   fCrateLoc.push_back(new PdataLoc("loStdcL", 4, 13, 83));  
   fCrateLoc.push_back(new PdataLoc("loStdcR", 2, 14, 83));  

   fCrateLoc.push_back(new PdataLoc("s0La", 3, 25, 16));
   fCrateLoc.push_back(new PdataLoc("s0Lb", 3, 25, 22));
   fCrateLoc.push_back(new PdataLoc("s0Ra", 1, 25, 0));
   fCrateLoc.push_back(new PdataLoc("s2Rb", 1, 25, 8));

   fCrateLoc.push_back(new PdataLoc("s2La", 3, 25, 16));
   fCrateLoc.push_back(new PdataLoc("s2Lb", 3, 25, 22));
   fCrateLoc.push_back(new PdataLoc("s2Ra", 1, 25, 0));
   fCrateLoc.push_back(new PdataLoc("s2Rb", 1, 25, 8));


   return 0;
}


//_____________________________________________________________________________
Int_t ParityData::Decode(const THaEvData& evdata)
{

  Int_t i;
  //  Int_t ldebug=1;

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

  cout << "Into Decode "<<evtype<<"  "<<fCrateLoc.size()<<endl;


  for( Iter_t p = fCrateLoc.begin(); p != fCrateLoc.end(); p++) {
    PdataLoc *dataloc = *p;

// bit pattern of triggers
    for (UInt_t i = 0; i < bits.GetNbits(); i++) {
      if ( dataloc->ThisIs(Form("bit%d",i+1)) ) TrigBits(i+1,dataloc);
    }


    if ( dataloc->ThisIs("upQadcL") ) upQadcL  = dataloc->Get();  
    if ( dataloc->ThisIs("upQadcR") ) upQadcR  = dataloc->Get();  
    if ( dataloc->ThisIs("loQadcL") ) loQadcL  = dataloc->Get();  
    if ( dataloc->ThisIs("loQadcR") ) loQadcR  = dataloc->Get();  

    if ( dataloc->ThisIs("upQtdcL") ) upQtdcL  = dataloc->Get();  
    if ( dataloc->ThisIs("upQtdcR") ) upQtdcR  = dataloc->Get();  
    if ( dataloc->ThisIs("loQtdcL") ) loQtdcL  = dataloc->Get();  
    if ( dataloc->ThisIs("loQtdcR") ) loQtdcR  = dataloc->Get();  

    if ( dataloc->ThisIs("atQadcL") ) atQadcL  = dataloc->Get();  
    if ( dataloc->ThisIs("atQadcR") ) atQadcR  = dataloc->Get();  
    if ( dataloc->ThisIs("atQtdcL") ) atQtdcL  = dataloc->Get();  
    if ( dataloc->ThisIs("atQtdcR") ) atQtdcR  = dataloc->Get();  

    if ( dataloc->ThisIs("upSadcL") ) upSadcL  = dataloc->Get();  
    if ( dataloc->ThisIs("upSadcR") ) upSadcR  = dataloc->Get();  
    if ( dataloc->ThisIs("upStdcL") ) upStdcL  = dataloc->Get();  
    if ( dataloc->ThisIs("upStdcR") ) upStdcR  = dataloc->Get();  

    if ( dataloc->ThisIs("loSadcL") ) loSadcL  = dataloc->Get();  
    if ( dataloc->ThisIs("loSadcR") ) loSadcR  = dataloc->Get();  
    if ( dataloc->ThisIs("loStdcL") ) loStdcL  = dataloc->Get();  
    if ( dataloc->ThisIs("loStdcR") ) loStdcR  = dataloc->Get();  

    if ( dataloc->ThisIs("s0La") ) s0La  = dataloc->Get();
    if ( dataloc->ThisIs("s0Lb") ) s0Lb  = dataloc->Get();
    if ( dataloc->ThisIs("s0Ra") ) s0Ra  = dataloc->Get();
    if ( dataloc->ThisIs("s0Rb") ) s0Rb  = dataloc->Get();

    if ( dataloc->ThisIs("s2La") ) s2La  = dataloc->Get();
    if ( dataloc->ThisIs("s2Lb") ) s2Lb  = dataloc->Get();
    if ( dataloc->ThisIs("s2Ra") ) s2Ra  = dataloc->Get();
    if ( dataloc->ThisIs("s2Rb") ) s2Rb  = dataloc->Get();

  }

  DoBpm();

  DoKine();

  if (PARDATA_PRINT) Print();

  return 0;
}


//_____________________________________________________________________________
Int_t ParityData::DoBpm( ) {


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

  Double_t ebeam = 1.05;     // E_beam - < dE/dx >

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

  if (pL > -999 && thL > -999 && phiL > -999 &&
      pL < 1e32 && thL < 1e32 && phiL < 1e32) {
    q2L = 2*ebeam*(pscale*pL)*(1-((costhL-sinthL*phiL)/TMath::Sqrt(1+thL*thL+phiL*phiL)));
  }

  q2R = -999;
  if (pR > -999 && thR > -999 && phiR > -999 &&
      pR < 1e32 && thR < 1e32 && phiR < 1e32) {
// note : sign convention requires costhR + sinthR*phiR 
//                                    not -
    q2R = 2*ebeam*(pscale*pR)*(1-((costhR+sinthR*phiR)/TMath::Sqrt(1+thR*thR+phiR*phiR)));
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
  cout << "PREX Upper Quartz ADC L, R "<<upQadcL<<"  "<<upQadcR<<endl;  
  cout << "PREX Lower Quartz ADC L, R "<<loQadcL<<"  "<<loQadcR<<endl;  
  cout << "PREX Upper Quartz TDC L, R "<<upQtdcL<<"  "<<upQtdcR<<endl;  
  cout << "PREX Lower Quartz ADC L, R "<<loQadcL<<"  "<<loQadcR<<endl;  
  cout << "PREX A_T Quartz ADC L, R "<<atQadcL<<"  "<<atQadcR<<endl;    
  cout << "PREX A_T Quartz TDC L, R "<<atQtdcL<<"  "<<atQtdcR<<endl;    
  cout << "PREX Upper Scint ADC L, R "<<upSadcL<<"  "<<upSadcR<<endl;   
  cout << "PREX Upper Scint TDC L, R "<<upStdcL<<"  "<<upStdcR<<endl;   
  cout << "PREX Lower Scint ADC L, R "<<loSadcL<<"  "<<loSadcR<<endl;   
  cout << "PREX Lower Scint TDC L, R "<<loStdcL<<"  "<<loStdcR<<endl;  
  cout << "event types,  CODA = "<<evtype;
  cout << "   bit pattern = "<<evtypebits<<endl;
  cout << "S0 detectors Left HRS  "<<s0La<<"  "<<s0Lb<<endl;
  cout << "S0 detectors Right HRS  "<<s0Ra<<"  "<<s0Rb<<endl;
  cout << "S2 detectors Left HRS  "<<s2La<<"  "<<s2Lb<<endl;
  cout << "S2 detectors Right HRS  "<<s2Ra<<"  "<<s2Rb<<endl;
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

