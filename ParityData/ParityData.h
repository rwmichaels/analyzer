#ifndef ROOT_ParityData
#define ROOT_ParityData

//////////////////////////////////////////////////////////////////////////
//
// ParityData
//
//////////////////////////////////////////////////////////////////////////

#define ILEFT   1
#define IRIGHT  2

#include "THaApparatus.h"
#include "TBits.h"
#include "TString.h"
#include "THaVar.h"
//#include "TH1.h"
//#include "TH2.h"
#include <vector>
#include <iostream>
#include <string>
#include <cstring>

using namespace std;

class TH1;
class THaScaler;

class PdataLoc {
// Utility class used by ParityData.
// Data location, either in (crates, slots, channel), or
// relative to a unique header in a crate or in an event.
   static const Int_t MxHits=16;
 public:
   // c'tor for (crate,slot,channel) selection
   PdataLoc ( const char* nm, Int_t cra, Int_t slo, Int_t cha ) :
     crate(cra), slot(slo), chan(cha), header(0), ntoskip(0), 
     name(nm), search_choice(0) {}
   // c'tor for header search (note, the only diff to above is 3rd arg is UInt_t)
   PdataLoc ( const char* nm, Int_t cra, UInt_t head, Int_t skip ) :
     crate(cra), slot(0), chan(0), header(head), ntoskip(skip),
     name(nm), search_choice(1) {}
   Bool_t IsSlot() { return (search_choice == 0); }
   void Clear() { ndata=0;  loaded_once = kFALSE; }
   void Print() { 
     cout << "DataLoc "<<name<<"   crate "<<crate;
     if (IsSlot()) {
       cout << "  slot "<<slot<<"   chan "<<chan<<endl;
     } else {
       cout << "  header "<<header<<"  ntoskip "<<ntoskip<<endl;
     }
   } 
   void Load(UInt_t data) {
     if (ndata<MxHits) rdata[ndata++]=data;
     loaded_once = kTRUE;
   }
   Bool_t DidLoad() { return loaded_once; }
   Int_t NumHits() { return ndata; }
   UInt_t Get(Int_t i=0) { 
     return (i >= 0 && ndata > i) ? rdata[i] : 0; }
   Bool_t ThisIs(const char* aname) { return (strstr(name.c_str(),aname) != 0);}
   ~PdataLoc() {}
   
   Int_t  crate, slot, chan;   // where to look in crates
   UInt_t header;              // header (unique either in data or in crate)
   Int_t ntoskip;              // how far to skip beyond header
   const std::string name;     // name of the variable in global list.
   
   UInt_t rdata[MxHits];       //[ndata] raw data (to accom. multihit chanl)
   Int_t  ndata;               // number of relevant entries
 private:
   Int_t  search_choice;       // whether to search in crates or rel. to header
   Bool_t loaded_once;
   PdataLoc();
   PdataLoc(const PdataLoc& dataloc);
   PdataLoc& operator=(const PdataLoc& dataloc);
};


class ParityData : public THaApparatus {
  
public:
   ParityData( const char* name="D", const char* description="" );
   virtual ~ParityData();

   virtual EStatus Init( const TDatime& run_time );
   virtual Int_t   End(THaRunBase* r=0);
   virtual void    WriteHist(); 
   virtual Int_t   Reconstruct() { return 0; }
   virtual Int_t   Decode( const THaEvData& );

protected:
   virtual PdataLoc* DefineChannel(PdataLoc*, EMode, const char* desc="automatic");

private:
   TBits  bits;
   Int_t  *trigcnt;

   UInt_t evtypebits, evtype;
   Int_t WindowSize;   // this can be at most 500, per design of FADC

   std::ofstream *fDebugFile;  // debug output

   // Global variable definitions for this class
   vector<string> GloVars;

   Int_t Nvars=0;
   Double_t *dvars;
   
   // Q^2 on left (L) and right (R) HRS.  (can also use Kinematics class)
   Double_t q2L,q2R;

   std::vector < PdataLoc* > fCrateLoc;   // Raw Data locations by crate, slot, channel
   std::vector < PdataLoc* > fWordLoc;    // Raw Data locations relative to header word
   std::vector<PdataLoc*> fDataLocs;  // these are all the WordLocs and CrateLocs combined

   Int_t IptrFadcL, IptrFadcR;   // pointer to fCrateLoc where FADC is

   // Quartz FADC data arrays -- for raw mode
   Double_t *QLfadclo, *QLfadcup;  // L-HRS lower, upper 
   Double_t *QRfadclo, *QRfadcup;  // R-HRS lower, upper 
   // A_T detector FADC data arrays -- for raw mode
   Double_t *ATLfadc1, *ATLfadc2;  // L-HRS detectors 1 and 2
   Double_t *ATRfadc1, *ATRfadc2;  // R-HRS detectors 1 and 2
   // Quartz FADC data -- for pulse integral mode ("I")
   Double_t QLIfadclo,  QLIfadcup;    // L-HRS lower, upper 
   Double_t QRIfadclo,  QRIfadcup;    // R-HRS lower, upper 
   Double_t QLIftimelo,  QLIftimeup;  // L-HRS time over thr, L-HRS lower, upper
   Double_t QRIftimelo,  QRIftimeup;  // R-HRS time over thr, L-HRS lower, upper
   // A_T detector FADC data -- for pulse integral mode ("I")
   Double_t ATLIfadc1,  ATLIfadc2;    // integrated adc L-HRS AT#1 & 2
   Double_t ATRIfadc1,  ATRIfadc2;    // integrated adc R-HRS AT#1 & 2
   Double_t ATLIftime1, ATLIftime2;   // time over thr, L-HRS AT#1 & 2
   Double_t ATRIftime1, ATRIftime2;   // time over thr, R-HRS AT#1 & 2
   
   virtual void Clear( Option_t* opt="" );
   virtual void Print( Option_t* opt="" ) const;
   std::vector<TH1* > hist;
   Int_t DefaultMap();
   void  PrintMap();
   Int_t DoBpm();
   Int_t DoKine();
   void TrigBits(UInt_t ibit, PdataLoc *dataloc);
   Int_t   DecodeFadc(Int_t iarm, Int_t iptr, const THaEvData& );
   static std::vector<std::string> vsplit(const std::string& s);
   Int_t SetupParData( const TDatime* runTime = NULL, EMode mode = kDefine );
   virtual void BookHist(); 

   static UInt_t header_str_to_base16(const char* hdr);

   static const int PARDATA_VERBOSE = 1;

   ClassDef(ParityData,0)  
};

#endif








