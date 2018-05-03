
#ifndef AdSimple_h
#define AdSimple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
using namespace std;

const Int_t kMaxinputHeaders = 2;

class AdSimple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   Short_t         site;
   Short_t         detector;
   UInt_t          triggerType;
   Int_t           triggerTimeSec;
   Int_t           triggerTimeNanoSec;
   Float_t         energy;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         t;

   TBranch        *b_Rec_AdSimple_site;   //!
   TBranch        *b_Rec_AdSimple_detector;   //!
   TBranch        *b_Rec_AdSimple_triggerType;   //!
   TBranch        *b_Rec_AdSimple_triggerTimeSec;   //!
   TBranch        *b_Rec_AdSimple_triggerTimeNanoSec;   //!
   TBranch        *b_Rec_AdSimple_energy;   //!
   TBranch        *b_Rec_AdSimple_x;   //!
   TBranch        *b_Rec_AdSimple_y;   //!
   TBranch        *b_Rec_AdSimple_z;   //!
   TBranch        *b_Rec_AdSimple_t;   //!

   AdSimple(TTree *tr);
   virtual ~AdSimple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AdSimple_cxx
AdSimple::AdSimple(TTree *tr) 
{
   Init(tr);
}

AdSimple::~AdSimple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AdSimple::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AdSimple::LoadTree(Long64_t entry)
{
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AdSimple::Init(TTree *tree)
{
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("site",kTRUE);
   fChain->SetBranchAddress("site", &site, &b_Rec_AdSimple_site);
   fChain->SetBranchStatus("detector",kTRUE);
   fChain->SetBranchAddress("detector", &detector, &b_Rec_AdSimple_detector);
   fChain->SetBranchStatus("triggerType",kTRUE);
   fChain->SetBranchAddress("triggerType", &triggerType, &b_Rec_AdSimple_triggerType);
   fChain->SetBranchStatus("triggerTimeSec",kTRUE);
   fChain->SetBranchAddress("triggerTimeSec", &triggerTimeSec, &b_Rec_AdSimple_triggerTimeSec);
   fChain->SetBranchStatus("triggerTimeNanoSec",kTRUE);
   fChain->SetBranchAddress("triggerTimeNanoSec", &triggerTimeNanoSec, &b_Rec_AdSimple_triggerTimeNanoSec);
   fChain->SetBranchStatus("energy",kTRUE);
   fChain->SetBranchAddress("energy", &energy, &b_Rec_AdSimple_energy);
   fChain->SetBranchStatus("x",kTRUE);
   fChain->SetBranchAddress("x", &x, &b_Rec_AdSimple_x);
   fChain->SetBranchStatus("y",kTRUE);
   fChain->SetBranchAddress("y", &y, &b_Rec_AdSimple_y);
   fChain->SetBranchStatus("z",kTRUE);
   fChain->SetBranchAddress("z", &z, &b_Rec_AdSimple_z);
   fChain->SetBranchStatus("t",kTRUE);
   fChain->SetBranchAddress("t", &t, &b_Rec_AdSimple_t);
   Notify();
}

Bool_t AdSimple::Notify()
{

   return kTRUE;
}

void AdSimple::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AdSimple::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifdef AdSimple_cxx
