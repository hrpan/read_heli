#ifndef CalibStats_h
#define CalibStats_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <vector>

const Int_t kMaxinputHeaders2 = 1;

class CalibStats {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   Float_t         MaxQ;
   Float_t         MaxQ_2inchPMT;
   Float_t         NominalCharge;
   Float_t         Quadrant;
   Float_t         integralLiveTime_blocked_trigger_ms;
   Float_t         integralLiveTime_buffer_full_ms;
   Float_t         integralRunTime_ms;
   Int_t           nHit;

   TBranch        *b_MaxQ;   //!
   TBranch        *b_MaxQ_2inchPMT;   //!
   TBranch        *b_NominalCharge;   //!
   TBranch        *b_Quadrant;   //!
   TBranch        *b_integralLiveTime_blocked_trigger_ms;   //!
   TBranch        *b_integralLiveTime_buffer_full_ms;   //!
   TBranch        *b_integralRunTime_ms;   //!
   TBranch        *b_nHit;   //!

   CalibStats(TTree *tr);
   virtual ~CalibStats();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef CalibStats_cxx
CalibStats::CalibStats(TTree *tr)
{
   Init(tr);
}

CalibStats::~CalibStats()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CalibStats::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CalibStats::LoadTree(Long64_t entry)
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

void CalibStats::Init(TTree *tree)
{
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("MaxQ",kTRUE);
   fChain->SetBranchAddress("MaxQ", &MaxQ, &b_MaxQ);
   fChain->SetBranchStatus("MaxQ_2inchPMT",kTRUE);
   fChain->SetBranchAddress("MaxQ_2inchPMT", &MaxQ_2inchPMT, &b_MaxQ_2inchPMT);
   fChain->SetBranchStatus("NominalCharge",kTRUE);
   fChain->SetBranchAddress("NominalCharge", &NominalCharge, &b_NominalCharge);
   fChain->SetBranchStatus("Quadrant",kTRUE);
   fChain->SetBranchAddress("Quadrant", &Quadrant, &b_Quadrant);
   fChain->SetBranchStatus("integralLiveTime_blocked_trigger_ms",kTRUE);
   fChain->SetBranchAddress("integralLiveTime_blocked_trigger_ms", &integralLiveTime_blocked_trigger_ms, &b_integralLiveTime_blocked_trigger_ms);
   fChain->SetBranchStatus("integralLiveTime_buffer_full_ms",kTRUE);
   fChain->SetBranchAddress("integralLiveTime_buffer_full_ms", &integralLiveTime_buffer_full_ms, &b_integralLiveTime_buffer_full_ms);
   fChain->SetBranchStatus("integralRunTime_ms",kTRUE);
   fChain->SetBranchAddress("integralRunTime_ms", &integralRunTime_ms, &b_integralRunTime_ms);
   fChain->SetBranchStatus("nHit",kTRUE);
   fChain->SetBranchAddress("nHit", &nHit, &b_nHit);
   Notify();
}

Bool_t CalibStats::Notify()
{
   return kTRUE;
}

void CalibStats::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CalibStats::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifdef CalibStats_cxx
