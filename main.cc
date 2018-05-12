#include <TROOT.h>
#include "AdSimple.C"
#include "CalibStats.C"
#include <iostream>
#include "TH1F.h"
#include "PhyEvent.h"
#include <deque>
#include <TFile.h>
#include <algorithm>
#include <fstream>
//#include "Hists.h"
#include "TChain.h"
#include "outfile.h"
#include <cstdlib>
#include <cmath>
#define startT 1324339200 
#define maxT 5000000000
using namespace std;

const long long vSH=400400000;
const long long vAD=1400000;
const long long vWS=600000;
const long long vpre=2000;
const long long capT=400000;
const long long capT2=1000;

const double distCut = 1500;

const double nTag_max = 800;

void printProgress(size_t current, size_t total);
void selectIBD(deque<PhyEvent> &candid,deque<PhyEvent> &muons,outfile *out);

int main(int argc,char** argv) {
	gROOT->ProcessLine("#include <vector>");
	if( argc != 3 ) {
    	cerr << "# of arguments error." << endl;
        return 1;
	}

	ifstream list(argv[1]);
	string fname;

	TChain *adch = new TChain("/Event/Rec/AdSimple");
	TChain *csch = new TChain("/Event/Data/CalibStats");

	while(getline(list,fname)){
		cout << fname << endl;
		adch->Add(fname.c_str());
		csch->Add(fname.c_str());
	}

	AdSimple adrec(adch);
	CalibStats calst(csch);

	TFile *outf = new TFile(argv[2],"RECREATE");
	outf->cd();
	outfile *out;

	size_t nevent = adrec.fChain->GetEntries();

	deque<PhyEvent> evtbuf[4],candid[4],muon[4];

	long long mu[4]={0,0,0,0};
	int det,site,ads;

	cout << "nevent: " << nevent << endl;

	for( size_t jentry = 0; jentry < nevent; ++jentry ) {

		if(jentry % (nevent/1000) == 0)
			printProgress(jentry, nevent);
		
		adrec.GetEntry(jentry);
		calst.GetEntry(jentry);

		if(jentry==0){
			site = adrec.site;
			if(site==1||site==2)
				ads=2;
			else
				ads=4;
			out = new outfile();
		}

		PhyEvent evt(calst,adrec);
		det = evt.detector-1;
		
		if(det==6) continue;

		if(!evt.isGood) continue;

		if(evt.isADMu){
			mu[det]=max(mu[det],evt.t+vAD);
			muon[det].push_back(evt);
			continue;
		}else if(evt.isWPMu){
			for(int i=0;i<ads;++i){
				bool flag = false;
				while(!evtbuf[i].empty()){
					bool pass_n = evt.t-evtbuf[i][0].t > vpre;
					flag = flag || pass_n;
					if(pass_n)
						candid[i].push_back(evtbuf[i][0]);
					
					evtbuf[i].pop_front();
				}
				if(flag)
					selectIBD(candid[i],muon[i],out);
//				for(int j=evtbuf[i].size()-1;j>=0&&evtbuf[i][j].dtnWP==-1;--j)	
//					evtbuf[i][j].dtnWP=evt.t-evtbuf[i][j].t;
				mu[i]=max(mu[i],evt.t+vWS);
			}
			continue;
		}

		if(det>=ads) 
			continue;

		if(evt.isPrompt){

			evt.passMu = evt.t > mu[det];

			while(!muon[det].empty()&&evt.t-muon[det][0].t>maxT)
				muon[det].pop_front();

			if(evt.isDelay){
				for(size_t i=muon[det].size();i>0;--i){
					float nTag_dt = float(evt.t-muon[det][i-1].t) / 1e3;
					if( nTag_dt > nTag_max )
						break;
					muon[det][i-1].nTag_e.push_back(evt.e);	
					muon[det][i-1].nTag_dt.push_back(nTag_dt);	
				}
			}

/*
			for(size_t i=0;i<muon[det].size();++i){
				evt.nPE.push_back(muon[det][i].nPESum);
				evt.dtlSH.push_back(double(evt.t-muon[det][i].t)/1000000.0);
			}
*/

			if(evt.passMu) 	evtbuf[det].push_back(evt);
		}
	} // end of event loop

	outf->Write();

	return 0;
}

void printProgress(size_t current, size_t total){

	printf("\r%zu/%zu ( %.1f%% )",current+1,total,100.0*(current+1)/total);
	fflush(stdout);

}

void selectIBD(deque<PhyEvent> &candid,deque<PhyEvent> &muons,outfile *out){

	size_t c_size = candid.size();
	
	while(c_size>3){

		PhyEvent &evt_p = candid[1];
		PhyEvent &evt_d = candid[2];

		if(candid[c_size-1].t - evt_d.t < capT) break;

		bool delay = evt_d.isDelay;
		bool cap = (evt_d.t - evt_p.t < capT && evt_d.t - evt_p.t > capT2);
		bool multP = evt_d.t - candid[0].t > 2 * capT;
		bool multD = true;
		for(size_t dIDX=3;dIDX<c_size&&multD;++dIDX)
			multD = !(candid[dIDX].isDelay && (candid[dIDX].t - evt_d.t) < capT);
		
		float dx = evt_d.x - evt_p.x;
		float dy = evt_d.y - evt_p.y;
		float dz = evt_d.z - evt_p.z;

		bool dist = sqrt( dx * dx + dy * dy + dz * dz ) < distCut;
		
		if( delay && cap && multP && multD && dist )
			out->FillIBD(evt_p,evt_d,muons);

		candid.pop_front();
		--c_size;

		
	}

}
