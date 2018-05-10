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
#include "math.h"
#include "outfile.h"
#include <cstdlib>
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

	deque<deque<PhyEvent> > evtbuf(4),candid(4),muon(4);

	long long mu[4]={0,0,0,0};
	long long muSH[4]={0,0,0,0};
	int det,site,ads;

	cout << "nevent: " << nevent << endl;

	for( size_t jentry = 0; jentry < nevent; ++jentry ) {
		if(jentry % (nevent/1000) == 0){
			printf("\r%zu/%zu ( %.1f%% )",jentry+1,nevent,100.0*(jentry+1)/nevent);
			fflush(stdout);
		}
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
			for(int i=evtbuf[det].size()-1;i>=0&&evtbuf[det][i].dtnAD==-1;--i)
				evtbuf[det][i].dtnAD=evt.t-evtbuf[det][i].t;
			if(evt.isSHMu)
				muSH[det]=evt.t;
			mu[det]=max(mu[det],evt.t+vAD);
			muon[det].push_back(evt);
			continue;
		}else if(evt.isWPMu){
			for(int i=0;i<ads;++i){
				for(int j=evtbuf[i].size()-1;j>=0&&evtbuf[i][j].dtnWP==-1;--j)	
					evtbuf[i][j].dtnWP=evt.t-evtbuf[i][j].t;
				mu[i]=max(mu[i],evt.t+vWS);
			}
			continue;
		}

		for(int i=0;i<ads;++i){
//			while(!evtbuf[i].empty() && evtbuf[i][0].dtnAD!=-1 && evtbuf[i][0].dtnWP!=-1){
			while(!evtbuf[i].empty() && evtbuf[i][0].dtnWP!=-1){
				//bool passnMu = evtbuf[i][0].dtnAD>vpre && evtbuf[i][0].dtnWP>vpre;
				bool passnMu = evtbuf[i][0].dtnWP>vpre;
				evtbuf[i][0].passMu = passnMu && evtbuf[i][0].passMu;
				if(evtbuf[i][0].passMu)
					candid[i].push_back(evtbuf[i][0]);

//---------------------------------------------------------------------------------
//					He/Li SPECTRUM
//---------------------------------------------------------------------------------

				size_t c_size = candid[i].size();
				
				while(c_size>3 && candid[i][c_size-1].t-candid[i][2].t>capT){
					bool delay = candid[i][2].isDelay;
					bool cap = (candid[i][2].t - candid[i][1].t < capT &&
								candid[i][2].t - candid[i][1].t > capT2);
					bool multP = candid[i][2].t - candid[i][0].t > 2 * capT;
					bool multD = true;

					for(size_t dIDX=3;dIDX<c_size&&multD;++dIDX)
						multD = !(candid[i][dIDX].isDelay && (candid[i][dIDX].t - candid[i][2].t) < capT);
					
					float dx = candid[i][2].x - candid[i][1].x;
					float dy = candid[i][2].y - candid[i][1].y;
					float dz = candid[i][2].z - candid[i][1].z;

					bool dist = sqrt( dx * dx + dy * dy + dz * dz ) < distCut;
					
					if( delay && cap && multP && multD && dist )
						out->FillIBD(candid[i][1],candid[i][2],muon[i]);

					candid[i].pop_front();
					--c_size;
				}
			
				evtbuf[i].pop_front();
			}
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

			if(!evt.passMu) continue;

			evtbuf[det].push_back(evt);
		}
	} // end of event loop

	outf->Write();

	return 0;
}

