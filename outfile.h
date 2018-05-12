#include <TFile.h>
#ifndef outfile_h
#define outfile_h
#include "PhyEvent.h"
#include <cmath>
class outfile {

public : 
//	TFile *file;// = new TFile("out.root","RECREATE");
	TTree *tr;// = new TTree("tr","");

	outfile();
	int site;
	int detector;
	int date;
	float ed,ep;
	float xp,yp,zp,xd,yd,zd;
	float dt;
	float dist;
	vector<double> dtlSH,nPESum;
	vector<int> nTag_n;
	vector<float>  nTag_e, nTag_dt;
	void FillIBD(PhyEvent &evtp,PhyEvent &evtd);
	void FillIBD(PhyEvent &evtp,PhyEvent &evtd, deque<PhyEvent> &mu);

};

outfile::outfile() {

//	file = new TFile(filename,"RECREATE");
	tr = new TTree("Heli",""); 


	tr->Branch("site",&site);	
	tr->Branch("detector",&detector);
	tr->Branch("date",&date);
	tr->Branch("ep",&ep);
	tr->Branch("ed",&ed);
	tr->Branch("dt",&dt);
	tr->Branch("xp",&xp);
	tr->Branch("yp",&yp);
	tr->Branch("zp",&zp);
	tr->Branch("xd",&xd);
	tr->Branch("yd",&yd);
	tr->Branch("zd",&zd);
	tr->Branch("dist",&dist);
	tr->Branch("dtlSH",&dtlSH);
	tr->Branch("nPESum",&nPESum);
	tr->Branch("nTag_n",&nTag_n);
	tr->Branch("nTag_e",&nTag_e);
	tr->Branch("nTag_dt",&nTag_dt);
}


void outfile::FillIBD(PhyEvent &evtp,PhyEvent &evtd) {
	site = evtp.site;
	detector = evtp.detector;
	date = evtp.date;	
	ep = evtp.e;
	ed = evtd.e;
	dt = double(evtd.t-evtp.t)/1000.0;	
	xp = evtp.x;
	yp = evtp.y;
	zp = evtp.z;
	xd = evtd.x;
	yd = evtd.y;
	zd = evtd.z;
	dist = sqrt( (xp-xd) * (xp-xd) + (yp-yd) * (yp-yd) + (zp-zd) * (zp-zd) );
	dtlSH = evtp.dtlSH;
	nPESum = evtp.nPE;
	tr->Fill();
}

void outfile::FillIBD(PhyEvent &evtp, PhyEvent &evtd, deque<PhyEvent> &mu){
	site = evtp.site;
	detector = evtp.detector;
	date = evtp.date;	
	ep = evtp.e;
	ed = evtd.e;
	dt = double(evtd.t-evtp.t)/1000.0;	
	xp = evtp.x;
	yp = evtp.y;
	zp = evtp.z;
	xd = evtd.x;
	yd = evtd.y;
	zd = evtd.z;
	dist = sqrt( (xp-xd) * (xp-xd) + (yp-yd) * (yp-yd) + (zp-zd) * (zp-zd) );

	dtlSH.clear();
	nPESum.clear();
	nTag_e.clear();
	nTag_dt.clear();
	nTag_n.clear();

	for(size_t i=0;i<mu.size();++i){
		dtlSH.push_back(double(evtp.t-mu[i].t)/1e6);
		nPESum.push_back(mu[i].nPESum);
		nTag_n.push_back(mu[i].nTag_e.size());
		nTag_dt.insert(nTag_dt.end(),mu[i].nTag_dt.begin(),mu[i].nTag_dt.end());
		nTag_e.insert(nTag_e.end(),mu[i].nTag_e.begin(),mu[i].nTag_e.end());
	}

	tr->Fill();

}
#endif
