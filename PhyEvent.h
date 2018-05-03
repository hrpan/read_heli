#ifndef PHYEVENT_H
#define PHYEVENT_H

#include "TTimeStamp.h"
//#include "CalibReadoutHeader.C"
//#include "AdScaled.C"

const double startT=1324339200;

class PhyEvent{
    public :

		// event context
		int         site;
		int         detector;
		
		// most frequently used info (from RecData)
		float       e;
		float		x,y,z;
		float       nPESum;
		long long	t;
		long long	trigTime_s;
		long long	trigTime_ns;
		double	 	dtnAD;
		double	 	dtnWP;
		bool		passMu;
		int 		nHit;
		int			date;

		vector<double> nPE;
		vector<double> dtlSH;

		vector<float> nTag_e;	
		vector<float> nTag_dt;	
		
		//Tags
		bool		isGood;
		bool		isDelay;
		bool		isPrompt;
		bool		isWPMu;
		bool		isADMu;
		bool		isSHMu;

		PhyEvent(CalibStats &calst, AdSimple &adrec);

};

PhyEvent::PhyEvent(CalibStats &calst, AdSimple &adrec) {

	site = adrec.site;
	detector = adrec.detector; 
	e = adrec.energy;
	x = adrec.x;
	y = adrec.y;
	z = adrec.z;
	
	trigTime_s = adrec.triggerTimeSec;
	trigTime_ns = adrec.triggerTimeNanoSec;

	t = trigTime_s*1000000000+trigTime_ns;

	date = int((trigTime_s-startT)/86400);

	nHit = calst.nHit;
	nPESum = calst.NominalCharge;

	dtnAD=-1;
	dtnWP=-1;

	bool isFlash = calst.Quadrant*calst.Quadrant+(calst.MaxQ/0.45)*(calst.MaxQ/0.45)>1|| (calst.MaxQ_2inchPMT>100&&e<20.); 		
	bool isCross = (adrec.triggerType==0x10000002);
	bool isForce = (adrec.triggerType==0x10000020);
	bool isRPCNS = (adrec.triggerType==0x10020000);
	isGood = (!isFlash)&&(!isCross)&&(!isForce)&&(!isRPCNS);
	isDelay = (e<12.0)&&(e>1.5);
	isPrompt = (e<12.0)&&(e>0.7);
	isWPMu = (detector==5||detector==6)&&calst.nHit>12;
	isADMu = (detector<=4)&&(calst.NominalCharge>3000);
	isSHMu = (detector<=4)&&(calst.NominalCharge>300000);
}

#endif
