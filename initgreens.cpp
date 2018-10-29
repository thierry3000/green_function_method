/****************************************************************************
initgreens - initial tissue source strengths, given uniform solute field, initial g0
initial vessel source strengths based on uniform efflux rate from all vessels
use this version only if values are not available from previous call to greens
TWS Jan 08
Version 2.0, May 1, 2010.
Version 3.0, May 17, 2011.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void tissrate(int nsp, float *p, float *mtiss, float *mptiss);

void initgreens()
{
	extern int nnt,nnv,nsp;
	extern int *lowflow,*mainseg,*permsolute;
	extern int *oxygen,*diffsolute; //added April 2010

	extern float vol,errfac,tlength;
	extern float tlengthq,tlengthqhd;//added 8/09

	extern float *mtiss,*mptiss,*epsvessel,*epstissue,*eps,*errvessel,*errtissue,*pinit,*p;
	extern float *g0,*qtsum,*pref,*ds,*qq,*hd;
	extern float **qt,**qv,**pt,**pv,**tissparam,*mptissref;

	int isp,i,itp;

	for(isp=1; isp<=nsp; isp++)	pinit[isp] = g0[isp];
	tissrate(nsp,pinit,mtiss,mptiss);
	for(isp=1; isp<=nsp; isp++){
		mptissref[isp] = mptiss[isp];
		qtsum[isp] = 0.;
		for(itp=1; itp<=nnt; itp++){
			qt[itp][isp] = mtiss[isp]*vol;
			pt[itp][isp] = pinit[isp];
			qtsum[isp] += qt[itp][isp];
		}
		for(i=1; i<=nnv; i++){
			qv[i][isp] = 0.;
			pv[i][isp] = 0.;
			if(permsolute[isp] == 1){
				if(oxygen[isp] == 1){
					qv[i][isp] = -qtsum[isp]*ds[mainseg[i]]*qq[mainseg[i]]*hd[mainseg[i]]/tlengthqhd;//modified 8/10
					if(lowflow[mainseg[i]] == 1) qv[i][isp] = 0.;  //low q*hd
				}
				else qv[i][isp] = -qtsum[isp]*ds[mainseg[i]]*qq[mainseg[i]]/tlengthq;//modified 8/09
				pv[i][isp] = pinit[isp];
			}
		}
	}
//set error bounds, proportional to errfac, relax critera for vessel and tissue source strengths - April 2010
	for(isp=1; isp<=nsp; isp++){
		epsvessel[isp] = fabs(qtsum[isp])/nnv*errfac*10.;
		epstissue[isp] = tissparam[1][isp]*vol*errfac*10.;  //June 2009 - requires tissparam[1][isp] to be max rate
		eps[isp] = pref[isp]*errfac;
	}
}