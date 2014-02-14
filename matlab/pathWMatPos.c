#include "math.h"
#include "mex.h"

#define MATIND(i,j,nR) ((j)*(nR)+(i))

void wmat(double *w,double *dwd,double *dwhv,int n,double cdi,double cax) {
    double v,t;
    int i,j;
    
    w[0]=1;
    dwd[0]=0;
    dwhv[0]=0;
    t=1;
    for(i=1;i<n;i++) {
        v=w[MATIND(i-1,0,n)]*cax;
        w[MATIND(i,0,n)]=v;
        w[MATIND(0,i,n)]=v;

        dwd[MATIND(i,0,n)]=0;
        dwd[MATIND(0,i,n)]=0;

        v=i*t;
        dwhv[MATIND(i,0,n)]=v;
        dwhv[MATIND(0,i,n)]=v;
        t=t*cax;
    }

    for(i=1;i<n;i++) {
        w[MATIND(i,i,n)]=cdi*cdi*w[MATIND(i-1,i-1,n)]+2*cax*w[MATIND(i-1,i,n)];
        dwd[MATIND(i,i,n)]=cdi*cdi*dwd[MATIND(i-1,i-1,n)]+2*cax*dwd[MATIND(i-1,i,n)]+2*cdi*w[MATIND(i-1,i-1,n)];
        dwhv[MATIND(i,i,n)]=2*(w[MATIND(i-1,i,n)]+cax*dwhv[MATIND(i-1,i,n)])+cdi*cdi*dwhv[MATIND(i-1,i-1,n)];
        for(j=i+1;j<n;j++) {
	    v=cdi*cdi*w[MATIND(i-1,j-1,n)]+cax*(w[MATIND(i-1,j,n)]+w[MATIND(i,j-1,n)]);
            w[MATIND(i,j,n)]=v;
            w[MATIND(j,i,n)]=v;

            v=cdi*cdi*dwd[MATIND(i-1,j-1,n)]+cax*(dwd[MATIND(i-1,j,n)]+dwd[MATIND(i,j-1,n)])+2*cdi*w[MATIND(i-1,j-1,n)];
            dwd[MATIND(i,j,n)]=v;
            dwd[MATIND(j,i,n)]=v;

            v=w[MATIND(i-1,j,n)]+w[MATIND(i,j-1,n)]+cax*(dwhv[MATIND(i-1,j,n)]+dwhv[MATIND(i,j-1,n)])+cdi*cdi*dwhv[MATIND(i-1,j-1,n)];
            dwhv[MATIND(i,j,n)]=v;
            dwhv[MATIND(j,i,n)]=v;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs<3)
        mexErrMsgTxt("Too few input parameters.");
    if(nrhs>3)
        mexErrMsgTxt("Too many input arguments.");
    if(nlhs>3)
        mexErrMsgTxt("Too many output arguments.");

    int n;
    double cdi,cax;
    
    n=(int)mxGetScalar(prhs[0]);
    cdi=sqrt((double)mxGetScalar(prhs[1]));
    cax=(double)mxGetScalar(prhs[2]);
    
    plhs[0]=mxCreateDoubleMatrix(n,n,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n,n,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n,n,mxREAL);
    double *w=mxGetPr(plhs[0]);
    double *dwd=mxGetPr(plhs[1]);
    double *dwhv=mxGetPr(plhs[2]);
    
    wmat(w,dwd,dwhv,n,cdi,cax);
    return;
}

