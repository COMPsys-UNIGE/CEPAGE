#include "mex.h"
#include "math.h"
#include "vectorField.hpp"


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    mwSize nStep;   /* Number of step */
    double dt;      /* Integration step */
    mwSize Nstati;       /* Number of state variables*/
    mwSize i,j;
    double *ptr;
    
    double *x0;     /* Inital condition */
    double *xOutMatrix; /* Out vector */
    double *dx;
    
    
    
    double t;
    
    dynSys *vectorField;
    double *I0;
    
    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 4 Input required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 1 Output required.");
    }
    
    
    initVectorField(&vectorField);
    
    Nstati =  mxGetScalar(prhs[0]);
    nStep = mxGetScalar(prhs[1]);
    dt = mxGetScalar(prhs[2]);
    x0 = mxGetPr(prhs[3]);
    
    dx = (double *)mxMalloc(Nstati*sizeof(double));
    
    I0 = (double *)mxMalloc(Nstati*sizeof(double));
    
    for(i=0;i<Nstati;i++)
        I0[i] = 0;
    
    plhs[0] = mxCreateDoubleMatrix(Nstati,nStep,mxREAL);
    
    
    xOutMatrix  = mxGetPr(plhs[0]);
    for(i=0;i<Nstati;i++)
        xOutMatrix[i] = x0[i];
    
    
    t = 0;
    
    for(i=1;i<nStep;i++)
    {
        
        ptr = xOutMatrix+i*Nstati;
        vectorField->getXdot(t,ptr-Nstati,dx,I0);
        
        for(j=0;j<Nstati;j++)
        {
            ptr[j] = ptr[j-Nstati]+dt*dx[j];
        }
        
        vectorField->resetStates(ptr);
        
        
        t+=dt;
        
    }
    
    mxFree(dx);
    
    delete(vectorField);
    mxFree(I0);
}
