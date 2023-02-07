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
    mwSize i,j,k;
    double *ptr;
    
    double *x0;     /* Inital condition */
    double *xOutMatrix; /* Out vector */
    double *dx;


    double t;
   
    
    
    /***** ADD FOR DELAYED SYSTEM ******/
    int Ndelay;
    double *delays;
            
    int *delayIndex;

    double *x0Del;
    
    double **xDel; // xDel matrix is Ndelay x Nstates
    
    /**********************************/
    
    dynSys *vectorField;
    
        double *I0;  
    

    
    /* check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 5 Input required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 1 Output required.");
    }
    

    
    initVectorField(&vectorField);
   
    
    Ndelay = vectorField->getNdelay();
    delays = (double *)mxMalloc(Ndelay*sizeof(double));
    
    
    vectorField->getDelays(delays);
    
    delayIndex = (int *)mxMalloc(Ndelay*sizeof(int));

    
    
    Nstati =  mxGetScalar(prhs[0]);   
    nStep = mxGetScalar(prhs[1]);    
    dt = mxGetScalar(prhs[2]);
    x0 = mxGetPr(prhs[3]);
    x0Del = mxGetPr(prhs[4]);
    
        I0 = (double *)mxMalloc(Nstati*sizeof(double));
    
    for(i=0;i<Nstati;i++)
        I0[i] = 0;
    
    dx = (double *)mxMalloc(Nstati*sizeof(double));
        
    plhs[0] = mxCreateDoubleMatrix(Nstati,nStep,mxREAL);
    
    xDel = (double **)mxMalloc(Ndelay*sizeof(double *));

    for(i=0;i<Ndelay;i++)
    {
        xDel[i] = (double *)mxMalloc(Nstati*sizeof(double));
        for(j=0;j<Nstati;j++)
        {
            xDel[i][j] = x0Del[i*Nstati+j];
        }
    }

    /* inizializzare indici delay a valori negativi */
    for(i=0;i<Ndelay;i++)
        delayIndex[i] = -(int)((delays[i]/dt)-0.5)+1; 

    xOutMatrix  = mxGetPr(plhs[0]);
    for(i=0;i<Nstati;i++)
      xOutMatrix[i] = x0[i];

    t = 0;
    
    
    for(i=1;i<nStep;i++)
    {
        

        ptr = xOutMatrix+i*Nstati;
        
        vectorField->getXdot(t,ptr-Nstati,dx,I0,xDel);
 
        for(j=0;j<Nstati;j++)
        {
            ptr[j] = ptr[j-Nstati]+dt*dx[j];
        }

        vectorField->resetStates(ptr);
        
        /* Update delayIndeces */
        for(k=0;k<Ndelay;k++)
            delayIndex[k]++;
        
        /* Get the right delayed states */
        for(k=0;k<Ndelay;k++)
            if(delayIndex[k] >= 0)
                for(j=0;j<Nstati;j++)
                {
                    xDel[k][j] = xOutMatrix[delayIndex[k]*Nstati+j];
                }
        
        t+=dt;
                
    }
        
    for(i=0;i<Ndelay;i++)
            mxFree(xDel[i]);

    mxFree(xDel);

    mxFree(dx);
    mxFree(delays);
    mxFree(delayIndex);
    
    delete(vectorField);
    mxFree(I0);

}
