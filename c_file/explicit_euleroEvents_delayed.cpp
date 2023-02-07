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
    mwSize i,j,ii;
    
    double *dx;
    
    double *oldState;
    double *currentState;
    double currentT;
    
    double **eventMatrix;
    
    double *phiOut;
    double *tOut;
    double *xFinale;
    
    double Vth;
    
    double period;
    
    double stopThreshold;
    
    double diff,tmp;
    
    mwSize *nEvent;
    
    mwSize minEventNumber;
    
    
    mwSize N;
    
    int *firstIndex;
    
    dynSys *vectorField;
    
    mwSize tmpIndex;
    
    
    /***** ADD FOR DELAYED SYSTEM ******/
    int Ndelay;
    double *delays;
    int maxDelays;
    
    int maxMemory; // = maxDelay x dt
    
    int *delayIndex;
    
    double *x0Del;
    
    double **xDelTot; // xDel matrix is (maxDelay x dt) x Nstates
    
    double **xDel; // xDel matrix is Ndelay x Nstates
    
    int currentPointer;
    
    /**********************************/
    
    double *I0;
    
    
    /* check for proper number of arguments */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 8 Input required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 3 Output required.");
    }
    
    initVectorField(&vectorField);
    
    
    Nstati =  mxGetScalar(prhs[0]);
    nStep = mxGetScalar(prhs[1]);
    dt = mxGetScalar(prhs[2]);
    oldState = mxGetPr(prhs[3]);
    Vth = mxGetScalar(prhs[4]);
    stopThreshold = mxGetScalar(prhs[5]);
    N =  mxGetScalar(prhs[6]);
    x0Del = mxGetPr(prhs[7]);
    
    Ndelay = vectorField->getNdelay();
    delays = (double *)mxMalloc(Ndelay*sizeof(double));
    
    vectorField->getDelays(delays);
    maxDelays = delays[0];
    for(i=1;i<Ndelay;i++)
        if(maxDelays < delays[i])
            maxDelays = delays[i];
    
    I0 = (double *)mxMalloc(Nstati*sizeof(double));
    
    for(i=0;i<Nstati;i++)
        I0[i] = 0;
    
    maxMemory = (int)(maxDelays/dt);
    
    delayIndex = (int *)mxMalloc(Ndelay*sizeof(int));
    
    xDel = (double **)mxMalloc(Ndelay*sizeof(double *));
    
    xDelTot = (double **)mxMalloc(maxMemory*sizeof(double *));
    
    for(i=0;i<Ndelay;i++)
        delayIndex[i] = (int)(((maxDelays-delays[i])/dt)-0.5);
    
    ii = 0;
    for(i=0;i<maxMemory;i++)
    {
        xDelTot[i] = (double *)mxMalloc(Nstati*sizeof(double));
        
        if(ii < Ndelay-1)
        {
            if(i >= delayIndex[ii+1])
                ii++;
        }
        
        for(j=0;j<Nstati;j++)
        {
            xDelTot[i][j] = x0Del[ii*Nstati+j];
        }
    }
    
    
    currentPointer = 0;
    
    currentState = (double *)mxMalloc(Nstati*sizeof(double));
    dx = (double *)mxMalloc(Nstati*sizeof(double));
    nEvent = (mwSize *)mxMalloc(N*sizeof(mwSize));
    eventMatrix = (double **)mxMalloc(N*sizeof(double *));
    
    firstIndex = (int *)mxMalloc((N+2*N*N)*sizeof(int));
    
    vectorField->getFirstIndex(firstIndex);
    
    
    for(j=0;j<N;j++)
    {
        eventMatrix[j] = (double *)mxMalloc(0);
        nEvent[j] = 0;
    }
    
    for(j=0;j<Nstati;j++)
    {
        currentState[j] = 0;
    }
    
    currentT = 0;
    
    
    for(i=1;i<nStep;i++)
    {
        for(j=0;j<Ndelay;j++)
            xDel[j] = xDelTot[delayIndex[j]];
        
        vectorField->getXdot(currentT,oldState,dx,I0,xDel);
        currentT += dt;
        
        for(j=0;j<Nstati;j++)
            currentState[j] = oldState[j]+dt*dx[j];
        
        
        vectorField->resetStates(currentState);
        
        for(j=0;j<N;j++)
        {
            
            tmpIndex = firstIndex[j];
            
            if((oldState[tmpIndex] < Vth) && (currentState[tmpIndex] > Vth))
            {
                nEvent[j]++;
                eventMatrix[j] = (double *)mxRealloc(eventMatrix[j],nEvent[j]*sizeof(double));
                eventMatrix[j][nEvent[j]-1] = (Vth-oldState[tmpIndex])/(currentState[tmpIndex]-oldState[tmpIndex])*(dt)+(currentT-dt);
            }
        }
        
        
        for(j=0;j<Nstati;j++)
        {
            oldState[j] = currentState[j];
        }
        
        for(j=0;j<Nstati;j++)
            xDelTot[currentPointer][j] = currentState[j];
        
        
        for(j=0;j<Ndelay;j++)
        {
            delayIndex[j]++;
            if(delayIndex[j] > maxMemory-1)
                delayIndex[j] = 0;
        }
        
        currentPointer++;
        if(currentPointer > maxMemory-1)
            currentPointer = 0;
        
    }
    
    
    minEventNumber = nEvent[0];
    
    for(i=1;i<N;i++)
        if(minEventNumber > nEvent[i])
            minEventNumber = nEvent[i];
    
    plhs[0] = mxCreateDoubleMatrix(minEventNumber,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(minEventNumber,N-1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(Nstati,1,mxREAL);
    
    tOut = mxGetPr(plhs[0]);
    phiOut = mxGetPr(plhs[1]);
    xFinale = mxGetPr(plhs[2]);
    
    if(minEventNumber > 1)
    {
        
        for(j=0;j<minEventNumber;j++)
        {   
            tOut[j] = eventMatrix[0][j];
            if(j == minEventNumber-1)
                period = eventMatrix[0][j] - eventMatrix[0][j-1];
            else
                period = eventMatrix[0][j+1] - eventMatrix[0][j];
            for(i=0;i<N-1;i++)
            {
                phiOut[i*minEventNumber+j] =  ((eventMatrix[i+1][j]-eventMatrix[0][j])/period);
            }
        }
    }
    
    
    for(j=0;j<Nstati;j++)
        xFinale[j] = currentState[j];
    
    mxFree(dx);
    mxFree(currentState);
    mxFree(nEvent);
    for(j=0;j<N;j++)
    {
        mxFree(eventMatrix[j]);
    }
    mxFree(eventMatrix);
    mxFree(firstIndex);
    
    
    mxFree(delayIndex);
    
    for(j=0;j<maxMemory;j++)
        mxFree(xDelTot[j]);
    
    mxFree(xDelTot);
    mxFree(xDel);
    mxFree(delays);
    
    
    delete(vectorField);
    mxFree(I0);
    
}
