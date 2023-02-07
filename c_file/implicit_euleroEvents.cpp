#include "mex.h"
#include "math.h"
#include "vectorField.hpp"
#include <stdio.h>

#define MAX_STEP 10000

void swap_row(double **mat, int p, int q, int R);
int gaussianElimination(double **mat,int R, double *res);
int forwardElim(double **mat,int R);
void solveLinearSystem(double **mat, int R, double *res);

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    mwSize nStep;   /* Number of step */
    double dt;      /* Integration step */
    mwSize Nstati;       /* Number of state variables*/
    mwSize i,j,k,ii;
    
    double *dx;
    
    double *oldState;
    double *currentState;
    double currentT;
    
    double **eventMatrix;
    
    double *phiOut;
    double *tOut;
    
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
    
	 double zeroVal; /* Zero for newton method */
    double err;
    
    double *DeltaX;
    
    double **J;
        
    int stepsCnt;
	
    double *I0;
    
    /* check for proper number of arguments */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 7 Input required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 2 Output required.");
    }
    
    initVectorField(&vectorField);
    
    
    Nstati =  mxGetScalar(prhs[0]);
    nStep = mxGetScalar(prhs[1]);
    dt = mxGetScalar(prhs[2]);
    oldState = mxGetPr(prhs[3]);
    Vth = mxGetScalar(prhs[4]);
    stopThreshold = mxGetScalar(prhs[5]);
    N =  mxGetScalar(prhs[6]);
    zeroVal =  mxGetScalar(prhs[7]);
    
    I0 = (double *)mxMalloc(Nstati*sizeof(double));
    
	DeltaX = (double *)mxMalloc(Nstati*sizeof(double));
    
    J = (double **)mxMalloc(sizeof(double *)*Nstati);
    for(i=0;i<Nstati;i++)
    {
        J[i] = (double *)mxMalloc(sizeof(double)*(Nstati+1));
        for(j=0;j<Nstati+1;j++)
            J[i][j] = 0;
    }
	
    for(i=0;i<Nstati;i++)
        I0[i] = 0;
    
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
        currentState[j] = oldState[j];
    }
    
    currentT = 0;
    
    for(i=1;i<nStep;i++)
    {
        
        
		
		stepsCnt = 0;
        do
        {
            for(j=0;j<Nstati;j++)
            {
                for(k=0;k<Nstati+1;k++)
                {
                    J[j][k] = 0;
                }
            }
            
            vectorField->getXdot(currentT,currentState,dx,I0);
            vectorField->getJacobian(currentT,currentState,J);
            err = 0;
		
            
            // Preparo le matrici per risolvere il sistema
            for(j=0;j<Nstati;j++)
            {
                for(k=0;k<Nstati;k++)
                {
                    J[j][k] *= -dt;
                    if(j == k)
                        J[j][k] = 1.0+J[j][k];
                }
                J[j][Nstati] = -currentState[j]+dx[j]*dt+oldState[j];
                err += J[j][Nstati]*J[j][Nstati];
                
            }
            
            
                // Risolvo il sistema
                ii = gaussianElimination(J,Nstati,DeltaX);
                
                if(ii != -1)
                    break;
                else
                {
                    for(j=0;j<Nstati;j++)
                    {
                        currentState[j] += DeltaX[j];
                    }
                }
            
            stepsCnt++;
            
            if(stepsCnt > MAX_STEP)
            {
                mexPrintf("Maximum step size reached! Please use smaller time step!\n");
                break;
            }
            
        } while(err > zeroVal);
		
        currentT += dt;      

        vectorField->resetStates(currentState);
        
		 if(stepsCnt > MAX_STEP)
            break;
        
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
    }
    
    
    minEventNumber = nEvent[0];
    
    for(i=1;i<N;i++)
        if(minEventNumber > nEvent[i])
            minEventNumber = nEvent[i];
    
    plhs[0] = mxCreateDoubleMatrix(minEventNumber,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(minEventNumber,N-1,mxREAL);
    
    tOut = mxGetPr(plhs[0]);
    phiOut = mxGetPr(plhs[1]);
    
    if(minEventNumber > 1)
    {
        
        for(j=0;j<minEventNumber;j++)
        {   tOut[j] = eventMatrix[0][j];
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
    
    mxFree(dx);
    mxFree(currentState);
    mxFree(nEvent);
    for(j=0;j<N;j++)
    {
        mxFree(eventMatrix[j]);
    }
    mxFree(eventMatrix);
    mxFree(firstIndex);
    
	 mxFree(DeltaX);
    
    for(i=0;i<Nstati;i++)
        mxFree(J[i]);
    mxFree(J);
	
    delete(vectorField);
    mxFree(I0);
}

int gaussianElimination(double **mat,int R, double *res)
{
    /* reduction into r.e.f. */
    int singular_flag = forwardElim(mat,R);
    
    if (singular_flag != -1)
        return singular_flag;
    
    /* get solution to system and print it using
     * backward substitution */
    solveLinearSystem(mat,R,res);
    return -1;
}

// function for elemntary operation of swapping two rows
void swap_row(double **mat, int p, int q, int R)
{
    int s;
    double temp;
    //printf("Swapped rows %d and %d\n", p, q);
    
    for (s = 0; s <= R; s++)
    {
        temp = mat[p][s];
        mat[p][s] = mat[q][s];
        mat[q][s] = temp;
    }
}


int forwardElim(double **mat,int R)
{
    int i, j, k, p;
    double tmp;
    
    for (i = 0; i < R; i++)
    {
        
        p=i;
        //comparison to select the pivot
        for(j=i+1;j<R;j++)
        {
            if(fabs(mat[j][i]) > fabs(mat[i][i]))
                swap_row(mat,i,j,R);
        }
        
        //cheking for nullity of the pivots
        for(p=i;p<R;p++)
            if(mat[p][i] != 0)
                break;
        
        // Check for singularity
        if(p == R)
        {
            mexPrintf("SINGOLARE!\n");
            mexPrintf("\n");
            for(j=0;j<R;j++)
            {
                for(k=0;k<R+1;k++)
                {
                    mexPrintf("%f\t",mat[j][k]);
                }
                mexPrintf("\n");
            }
            return i;
        }
        else if(p != i)
            swap_row(mat,i,p,R);
        
        
        for(j=i+1;j<R;j++)
        {
            tmp = mat[j][i]/mat[i][i];
            for(k=i+1;k<R+1;k++)
                mat[j][k] -= tmp*mat[i][k];
        }
        
    }

    return -1;
}



void solveLinearSystem(double **mat, int R, double *res) // solve linear system
{
    int i, j;
    /* Start calculating from last equation up to the
     * first */
    for (i = R - 1; i >= 0; i--)
    {
        /* start with the RHS of the equation */
        res[i] = mat[i][R];
        
        /* Initialize j to i+1 since matrix is upper
         * triangular*/
        for (j = i + 1; j < R; j++)
        {
            /* subtract all the lhs values
             * except the coefficient of the variable
             * whose value is being calculated */
            res[i] -= mat[i][j] * res[j];
        }
        
        /* divide the RHS by the coefficient of the
         * unknown being calculated */
        res[i] = res[i] / mat[i][i];
    }
}