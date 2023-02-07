#include "mex.h"
#include "math.h"
#include "vectorField.hpp"
//#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/odeint.hpp>
#include <vector>

typedef std::vector< double > state_type;

mwSize Nstati;       /* Number of state variables*/
double *xC;
double *dxC;

double Vth;

int *firstIndex;

mwSize *nEvent;

mwSize minEventNumber;

mwSize N;

state_type oldState;

double oldT;

dynSys *vectorField;

double *I0;

/* The rhs of x' = f(x) */
void cppVectorialField( const state_type &x , state_type &dxdt , const double  t  )
{
    mwSize i;
    for(i=0;i<Nstati;i++)
        xC[i] = x[i];
    
    vectorField->getXdot(t,xC,dxC,I0);
    
    for(i=0;i<Nstati;i++)
        dxdt[i] = dxC[i];
}

struct push_back_events
{
    double **eventMatrix;
    int tmpIndex;
    push_back_events(double **eventMatrix_) : eventMatrix(eventMatrix_) { }
    
    void operator()( state_type &x , double t )
    {
        mwSize i,j;
        for(i=0;i<Nstati;i++)
            xC[i] = x[i];
        
        

        vectorField->resetStates(xC);

        
        for(i=0;i<Nstati;i++)
            x[i] = xC[i];    
        
       
        
        for(j=0;j<N;j++)
        {

            tmpIndex = firstIndex[j];
           
            if((oldState[tmpIndex] < Vth) && (x[tmpIndex] > Vth))
            {
                nEvent[j]++;
                eventMatrix[j] = (double *)mxRealloc(eventMatrix[j],nEvent[j]*sizeof(double));
                eventMatrix[j][nEvent[j]-1] = (Vth-oldState[tmpIndex])/(x[tmpIndex]-oldState[tmpIndex])*(t-oldT)+oldT;
            }
        }
        
        copy(x.begin(), x.end(), oldState.begin());
        oldT = t;
    }
};

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    using namespace std;
    using namespace boost::numeric::odeint;
    
    mwSize nStep;   /* Number of step */
    double dt;      /* Integration step */
    mwSize i,j;
    double *ptr;
    double Tfinale;
    
    
    double *x0;     /* Inital condition */
    double *xOutMatrix; /* Out vector */
    double *dx;
    
    /* Integrator tolerances */
    double abs_err = 1.0e-10;
    double rel_err = 1.0e-6;
    const double a_x = 1.0;
    const double a_dxdt = 1.0;
    
    double *phiOut;
    double *tOut;
    double period;
    double **eventMatrix;
    
    /* check for proper number of arguments */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 8 Input required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 2 Output required.");
    }
    
    Nstati =  mxGetScalar(prhs[0]);
    Tfinale = mxGetScalar(prhs[1]);
    dt = mxGetScalar(prhs[2]);
    x0 = mxGetPr(prhs[3]);
    Vth = mxGetScalar(prhs[4]);
    N =  mxGetScalar(prhs[5]);
    abs_err =  mxGetScalar(prhs[6]);
    rel_err =  mxGetScalar(prhs[7]);
    
    
    state_type xinit(Nstati);
    state_type x(Nstati);

    I0 = (double *)(mxMalloc(Nstati*sizeof(double)));
    for(i=0;i<Nstati;i++)
        I0[i] = 0;
    
    xC = (double *)(mxMalloc(Nstati*sizeof(double)));
    dxC = (double *)(mxMalloc(Nstati*sizeof(double)));
    
    firstIndex = (int *)mxMalloc((N+2*N*N)*sizeof(int));
    

    oldT = 0;
    nEvent = (mwSize *)mxMalloc(N*sizeof(mwSize));
    eventMatrix = (double **)mxMalloc(N*sizeof(double *));

    initVectorField(&vectorField);

    vectorField->getFirstIndex(firstIndex);

    
    for(i=0;i<Nstati;i++)
    {
        xinit[i] = x0[i];
        oldState.push_back(x0[i]);
    }
    
    
    for(j=0;j<N;j++)
    {
        eventMatrix[j] = (double *)mxMalloc(0);
        nEvent[j] = 0;
    }
    
    
    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
    
    
    typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
    controlled_stepper_type controlled_stepper(default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );
    
    size_t steps = integrate_adaptive( controlled_stepper , cppVectorialField , xinit , 0.0 , Tfinale, dt ,
            push_back_events(eventMatrix));
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
    
    for(j=0;j<N;j++)
    {
        mxFree(eventMatrix[j]);
    }
    mxFree(eventMatrix);
    mxFree(xC);
    mxFree(dxC);
    
    mxFree(firstIndex);
    
    delete(vectorField);
    mxFree(I0);
}
