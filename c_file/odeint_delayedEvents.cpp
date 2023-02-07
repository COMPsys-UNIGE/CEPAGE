#include "mex.h"
#include "math.h"
#include "vectorField.hpp"

#include <boost/numeric/odeint.hpp>

#include <array>
#include <deque>

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



/***** ADD FOR DELAYED SYSTEM ******/
int Ndelay;
double *delays;

int *delayIndex;
double *x0Del;

double **xDel; /* Ndelay x Nx  */

double longestDelay;

/*************************************/

/* The rhs of x' = f(x) */
void cppVectorialField( const state_type &x , state_type &dxdt , const double  t  )
{
    mwSize i;
    for(i=0;i<Nstati;i++)
        xC[i] = x[i];
    
    vectorField->getXdot(t,xC,dxC,I0,xDel);
    
    for(i=0;i<Nstati;i++)
        dxdt[i] = dxC[i];
}

struct push_back_events
{
    double **eventMatrix;

    std::deque< state_type >& m_states;
    std::deque< double >& m_times;
    
    push_back_events(double **eventMatrix_, std::deque< state_type > &states , std::deque< double > &times ) 
    : eventMatrix(eventMatrix_), m_states( states ) , m_times( times ) { }
    
    void operator()( state_type &x , double t )
    {
        mwSize i,j,indDel,d;
       int tmpIndex;     
        for(i=0;i<Nstati;i++)
            xC[i] = x[i];
        
        

        vectorField->resetStates(xC);

        
        for(i=0;i<Nstati;i++)
            x[i] = xC[i];    
        
        m_states.push_back( x );
        m_times.push_back( t );
        
        
        // Controlla se devo togliere elementi da storia
        while(t-longestDelay > m_times[0])
        {
            m_times.pop_front();
            m_states.pop_front();
            
             for(i=0;i<Ndelay;i++)
                delayIndex[i]--;
        }
        
        for(i=0;i<Ndelay;i++)
            if(delayIndex[i] < 0)
                delayIndex[i] = 0;
        
        
        /* Update delay index*/
        for(i=0;i<Ndelay;i++)
        {
            while(m_times[delayIndex[i]] < t-delays[i])
            {
                delayIndex[i]++;
            }
        }
        
        // Creo Xdel
        for(d=0;d<Ndelay;d++)
        {
            if(t-delays[d] < 0)
            {
                for(i=0;i<Nstati;i++)
                    xDel[d][i] = x0Del[i];
            }
            else
            {
                indDel = delayIndex[d];
                for(i=0;i<Nstati;i++)
                    xDel[d][i] = m_states[indDel][i];
            }
        }
        
        // Controllo eventi per phi
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
    double *Xfinale;
    double period;
    double **eventMatrix;
    
    deque<state_type> x_vec;
    deque<double> times;
    
    /* check for proper number of arguments */
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 9 Input required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 3 Output required.");
    }
    
    Nstati =  mxGetScalar(prhs[0]);
    Tfinale = mxGetScalar(prhs[1]);
    dt = mxGetScalar(prhs[2]);
    x0 = mxGetPr(prhs[3]);
    Vth = mxGetScalar(prhs[4]);
    N =  mxGetScalar(prhs[5]);
    x0Del = mxGetPr(prhs[6]);
    abs_err =  mxGetScalar(prhs[7]);
    rel_err =  mxGetScalar(prhs[8]);
    
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

    
    Ndelay = vectorField->getNdelay();
    delays = (double *)mxMalloc(Ndelay*sizeof(double));
    
    
    vectorField->getDelays(delays);

    longestDelay = -1;
    for(i=0;i<Ndelay;i++)
    {
        if(delays[i] > longestDelay)
            longestDelay = delays[i];
    }

    delayIndex = (int *)mxMalloc(Ndelay*sizeof(int));
    
    xDel = (double **)mxMalloc(Ndelay*sizeof(double *));
    
    for(i=0;i<Ndelay;i++)
    {
        xDel[i] = (double *)mxMalloc(Nstati*sizeof(double));
        for(j=0;j<Nstati;j++)
        {
            xDel[i][j] = x0Del[i*Nstati+j];
        }
    }

    for(i=0;i<Ndelay;i++)
    {
        delayIndex[i] = 0;
    }

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
            push_back_events(eventMatrix,  x_vec , times));
    minEventNumber = nEvent[0];
    
    for(i=1;i<N;i++)
        if(minEventNumber > nEvent[i])
            minEventNumber = nEvent[i];
    

    
    plhs[0] = mxCreateDoubleMatrix(minEventNumber,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(minEventNumber,N-1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(Nstati,1,mxREAL);
    
    tOut = mxGetPr(plhs[0]);
    phiOut = mxGetPr(plhs[1]);
    Xfinale = mxGetPr(plhs[2]);
    
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

    oldState = x_vec.back();
     for(j=0;j<Nstati;j++)
         Xfinale[j] = oldState[j];
    
    for(j=0;j<N;j++)
    {
        mxFree(eventMatrix[j]);
    }
    mxFree(eventMatrix);
    
    mxFree(nEvent);
    
    mxFree(xC);
    mxFree(dxC);
    
    mxFree(firstIndex);
    
    for(i=0;i<Ndelay;i++)
        mxFree(xDel[i]);
    mxFree(xDel);
    mxFree(delays);
    mxFree(delayIndex);
    
    delete(vectorField);
    mxFree(I0);
}
