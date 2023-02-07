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

state_type oldState;

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

struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;
    
    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }
    
    void operator()( state_type &x , double t )
    {
        mwSize i;
        for(i=0;i<Nstati;i++)
            xC[i] = x[i];
        
        
        vectorField->resetStates(xC);
        
        
        for(i=0;i<Nstati;i++)
            x[i] = xC[i];
        
        
        
        m_states.push_back( x );
        m_times.push_back( t );
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
    
    
    vector<state_type> x_vec;
    vector<double> times;
    
    /* check for proper number of arguments */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 6 Input required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 1 Output required.");
    }
    
    initVectorField(&vectorField);
    
    Nstati =  mxGetScalar(prhs[0]);
    Tfinale = mxGetScalar(prhs[1]);
    dt = mxGetScalar(prhs[2]);
    x0 = mxGetPr(prhs[3]);
    abs_err = mxGetScalar(prhs[4]);
    rel_err = mxGetScalar(prhs[5]);
    
    state_type xinit(Nstati);
    state_type x(Nstati);
    
    I0 = (double *)mxMalloc(Nstati*sizeof(double));
    
    for(i=0;i<Nstati;i++)
        I0[i] = 0;
    
    xC = (double *)(mxMalloc(Nstati*sizeof(double)));
    dxC = (double *)(mxMalloc(Nstati*sizeof(double)));
    
    
    for(i=0;i<Nstati;i++)
        xinit[i] = x0[i];
    
    
    
    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
    typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
    
    controlled_stepper_type controlled_stepper(default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );
    
    size_t steps = integrate_adaptive( controlled_stepper , cppVectorialField , xinit , 0.0 , Tfinale, dt ,
            push_back_state_and_time( x_vec , times ));
    
    plhs[0] = mxCreateDoubleMatrix(1,(int)(times.size()-1), mxREAL);
    times.pop_back();
    std::copy(times.begin(), times.end(), mxGetPr(plhs[0]));
    
    plhs[1] = mxCreateDoubleMatrix((int)(steps),Nstati, mxREAL);
    
    ptr = mxGetPr(plhs[1]);
    
    mxFree(xC);
    mxFree(dxC);
    
    //printf("%f %f %f\n",x_vec[i][0],x_vec[i][1],x_vec[i][2]);
    for(j=0;j<Nstati;j++)
        for(i=0;i<steps;i++)
            ptr[i+j*steps] = x_vec[i][j];
    
    x_vec.clear();
    
    
    delete(vectorField);
    mxFree(I0);
    
}
