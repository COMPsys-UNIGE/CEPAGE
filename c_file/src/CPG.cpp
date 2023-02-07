/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   CPG.cpp
 * Author: picio
 *
 * Created on 16 maggio 2017, 18.41
 */
#include "mex.h"
#include "../inc/CPG.hpp"

CPG::CPG()
{
    this->N = 0;
    this->neuroni = new neuron_model*[0];
    
    
    this->firstState = new int[0];
    
    
    NneuronDelaysIndex = new int[0];
    NinhSynDelaysIndex = new int[0];
    NexcSynDelaysIndex = new int[0];
    
    this->Ninh = 0;
    this->Nexc = 0;
    this->Nel = 0;
    this-> Nstati = 0;
    
    firstState = new int[0];
    
}

CPG::CPG(int N,neuron_model **neuroni, int Ninh, int Nexc, int Nel, synStruct_t **inhSyn, synStruct_t **excSyn , synStruct_t **elSyn, int Ndelay, double networkDelays[])
{
    int i,index,j,k;
    
    double *tmpVect;
    int Ndelaytmp;
    
    double diffMin;
    double iMin;
    
    synapse_model *activation;
    
    this->N = N;
    
    this->neuroni = new neuron_model*[N];
    
    Nstati = 0;
    
    for(i=0;i<N;i++)
    {
        this->neuroni[i] = neuroni[i]->clone();
        Nstati += neuroni[i]->getnx();  
    }
    
    this->Ninh = Ninh;
    this->Nexc = Nexc;
    this->Nel = Nel;
    
    
    this->inhSyn = new synStruct*[Ninh];
    this->excSyn = new synStruct*[Nexc];
    this->elSyn = new synStruct*[Nel];
    
    for(i=0;i<Ninh;i++)
    {
        this->inhSyn[i] = inhSyn[i]->clone();
        Nstati += inhSyn[i]->activation->getnx();;
    }
    for(i=0;i<Nexc;i++)
    {
        this->excSyn[i] = excSyn[i]->clone();
        Nstati += excSyn[i]->activation->getnx();;
    }
    for(i=0;i<Nel;i++)
        this->elSyn[i] = elSyn[i]->clone();
    
    
    this->Ndelay = Ndelay;//((int)sizeof(networkDelays))/sizeof(double);
    
    this->delays = new double[Ndelay];
    
    
    for(i=0;i<Ndelay;i++)
        delays[i] = networkDelays[i];
    
    this->firstState = new int[N+Ninh+Nexc];
    
    
    tmpVect = new double[Ndelay];
    
    this->neuronDelaysIndex = new vector<int>*[N];
    this->NneuronDelaysIndex = new int[N];
    
    index = 0;
    
    
    for(i=0;i<N;i++)
    {
        this->firstState[i] = index;
        index += (this->neuroni[i])->getnx();
        
        
        (this->neuroni[i])->getDelays(tmpVect);
        Ndelaytmp = (this->neuroni[i])->getNdelay();
        
        neuronDelaysIndex[i] = new vector<int>();
        for(j=0;j<Ndelaytmp;j++)
        {
            iMin = 0;
            diffMin = abs(networkDelays[0] - tmpVect[j]);
            for(k=1;k<Ndelay;k++)
            {
                if(diffMin > abs(networkDelays[k] - tmpVect[j]))
                {
                    diffMin = abs(networkDelays[k] - tmpVect[j]);
                    iMin = k;
                }
            }
            
            neuronDelaysIndex[i]->push_back(iMin);
        }
        
        
        NneuronDelaysIndex[i] = Ndelaytmp;
        
        
    }
 
    
    this->inhSynDelaysIndex = new vector<int>*[Ninh];
    this->NinhSynDelaysIndex = new int[Ninh];
    
    for(i=0;i<Ninh;i++)
    {
        this->firstState[N+i] = index;
        
        activation = this->inhSyn[i]->activation;
        
        index += activation->getnx();
        
        activation->getDelays(tmpVect);
        Ndelaytmp = activation->getNdelay();
        inhSynDelaysIndex[i] = new vector<int>();
        for(j=0;j<Ndelaytmp;j++)
        {
            iMin = 0;
            diffMin = abs(networkDelays[0] - tmpVect[j]);
            for(k=1;k<Ndelay;k++)
            {
                if(diffMin > abs(networkDelays[k] - tmpVect[j]))
                {
                    diffMin = abs(networkDelays[k] - tmpVect[j]);
                    iMin = k;
                }
            }
            
            inhSynDelaysIndex[i]->push_back(iMin);
            
        }
        
        NinhSynDelaysIndex[i] = Ndelaytmp;
        
    }
    

    
    this->excSynDelaysIndex = new vector<int>*[Nexc];
    this->NexcSynDelaysIndex = new int[Nexc];
    
    for(i=0;i<Nexc;i++)
    {
        this->firstState[N+Ninh+i] = index;
        
        activation = this->excSyn[i]->activation;
        
        index += activation->getnx();
        
        
        activation->getDelays(tmpVect);
        Ndelaytmp = activation->getNdelay();
        
        excSynDelaysIndex[i] = new vector<int>();
        for(j=0;j<Ndelaytmp;j++)
        {
            iMin = 0;
            diffMin = abs(networkDelays[0] - tmpVect[j]);
            for(k=1;k<Ndelay;k++)
            {
                if(diffMin > abs(networkDelays[k] - tmpVect[j]))
                {
                    diffMin = abs(networkDelays[k] - tmpVect[j]);
                    iMin = k;
                }
            }
            
            excSynDelaysIndex[i]->push_back(iMin);
        }
        
        NexcSynDelaysIndex[i] = Ndelaytmp;
        
    }
    
    
    
    
    delete(tmpVect);
    
}


void CPG::getXdot(double t, double *x, double *xdot,double *Iext,double **Xprec)
{
    int N = this->N;
    int i,j,k;
    neuron_model **neuroni = this->neuroni;
    
    int Ninh = this->Ninh;
    int Nexc = this->Nexc;
    int Nel = this->Nel;
    
    synStruct **inhSyn = this->inhSyn;
    synStruct **excSyn = this->excSyn;
    synStruct **elSyn = this->elSyn;
    
    int *Vindex = this->firstState;
    
    double Isyn,Vi,Vj,g,Esyn;
    
    synapse_model *a;
    
    int maxNdelay = this->Ndelay;
    int *delays_;
    
    
    double *Vjprec;
    double **Xprec_;
    
    double *totI;
    
    int iInh,iExc,iEl;
    
    
    totI = new double[1];
    
    Xprec_ = new double *[maxNdelay];//(double**)malloc(maxNdelay*sizeof(double *));
    
    Vjprec = new double[maxNdelay];//(double *)malloc(maxNdelay*sizeof(double));
    
    iInh = 0;
    iExc = 0;
    iEl = 0;
    
    /* Compute neurons differentials */
    for(i=0;i<N;i++)
    {
        Isyn = 0;
        Vi = x[Vindex[i]];
        
        while(iInh < Ninh)
        {
            if (inhSyn[iInh]->i > i)
            {
                break;
            }
            else
            {                

                j = inhSyn[iInh]->j;
                g = inhSyn[iInh]->g;
                Esyn = inhSyn[iInh]->Esyn;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                
                
                for(k=0;k<NinhSynDelaysIndex[iInh];k++)
                {
                    Xprec_[k] = Xprec[inhSynDelaysIndex[iInh]->at(k)]+Vindex[N+iInh];
                    Vjprec[k] = Xprec[inhSynDelaysIndex[iInh]->at(k)][Vindex[j]];
                }
                
                Isyn += g*(Esyn-Vi)*(a->getActivation(x+Vindex[N+iInh],Vj,Vjprec));
                a->getXdot(t,x+Vindex[N+iInh],xdot+Vindex[N+iInh],Vj,Xprec_,Vjprec);
                
                iInh++;
            }
        }
        
        while(iExc < Nexc)
        {
            if (excSyn[iExc]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                g = excSyn[iExc]->g;
                Esyn = excSyn[iExc]->Esyn;
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
         
                for(k=0;k<NexcSynDelaysIndex[iExc];k++)
                {
                    Xprec_[k] = Xprec[excSynDelaysIndex[iExc]->at(k)]+Vindex[N+Ninh+iExc];
                    Vjprec[k] = Xprec[excSynDelaysIndex[iExc]->at(k)][Vindex[j]];
                }

                a->getXdot(t,x+Vindex[N+Ninh+iExc],xdot+Vindex[N+Ninh+iExc],Vj,Xprec_,Vjprec);
                Isyn += g*(Esyn-Vi)*(a->getActivation(x+Vindex[N+Ninh+iExc],Vj,Vjprec));
                
                iExc++;
            }
        }
        
        while(iEl < Nel)
        {
            if (elSyn[iEl]->i > i)
            {
                break;
            }
            else
            {
                j = elSyn[iEl]->j;
                g = elSyn[iEl]->g;
                
                Vj = x[Vindex[j]];
                Isyn += g*(Vj-Vi);
                
                iEl++;
            }
        }
        
        for(k=0;k<NneuronDelaysIndex[i];k++)
        {
            Xprec_[k] = Xprec[neuronDelaysIndex[i]->at(k)]+Vindex[i];
        }
        
        totI[0] = Isyn + Iext[i];
        
        neuroni[i]->getXdot(t,x+Vindex[i],xdot+Vindex[i],totI,Xprec_);
        
        
    }
    delete[](Xprec_);
    delete(Vjprec);
    delete(totI);
    
}



void CPG::getXdot(double t, double *x, double *xdot,double *Iext)
{
    int N = this->N;
    int i,j;
    neuron_model **neuroni = this->neuroni;
    
    int Ninh = this->Ninh;
    int Nexc = this->Nexc;
    int Nel = this->Nel;
    
    synStruct **inhSyn = this->inhSyn;
    synStruct **excSyn = this->excSyn;
    synStruct **elSyn = this->elSyn;
    
    int *Vindex = this->firstState;
    
    double Isyn,Vi,Vj,g,Esyn;
    
    double *totI;
    
    synapse_model *a;
    
    int iInh,iExc,iEl;
    
    totI = new double[1];

    
    /* Compute neurons differentials*/
    
    iInh = 0;
    iExc = 0;
    iEl = 0;
    
    for(i=0;i<N;i++)
    {
        Isyn = 0;
        Vi = x[Vindex[i]];
        
        while(iInh < Ninh)
        {
            if (inhSyn[iInh]->i > i)
            {
                break;
            }
            else
            {
                j = inhSyn[iInh]->j;
                g = inhSyn[iInh]->g;
                Esyn = inhSyn[iInh]->Esyn;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                
                Isyn += g*(Esyn-Vi)*(a->getActivation(x+Vindex[N+iInh],Vj));
                a->getXdot(t,x+Vindex[N+iInh],xdot+Vindex[N+iInh],Vj);
                
                iInh++;
            }
        }
        
        while(iExc < Nexc)
        {
            if (excSyn[iExc]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                g = excSyn[iExc]->g;
                Esyn = excSyn[iExc]->Esyn;
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
                
                a->getXdot(t,x+Vindex[N+Ninh+iExc],xdot+Vindex[N+Ninh+iExc],Vj);
                Isyn += g*(Esyn-Vi)*(a->getActivation(x+Vindex[N+Ninh+iExc],Vj));
                
                iExc++;
            }
        }
        
        while(iEl < Nel)
        {
            if (elSyn[iEl]->i > i)
            {
                break;
            }
            else
            {
                j = elSyn[iEl]->j;
                g = elSyn[iEl]->g;
                
                Vj = x[Vindex[j]];
                Isyn += g*(Vj-Vi);
                
                iEl++;
            }
        }
        
        totI[0] = Isyn + Iext[i];
        
        neuroni[i]->getXdot(t,x+Vindex[i],xdot+Vindex[i],totI);
    }
    
     delete(totI);
}

bool CPG::getResetConditions(double *x)
{
    int N = this->N;
    int i,j;
    neuron_model **neuroni = this->neuroni;
    int *Vindex = this->firstState;
    double Vj;
    
    synStruct **inhSyn = this->inhSyn;
    synStruct **excSyn = this->excSyn;
    synStruct **elSyn = this->elSyn;
    
    synapse_model *a;
    
    int iInh,iExc;
    
    iInh = 0;
    iExc = 0;
    
    
    
    for(i=0;i<N;i++)
    {
        
        if((neuroni[i])->getResetConditions(x+Vindex[i]))
        {
            return true;
            break;
        }
        
        while(iInh < Ninh)
        {
            if (inhSyn[iInh]->i > i)
            {
                break;
            }
            else
            {
                j = inhSyn[iInh]->j;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                
                if (a->getResetConditions(x+Vindex[N+iInh],Vj))
                    return true;
                
                iInh++;
            }
        }
        
        while(iExc < Nexc)
        {
            if (excSyn[iExc]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
                if (a->getResetConditions(x+Vindex[N+Ninh+iExc],Vj))
                    return true;
                
                iExc++;
            }
        }
    }
    return false;
}


void CPG::resetStates(double *x)
{
    int N = this->N;
    int i,j;
    neuron_model **neuroni = this->neuroni;
    int *Vindex = this->firstState;
    double Vj;
    
    synStruct **inhSyn = this->inhSyn;
    synStruct **excSyn = this->excSyn;
    synStruct **elSyn = this->elSyn;
    
    synapse_model *a;
    
    int iInh,iExc;
    
    iInh = 0;
    iExc = 0;
    
    for(i=0;i<N;i++)
    {
        if((neuroni[i])->getResetConditions(x+Vindex[i]))
        {
            (neuroni[i])->resetStates(x+Vindex[i]);
        }
        
        while(iInh < Ninh)
        {
            if (inhSyn[iInh]->i > i)
            {
                break;
            }
            else
            {
                j = inhSyn[iInh]->j;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                if (a->getResetConditions(x+Vindex[N+iInh],Vj))
                    a->resetStates(x+Vindex[N+iInh],Vj);
                
                iInh++;
            }
        }
        
        while(iExc < Nexc)
        {
            if (excSyn[iExc]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
                if (a->getResetConditions(x+Vindex[N+Ninh+iExc],Vj))
                    a->resetStates(x+Vindex[N+Ninh+iExc],Vj);
                
                iExc++;
            }
        }
        
        
    }
    
    return;
}


void CPG::getFirstIndex(int *firstIndex)
{
    int i;
    for(i=0;i<N+Ninh+Nexc;i++)
        firstIndex[i] = this->firstState[i];
}

CPG::~CPG()
{
    int i;

    delete(this->firstState);
        
    for(i=0;i<N;i++)
    {
        delete(this->neuroni[i]);
        delete(this->neuronDelaysIndex[i]);
    }
    delete[](neuroni);
    delete[](this->neuronDelaysIndex);
  
    for(i=0;i<Ninh;i++)
    {
        delete(this->inhSyn[i]);
        delete(this->inhSynDelaysIndex[i]);
    }
    delete[](inhSyn);
    delete[](this->inhSynDelaysIndex);
    
    for(i=0;i<Nexc;i++)
    {
        delete(this->excSyn[i]);
        delete(this->excSynDelaysIndex[i]);
    }
    delete[](excSyn);
    delete[](this->excSynDelaysIndex);
    
    for(i=0;i<Nel;i++)
        delete(this->elSyn[i]);
    delete[](elSyn);
    
    delete(this->delays);
    
    delete(this->NneuronDelaysIndex);
    delete(this->NinhSynDelaysIndex);
    delete(this->NexcSynDelaysIndex);
    


}

void CPG::getJacobian(double t, double *x, double **J)
{
    int N = this->N;
    int i,j,k;
    neuron_model **neuroni = this->neuroni;
    
    int Ninh = this->Ninh;
    int Nexc = this->Nexc;
    int Nel = this->Nel;
    
    synStruct **inhSyn = this->inhSyn;
    synStruct **excSyn = this->excSyn;
    synStruct **elSyn = this->elSyn;
    
    int *Vindex = this->firstState;
    
    double Isyn,Vi,Vj,g,Esyn;
    
    double *totI;
    
    synapse_model *a;
    
    int iInh,iExc,iEl,ii;
    
    int totStati = this->Nstati;
    
    totI = new double[1];

    double **tmp,**tmp2;
    
    tmp = (double **)mxMalloc(Nstati*sizeof(double *));
    
    tmp2 = (double **)mxMalloc(Nstati*sizeof(double *));
    for(i=0;i<Nstati;i++)
        tmp2[i] = (double *)mxMalloc(Nstati*sizeof(double));
    /* Compute neurons jacobians*/
    
    iInh = 0;
    iExc = 0;
    iEl = 0;
    
    for(i=0;i<N;i++)
    {
        Isyn = 0;
        ii = Vindex[i];
        
        totI[0] = 0;

        /* First neuron jacobian */
        for(j=0;j<this->neuroni[i]->getnx();j++)
            tmp[j] = J[ii+j]+ii;
        neuroni[i]->getJacobian(t,x+Vindex[i],tmp);
        
        Vi = x[Vindex[i]];
        
        
        
        while(iInh < Ninh)
        {
            if (inhSyn[iInh]->i > i)
            {
                break;
            }
            else
            {
                j = inhSyn[iInh]->j;
                g = inhSyn[iInh]->g;
                Esyn = inhSyn[iInh]->Esyn;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                
                a->getDA(x+Vindex[N+iInh],tmp2,Vj);
                J[ii][ii]  -= g*(a->getActivation(x+Vindex[N+iInh],Vj))*(neuroni[i]->dfdi(t,x+Vindex[i]));
                J[ii][Vindex[j]] += g*(Esyn-Vi)*tmp2[0][0]*neuroni[i]->dfdi(t,x+Vindex[i]);
                
                for(j=0;j<a->getnx();j++)
                    J[ii][Vindex[N+iInh]+j] = g*(Esyn-Vi)*tmp2[0][j+1]*neuroni[i]->dfdi(t,x+Vindex[i]);

                iInh++;
            }
        }
        
        while(iExc < Nexc)
        {
            if (excSyn[iExc]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                g = excSyn[iExc]->g;
                Esyn = excSyn[iExc]->Esyn;
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
                
                a->getDA(x+Vindex[N+Ninh+iExc],tmp2,Vj);
                J[ii][ii]  -= g*(a->getActivation(x+Vindex[N+Ninh+iExc],Vj))*(neuroni[i]->dfdi(t,x+Vindex[i]));
                J[ii][Vindex[j]] += g*(Esyn-Vi)*tmp2[0][0]*neuroni[i]->dfdi(t,x+Vindex[i]);
                
                for(j=0;j<a->getnx();j++)
                    J[ii][Vindex[N+Ninh+iExc]+j] = g*(Esyn-Vi)*tmp2[0][j+1]*neuroni[i]->dfdi(t,x+Vindex[i]);
                
                iExc++;
            }
        }
        
        
        
        while(iEl < Nel)
        {
            if (elSyn[iEl]->i > i)
            {
                break;
            }
            else
            {
                j = elSyn[iEl]->j;
                g = elSyn[iEl]->g;
                
                Vj = x[Vindex[j]];
                J[ii][ii]  -= g*(neuroni[i]->dfdi(t,x+Vindex[i]));
                
                iEl++;
            }
        }
        
    }
    
    
    /* Synapses state jac */
    for(iInh=0;iInh < Ninh;iInh++)
        {
                j = inhSyn[iInh]->j;
                g = inhSyn[iInh]->g;
                Esyn = inhSyn[iInh]->Esyn;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                
                ii = Vindex[N+iInh];
                
                a->getJacobian(t,x+ii,tmp2,Vj);

                for(i=0;i<a->getnx();i++)
                {
                    for(k=0;k<a->getnx();k++)
                    {
                        J[ii+i][ii+k] = tmp2[i][k+1];
                    }
                    J[ii+i][Vindex[j]] = tmp2[i][0];
                }
        }
        
        for(iExc=0;iExc< Nexc;iExc++)
        {
                j = excSyn[iExc]->j;
                g = excSyn[iExc]->g;
                Esyn = excSyn[iExc]->Esyn;
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
                ii = Vindex[N+Ninh+iExc];
                a->getJacobian(t,x+ii,tmp2,Vj);
                
                for(i=0;i<a->getnx();i++)
                {
                    for(j=0;j<a->getnx();j++)
                        J[ii+i][ii+j] = tmp2[i][j+1];
                    J[ii+i][Vindex[j]] = tmp2[i][0];
                }
                
            
        }
    
    
    for(i=0;i<Nstati;i++)
        mxFree(tmp2[i]);
    mxFree(tmp2);
    
    mxFree(tmp);
    
     delete(totI);
}

void CPG::getJacobian(double t, double *x, double **J,double **Xprec)
{

}
