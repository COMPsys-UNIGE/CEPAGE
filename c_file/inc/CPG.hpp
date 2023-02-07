/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CPG.h
 * Author: picio
 *
 * Created on 16 maggio 2017, 18.41
 */

#ifndef CPG_H
#define CPG_H


#include <vector>


#include "dynSys.hpp"
#include "neuron_model.hpp"
#include "synapse_model.hpp"

#include <cmath>

using namespace std;

class CPG : public dynSys {
    
public: struct synStruct
{
    int i;
    int j;
    double g;
    double Esyn;
    synapse_model *activation;
    
    synStruct(const synStruct &s)
    {
        this->i = s.i;
        this->j = s.j;
        this->g = s.g;
        this->Esyn = s.Esyn;
        this->activation = (s.activation)->clone();
    };
    
    synStruct(int i, int j, double g, double Esyn, synapse_model *activation)
    {
        this->i = i;
        this->j = j;
        this->g = g;
        this->Esyn = Esyn;
        this->activation = activation->clone();
    };
    
    synStruct *clone() const
    {
        return new synStruct(*this);
    }
    
    ~synStruct()
    {
        delete(activation);
    }
};

typedef synStruct synStruct_t; 
    
    
private:
    int N;
    neuron_model **neuroni;
         
    vector<int> **neuronDelaysIndex;
    vector<int> **inhSynDelaysIndex;
    vector<int> **excSynDelaysIndex;
 
    int *NneuronDelaysIndex;
    int *NinhSynDelaysIndex;
    int *NexcSynDelaysIndex;
    
    int Ninh;
    int Nexc;
    int Nel;
    
    synStruct_t **inhSyn;
    synStruct_t **excSyn;
    synStruct_t **elSyn;
    
    
    int *firstState;
    
    int Nstati;
    
public:
    CPG();
	CPG(int N,neuron_model **neuroni, int Ninh, int Nexc, int Nel, synStruct_t **inhSyn, synStruct_t **excSyn , synStruct_t **elSyn, int Ndelay, double networkDelays[]);
    
    /* Xprec is a Ndelays x Nstates matrix */
    void getXdot(double t, double *x, double *xdot,double *Iext,double **Xprec);
    void getXdot(double t, double *x, double *xdot,double *Iext);
    
    void getJacobian(double t, double *x, double **J);
    void getJacobian(double t, double *x, double **J,double **Xprec);
    
    bool getResetConditions(double *x);
    
    void resetStates(double *x);
    
    void getFirstIndex(int *firstIndex);
    
    virtual ~CPG();


};

#endif /* CPG_H */

