/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dynSys.h
 * Author: teo
 *
 * Created on July 12, 2017, 11:27 AM
 */

#ifndef DYNSYS_H
#define DYNSYS_H


class dynSys 
{
    
protected:
        int Ndelay = 0;
        double *delays;
    
public:
   
    
    virtual void getXdot(double t, double *x, double *xdot,double *Iext) = 0;
    virtual void getXdot(double t, double *x, double *xdot,double *Iext,double **Xprec);
    
    virtual void getJacobian(double t, double *x, double **J) = 0;
    virtual void getJacobian(double t, double *x, double **J,double **Xprec) = 0;
    virtual bool getResetConditions(double *x) = 0;
    virtual bool getResetConditions(double *x,double **Xprec);
    
    virtual void resetStates(double *x) = 0;
    virtual void resetStates(double *x,double **Xprec);
    
    virtual void getFirstIndex(int *firstIndex);
    
    int getNdelay();
    void getDelays(double *del);
    
    virtual ~dynSys()
    {}

};

#endif /* DYNSYS_H */

