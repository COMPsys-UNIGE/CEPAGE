/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dynSys.cpp
 * Author: teo
 * 
 * Created on July 12, 2017, 11:27 AM
 */

#include "../inc/dynSys.hpp"

void dynSys::getFirstIndex(int *firstIndex)
{
    firstIndex = new int[1];
    firstIndex[0] = 0;
}

int dynSys::getNdelay()
{
    return Ndelay;
}

void dynSys::getDelays(double *del)
{
    int i;
    for(i=0;i<Ndelay;i++)
        del[i] = this->delays[i];
}


void dynSys::getXdot(double t, double *x, double *xdot,double *Iext,double **Xprec)
{
    this->getXdot(t,x,xdot,Iext);
}


bool dynSys::getResetConditions(double *x,double **Xprec)
{
    return this->getResetConditions(x);
}
    
void dynSys::resetStates(double *x,double **Xprec)
{
    this->resetStates(x);
}
