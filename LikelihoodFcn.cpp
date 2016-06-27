//
//  main.cpp
//  Lkhd_For_Ingyin - XCode
//
//  Created by Joshua Rosaler on 17/08/2015.
//  Copyright (c) 2015 Joshua Rosaler. All rights reserved.
//

// Lkhd_for_Ingyin.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
using namespace std;

/*Overview:
 The objective is to write a class with a method that accepts measurements of a track and outputs the value in ADC
 counts corresponding to the MPV of the Landau distribution into which the 6 (or however many) track hits most probably
 fall. This is done by maximizing the Likelihood function with respect to the free parameter u1_ag, which is the MPV of
 the Landau (explained in greater detail below).
 
 We fit our approximation of the Landau, a bifurcated Gaussian with an exponential tail possessing 4 free
 parameters (lower, upper widths and MPV for the bifurcated gaussian part, and the place where tail begins) to each of 6
 single hit histograms taken from our detector. The last histogram,
 the aggregate, is the distribution over all single hits for several thousand tracks, and corresponds to the true Landau.
 Each of the 5 other histograms is the distribution over all single hits reading out over a particular number of strips.
 That is, we have a distribution for single hits reading out over 1 strip, another for those reading out over 2, etc.
 After fitting each of these distributions, and collecting the corresponding parameters for each, we expressed,
 by "brute force", each of the parameters as a linear function of u1_ag (the MPV of the agg distribution).
 This allowed us to reduce the number of parameters in each of the five single hit distributions, which we determined
 from the fits, to one. Since the likelihood function for a single track is simply the product of the Landau distributions
 for each hit, it goes from being a function of 4 variables, or parameters, to being a function of 1, which makes the
 process of maximization much easier.
 
 
 
 Design:
 (explanation of most important functions)
 Find_Max_Lkhd:
 The entire class revolves around the function Find_Max_Lhkd. This function
 accepts 3 arrays of data from a single track, called ADC, ns, and side, and loops through a range of values for u1_ag
 until it finds one that maximizes the Likelihood function, or that,equivalently, minimizes -2ln(L).
 ADC lists the depositions for the six hits, ns lists strip numbers, and side lists whether each hit is phi, z or
 the avg of those two.
 
 fLikelihood:
 This is effectively the Likelihood function, or rather -2 times the natural log of it.
 It accepts the three arrays ADC,ns, and side, as well as a value of the parameter u1_ag, and returns -2ln(L).
 It calculates this by multiplying the single hit distributions for all the hits and taking -2ln of this product.
 It actually takes the log first, sums these and multiplies by -2, which is equivalent.
 
 fLan:
 This function incapsulates all 6 single hit distributions.
 It inputs values of the deposition (any positive real number), strip number (1 to 5; 0 is used to refer to the agg
 distribution) and side (0 for phi, 1 for z, and 2 for the avg) for a single hit, as well as a value for u1_ag, to which
 the parameters for all 6 distributions are linearly related, and returns a value for the distribution into which that
 particular hit falls. For example, if we have a phi-side hit with an ADC readout of 40 over 3 strips, the function will
 return the value of the three-strip distribution at 40. If we have a hit with an ADC count also of 40, but reading out
 over 5 strips, the return value of the function will be different because it is instead evaluating the 5-strip
 distribution at 40. The side (phi, z or avg) matters because if the hit is z-side, then we cannot use individual strip
 distributions and instead must evaluate the probability of the hit within the aggregate distribution. If the hit is
 phi-side or avg, though, it is ok to go to the individual strip distributions.
 
 Landau and Gauss:
 These simply evaluate the Landau (or rather, our approximation to it) and the bifurcated Gaussian functions,
 given a set of parameter values and a real number x at which to evaluate the distribution.
 
 Set_Norm/Norm:
 In order for the Likelihood function to work, we need to normalize each of the single hit distributions.
 But these normalizations depend on the parameters for the distributions, which in turn depend on u1_ag,
 the value of which changes in the loop in Find_Max_Lkhd. In order to avoid integrating every
 distribution every time we use input a new value
 of u1_ag into fLikelihood, we tabulate the values of the normalization as a function of u1_ag for each of the
 six single hit distributions at 1000 values of u1_ag ranging from 0 to 100 ADC counts. So, we end up with a 6 by 1000
 array (normalizations[6][1000]) that holds all the normalizations for the.six distributions Then, when we want to
 find the value of the normalization for a single hit distribution for any arbitrary value of u1_ag (between 1 and 100),
 we can interpolate between the tabulated values, which is what the function Norm does. More specifically, Norm
 accepts a value for ns from 1 to 6, which tells it which distribution to get the norm of, and a value of u1_ag, which
 indicates which two values in the tabulation to interpolate between
 (note: ns=6 refers not to a strip number but the aggregate).
 
 
 */





class Likelihood
{
public:
    Likelihood();
    ~Likelihood();
    
    //these functions retrieve and assign the elements of the arrays declared below.
    double Get_sigma1(int n);
    void Set_sigma1(double * ps1array);
    double Get_sigma2(int n);
    void Set_sigma2(double * ps2array);
    double Get_mu1(int n);
    void Set_mu1(double * pu1array);
    double Get_mu2(int n);
    void Set_mu2(double * pu2array);
    
    double Get_norms(int n, int m);//retrieves array elements
    void Set_norms_tab();//tabulate values for normalization
    double Norm(int ns, double u1);//use interpolation to retrieve norm values over continuous range
    
    double Gauss(double sigma_1, double sigma_2, double u, double x);//Gaussian function, takes 3 parameters (
    // sigma1=lower width, sigma2=upper width, u=MPV) and argument (x)
    double Landau(double sigma_1, double sigma_2, double u1, double u2, double x);//approx. to Landau function.
    //is a bifurcated Gaussian with an exponential tail. takes same 3 parameters as Gauss, plus u2, where
    //the tail starts.
    double fLan(double u1, int ns, double adc, int side); //this takes data from a single hit
    //(side = phi,z, or avg; ns = number of strips on readout; adc = deposition in adc counts)
    //and a specified value for u1_ag (u1), on which all 6 ns distributions depend, and returns the value of the
    //single hit distribution corresponding to ns evaluated at adc. the values of side are 0 for phi, 1 for z, 2 for avg.
    double fLikelihood(double u1,double * pADC, int * pns, int * pside);//this is the Likelihood function.
    //it accepts 3 arrays for a single track  and multiplies the values of fLan corresponding to all the hits.
    //it depends on u1_ag (u1) since fLan does.
    double Find_MaxLkhd(double * pADC, int * pns, int * pside);//this finds the value of u1_ag that maximizes the
    //Likelihood function.
    
private:
    double sigma1[6];//holds coefficients for sigma1 expressed as a linear function of u1_ag.
    double sigma2[6];//" " " sigma2 " " " " " ".
    double mu1[6];//" " " mu1 " " " " " ".
    double mu2[6];//" " " mu2 " " " " " ".
    double normalizations[6][1000];//these are the tabulated normalizations for each of the 6 single hit distributions.
};




Likelihood::Likelihood() //set values for coefficient arrays (coefficients determined experimentally)
{
    //These arrays give parameters of single hit dist (for agg and for all ns) as linear func's
    //of u1_ag.
    double s1array[6] = {0.230004,0.188371,0.190532,0.213832,0.355758,0.247916};
    double * ps1array = s1array;
    Set_sigma1(ps1array);
    
    double s2array[6] = {0.178224,0.232363,0.285314,0.631839,1.558563,0.328999};
    double * ps2array = s2array;
    Set_sigma2(ps2array);
    
    double u1array[6] = {0.829174,0.962608,1.143752,1.298738,1.726436,1.000000};
    double * pu1array = u1array;
    Set_mu1(pu1array);
    
    double u2array[6] = {1.047041,1.220573,1.390308,1.935205,3.282429,1.346366};
    double * pu2array = u2array;
    Set_mu2(pu2array);
}



Likelihood::~Likelihood()
{
}




double Likelihood::Get_sigma1(int n)
{
    if (n<0||n>=6)
        return 999999;
    return sigma1[n];
}

void Likelihood::Set_sigma1(double * parray)
{
    for (int i=0; i<6; i++)
    {
        sigma1[i]=parray[i];
    }
}



double Likelihood::Get_sigma2(int n)
{
    if (n<0||n>=6)
        return 999999;
    return sigma2[n];
}

void Likelihood::Set_sigma2(double * parray)
{
    for (int i=0; i<6; i++)
    {
        sigma2[i]=parray[i];
    }
}



double Likelihood::Get_mu1(int n)
{
    if (n<0||n>=6)
        return 999999;
    return mu1[n];
}

void Likelihood::Set_mu1(double * parray)
{
    for (int i=0; i<6; i++)
    {
        mu1[i]=parray[i];
    }
}



double Likelihood::Get_mu2(int n)
{
    if (n<0||n>=6)
        return 999999;
    return mu2[n];
}

void Likelihood::Set_mu2(double * parray)
{
    for (int i=0; i<6; i++)
    {
        mu2[i]=parray[i];
    }
}





double Likelihood::Get_norms(int n,int m)//retrieves specified element of tabulated normalizations.
{
    if (n<0||n>5||m<0||m>999)//if input out of range, return some recognizeably incorrect number.
        return 999999;
    return normalizations[n][m];
}

void Likelihood::Set_norms_tab()//sets normalizations for single hit distribution (fLan), for each strip.
{
    for(int ns=1; ns<7; ns++)//ns=6 corresponds to aggregate ditribution.
    {
        normalizations[ns-1][0]=0;//set first element,which we will not refer to later on, as 0.
        //split integral into two parts, bifurcated Gaussian and tail,
        for(int i=1; i<1000; i++)
        {
            cout<<"In Set_norms_tab; ns: "<<ns<<", i: "<<i<<endl;
            double u1=.1*i;//loop through values of u1_ag
            double dy=0.001;//set interval for summing
            double y=0.0;//start off at 0.
            double sumG=0.0;//start sum off at 0.
            
            while (y>=0 && y<=mu2[ns-1]*u1)//integrate Gaussian part
            {
                y=y+dy;//increment variable
                sumG=sumG+dy*Gauss(sigma1[ns-1]*u1, sigma2[ns-1]*u1, mu1[ns-1]*u1, y);//perform correpsonding increment
                //on sum.
            }
            
            //Can calculate integral of tail exactly, w/o loop.
            double k=exp(-(mu2[ns-1]*u1 - mu1[ns-1]*u1)*(mu2[ns-1]*u1 - mu1[ns-1]*u1)/(2.0*sigma2[ns-1]*u1*sigma2[ns-1]*u1));
            //k was calculated as a function of the parameters u1,u2,s1,s2, by setting the functions and their derivatives
            //equal at u2.
            double t=sigma2[ns-1]*u1*sigma2[ns-1]*u1/(mu2[ns-1]*u1-mu1[ns-1]*u1);
            //same for t.
            double norm = sumG+k*t;//kt is integral of exponential tail, total norm is sum of this and gauss integral.
            normalizations[ns-1][i] = norm;//set array element equal to norm.
        }
    }
}

double Likelihood::Norm(int ns, double u1)//Interpolates between tabulated norms. Before calling, must call Set_norms_tab.
//accepts ns to specify which distribution to get the norm of, and u1_ag, which fixes param values.
{
    //to perform linear interpolation, we will want to know the slope dy/dx between successive values of the norm w.r.t.
    //the variable u1, which represents u1_ag. Since our interpolation increments u1_ag by .1, our dx=.1, as well.
    //given u1_ag, we want to find the two entries in the array normalizations[6][1000] corresponding to the values
    //of u1_ag immediately below and above it.
    if(ns<0||ns>5||u1<.1||u1>100) //Check that inputs are within range.
    {
        cout<<"ns or u1 out of range"<<endl;
        return 0;
    }
    double dx = .1;//set dx=.1
    int n = int(10*u1);//in order to find integer place of u1 value within norm array, multiply value of u1 by 10 and
    //chop off decimal part using int funciton.
    double y1 = normalizations[ns-1][n];//get lower norm
    double y2 = normalizations[ns-1][n+1];//get upper norm
    double slope = (y2 - y1)/dx;//divide difference by increment in u1_ag between two norms.
    double ans = y1 + slope*(u1-.1*n);//to find norm at u1, add slope*(difference between u1 and closest value of u1_ag
    //that is < u1 and for which norm is tabulated = u1-.1*n) to y1.
    return ans;
}











double Likelihood::Gauss(double sigma_1, double sigma_2, double u, double x)//this is just a bifurcated Gaussian
{
    double ans;
    if (sigma_1==0||sigma_2==0)//just in case one of the widths is 0.
    {
        ans=0;
    }
    else
    {
        if(x<=u)
        {
            ans=exp(-(x-u)*(x-u)/(2*sigma_1*sigma_1));//lower half of bif gaus
        }
        else
        {
            ans=exp(-(x-u)*(x-u)/(2*sigma_2*sigma_2));//upper half of bif gaus
        }
    }
    return ans;
}


double Likelihood::Landau(double sigma_1, double sigma_2, double u1, double u2, double x)
{
    double ans;
    double k=exp(-(u2-u1)*(u2-u1)/(2.0*sigma_2*sigma_2));//k is coefficient of exp tail, calculated in terms of the params
    //by setting tail and bif gaus equal at u2, and also by setting their first derivatives equal there.
    double t=sigma_2*sigma_2/(u2-u1);//term in exp tail determined in same way as k.
    if(sigma_1==0||sigma_2==0||u2<=u1)//just in case one of param's is 0.
    {
        cout<<"In Landau (sigma_1 or sigma_2=0 or u1<=u2)"<<endl;
        if (sigma_1==0)
            cout<<"In Landau (sigma_1==0)"<<endl;
        if (sigma_2==0)
            cout<<"In Landau (sigma_2==0)"<<endl;
        if (u2<=u1)
            cout<<"In Landau (u1<=u2); u1: "<<u1<<", u2: "<<u2<<endl;
        ans=0;
    }
    else//Landau function
    {
        if (x<=u2)
        {
            ans=Gauss(sigma_1,sigma_2,u1,x);//bif gaus part
        }
        else
        {
            ans=k*exp(-(x-u2)/t);//exp tail part
        }
    }
    return ans;
}



double Likelihood::fLan(double u1, int ns, double x, int side)//given
//u1_ag, ns (number of strips over which a hit reads out), and side (phi, z, or avg), this function evaluates the
//corresponding single hit distribution at x, which represents the charge deposition readout in ADC counts for a single
//hit.
{
    if (u1==0)//can't have u1_ag = 0, because then all param's for all dist's are 0.
    {
        cout<<"In fLan, u1==0"<<endl;
        return 1; //so it does not effect likelihood when log is taken.
    }
    if(side==0||side==2)  //if side is phi or avg, assuming ns for an avg hit is taken from phi, we may refer to the
        //strip number to determine which distribution to use.
    {
        if(ns<1||ns>5) //make sure ns within range
        {
            cout<<"Strip # out of range"<<endl;
            return 1;
        }
        
        double ans=Landau(sigma1[ns-1]*u1, sigma2[ns-1]*u1, mu1[ns-1]*u1, mu2[ns-1]*u1, x);//use parameter coeff arrays
        //sigma1[6], sigma2[6], mu1[6], mu2[6] to get appropriate param's for Landau.
        ans=ans/(Norm(ns,u1));//normalize appropriately
        return ans;
    }
    else//if side is z, we cannot refer to strip number, so we use agg. here ns=6 corresponds to agg.
    {
        ns=6;
        double ans=Landau(sigma1[ns-1]*u1, sigma2[ns-1]*u1, mu1[ns-1]*u1, mu2[ns-1]*u1, x);//use coeff arrays
        ans=ans/(Norm(ns,u1));//normalize appropriately
        return ans;
    }
}


double Likelihood::fLikelihood(double u1, double * pADC, int * pns, int * pside)//takes hit arrays (ADC, ns, side) and
//evaluates -2 ln of the likelihood function for specified value of u1_ag parameter.
{
    double L=0;
    if (u1==0)//in case u1=0
    {
        cout<<"In fLikelihood, u1==0"<<endl;
        return 0;
    }
    for (int i=0;i<8;i++)//loop through hits on a given track
    {
        if (pns[i]==0||pADC[i]==0)//skip hits that have 0 for ns or ADC
        {;}
        else
        {
            L=L+log(fLan(u1,pns[i],pADC[i],pside[i]));//summing logs of individual distribution points is equivalent
            //to taking ln of product of distribution values.
        }
    }
    double ans=-2*L;//this is the function we ultimately would like to minimize.
    return ans;
}


double Likelihood::Find_MaxLkhd(double * pADC, int * pns, int * pside)//Accepts track arrays and finds u1_ag for which
{																	  //-2*ln(L) is minimum.
    for (int i=2;i<1000;i++)//start at i=2 to avoid u1=0. loop from u1_ag = .2 to u1_ag = 100.
    {
        cout<<"In Find_MaxLkhd, i: "<<i<<endl;
        if (fLikelihood(.1*i,pADC,pns,pside)<fLikelihood(.1*(i-1),pADC,pns,pside) && fLikelihood(.1*i,pADC,pns,pside)<fLikelihood(.1*(i+1),pADC,pns,pside))
            //if -2*ln(L) at a particular value of u1_ag is less than -2*ln(L) for the values of u1_ag immediately above and
            //below it, then that value maximizes the likelihood function.
        {	
            double ans = .1*i;
            return ans;
        }
    }
    return 0;//if nothing turns up,which shouldnt happen, return 0
}











int main(int argc, char* argv[])
{
    Likelihood function1;//declare an instance of the class
    function1.Set_norms_tab();//Set tabulation for norms for Landau. Will take a long time, but only has to be done
    //once at beginning. May do many tracks without calling this again. 
    
    //Set up input arrays for Find_MaxLkhd
    double ADC[8] = {26.4,28.2,42.0,25.8,48.7,36.3,0,0};//these are phi-side hits and strip #'s generated randomly by 
    //other code.
    double * pADC = ADC;
    int ns[8] = {3,2,3,1,2,4,0,0};
    int * pns = ns;
    int side[8] = {0,0,0,0,0,0,0,0};
    int * pside = side;
    
    cout<<"Max Likelihood at u1_ag= "<<function1.Find_MaxLkhd(pADC,pns,pside)<<endl;
    
    return 0; 
}
