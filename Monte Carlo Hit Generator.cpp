// hit_generator.cpp : Defines the entry point for the console application.
//




// MaxLikelihoodFunc.cpp : Defines the entry point for the console application.
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
using std::endl;
using namespace std;



class Likelihood
{
public:
    Likelihood();
    //Likelihood(double * psigma1, double * psigma2, double * pmu1, double * pmu2, double * pphi_ns, double  * pphi_ADC, double * pz_ns, double * pz_ADC);
    ~Likelihood();
    
    double Get_sigma1(int n);
    void Set_sigma1(double * parray1);
    double Get_sigma2(int n);
    void Set_sigma2(double * parray2);
    double Get_mu1(int n);
    void Set_mu1(double * parray3);
    double Get_mu2(int n);
    void Set_mu2(double * parray4);
    int Get_phi_ns(int n);//n is place in array
    void Set_phi_ns(int * parray5);
    double Get_phi_ADC(int n);
    void Set_phi_ADC(double * parray6);
    int Get_z_ns(int n);
    void Set_z_ns(int * parray7);
    double Get_z_ADC(int n);
    void Set_z_ADC(double * parray8);
    
    double Gauss(double sigma_1, double sigma_2, double u, double x);
    double Landau(double sigma_1, double sigma_2, double u1, double u2, double x);
    double Get_norms(int n, int m);//retrieves array elements
    void Set_norms_tab();//tabulate values for normalization
    double Norm(int ns, double u1);//use interpolation to retrieve norm values over continuous range
    double fLan(double u1, int ns, double adc, bool is_z);
    double fLikelihood(double u1);
     double Find_MaxLkhd(double * pADC, int * pns, int * pside);
    double Cum_Lan(double u1, int ns, double x, bool is_z);
    void pick_ns();
    void hit_generator(double u1);
    
private:
    double sigma1[6];//;={0.230004,0.188371,0.190532,0.213832,0.355758,0.247916};
    double sigma2[6];//={0.178224,0.232363,0.285314,0.631839,1.558563,0.328999};
    double mu1[6];//={0.829174,0.962608,1.143752,1.298738,1.726436,1.000000};
    double mu2[6];//={1.047041,1.220573,1.390308,1.935205,3.282429,1.346366};
    //these come from curve fitting.
    
    int phi_ns[8];//={1,3,2,4,3,3,5,4};
    double phi_ADC[8];//={30,37,34,40,37,37,43,40};
    int z_ns[8];//={3,3,4,1,5,3,2,3};
    double z_ADC[8];//={37,37,40,30,43,37,34,37};
    double normalizations[6][1000];//these are the tabulated normalizations for each of the 6 single hit distributions.
    //these are made up, and ultimately should come from a measurements of a particle track.
};




Likelihood::Likelihood() //set default values for arrays
{
    /*cout<<"CONSTRUCTOR..."<<endl;
     sigma1[6]={0.230004,0.188371,0.190532,0.213832,0.355758,0.247916};
     for(int j=0;j<6;j++)
     {
     cout<<sigma1[j]<<endl;
     }
     sigma2[6]={0.178224,0.232363,0.285314,0.631839,1.558563,0.328999};
     mu1[6]={0.829174,0.962608,1.143752,1.298738,1.726436,1.000000};
     mu2[6]={1.047041,1.220573,1.390308,1.935205,3.282429,1.346366};
     
     phi_ns[8]={1,3,2,4,3,3,5,4};
     for(int i=0;i<8;i++)
     {
     cout<<phi_ns[i]<<endl;
     }
     phi_ADC[8]={30.0,37.0,34.0,40.0,37.0,37.0,43.0,40.0};
     z_ns[8]={3,3,4,1,5,3,2,3};
     z_ADC[8]={37.0,37.0,40.0,30.0,43.0,37.0,34.0,37.0};*/
}

/*Likelihood::Likelihood(double * psigma1, double * psigma2, double * pmu1, double * pmu2, double * pphi_ns, double  * pphi_ADC, double * pz_ns, double * pz_ADC)
 {
	for (int i=0; i<6; i++)
	{
 sigma1[i]=psigma1[i];
 sigma2[i]=psigma2[i];
 mu1[i]=pmu1[i];
 mu2[i]=pmu2[i];
	}
 
	for (int j=0; j<8; j++)
	{
 phi_ns[j]=pphi_ns[j];
 phi_ADC[j]=pphi_ADC[j];
 z_ns[j]=pz_ns[j];
 z_ADC[j]=pz_ADC[j];
	}
 
 }*/

Likelihood::~Likelihood()
{
    //	cout<<"DESTRUCTOR"<<endl;
    //	for(int i=0;i<8;i++)
    //	{
    //		cout<<phi_ns[i]<<endl;
    //	}
    //	for(int j=0;j<6;j++)
    //	{
    //		cout<<sigma1[i]<<endl;
    //	}
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





int Likelihood::Get_phi_ns(int n)
{
    //cout<<"phi_ns[n]: "<<phi_ns[n]<<endl;
    if (n<0||n>=8)
        return 999999;
    else
        return phi_ns[n];
}

void Likelihood::Set_phi_ns(int * parray)
{
    for (int i=0; i<8; i++)
    {
        phi_ns[i]=parray[i];
    }
}



double Likelihood::Get_phi_ADC(int n)
{
    if (n<0||n>=8)
        return 999999;
    return phi_ADC[n];
}

void Likelihood::Set_phi_ADC(double * parray)
{
    for (int i=0; i<8; i++)
    {
        phi_ADC[i]=parray[i];
    }
}



int Likelihood::Get_z_ns(int n)
{
    if (n<0||n>=8)
        return 999999;
    return z_ns[n];
}

void Likelihood::Set_z_ns(int * parray)
{
    for (int i=0; i<8; i++)
    {
        z_ns[i]=parray[i];
    }
}



double Likelihood::Get_z_ADC(int n)
{
    if (n<0||n>=8)
        return 999999;
    return z_ADC[n];
}

void Likelihood::Set_z_ADC(double * parray)
{
    for (int i=0; i<8; i++)
    {
        z_ADC[i]=parray[i];
    }
}




double Likelihood::Gauss(double sigma_1, double sigma_2, double u, double x)
{
    double ans;
    if (sigma_1==0||sigma_2==0)
    {
        ans=0;
    }
    else
    {
        if(x<=u)
        {
            ans=exp(-(x-u)*(x-u)/(2*sigma_1*sigma_1));
        }
        else
        {
            ans=exp(-(x-u)*(x-u)/(2*sigma_2*sigma_2));
        }
    }
    return ans;
}


double Likelihood::Landau(double sigma_1, double sigma_2, double u1, double u2, double x)
{
    double ans;
    double k=exp(-(u2-u1)*(u2-u1)/(2.0*sigma_2*sigma_2));
    double t=sigma_2*sigma_2/(u2-u1);
    if(sigma_1==0||sigma_2==0||u2<=u1)
    {
        //double t=sigma_2*sigma_2/(u2-u1);
        cout<<"In Landau (sigma_1 or sigma_2=0 or u1<=u2)"<<endl;
        ans=0;
    }
    else
    {
        if (x<=u2)
        {
            ans=Gauss(sigma_1,sigma_2,u1,x);
            
        }
        else
        {
            ans=k*exp(-(x-u2)/t);
        }
    }
    return ans;
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




double Likelihood::fLan(double u1, int ns, double adc, bool is_z)
{
    double normalizations[6]={15.9344, 16.8785, 20.5540, 35.6729, 81.7894, 23.4938}; //ONLY WORKS FOR U1=30! fix!
    
    if(!is_z)  //if z-side, take last members of arrays sigma1, mu1, etc. ns=6 no longer corresponds to a strip # here.
    {
        if(ns<1||ns>5) //ns=6, is_z=false is the same as is_z=true
        {
            cout<<"Strip # out of range"<<endl;
            return 999999;
        }
        
        //cout<<"ns: "<<ns<<endl;
        //cout<<"COEFFICIENTS"<<endl;
        //cout<<"sigma1: "<<sigma1[ns-1]*u1<<", sigma2: "<<sigma2[ns-1]*u1<<", mu1: "<<mu1[ns-1]*u1<<", mu2:"<<mu2[ns-1]*u1<<endl;
        double ans=Landau(sigma1[ns-1]*u1, sigma2[ns-1]*u1, mu1[ns-1]*u1, mu2[ns-1]*u1, adc);
        ans=ans/(normalizations[ns-1]);
        return ans;
    }
    else
    {
        ns=6;
        //cout<<"ns: "<<ns<<endl;
        //cout<<"COEFFICIENTS"<<endl;
        //cout<<"sigma1: "<<sigma1[ns-1]<<", sigma2: "<<sigma2[ns-1]<<", mu1: "<<mu1[ns-1]<<", mu2:"<<mu2[ns-1]<<endl;
        double ans=Landau(sigma1[ns-1]*u1, sigma2[ns-1]*u1, mu1[ns-1]*u1, mu2[ns-1]*u1, adc);
        ans=ans/(normalizations[ns-1]);
        //cout<<"In fLan, agg; ans: "<<ans<<endl;
        return ans;
    }
}


double Likelihood::fLikelihood(double u1)
{
    double L=0;
    
    for (int i=0;i<8;i++)
    {
        if (phi_ns[i]==0||phi_ADC[i]==0)
        {;}
        else
        {
            L=L+log(fLan(u1,phi_ns[i],phi_ADC[i],false));
        }
        if (z_ADC[i]==0)
        {;}
        else
        {
            L=L+log(fLan(u1,z_ns[i],z_ADC[i],true));
        }
    }
    //	cout<<"L: "<<L<<endl;
    double ans=-2*L;
    //	cout<<"ans: "<<ans<<endl;
    return ans;
}


double Likelihood::Find_MaxLkhd(double * pADC, int * pns, int * pside)//Accepts track arrays and finds u1_ag for which
{																	  //-2*ln(L) is minimum.
    for (int i=2;i<1000;i++)//start at i=2 to avoid u1=0. loop from u1_ag = .2 to u1_ag = 100.
    {
        cout<<"In Find_MaxLkhd, i: "<<i<<endl;
        if (fLikelihood(.1*i)<fLikelihood(.1*(i-1)) && fLikelihood(.1*i)<fLikelihood(.1*(i+1)))
            //if -2*ln(L) at a particular value of u1_ag is less than -2*ln(L) for the values of u1_ag immediately above and
            //below it, then that value maximizes the likelihood function.
        {
            double ans = .1*i;
            return ans;
        }
    }
    return 0;//if nothing turns up,which shouldnt happen, return 0
}





void Likelihood::pick_ns() //sets ns array for phi side
{
    float x;
    
    // Set evil seed (initial seed)
    
    int ns_array[8]={0,0,0,0,0,0,0,0};//we want to fill this array with random strip numbers
    double array[6]={0.00,0.12,0.48,0.85,0.97,1.00}; //gives percentages for strip numbers.
    for (int i = 0; i < 6; i++)
    {
        x = (float) rand()/(RAND_MAX + 1.0);
        cout<<"In pick_ns, x: "<<x<<endl;
        for(int j = 1; j < 8 ; j++)
        {
            if (array[j-1]<=x && x<array[j])
            {
                cout<<"strip number: "<<j<<endl;
                ns_array[i]=j;
            }
        }
    }
    /*for (int j=0;j<8;j++)//Check
     {
     cout<<"In pick_ns, ns_array: "<<ns_array[j]<<endl;
     }
     */
    int *pns_array = ns_array;
    Set_phi_ns(pns_array);
}


double Likelihood::Cum_Lan(double u1, int ns, double x, bool is_z) //Landau function integrated out to x.
{
    
    double dy=.01;
    double y=0.0;
    double sumL=0.0;
    
    while (y>=0 && y<=x)
    {
        //cout<<"In Cum_Lan; y, fLan, sumL: "<<y<<", "<<fLan(u1,ns,y,is_z)<<", "<<sumL<<endl;
        y=y+dy;
        sumL=sumL+dy*fLan(u1,ns,y,is_z);
    }
    return sumL;
}


void Likelihood::hit_generator(double u1)//sets ADC arrays for both phi and z sides.
{
    pick_ns();//set phi_ns array
    
    float x;
    float y;
    //Set Phi-side arrays
    for (int i = 0; i < 6; i++)
    {
        x = (float) rand()/(RAND_MAX + 1.0);
        double ns = phi_ns[i];
        
        if (x>Cum_Lan(u1,ns,150,false))//provide for unlikely circumstance that x>.9999, or whatever the upper limit is on Cum_Lan here.
        {
            cout<<"random number x too close to 1"<<endl;
            phi_ADC[i]=0;
        }
        else
        {
            for (int j=0; j<1501; j++)
            {
                cout<<"phi-side inner loop, i,j: "<<i<<", "<<j<<endl;
                double e_plus = (Cum_Lan(u1, ns, .1*(j+1), false) - Cum_Lan(u1, ns, .1*j, false))/2;//set upper uncertainty
                double e_minus = (Cum_Lan(u1, ns, .1*j, false) - Cum_Lan(u1, ns, .1*(j-1), false))/2;//set lower uncertainty
                if (j==0)//for j-1 here, e_minus is not defined, so we set it to 0.
                {
                    e_minus=0;
                }
                
                if (Cum_Lan(u1, ns,.1*j, false)-e_minus<=x && x<Cum_Lan(u1, ns, .1*j, false)+e_plus)
                {
                    phi_ADC[i]=.1*j;
                }
            }
        }
    }
    
    //Set z-side arrays
    for (int k = 0; k < 6; k++)
    {
        //cout<<"z-side, k: "<<k<<endl;
        y = (float) rand()/(RAND_MAX + 1.0);
        double ns = 1; //Set dummy ns, since strip # does not matter on z-side.
        if (y>Cum_Lan(u1,ns,150,true))//provide for unlikely circumstance that y>.9999, or whatever the upper limit is on Cum_Lan here.
        {
            cout<<"random number y too close to 1"<<endl;
            z_ADC[k]=0;
        }
        else
        {
            for (int m=0; m<1501; m++)
            {
                cout<<"z-side inner loop, k,m: "<<k<<", "<<m<<endl;
                //cout<<"(Cum_Lan(u1, ns, .1*m, true): "<<Cum_Lan(u1, ns, .1*m, true)<<endl;
                double e_plus = (Cum_Lan(u1, ns, .1*(m+1), true) - Cum_Lan(u1, ns, .1*m, true))/2;//set upper uncertainty
                double e_minus = (Cum_Lan(u1, ns, .1*m, true) - Cum_Lan(u1, ns, .1*(m-1), true))/2;//set lower uncertainty
                if (m==0)//for j-1 here, e_minus is not defined, so we set it to 0.
                {
                    e_minus=0;
                }
                
                if (Cum_Lan(u1, ns, .1*m, true) - e_minus <= y  &&  y < Cum_Lan(u1, ns, .1*m, true) + e_plus)
                {
                    z_ADC[k]=.1*m;
                }
            }
        }
    }
}




int main(int argc, char* argv[])
{
    //cout<<"check"<<endl;
    srand( (unsigned)time( NULL ) );
    
    Likelihood function1;
    
    double sigma1_a[6]={0.230004,0.188371,0.190532,0.213832,0.355758,0.247916};
    double *psigma1_a=sigma1_a;
    function1.Set_sigma1(psigma1_a);
    
    double sigma2_a[6]={0.178224,0.232363,0.285314,0.631839,1.558563,0.328999};
    double *psigma2_a=sigma2_a;
    function1.Set_sigma2(psigma2_a);
    
    double mu1_a[6]={0.829174,0.962608,1.143752,1.298738,1.726436,1.000000};
    double *pmu1_a=mu1_a;
    function1.Set_mu1(pmu1_a);
    
    double mu2_a[6]={1.047041,1.220573,1.390308,1.935205,3.282429,1.346366};
    double *pmu2_a=mu2_a;
    function1.Set_mu2(pmu2_a);
    
    int phi_ns_a[8]={0,0,0,0,0,0,0,0};
    int *pphi_ns_a=phi_ns_a;
    function1.Set_phi_ns(pphi_ns_a);
    
    double phi_ADC_a[8]={0,0,0,0,0,0,0,0};
    double *pphi_ADC_a=phi_ADC_a;
    function1.Set_phi_ADC(pphi_ADC_a);
    
    int z_ns_a[8]={0,0,0,0,0,0,0,0};
    int *pz_ns_a=z_ns_a;
    function1.Set_z_ns(pz_ns_a);
    
    double z_ADC_a[8]={0,0,0,0,0,0,0,0};
    double *pz_ADC_a=z_ADC_a;
    function1.Set_z_ADC(pz_ADC_a);
    
    
    function1.hit_generator(30); //generate random hit for a given input for u1
    
    //output numbers for random hit, as well as fixed parameters sigma1, sigma2, mu1, mu2
    for (int i=0;i<8;i++)
    {
        cout<<"phi_ns: "<<function1.Get_phi_ns(i)<<endl;
    }
    for (int j=0;j<8;j++)
    {
        cout<<"phi_ADC: "<<function1.Get_phi_ADC(j)<<endl;
    }
    for (int k=0;k<8;k++)
    {
        cout<<"z_ns: "<<function1.Get_z_ns(k)<<endl;
    }
    for (int l=0;l<8;l++)
    {
        cout<<"z_ADC: "<<function1.Get_z_ADC(l)<<endl;
    }
    for (int m=0;m<6;m++)
    {
        cout<<"sigma1: "<<function1.Get_sigma1(m)<<endl;
    }
    for (int n=0;n<6;n++)
    {
        cout<<"sigma2: "<<function1.Get_sigma2(n)<<endl;
    }
    for (int o=0;o<6;o++)
    {
        cout<<"mu1: "<<function1.Get_mu1(o)<<endl;
    }
    for (int p=0;p<6;p++)
    {
        cout<<"mu2: "<<function1.Get_mu2(p)<<endl;
    }
    
    /*
    //find min of likelihood function and see if it agrees well with input u1.
    for (int q=400; q<900;q++)
    {
        cout<<"Likelihood function("<<.05*q<<"): "<<function1.fLikelihood(.05*q)<<endl;
    }
    
    */
    
    
    
    /*//Check for pick_ns
     Likelihood random1;
     random1.pick_ns();
     
     for (int q=0; q<8; q++)
     {
     cout<<"random1: "<<random1.Get_phi_ns(q)<<endl;
     }
     
     
     Likelihood random2;
     random2.pick_ns();
     
     for (int r=0; r<8; r++)
     {
     cout<<"random2: "<<random2.Get_phi_ns(r)<<endl;
     }
     
     Likelihood random3;
     random3.pick_ns();
     
     for (int s=0; s<8; s++)
     {
     cout<<"random3: "<<random3.Get_phi_ns(s)<<endl;
     }
     
     cout<<"function1.fLan(30,3,25,false): " <<function1.fLan(30,3,25,false)<<endl;
     cout<<"Cum_Lan: "<<function1.Cum_Lan(30,3,100,false)<<endl;
     */
    
    return 0;
}


