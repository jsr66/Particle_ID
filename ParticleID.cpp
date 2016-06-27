// ibg_class_lkhd_revised UC final.cpp : Defines the entry point for the console application.
//

/*PARTICLE ID IN THE SILICON VERTEX DETECTOR AT CDF
The purpose of this program is to perform particle identification for an unknown track in the silicon 
vertex detector. Given two
track measurements - particle momentum and measured charge deposition - the program outputs the 
relative probability, in terms of normed differences between expected ionization and measured 
ionization, that the particle is either an electron, muon, pion, kaon, or proton.

Using large numbers of electron, muon, pion, kaon, and proton tracks from our detector, and applying
a particular method of combining the values of charge deposition for the different hits
along a track (e.g. trunc, trunc3, likelihood), we obtain 
data that allows us to determine parameters for the universal curve through fitting, and also to 
encapsulate the functional variation of the upper and lower widths of the distribution with mean 
ionization 
(note: the terms ionization and charge deposition are used here interchangeably). 

The variation of upper and lower widths with mean ionization is determined by grouping proton tracks,
which span the widest range of charge deposition, into momentum 
bins. For each bin, we fit a bifurcated gaussian to the measured distribution of tracks over ionization 
and record the mpv (which we also refer to loosely as the mean) and upper and lower widths. Using this
information,we fit the function s = sqrt(5.2+k*k*I*I), where s is either the upper or lower width,
I is the mean ionization (which varies with betagamma according to the specifications of the universal 
curve) and k is a fitting parameter, to determine how the upper and  lower widths vary with 
ionization, and in turn with betagamma. We find values of k for both the upper and lower widths. 

Having encapsulated into functional form the variation of the distribution widths with mean 
ionization, it is possible, given measurements I_meas and p of the deposition and momentum for
an unknown track, it is possible to get a number of pieces of information. If, for, each pair 
(I_meas, p) that the detector gives us, we guess arbitrarily that the measurements correspond to a 
particular particle, we can determine betagamma from the mass, from beta gamma the expected ionization
(from the universal curve), and from the expected ionization, the upper and lower widths. This
information allows us to determine how many standard deviations the value I_meas is from what we 
expect (known as the normed difference) given our assumption about the kind of particle. We do this 
for each of 5 the particle species  mentioned above, and determine that the species yielding the 
smallest normed difference is the one most likely to have produced the observed track.     
*/



// * don't forget to revise Universal Curve parameters (I_0, k_1,...,k_4) and parameters for s_exp 
//(k_1sis, k_2sis) depending on the method used to combine the hits along a track (e.g. trunc, trunc3,
//likelihood). 



// ibg_class.cpp : Defines the entry point for the console application.
//
// #include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <string>
#include <stdlib.h>
#include <string.h>

const double m_el=0.000511;//electron mass
const double m_mu=0.106;//muon mass
const double m_pi=0.1396;//pion mass
const double m_ka=0.495;//kaon mass
const double m_pr=0.938;//proton mass

class I_pred
{
public:
	I_pred(double I_0_init,double k_1_init,double k_2_init,double k_3_init,double k_4_init,double k_5_init,double k_1sis_init,double k_2sis_init);//constructor function
	I_pred();//constructor function
	~I_pred();//destructor function
	//functions for setting and retrieving various parameters
	double GetI_0() {return I_0;}
	void SetI_0(double value1) {I_0=value1;}
	double Getk_1(){return k_1;}
	void Setk_1(double value2) {k_1=value2;}
	double Getk_2(){return k_2;}
	void Setk_2(double value3) {k_2=value3;}
	double Getk_3(){return k_3;}
	void Setk_3(double value4) {k_3=value4;}
	double Getk_4(){return k_4;}
	void Setk_4(double value5) {k_4=value5;}
	double Getk_5(){return k_5;}
	void Setk_5(double value8) {k_4=value8;}
	double Getk_1sis(){return k_1sis;}
	void Setk_1sis(double value6) {k_1sis=value6;}
	double Getk_2sis(){return k_2sis;}
	void Setk_2sis(double value7) {k_2sis=value7;}
	
	
	double I_exp(double bg);//universal curve parametrization as a function of betagamma
	double s_exp(double bg, bool above);//expected bifurcated gaussian width as a function of bg
	
	//expected charge deposition as a function of momentum for each particle species
	double I_expected_el(double p);
	double I_expected_mu(double p);
	double I_expected_pi(double p);
	double I_expected_ka(double p);
	double I_expected_pr(double p);

	bool is_above(double I_meas,double bg);//given charge deposition and bg, 
	//determines whether it's above or below universal curve

	//difference between measured charge deposition and what we expect from universal curve
	double uncertainty_el(double p,double I_meas);
	double uncertainty_mu(double p,double I_meas);
	double uncertainty_pi(double p,double I_meas);
	double uncertainty_ka(double p,double I_meas);
	double uncertainty_pr(double p,double I_meas);

	//number of sigma difference between measured charge deposition and charge deposition predicted
	//by universal curve
	double normdiff_el(double p,double I_meas);
	double normdiff_mu(double p,double I_meas);
	double normdiff_pi(double p,double I_meas);
	double normdiff_ka(double p,double I_meas);
	double normdiff_pr(double p,double I_meas);

	//tells us particle that gives smallest normed difference, given measurements of particle momentum
	//and charge deposition
	double particle_ID(double p, double I_meas, int particle);

private://set UC parameters and upper and lower width parameters as private variables 
	double I_0;
	double k_1;
	double k_2;
	double k_3;
	double k_4;
	double k_5;
	double k_1sis;
	double k_2sis;
	
};

I_pred::I_pred()//Initialize variables
{
	I_0=0.8782; //lkhd, COT param's, from file: UC_likelihood_bifgaus_fit_COTfunc_newpar_josh.gif
	k_1=149.3;
	k_2=34.57;
	k_3=-11.79;
	k_4=0.3863;
	k_5=88.21;
	k_1sis=0.081;//lkhd values for param of upper, lower widths
	k_2sis=0.157;	
}

I_pred::I_pred(double I_0_init,double k_1_init,double k_2_init,double k_3_init,double k_4_init,double k_5_init, double k_1sis_init,double k_2sis_init)//for resetting param's
{
	I_0=I_0_init;
	k_1=k_1_init;
	k_2=k_2_init;
	k_3=k_3_init;
	k_4=k_4_init;
	k_5=k_5_init;
	k_1sis=k_1sis_init;
	k_2sis=k_2sis_init;
}



I_pred::~I_pred()//destructor
{
}

double I_pred::I_exp(double x)//Universal Curve parametrization
{
	double answer;

	answer= I_0*pow((1+x*x)/(x*x),k_4)*(log(x/(x+k_1))+k_2) + k_3*((x/sqrt(1+x*x))-1) + k_5*((x/sqrt(1+x*x))-1)*((x/sqrt(1+x*x))-1);
	
	//revised, COT parametrization of the Universal Curve

	return answer;
}

double I_pred::s_exp(double bg, bool above)//distribution 
//width as a function of bg, k_1 and k_2 are free param's
{
	if (above)
		return sqrt(5.2+k_2sis*k_2sis*I_exp(bg)*I_exp(bg));
	else
		return sqrt(5.2+k_1sis*k_1sis*I_exp(bg)*I_exp(bg));
}




//expected deposition as a function of momentum p. we use the relativisitic relation p = m (betagamma).
//this depends on the mass of the particle being considered, so we have a separate function for each
//particle species
double I_pred::I_expected_el(double p)
{	
	double betagamma=p/m_el;
	return I_exp(betagamma);
}


double I_pred::I_expected_mu(double p)
{	
	double betagamma=p/m_mu;
	return I_exp(betagamma);
}

double I_pred::I_expected_pi(double p)
{
	double betagamma=p/m_pi;
	return I_exp(betagamma);
}

double I_pred::I_expected_ka(double p)
{	
	double betagamma=p/m_ka;
	return I_exp(betagamma);
}

double I_pred::I_expected_pr(double p)
{	
	double betagamma=p/m_pr;
	return I_exp(betagamma);
}




//this tells us, given I_measured and betagamma, whether the meeasurement of charge deposition is above
//or below what Universal Curve predicts
bool I_pred::is_above(double I_meas,double bg)
{
	bool ans;
	if (I_meas>I_exp(bg))
		ans=true;
	else
		ans=false;
	return ans;
}




//gives expected uncertainty for a particular particle species, given measurements of momentum and 
//charge deposition.
double I_pred::uncertainty_el(double p,double I_meas)
{	
	double gb=p/m_el;
	double ans=s_exp(gb,is_above(I_meas,gb));
	return ans;
}

double I_pred::uncertainty_mu(double p,double I_meas)
{	
	double gb=p/m_mu;
	double ans=s_exp(gb,is_above(I_meas,gb));
	return ans;
}

double I_pred::uncertainty_pi(double p,double I_meas)
{	
	double gb=p/m_pi;
	double ans=s_exp(gb,is_above(I_meas,gb));
	return ans;
}

double I_pred::uncertainty_ka(double p,double I_meas)
{	
	double gb=p/m_ka;
	double ans=s_exp(gb,is_above(I_meas,gb));
	return ans;
 }

double I_pred::uncertainty_pr(double p,double I_meas)
{	
	double gb=p/m_pr;
	double ans=s_exp(gb,is_above(I_meas,gb));
	return ans;
}

  
 	
//gives difference in number of sigma between measured deposition and deposition we expect from 
//universal curve.
double I_pred::normdiff_el(double p,double I_meas)
 {	
	double gb=p/m_el;
	double diff=(I_meas - I_expected_el(p));
	double normdiff=diff/uncertainty_el(p,I_meas);
	return normdiff;
 }
	
double I_pred::normdiff_mu(double p,double I_meas)
 {
	double gb=p/m_mu;
	double diff=(I_meas - I_expected_mu(p));
	double normdiff=diff/uncertainty_mu(p,I_meas);
	return normdiff;
}

double I_pred::normdiff_pi(double p,double I_meas)
 {	
	double gb=p/m_pi;
	double diff=(I_meas - I_expected_pi(p));
	double normdiff=diff/uncertainty_pi(p,I_meas);
	return normdiff;
  }

double I_pred::normdiff_ka(double p,double I_meas)
 {	
	double gb=p/m_ka;
	double diff=(I_meas - I_expected_ka(p));
	double normdiff=diff/uncertainty_ka(p,I_meas);
	return normdiff;
  }

double I_pred::normdiff_pr(double p,double I_meas)
 {	
	double gb=p/m_pr;
	double diff=(I_meas - I_expected_pr(p));
	double normdiff=diff/uncertainty_pr(p,I_meas);
	return normdiff;
  }

 
//values for particle variable in particle_ID function below.
//el: 1
//mu: 2
//pi: 3
//ka: 4
//pr: 5
double I_pred::particle_ID(double p, double I_meas, int particle)//lists normed difference between 
//measured deposition and deposition we expect from universal curve, given measured momentum and
//and charge deposition,as well as a guess (1,2,3,4,or 5) about what kind of particle we're looking at.  
{
	double NormDiff;
	
	switch(particle)//switch matches particle number we input with case number below, and executes
	//command corresponding to particular case number. 
	{
	case 1:
		NormDiff=normdiff_el(p,I_meas);
		break;
	case 2:
		NormDiff=normdiff_mu(p,I_meas);
		break;
	case 3:
		NormDiff=normdiff_pi(p,I_meas);
		break;
	case 4:
		NormDiff=normdiff_ka(p,I_meas);
		break;
	case 5:
		NormDiff=normdiff_pr(p,I_meas);
		break;
	default://in case entered particle integer is not a number from 1 to 5.
		cout<<"Error!"<<endl;
		NormDiff=999;
	}
	
	return NormDiff;
}
	


int main(int argc, char* argv[])
{
	I_pred Ionization1;//declare class member

	//particle labels
	cout<<"electron: EL"<<endl;
	cout<<"muon: MU"<<endl;
	cout<<"pion: PI"<<endl;
	cout<<"kaon: KA"<<endl;
	cout<<"proton: PR"<<endl;
	cout<<"ERROR: X!"<<endl;
	//define input variables for momentum and measured ionization
	double p;
	double I_meas;
	
	//allow input of measured values from user
	cout<<endl;
	cout<<"Enter particle momentum (in GeV/c): ";
	cin>>p;
	cout<<endl;
	cout<<"Enter measured charge deposition of particle (in ADC counts): ";
	cin>>I_meas;
	cout<<endl;
	
	//set array of particle labels
	char particle_array[6][3] ={"EL","MU","PI","KA","PR","X!"};
	
	//initialize array to hold normed differences for different particles
	double ndiff_array[5];

	//output normed differences between measured values and uc-predicted values for each particle
	cout<<"Separations Between Measured Deposition and UC-Predicted Deposition "<<endl; 
	cout<<"(+ sign indicates measured deposition value above predicted value"<<endl;
	cout<<"- sign indicates measured deposition value below predicted value): "<<endl;
	for(int i=1;i<6;i++)
	{	
		double Normdiff=Ionization1.particle_ID(p,I_meas,i);
		ndiff_array[i-1] = Normdiff;
		cout<<particle_array[i-1]<<": "<<Normdiff<<endl;
	}	
	
	//determine which normed difference is smallest and output label for corresponding particle
	double min = abs(ndiff_array[0]);//initialize to electron value
	int number=0;//gives position in particle array...initialize to error (X!), in case no particle is 
					//identified

	for(int j=1;j<6;j++)
	{
		double x = abs(ndiff_array[j-1]);
		if (x<min)//PROBLEM IS HERE. X CAN NEVER BE LESS THAN ITSELF FOR AN ELECTRON
		{	
			min = x;
			number = j-1;
		
		}  
	}
	
	cout<<endl;
	cout<<endl;
	cout<<endl;
	
	
	
	cout<<"particle I.D.:  "<<particle_array[number]<<endl;


	cout<<endl;
	cout<<endl;
	
	return 0;
}







