#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <time.h>
#include "Parton_Shower_Lib.h"
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_exp.h>


#define PI 3.14159265359
#define NewtonsInteration 6
#define DIPOLE_INITIAL_CAPACITY 200000
#define MAXforAreaCalc 100000
#define as 0.1187
#define MZ 91.187
#define Exp 2.7182818284
#define RescaleAlphaFactor 1.0
#define NPcorr 2.0
#define c0 1.0
#define c1 0.0
#define c2 1.0




////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


double LambertW(const double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p,e,t,w;
  //  const int dbgW=0;

  //  if (dbgW) fprintf(stderr,"LambertW: z=%g\n",z);
  if (z<-em1 || isinf(z) || isnan(z)) { 
    fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z); exit(1); 
  }
  if (0.0==z) return 0.0;
  if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
    double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
    return 
     -1.0
     +2.331643981597124203363536062168*r
     -1.812187885639363490240191647568*q
     +1.936631114492359755363277457668*r*q
     -2.353551201881614516821543561516*q2
     +3.066858901050631912893148922704*r*q2
     -4.175335600258177138854984177460*q3
     +5.858023729874774148815053846119*r*q3
     -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }
  /* initial approx for iteration... */
  if (z<1.0) { /* series near 0 */
    p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else 
    w=log(z); /* asymptotic */
  if (z>3.0) w-=log(w); /* useful? */
  for (i=0; i<10; i++) { /* Halley iteration */
    e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  fprintf(stderr,"LambertW: No convergence at z=%g, exiting.\n",z); 
  exit(1);
}

double ProductLog(const double z) {
  return LambertW(z);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////

double Alphas(double MU, double NF, double TheCA){
  double MUeval;
  MUeval = MU;

  return pow(pow(as,-1) + (((11*TheCA)/3. - (4*NF*(0.5))/3.)*log(MUeval*pow(MZ,-1))*pow(PI,-1))/2.,-1);
}


//NOTE:
//mu das/dmu=-2*as*((as^2*beta1)/(16*PI^2) + (as*beta0)/(4*PI) + (as^3*beta2)/(64*PI^3))

double AlphasNLO(double MU, double *params){
  double CA = *(params+1);
  double NF = *(params+2);
  double CF = *(params+3);
  double beta0 = 3.66667 * CA - 0.666667 * NF; 
  double beta1 = 11.333333333*pow(CA,2.0) - 3.3333333333*CA*NF - 2.*CF*NF;
  double LOGmuDIVmz = log(MU/MZ);
  //  if(MU<*(params+6)){
  //    LOGmuDIVmz = log((MU+*(params+6))/MZ);
  //  }
  
  return 1/(1/as + (beta0*LOGmuDIVmz)/(2*PI) + (beta1*log(1 + (as*beta0*LOGmuDIVmz)/(2*PI)))/
	    (4*beta0*PI));
}

//This function has a relative accuracy of 0.0001 or better throughout
//the range of as values defined from
//0.5 GeV to 2000GeV
//accounts for a shift from 5 to 4 flavors at the mb scale 4.18
//otherwise assumes 5 flavor running, QCD
double InvertASqcdNLO(double AS){
  double LogT = log(AS);

  if(as<.2227){
    return pow(Exp, 0.05608178720857118 + 3.1109033999914653*pow(LogT, 2) + 
   1.0782875327246744*pow(LogT, 3) + 0.26455572472979894*pow(LogT, 4))*
      pow(AS, 2.220161822416135);
  }else{
    return pow(Exp, -0.8895201538903502 + 0.39099829467755537*pow(LogT, 2) + 
	       0.000733363196797558*pow(LogT, 3) + 0.09057059567440975*pow(LogT, 4))*
      pow(AS, -0.6527710045891946);
  }
};


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//We solve for RF in the equation t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}]/PI, RF < RI
double ComputeSplitAngle(double Q, double RI, double t, double  NF, double CA){
  double TanRIdiv2 = tan(RI/2);
  double RescaleT = t/RescaleAlphaFactor;
    
  return 2*atan(TanRIdiv2*pow(Exp,6*PI*pow(as,-1)*(-1 + pow(Exp,((-11*CA + 2*NF)*RescaleT)/6.))*pow(11*CA - 2*NF,-1))*pow(Q*TanRIdiv2*pow(MZ,-1),-1 + pow(Exp,((-11*CA + 2*NF)*RescaleT)/6.)));

};

double NoRunComputeSplitAngle(double Q, double RI, double t, double  NF, double CA){
  double TanRIdiv2 = tan(RI/2);
  double RescaleT = t/RescaleAlphaFactor;
    
  return 2*atan(TanRIdiv2*pow(Exp,-PI*t/as));

};

//We solve for RF in the equation t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}]/PI, RF < RI, with alpha_s running at two loop
double ComputeSplitAngleNLO(double Q, double RI, double t, double *params){
  double TanRIdiv2 = tan(RI/2);
  double RescaleT = t/RescaleAlphaFactor;
  double CA = *(params+1);
  double NF = *(params+2);
  double CF = *(params+3);
  double beta0 = 3.66667 * CA - 0.666667 * NF; 
  double beta1 = 11.333333333*pow(CA,2.0) - 3.3333333333*CA*NF - 2.*CF*NF;
  double aini = AlphasNLO(Q*TanRIdiv2, params);
  
  return 2.0*atan(pow(Q,-1.0)*InvertASqcdNLO(-4*beta0*PI*pow(beta1, -1)* ProductLog(-(aini*beta1*pow(beta0, -1)*pow(Exp, -(aini*beta1)/(4*beta0*PI) - (beta0*t)/2)*pow(PI, -1))/4)));
};

//We solve for KT in the equation t(RI,RF) = - Int[as(m)/m,{m,Q,KT}]/PI, RF < RI
double ComputeKTRunningCoupling(double Q, double RI, double t, double  NF, double CA){
  double TanRIdiv2 = tan(RI/2);
  double RescaleT = t/RescaleAlphaFactor;
  
  return RescaleAlphaFactor*MZ*pow(Exp,6*PI*pow(as,-1)*(-1 + pow(Exp,((-11*CA + 2*NF)*RescaleT)/6.))*pow(11*CA - 2*NF,-1))*pow(Q*RI*pow(MZ,-1),pow(Exp,((-11*CA + 2*NF)*RescaleT)/6.));
};

////////////////////////////////////////////////////////////////
double AlphasIntegralLO(double Q, double RI, double RF, double NF, double TheCA){
  double TanRIdiv2 = tan(RI/2);
  double TanRFdiv2 = tan(RF/2);
 
  return 6*log((6*PI + as*(11*TheCA - 2*NF)*log(Q*TanRIdiv2*pow(MZ,-1)))*pow(6*PI + as*(11*TheCA - 2*NF)*log(Q*TanRFdiv2*pow(MZ,-1)),-1))*pow(11*TheCA - 2*NF,-1);

};
/*
*/

////////////////////////////////////////////////////////////////
// Int[as[mu]/mu,{mu,MUI,MUF}]
double AlphasIntegral(double MUF, double MUI, double *params){
  double CA = *(params+1);
  double NF = *(params+2);
  double CF = *(params+3);
  double beta0 = 3.66667 * CA - 0.666667 * NF; 
  double beta1 = 11.333333333*pow(CA,2.0) - 3.3333333333*CA*NF - 2.*CF*NF;
  double AsMUI = AlphasNLO(MUI, params);
  double AsMUF = AlphasNLO(MUF, params);
    
  
  return (AsMUF*beta1)/(2*beta0*beta0) - (AsMUI*beta1)/(2*beta0*beta0) + (2*PI*log(AsMUI/AsMUF))/beta0;
}


// Int[as[mu]^2/mu,{mu,MUI,MUF}]
double AlphasSquIntegral(double MUF, double MUI, double *params){
  double CA = *(params+1);
  double NF = *(params+2);
  double CF = *(params+3);
  double beta0 = 3.66667 * CA - 0.666667 * NF; 
  double beta1 = 11.333333333*pow(CA,2.0) - 3.3333333333*CA*NF - 2.*CF*NF;
  double AsMUI = AlphasNLO(MUI, params);
  double AsMUF = AlphasNLO(MUF, params);
    
return (pow(AsMUF,2.0)*beta1)/(4*beta0*beta0) - (pow(AsMUI,2.0)*beta1)/(4*beta0*beta0) - (2*AsMUF*PI)/beta0 + 
 (2*AsMUI*PI)/beta0;
}



// Int[as[mu]Log[mu/MUI]/mu,{mu,MUI,MUF}]
double AlphasLogIntegral(double MUF, double MUI, double *params){
  double CA = *(params+1);
  double NF = *(params+2);
  double CF = *(params+3);
  double beta0 = 3.66667 * CA - 0.666667 * NF; 
  double beta1 = 11.333333333*pow(CA,2.0) - 3.3333333333*CA*NF - 2.*CF*NF;
  double AsMUI = AlphasNLO(MUI, params);
  double AsMUF = AlphasNLO(MUF, params);
  double L = log(AsMUF/AsMUI);  

  return -(beta1*PI*pow(AsMUI, -1)*pow(beta0, -3)*(2*AsMUF + AsMUI*(-2 - 2*L + pow(L, 2))))/
    2 + 4*pow(AsMUF, -1)*pow(beta0, -2)*pow(PI, 2) - 
    4*pow(AsMUI, -1)*pow(beta0, -2)*pow(PI, 2) + 4*L*pow(AsMUI, -1)*pow(beta0, -2)*
    pow(PI, 2);
}


// Int[as[mu]^2 Log[mu/MUI]/mu,{mu,MUI,MUF}]
double AlphasSquLogIntegral(double MUF, double MUI, double *params){
  double CA = *(params+1);
  double NF = *(params+2);
  double CF = *(params+3);
  double beta0 = 3.66667 * CA - 0.666667 * NF; 
  double beta1 = 11.333333333*pow(CA,2.0) - 3.3333333333*CA*NF - 2.*CF*NF;
  double AsMUI = AlphasNLO(MUI, params);
  double AsMUF = AlphasNLO(MUF, params);
  double L = log(AsMUF/AsMUI);  

  return -(beta1*PI*pow(AsMUI, -1)*(2*AsMUF*AsMUI*(-2 + L) + pow(AsMUF, 2) + 3*pow(AsMUI, 2))*
    pow(beta0, -3))/2 - 4*pow(beta0, -2)*pow(PI, 2) - 
    4*L*pow(beta0, -2)*pow(PI, 2) + 4*AsMUF*pow(AsMUI, -1)*pow(beta0, -2)*pow(PI, 2);

}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

//Must multiply by CF or CA to get quark or gluon respectively,
//we have stripped that factor off
//Gamma[as] = (as/(4*PI) Gamma0 + (as/(4*PI))^2 Gamma1 + ....
//This follows the conventions of 1401.4460, eqn. 3.26 and 3.27
//we give
// Int[Gamma[as]log[mu/MUI],{mu,MUI,MUF}] + log[MUI/Q]Int[Gamma[as],{mu,MUI,MUF}]

double CuspIntegralLL(double Q, double MUF, double MUI, double *params){
  double CA = *(params+1);
  double NF = *(params+2);
  double CF = *(params+3);
  double Gamma0 = 4;
  double Gamma1 = 16.618305243*CA - 4.4444444*NF;
  
  return (Gamma0/(4*PI))*(AlphasLogIntegral(MUF, MUI, params)+log(MUI/Q)*AlphasIntegral(MUF, MUI, params) );
}

double CuspIntegralNLL(double Q, double MUF, double MUI, double *params){
  double CA = *(params+1);
  double NF = *(params+2);
  double CF = *(params+3);
  double Gamma0 = 4;
  double Gamma1 = 16.618305243*CA - 4.4444444*NF;
  
  return (Gamma0/(4*PI))*(AlphasLogIntegral(MUF, MUI, params)+log(MUI/Q)*AlphasIntegral(MUF, MUI, params) )+(Gamma1/pow(4*PI,2.0))*(AlphasSquLogIntegral(MUF, MUI, params)+log(MUI/Q)*AlphasSquIntegral(MUF, MUI, params) );
}


double NoLogCuspIntegralNLL(double MUF, double MUI, double *params){
  double CA = *(params+1);
  double NF = *(params+2);
  double CF = *(params+3);
  double Gamma0 = 4;
  double Gamma1 = 16.618305243*CA - 4.4444444*NF;
  
  return (Gamma0/(4*PI))*(AlphasIntegral(MUF, MUI, params) )+(Gamma1/pow(4*PI,2.0))*(AlphasSquIntegral(MUF, MUI, params) );
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void PSemissions_init(PSEmissionsList *emissions) {
  // initialize size and capacity
  emissions->size = 0;
  emissions->capacity = DIPOLE_INITIAL_CAPACITY;

  // allocate memory for dipole->data
  emissions->emissions_Label = malloc(sizeof(int ) * emissions->capacity * 1);
  emissions->Flavor = malloc(sizeof(int ) * emissions->capacity * 1);
  emissions->Dead = malloc(sizeof(int ) * emissions->capacity * 1);
  emissions->Polar = malloc(sizeof(double) * emissions->capacity * 1);
  emissions->Azimuth = malloc(sizeof(double) * emissions->capacity * 1);
  emissions->Efrac = malloc(sizeof(double) * emissions->capacity * 1);
  emissions->ParentSplitAngle = malloc(sizeof(double) * emissions->capacity * 1);
  emissions->ParentEfrac = malloc(sizeof(double) * emissions->capacity * 1);
}

void PSemissions_append(PSEmissionsList *emissions, double *value, int PartonFlavor, int PartonDead, int Label) {
  int AddedIndex;
  // make sure there's room to expand into
  PSemissions_double_capacity_if_full(emissions);

  // append the value and increment dipole->size
  AddedIndex=emissions->size++;
  emissions->Flavor[AddedIndex] = PartonFlavor;
  emissions->Dead[AddedIndex] = PartonDead;
  emissions->Polar[AddedIndex] = *(value); 
  emissions->Azimuth[AddedIndex] = *(value+1);
  emissions->Efrac[AddedIndex] = *(value+2);
  emissions->ParentSplitAngle[AddedIndex] = *(value+3);
  emissions->ParentEfrac[AddedIndex] = *(value+4);
  emissions->emissions_Label[AddedIndex] = Label;
}

//Do not think I need this routine, wrote it anyways
void PSemissions_edit(PSEmissionsList *emissions, int index, double *value, int PartonFlavor, int PartonDead, int Label) {
  if (index >= emissions->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, emissions->size);
    exit(1);
  }
  // assign value to data at index
  emissions->Flavor[index] = PartonFlavor;
  emissions->Dead[index] = PartonDead;
  emissions->Polar[index] = value[0]; 
  emissions->Azimuth[index] = value[1];
  emissions->Efrac[index] = value[2];
  emissions->ParentSplitAngle[index] = value[3];
  emissions->ParentEfrac[index] = *(value+4);
  emissions->emissions_Label[index] = Label;
}
//this is a kludge, put I do not know how else to do it.
// Make functions to fetch each piece of information about a dipole
int PSemissions_get_Label(PSEmissionsList *emissions, int index) {
  if (index >= emissions->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, emissions->size);
    exit(1);
  }
return emissions->emissions_Label[index];
}

double PSemissions_get_Polar(PSEmissionsList *emissions, int index) {
  if (index >= emissions->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, emissions->size);
    exit(1);
  }
return emissions->Polar[index];
}
double PSemissions_get_Azimuth(PSEmissionsList *emissions, int index) {
  if (index >= emissions->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, emissions->size);
    exit(1);
  }
return emissions->Azimuth[index];
}
double PSemissions_get_Efrac(PSEmissionsList *emissions, int index) {
  if (index >= emissions->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, emissions->size);
    exit(1);
  }
return emissions->Efrac[index];
}
double PSemissions_get_SplitAngle(PSEmissionsList *emissions, int index) {
  if (index >= emissions->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, emissions->size);
    exit(1);
  }
return emissions->ParentSplitAngle[index];
}

double PSemissions_get_ParentEfrac(PSEmissionsList *emissions, int index) {
  if (index >= emissions->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, emissions->size);
    exit(1);
  }
return emissions->ParentEfrac[index];
}

int PSemissions_get_Flavor(PSEmissionsList *emissions, int index) {
  if (index >= emissions->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, emissions->size);
    exit(1);
  }
return emissions->Flavor[index];
}

int PSemissions_get_Death(PSEmissionsList *emissions, int index) {
  if (index >= emissions->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, emissions->size);
    exit(1);
  }
return emissions->Dead[index];
}

///////////////////////////////////////////////////////

void RetrieveParton(PSEmissionsList *emissions, int TheChosen, double *Momentum, int *Flavor, int *Death, int *Label){
	  *(Momentum) = PSemissions_get_Polar(emissions, TheChosen);
	  *(Momentum+1) = PSemissions_get_Azimuth(emissions, TheChosen);
	  *(Momentum+2) = PSemissions_get_Efrac(emissions, TheChosen);
	  *(Momentum+3) = PSemissions_get_SplitAngle(emissions,TheChosen);
	  *(Momentum+4) = PSemissions_get_ParentEfrac(emissions,TheChosen);
	  *(Flavor) = PSemissions_get_Flavor(emissions, TheChosen);
	  *(Death) = PSemissions_get_Death(emissions, TheChosen);
	  *(Label) = PSemissions_get_Label(emissions, TheChosen);//This is zero or one. If one, this is the current WTA axis
	  return;
}


///////////////////////////////////////////////////////
int PSemissions_double_capacity_if_full(PSEmissionsList *emissions) {
  if (emissions->size >= emissions->capacity) {
    // double vector->capacity and resize the allocated memory accordingly
    printf("That is alot of emissions."); 
    printf(" Current Carrying Capacity: %d.\n", (emissions->capacity)*2);
    emissions->capacity *= 2;
 
    emissions->emissions_Label = realloc(emissions->emissions_Label,sizeof(int) * emissions->capacity * 1);
    emissions->Flavor = realloc(emissions->Flavor,sizeof(int) * emissions->capacity * 1);
    emissions->Dead = realloc(emissions->Flavor,sizeof(int) * emissions->capacity * 1);
    emissions->Polar = realloc(emissions->Polar,sizeof(double) * emissions->capacity * 1);
    emissions->Azimuth = realloc(emissions->Azimuth,sizeof(double) * emissions->capacity * 1);
    emissions->Efrac = realloc(emissions->Efrac,sizeof(double) * emissions->capacity * 1);
    emissions->ParentSplitAngle = realloc(emissions->ParentSplitAngle,sizeof(double) * emissions->capacity * 1);
    emissions->ParentEfrac = realloc(emissions->ParentSplitAngle,sizeof(double) * emissions->capacity * 1);
  }
  return 0;
}

void PSemissions_free(PSEmissionsList *emissions) {
  free(emissions->emissions_Label);
  free(emissions->Flavor);
  free(emissions->Dead);
  free(emissions->Polar);
  free(emissions->Azimuth);
  free(emissions->Efrac);
  free(emissions->ParentSplitAngle);
  free(emissions->ParentEfrac);
}

////////////////////////////////////////////////////////////////

int  WriteNewPSEmission(PSEmissionsList *emissions, double *Momentum, int PartonFlavor,
			int PartonDead, int *EmissionsLabel){
      //We count this gluon as a jet
      *(EmissionsLabel) = *(EmissionsLabel)+1;
      PSemissions_append(emissions, Momentum, PartonFlavor, PartonDead, *(EmissionsLabel) );

      return 0;
};



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


double IntToZmaxPgTOgg(double Zmax, void *params){
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  double CA = *(PointerToParams+1);
  double Zmin = CUTOFF;

  //return -2*CA*Zmax + 2*CA*Zmin - 2*CA*log(1 - Zmax) + 2*CA*log(1 - Zmin) + (CA*pow(Zmax,2))/2. - (CA*pow(Zmax,3))/3. - (CA*pow(Zmin,2))/2. + (CA*pow(Zmin,3))/3.;
  return 2.0*CA*log(Zmax/Zmin);

}

double IntToZmaxPgTOggSimplified(double Zmax, void *params){
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  double CA = *(PointerToParams+1);
  double Zmin = CUTOFF;

  return  -2.0*CA*log(Zmax/Zmin);

}

double PgTOgg(double Z, void *params){
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  double CA = *(PointerToParams+1);

  //return 2*CA*(((1 - Z)*Z)/2. + Z*pow(1 - Z,-1));
  return 2*CA/Z;

}

double PgTOggSimplified(double Z, void *params){
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  double CA = *(PointerToParams+1);

  return 2*CA/Z;

}

double IntToZmaxPgTOqqbar(double Zmax, void *params){
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  double NF = *(PointerToParams+2);
  double Zmin = CUTOFF;

  //return (NF*(Zmax - Zmin)*(3 - 3*Zmax - 3*Zmin + 2*Zmax*Zmin + 2*pow(Zmax,2) + 2*pow(Zmin,2)))/6.;
  return 0.0;

}

double PgTOqqbar(double Z, void *params){
  double *PointerToParams = (double *) params;
  double NF = *(PointerToParams+2);

  //return (NF*(1 - 2*Z + 2*pow(Z,2)))/2.;
  return 0.0;
}



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

double IntToZmaxPqTOqg(double Zmax, void *params){
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  double CF = *(PointerToParams+3);
  double Zmin = CUTOFF;
 
  //return (CF*(-4*Zmax + 4*Zmin + 4*log(Zmax) - 4*log(Zmin) + pow(Zmax,2) - pow(Zmin,2)))/2.;
  return 2.0*CF*gsl_sf_expint_Ei(Zmax) - 2.0*CF*gsl_sf_expint_Ei(Zmin); // e^z variant
  //return CF*gsl_sf_expint_Ei(Zmax*Zmax) - CF*gsl_sf_expint_Ei(Zmin*Zmin); // e^z2 variant

}


double PqTOqg(double Z, void *params){
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  double CF = *(PointerToParams+3);

  //return -2*CF + CF*Z + 2*CF*pow(Z,-1);
  return 2.0*CF*exp(Z)/Z; // e^z variant
  //return 2.0*CF*exp(Z*Z)/Z; // e^z2 variant

}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////



double UnitStep(const double x){
  if(x > 0.0){
    return 1.0;
  }else{
    return 0.0;
  };
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

double InvertPgToALLDistr(double RND, void *params){
  int i;
  double accuracy, f, df, CurrentZ, NewZ,A,B;
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  double CA = *(PointerToParams+1);
 
  //We iterate at least NewtonsInteration times
  //It often takes some time before we move away from the initial value
  //so we only increase the interation variable once we drop below
  //the accuracy threshold
  /*i = 0;
  CurrentZ = 1 - CUTOFF;
  double NORM = (IntToZmaxPgTOgg(1-CUTOFF, params)+IntToZmaxPgTOqqbar(1-CUTOFF, params));

  while(i<=NewtonsInteration){
    f = RND -(IntToZmaxPgTOgg(CurrentZ, params)+IntToZmaxPgTOqqbar(CurrentZ, params))/NORM;
    df = -1*(PgTOgg(CurrentZ, params) +  PgTOqqbar(CurrentZ, params) )/NORM;
    NewZ = CurrentZ - f/df;
    accuracy = sqrt(pow(NewZ-CurrentZ,2))/NewZ;
    if(accuracy < 0.0001){i = i + 1;};
    CurrentZ = NewZ;
  }*/

  //replace with simplified version for CA/z splitting
  i = 0;
  CurrentZ = 1 - CUTOFF;
  double NORM = (IntToZmaxPgTOgg(1-CUTOFF, params));//(IntToZmaxPgTOgg(1-CUTOFF, params)+IntToZmaxPgTOqqbar(1-CUTOFF, params))
  double LowerBoundPos = CUTOFF;
  double UpperBoundPos = 1.0-CUTOFF;

  while(i<=NewtonsInteration){
    CurrentZ = 0.5*(LowerBoundPos+UpperBoundPos);
    f = RND - (IntToZmaxPgTOgg(CurrentZ, params))/NORM;//RND -(IntToZmaxPgTOgg(CurrentZ, params)+IntToZmaxPgTOqqbar(CurrentZ, params))/NORM
    accuracy = sqrt(pow(LowerBoundPos-CurrentZ,2))/CurrentZ;
    if(accuracy < 0.0001){i = i + 1;};
    if( f > 0 ){
      LowerBoundPos = CurrentZ;
    }else{
      UpperBoundPos = CurrentZ;
    };
  }
    
  //  if(CurrentZ>0.9){
  //    CurrentZ = BisectionForGluon(RND, params);
  //  }

  return CurrentZ;

}

double InvertPgToALLDistrSimplified(double RND, void *params){
  int i;
  double accuracy, f, df, CurrentZ, NewZ,A,B;
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  double CA = *(PointerToParams+1);
 
  //We iterate at least NewtonsInteration times
  //It often takes some time before we move away from the initial value
  //so we only increase the interation variable once we drop below
  //the accuracy threshold
  i = 0;
  CurrentZ = 1 - CUTOFF;
  double NORM = (IntToZmaxPgTOggSimplified(1-CUTOFF, params));//(IntToZmaxPgTOgg(1-CUTOFF, params)+IntToZmaxPgTOqqbar(1-CUTOFF, params))
  double LowerBoundPos = CUTOFF;
  double UpperBoundPos = 1.0-CUTOFF;

  while(i<=NewtonsInteration){
    CurrentZ = 0.5*(LowerBoundPos+UpperBoundPos);
    f = RND - (IntToZmaxPgTOggSimplified(CurrentZ, params))/NORM;//RND -(IntToZmaxPgTOgg(CurrentZ, params)+IntToZmaxPgTOqqbar(CurrentZ, params))/NORM
    accuracy = sqrt(pow(LowerBoundPos-CurrentZ,2))/CurrentZ;
    if(accuracy < 0.0001){i = i + 1;};
    if( f > 0 ){
      LowerBoundPos = CurrentZ;
    }else{
      UpperBoundPos = CurrentZ;
    };
  }
    
  //  if(CurrentZ>0.9){
  //    CurrentZ = BisectionForGluon(RND, params);
  //  }

  return CurrentZ;


}

double InvertPqToALLDistr(double RND, void *params){
  int i;
  double accuracy, f, df, CurrentZ, NewZ;
  double *PointerToParams = (double *) params;
  double CUTOFF = *(PointerToParams);
  
  //We iterate at least NewtonsInteration times
  //It often takes some time before we move away from the initial value
  //so we only increase the interation variable once we drop below
  //the accuracy threshold
  /*i = 0;
  CurrentZ = CUTOFF;
  double NORM = IntToZmaxPqTOqg(1-CUTOFF, params);

  while(i<=NewtonsInteration){
    f = RND - IntToZmaxPqTOqg(CurrentZ, params)/NORM;
    df = -1*PqTOqg(CurrentZ, params)/NORM;

    NewZ = CurrentZ - f/df;
    accuracy = sqrt(pow(NewZ-CurrentZ,2))/NewZ;
    if(accuracy < 0.0001){i = i + 1;};
    CurrentZ = NewZ;
  }*/

  i = 0;
  CurrentZ = 1 - CUTOFF;
  double NORM = (IntToZmaxPqTOqg(1-CUTOFF, params));//(IntToZmaxPgTOgg(1-CUTOFF, params)+IntToZmaxPgTOqqbar(1-CUTOFF, params))
  double LowerBoundPos = CUTOFF;
  double UpperBoundPos = 1.0-CUTOFF;

  while(i<=NewtonsInteration){
    CurrentZ = 0.5*(LowerBoundPos+UpperBoundPos);
    f = RND - (IntToZmaxPqTOqg(CurrentZ, params))/NORM;//RND -(IntToZmaxPgTOgg(CurrentZ, params)+IntToZmaxPgTOqqbar(CurrentZ, params))/NORM
    accuracy = sqrt(pow(LowerBoundPos-CurrentZ,2))/CurrentZ;
    if(accuracy < 0.0001){i = i + 1;};
    if( f > 0 ){
      LowerBoundPos = CurrentZ;
    }else{
      UpperBoundPos = CurrentZ;
    };
  }

  //  if(CurrentZ<0.1){
  //    CurrentZ = BisectionForQuark(RND, params);
  //  }

  return CurrentZ;

  //The return Z is always the Z of the gluon
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*

//This guarantees that the local transverse momentum,
//that is, the transverse momentum in the plane perpendicular to the parent
//sums to zero
void TransverseMomentumAngles(double SplitAng, double Z, double *TheAngles){
  double CurrentTheta = SplitAng/2;
  double NewTheta;
  double acc;
  int i;
  i = 0;
  while( (i<=4) ){
    NewTheta = CurrentTheta -(Z*sin (SplitAng - CurrentTheta) - (1-Z)* sin(CurrentTheta)  )/(-Z*cos (SplitAng - CurrentTheta) - (1-Z)* cos(CurrentTheta)  );
    acc = sqrt(pow(NewTheta-CurrentTheta, 2.0)) /  sqrt(pow(NewTheta, 2.0));
    if(acc <  0.001 ){i=i+1;}
    CurrentTheta = NewTheta;
  }
  //  *(TheAngles) = SplitAng - NewTheta;//THIS IS THE ANGLE THAT THE DAUGHTER WITH MOMENTUM FRACTION Z MAKES TO THE PARENT
  //  *(TheAngles+1) = NewTheta;
  
  *(TheAngles) = SplitAng-CurrentTheta;//THIS IS THE ANGLE THAT THE DAUGHTER WITH MOMENTUM FRACTION Z MAKES TO THE PARENT
  *(TheAngles+1) = CurrentTheta;
  return;
}

*/

void TransverseMomentumAngles(double SplitAngle, double Z, double *TheAngles){


  double ZfracThetaToParent = acos((Z + (1.0 - Z)*cos(SplitAngle))/sqrt(1.0 - 2.0*(1 - Z)*Z*(1.0 - cos(SplitAngle))));
  double OnemZfracThetaToParent = SplitAngle - ZfracThetaToParent;
  
  *(TheAngles) = ZfracThetaToParent;//THIS IS THE ANGLE THAT THE DAUGHTER WITH MOMENTUM FRACTION Z MAKES TO THE PARENT
  *(TheAngles+1) = OnemZfracThetaToParent;
  return;
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void DetermineSplittingDirections(double *AngularPos, double *Daughter1, double *Daughter2,
				  double Z, double SplitAngle, const gsl_rng *r){

  double ParentDir[3], Daughter1dir[3], Daughter2dir[3];
  double NormalA[3], NormalB[3], RandomBasisVec[3];
  double NORMALIZEIT, DotToRandom, RND;
  double TransMomAngles[2]={0,0};
  
    TransverseMomentumAngles(SplitAngle,Z,TransMomAngles);

    ParentDir[0] = sin(AngularPos[0])*cos(AngularPos[1]);
    ParentDir[1] = sin(AngularPos[0])*sin(AngularPos[1]);
    ParentDir[2] = cos(AngularPos[0]);

    //We generate a random vector; this will most likely be linearly independent of the parent direction
    //then we orthonormalize to get the basis vectors in the transverse plane of the parent direction
    //Since this vector was random, we get a random distribution in the azimuthal of the transverse plane
    RandomBasisVec[0] = 2 * gsl_rng_uniform(r)-1;
    RandomBasisVec[1] = 2 * gsl_rng_uniform(r)-1;
    RandomBasisVec[2] = 2 * gsl_rng_uniform(r)-1;
      
    DotToRandom = RandomBasisVec[0]*ParentDir[0] + RandomBasisVec[1]*ParentDir[1] + RandomBasisVec[2]*ParentDir[2];
    NORMALIZEIT = sqrt( pow(DotToRandom*ParentDir[0]-RandomBasisVec[0],2.)+pow(DotToRandom*ParentDir[1]-RandomBasisVec[1],2.)+pow(DotToRandom*ParentDir[2]-RandomBasisVec[2],2.)  );
      
    NormalA[0] = (DotToRandom*ParentDir[0]-RandomBasisVec[0])/ NORMALIZEIT;
    NormalA[1] = (DotToRandom*ParentDir[1]-RandomBasisVec[1])/ NORMALIZEIT;
    NormalA[2] = (DotToRandom*ParentDir[2]-RandomBasisVec[2])/ NORMALIZEIT;

    //Compute cross product to get third orthonormal vector
    //WE DONT ACTUALLY NEED THIS
    //BUT THIS GUARANTEES azimuthal symmetry in transverse plan
    //once we apply the rotation
    NormalB[0] = ParentDir[2]*NormalA[1] - NormalA[2]*ParentDir[1];
    NormalB[1] = ParentDir[0]*NormalA[2] - NormalA[0]*ParentDir[2];
    NormalB[2] = ParentDir[1]*NormalA[0] - NormalA[1]*ParentDir[0];
  
    RND = 2 * PI * gsl_rng_uniform(r);

    //Sum of these two vectors have zero transverse momentum
    //compared to direction of the parent 
    //angle between them is the splitting angle determined by the evolution
    Daughter1dir[0] = cos( TransMomAngles[0] ) * ParentDir[0] + sin( TransMomAngles[0] ) *( cos(RND) * NormalA[0] + sin(RND) * NormalB[0]);
    Daughter1dir[1] = cos( TransMomAngles[0] ) * ParentDir[1] + sin( TransMomAngles[0] ) *( cos(RND) * NormalA[1] + sin(RND) * NormalB[1]);
    Daughter1dir[2] = cos( TransMomAngles[0] ) * ParentDir[2] + sin( TransMomAngles[0] ) *( cos(RND) * NormalA[2] + sin(RND) * NormalB[2]);

    Daughter2dir[0] = cos( TransMomAngles[1] ) * ParentDir[0] - sin( TransMomAngles[1] ) *( cos(RND) * NormalA[0] + sin(RND) * NormalB[0]);
    Daughter2dir[1] = cos( TransMomAngles[1] ) * ParentDir[1] - sin( TransMomAngles[1] ) *( cos(RND) * NormalA[1] + sin(RND) * NormalB[1]);
    Daughter2dir[2] = cos( TransMomAngles[1] ) * ParentDir[2] - sin( TransMomAngles[1] ) *( cos(RND) * NormalA[2] + sin(RND) * NormalB[2]);

    /*
    Include this snippet to check local momentum conservation.
    double Test = sqrt(pow(Z*sin( TransMomAngles[0] ) - (1-Z)*sin( TransMomAngles[1] ),2.0));
    if(Test>0.0001){printf("NO! Local momentum conservation violation %f\n",Test);}
    */
    
    Daughter1[0] = acos(Daughter1dir[2]);//Polar
    Daughter1[1] = PI+atan2(Daughter1dir[1],Daughter1dir[0]);//Azimuth

    Daughter2[0] = acos(Daughter2dir[2]);//Polar
    Daughter2[1] = PI+atan2(Daughter2dir[1],Daughter2dir[0]);//Azimuth

return;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void DetermineSplittingDirectionsNoRecoil(double *AngularPos, double *Daughter1, double *Daughter2,
					  double Z, double SplitAngle, const gsl_rng *r){

  double ParentDir[3], Daughter1dir[3], Daughter2dir[3];
  double NormalA[3], NormalB[3], RandomBasisVec[3];
  double NORMALIZEIT, DotToRandom, RND;
  double TransMomAngles[2]={0,0};
  
  ParentDir[0] = sin(AngularPos[0])*cos(AngularPos[1]);
  ParentDir[1] = sin(AngularPos[0])*sin(AngularPos[1]);
  ParentDir[2] = cos(AngularPos[0]);
  
  //We generate a random vector; this will most likely be linearly independent of the parent direction
  //then we orthonormalize to get the basis vectors in the transverse plane of the parent direction
  //Since this vector was random, we get a random distribution in the azimuthal of the transverse plane
  RandomBasisVec[0] = 2 * gsl_rng_uniform(r)-1;
  RandomBasisVec[1] = 2 * gsl_rng_uniform(r)-1;
  RandomBasisVec[2] = 2 * gsl_rng_uniform(r)-1;
      
  DotToRandom = RandomBasisVec[0]*ParentDir[0] + RandomBasisVec[1]*ParentDir[1] + RandomBasisVec[2]*ParentDir[2];
  NORMALIZEIT = sqrt( pow(DotToRandom*ParentDir[0]-RandomBasisVec[0],2.)+pow(DotToRandom*ParentDir[1]-RandomBasisVec[1],2.)+pow(DotToRandom*ParentDir[2]-RandomBasisVec[2],2.)  );
      
  NormalA[0] = (DotToRandom*ParentDir[0]-RandomBasisVec[0])/ NORMALIZEIT;
  NormalA[1] = (DotToRandom*ParentDir[1]-RandomBasisVec[1])/ NORMALIZEIT;
  NormalA[2] = (DotToRandom*ParentDir[2]-RandomBasisVec[2])/ NORMALIZEIT;

  //Compute cross product to get third orthonormal vector
  //WE DONT ACTUALLY NEED THIS
  //BUT THIS GUARANTEES azimuthal symmetry in transverse plan
  //once we apply the rotation
  NormalB[0] = ParentDir[2]*NormalA[1] - NormalA[2]*ParentDir[1];
  NormalB[1] = ParentDir[0]*NormalA[2] - NormalA[0]*ParentDir[2];
  NormalB[2] = ParentDir[1]*NormalA[0] - NormalA[1]*ParentDir[0];
  
  RND = 2 * PI * gsl_rng_uniform(r);

  //Sum of these two vectors have zero transverse momentum
  //compared to direction of the parent 
  //angle between them is the splitting angle determined by the evolution
  Daughter1dir[0] = cos( SplitAngle ) * ParentDir[0] + sin( SplitAngle ) *( cos(RND) * NormalA[0] + sin(RND) * NormalB[0]);
  Daughter1dir[1] = cos( SplitAngle ) * ParentDir[1] + sin( SplitAngle ) *( cos(RND) * NormalA[1] + sin(RND) * NormalB[1]);
  Daughter1dir[2] = cos( SplitAngle ) * ParentDir[2] + sin( SplitAngle ) *( cos(RND) * NormalA[2] + sin(RND) * NormalB[2]);

  Daughter2dir[0] =  ParentDir[0];
  Daughter2dir[1] =  ParentDir[1];
  Daughter2dir[2] =  ParentDir[2];

  Daughter1[0] = acos(Daughter1dir[2]);//Polar
  Daughter1[1] = PI+atan2(Daughter1dir[1],Daughter1dir[0]);//Azimuth

  Daughter2[0] = acos(Daughter2dir[2]);//Polar
  Daughter2[1] = PI+atan2(Daughter2dir[1],Daughter2dir[0]);//Azimuth

return;
};


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

double ComputeVirtuality(double Efrac, double Q, double Z, double SplitAngle){

  return  Efrac * Q * sqrt( (2.0) * Z * (1.0-Z) * ( 1.0 - cos( SplitAngle) )  );
  
}

double ComputeKT(double Efrac, double Q, double Z, double SplitAngle){

  return  (Q*Efrac*Z*(1.0-Z)*sin(SplitAngle));///sqrt( 1.0 - 2.0*Z*(1.0-Z)*( 1.0 - cos(SplitAngle) ) );
  
}

double MinZFunction(double Efrac, double Q, double PhysSplitAngle, double MinQ, double p, double q ){

  return MinQ/(  pow(Efrac, q) * Q * pow(1.0*sin(PhysSplitAngle ) ,p) );

}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
 
//We choose an emission, and compute its splitting
void BuildSplitFunctionInversion(double *InvertPg, double *InvertPq, int NumInverts, double *params, const gsl_rng *r){
  int i;
  double Z,RND;
  
  printf("For Gluon.....\n"  );
  for(i=0; i<NumInverts;i=i+1){
    RND = (double)( i )/(double)( NumInverts );
    Z =  InvertPgToALLDistr( RND, params);
    *(InvertPg+i) = Z;
  }

  printf("For Quark.....\n"  );
  for(i=0; i<NumInverts;i=i+1){
    RND = (double)( i )/(double)( NumInverts );
    Z =  InvertPqToALLDistr( RND, params);
    *(InvertPq+i) = Z;
  }
return;
}

 
//We choose an emission, and compute its splitting
void BuildSplitFunctionInversionSimplified(double *InvertPg, double *InvertPq, int NumInverts, double *params, const gsl_rng *r){
  int i;
  double Z,RND;
  
  printf("For Gluon.....\n"  );
  for(i=0; i<NumInverts;i=i+1){
    RND = (double)( i )/(double)( NumInverts );
    Z =  InvertPgToALLDistrSimplified( RND, params);
    *(InvertPg+i) = Z;
  }

  return;
}



////////////////////////////////////////////////////////////////

void SplitSelectedParton(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
			 int *TheChosen, double *params, const gsl_rng *r, int *EndEvent, int Mode,
			 double *InvertPg, double *InvertPq, int NumInverts){
  int i, CurrentFlavor, FoundaZ, WTA1, WTA2,TryToSplit, CanSplit, PartonADead, PartonBDead;
  int AssignFlavor[2];
  double AngularPos[2];
  double MinZ,MaxAngle, MinQ, CF, CA, NF, Q, Rmax,
    TotalP, Efrac, CoinFlipForEmission, DT, prob, Z,
    RND,SplitAngle;
  double ThePartonPhaseSpace[(emissions->size)+1];
  int TranslateActiveEmissions[(emissions->size)+1];
  double Daughter1[5],Daughter2[5];
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);

  double Pg = (1.0)*(IntToZmaxPgTOgg(1-MinZ, params)+IntToZmaxPgTOqqbar(1-MinZ, params));
  double Pq = (IntToZmaxPqTOqg(1-MinZ, params));

  TotalP = 0;
  //First we must establish the total phase space to radiate
  CanSplit = -1;
  
  for(i=0;i< (emissions->size);i++){    
    if( (PSemissions_get_Death(emissions, i ) == 0)  ){

      //This array is important, it allows us to find for the given ActiveEmission
      //What is the corresponding position in the Total emissions List of all partons in the event
      //That way we do not have to create a list of active emissions
      //We just have a list that tells us how to find the emissions which can split!
      //If it is not on that list, it can't split
      CanSplit = CanSplit + 1;
      TranslateActiveEmissions[CanSplit] = i;
      
      if( PSemissions_get_Flavor(emissions, i) == 0){
	ThePartonPhaseSpace[CanSplit] = Pg; //HAVE TO REMEMBER TO INDEX BY CanSplit here
      }else{
	ThePartonPhaseSpace[CanSplit] = Pq; //HAVE TO REMEMBER TO INDEX BY CanSplit here
      };//IF
      TotalP = TotalP + ThePartonPhaseSpace[CanSplit];
    }//IF
  };//for

  //(CanSplit+1) is the number of partons that can split currently
  //No parton will be able to split
  if( (CanSplit+1) <= 0 ){
    *(EndEvent) = 1;
    return;
  };

  DT = - log( gsl_rng_uniform_pos(r)) /  TotalP ;
  *(t) = *(t) + DT;

  CoinFlipForEmission = gsl_rng_uniform_pos(r);
  
  i=-1;
  while(i < (CanSplit+1) ){
    //Probability for i-th dipole
    i++;
    prob = ThePartonPhaseSpace[i]/TotalP;
    //If probability is great than the coin flip, the parton is chosen
    if(CoinFlipForEmission-prob<0){
      *(TheChosen) = TranslateActiveEmissions[i];//We have chosen the ith active parton, now we find its index where it sits in the master list
      i = (emissions->size)+10;//This will force exiting the while loop
    }else{
      //else, we chew off that bit of probability from coin flip and go to the next parton.
      CoinFlipForEmission = CoinFlipForEmission-prob;
    };//ELSEIF
  };//while

  //These are the quantum numbers of the chosen parton
  AngularPos[0] = PSemissions_get_Polar(emissions, *(TheChosen));
  AngularPos[1] = PSemissions_get_Azimuth(emissions, *(TheChosen));
  CurrentFlavor = PSemissions_get_Flavor(emissions, *(TheChosen));
  Efrac = PSemissions_get_Efrac(emissions, *(TheChosen) );

  //We now compute the splitting directions
  //We note that t is defined to be:
  //t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}], RF < RI
  //so now we solve for RF to get split angle
  SplitAngle = ComputeSplitAngle(Q, Rmax, *t, NF, CA);

  //Now we determine the energy fraction

  if(CurrentFlavor == 0){
    //We chose a gluon! Fun!

      //Now we decide whether to split to gg or qqbar
      //We split to gluons
      //We generate a random number, then invert the integral of the PgTOgg splitting kernel to find the
      //Correct energy fraction
      RND = gsl_rng_uniform_pos(r);
      Z = *(InvertPg + (int)floor(RND*(double)(NumInverts)) );

      Daughter1[2] = Z * Efrac;
      Daughter2[2] = (1-Z) * Efrac;

      if(Z >= 0.5){
	WTA1 = 1;
	WTA2 = 0;
      }else{
	WTA1 = 0;
	WTA2 = 1;     
      }
      
	RND = gsl_rng_uniform_pos(r);
	//Now we decide whether it splits to qqbar or g
	if( RND < ( PgTOgg(Z, params) )/( PgTOgg(Z, params) + PgTOqqbar(Z, params) ) ){
	  AssignFlavor[0] = 0;
	  AssignFlavor[1] = 0;
	}else{
	  AssignFlavor[0] = gsl_rng_uniform_int(r, (int)ceil(NF) ) + 1;
	  AssignFlavor[1] = -AssignFlavor[0];
	};
  } else {
    //We chose a quark! Fun!
    FoundaZ = 0;
    RND = gsl_rng_uniform_pos(r);
    Z = *(InvertPq + (int)floor(RND*(double)(NumInverts)) );

    Daughter1[2] = Z * Efrac;//THIS IS THE GLUON
    Daughter2[2] = (1-Z) * Efrac;//This is the quark or anti-quark

    if(Z >= 0.5){
      WTA1 = 1;
      WTA2 = 0;
    }else{
      WTA1 = 0;
      WTA2 = 1;     
    }

    AssignFlavor[0] = 0;
    AssignFlavor[1] = CurrentFlavor;

  };//IF

  *(Angle) = SplitAngle;
  DetermineSplittingDirections(AngularPos, Daughter1, Daughter2, Z, SplitAngle, r);
  PartonADead = 0;
  PartonBDead = 0;
  if(Daughter1[2]>MinZ){
    PartonADead = 0;
  }else{
    PartonADead = 1;
  }
  if(Daughter2[2]>MinZ){
    PartonBDead = 0;
  }else{
    PartonBDead = 1;
  }

  Daughter1[3] = SplitAngle;
  Daughter2[3] = SplitAngle;
  Daughter1[4] = Efrac;
  Daughter2[4] = Efrac;
  PSemissions_edit(DaughterEmissions, 0, Daughter1, AssignFlavor[0], PartonADead, WTA1); 
  PSemissions_edit(DaughterEmissions, 1, Daughter2, AssignFlavor[1], PartonBDead, WTA2);

  return;

}


////////////////////////////////////////////////////////////////

void SplitSelectedPartonNoConserve(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
				   int *TheChosen, double *params, const gsl_rng *r, int *EndEvent, int Mode,
				   double *InvertPg, double *InvertPq, int NumInverts){
  int i, CurrentFlavor, FoundaZ, WTA1, WTA2,TryToSplit, CanSplit, PartonADead, PartonBDead;
  int AssignFlavor[2];
  double AngularPos[2];
  double MinZ,MaxAngle, MinQ, CF, CA, NF, Q, Rmax,
    TotalP, Efrac, CoinFlipForEmission, DT, prob, Z,
    RND,SplitAngle;
  double ThePartonPhaseSpace[(emissions->size)+1];
  int TranslateActiveEmissions[(emissions->size)+1];
  double Daughter1[5],Daughter2[5];
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);

  double Pg = IntToZmaxPgTOggSimplified(1-MinZ, params);
  double Pq = 0;

  TotalP = 0;
  //First we must establish the total phase space to radiate
  CanSplit = -1;
  
  for(i=0;i< (emissions->size);i++){    
    if( (PSemissions_get_Death(emissions, i ) == 0)  ){

      //This array is important, it allows us to find for the given ActiveEmission
      //What is the corresponding position in the Total emissions List of all partons in the event
      //That way we do not have to create a list of active emissions
      //We just have a list that tells us how to find the emissions which can split!
      //If it is not on that list, it can't split
      CanSplit = CanSplit + 1;
      TranslateActiveEmissions[CanSplit] = i;
      
      if( PSemissions_get_Flavor(emissions, i) == 0){
	ThePartonPhaseSpace[CanSplit] = Pg; //HAVE TO REMEMBER TO INDEX BY CanSplit here
      }else{
	ThePartonPhaseSpace[CanSplit] = Pq; //HAVE TO REMEMBER TO INDEX BY CanSplit here
      };//IF
      TotalP = TotalP + ThePartonPhaseSpace[CanSplit];
    }//IF
  };//for

  //(CanSplit+1) is the number of partons that can split currently
  //No parton will be able to split
  if( (CanSplit+1) <= 0 ){
    *(EndEvent) = 1;
    return;
  };

  CoinFlipForEmission = gsl_rng_uniform_pos(r);
  
  i=-1;
  while(i < (CanSplit+1) ){
    //Probability for i-th dipole
    i++;
    prob = ThePartonPhaseSpace[i]/TotalP;
    //If probability is great than the coin flip, the parton is chosen
    if(CoinFlipForEmission-prob<0){
      *(TheChosen) = TranslateActiveEmissions[i];//We have chosen the ith active parton, now we find its index where it sits in the master list
      i = (emissions->size)+10;//This will force exiting the while loop
    }else{
      //else, we chew off that bit of probability from coin flip and go to the next parton.
      CoinFlipForEmission = CoinFlipForEmission-prob;
    };//ELSEIF
  };//while

  //These are the quantum numbers of the chosen parton
  AngularPos[0] = PSemissions_get_Polar(emissions, *(TheChosen));
  AngularPos[1] = PSemissions_get_Azimuth(emissions, *(TheChosen));
  CurrentFlavor = PSemissions_get_Flavor(emissions, *(TheChosen));
  Efrac = PSemissions_get_Efrac(emissions, *(TheChosen) );

  //We now compute the splitting directions
  //We note that t is defined to be:
  //t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}], RF < RI
  //so now we solve for RF to get split angle

  RND = gsl_rng_uniform_pos(r) ;

  SplitAngle = 2.0*atan(tan( *(Angle)/2  )*pow(RND,-PI/(as*TotalP)) );
  *(Angle) = SplitAngle;

  //Now we determine the energy fraction

  RND = gsl_rng_uniform_pos(r);
  Z = *(InvertPg + (int)floor(RND*(double)(NumInverts)) );

  Daughter1[2] = Z * Efrac;
  Daughter2[2] = Efrac;

  if(Z >= 0.5){
    WTA1 = 1;
    WTA2 = 0;
  }else{
    WTA1 = 0;
    WTA2 = 1;     
  }
   

  AssignFlavor[0] = 0;
  AssignFlavor[1] = CurrentFlavor;

  DetermineSplittingDirectionsNoRecoil(AngularPos, Daughter1, Daughter2, Z, SplitAngle, r);
  PartonADead = 0;
  PartonBDead = 0;
  if(Daughter1[2]>MinZ){
    PartonADead = 0;
  }else{
    PartonADead = 1;
  }
  if(Daughter2[2]>MinZ){
    PartonBDead = 0;
  }else{
    PartonBDead = 1;
  }

  Daughter1[3] = SplitAngle;
  Daughter2[3] = SplitAngle;
  Daughter1[4] = Efrac;
  Daughter2[4] = Efrac;
  PSemissions_edit(DaughterEmissions, 0, Daughter1, AssignFlavor[0], PartonADead, WTA1); 
  PSemissions_edit(DaughterEmissions, 1, Daughter2, AssignFlavor[1], PartonBDead, WTA2);

  return;

}



////////////////////////////////////////////////////////////////

void SplitDirectedPartonNoConserve(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
				   int *TheChosen, double *params, const gsl_rng *r, int *EndEvent, int Mode,
				   double *InvertPg, double *InvertPq, int NumInverts){
  int i, CurrentFlavor, FoundaZ, WTA1, WTA2,TryToSplit, CanSplit, PartonADead, PartonBDead;
  int AssignFlavor[2];
  double AngularPos[2];
  double MinZ,MaxAngle, MinQ, CF, CA, NF, Q, Rmax,
    TotalP, Efrac, CoinFlipForEmission, DT, prob, Z,
    RND,SplitAngle;
  double ThePartonPhaseSpace[(emissions->size)+1];
  int TranslateActiveEmissions[(emissions->size)+1];
  double Daughter1[5],Daughter2[5];
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);

  double Pg = IntToZmaxPgTOggSimplified(1-MinZ, params);
  double Pq = 0;
  TotalP = Pg;

  //These are the quantum numbers of the chosen parton
  AngularPos[0] = PSemissions_get_Polar(emissions, *(TheChosen));
  AngularPos[1] = PSemissions_get_Azimuth(emissions, *(TheChosen));
  CurrentFlavor = PSemissions_get_Flavor(emissions, *(TheChosen));
  Efrac = PSemissions_get_Efrac(emissions, *(TheChosen) );

  //We now compute the splitting directions
  //We note that t is defined to be:
  //t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}], RF < RI
  //so now we solve for RF to get split angle

  RND = gsl_rng_uniform_pos(r) ;

  SplitAngle = 2.0*atan(tan( *(Angle)/2  )*pow(RND,-PI/(as*TotalP)) );
  *(Angle) = SplitAngle;

  //Now we determine the energy fraction

  RND = gsl_rng_uniform_pos(r);
  Z = *(InvertPg + (int)floor(RND*(double)(NumInverts)) );

  Daughter2[2] = Z * Efrac;
  Daughter1[2] = Efrac;
  WTA1 = 1;
  WTA2 = 0;
   

  AssignFlavor[0] = 0;
  AssignFlavor[1] = 0;

  DetermineSplittingDirectionsNoRecoil(AngularPos, Daughter1, Daughter2, Z, SplitAngle, r);
  PartonADead = 0;
  PartonBDead = 0;
  if(Daughter1[2]>MinZ){
    PartonADead = 0;
  }else{
    PartonADead = 1;
  }
  if(Daughter2[2]>MinZ){
    PartonBDead = 0;
  }else{
    PartonBDead = 1;
  }

  Daughter1[3] = SplitAngle;
  Daughter2[3] = SplitAngle;
  Daughter1[4] = Efrac;
  Daughter2[4] = Efrac;
  PSemissions_edit(DaughterEmissions, 0, Daughter1, AssignFlavor[0], PartonADead, WTA1); 
  PSemissions_edit(DaughterEmissions, 1, Daughter2, AssignFlavor[1], PartonBDead, WTA2);

  return;

}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void SplitSelectedPartonNLOThresh(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
				  int *TheChosen, double *params, const gsl_rng *r, int *EndEvent, int Mode,
				  double *InvertPg, double *InvertPq, int NumInverts){
  int i, CurrentFlavor, FoundaZ, WTA1, WTA2,TryToSplit, CanSplit, PartonADead, PartonBDead;
  int AssignFlavor[2];
  double AngularPos[2];
  double MinZ,MaxAngle, MinQ, CF, CA, NF, Q, Rmax,
    TotalP, Efrac, CoinFlipForEmission, DT, prob, Z,
    RND,SplitAngle;
  double ThePartonPhaseSpace[(emissions->size)+1];
  int TranslateActiveEmissions[(emissions->size)+1];
  double Daughter1[5],Daughter2[5];
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);
  double Gamma1 = (16.618305243*CA - 4.4444444*NF)/16.0;
  double beta0 = 3.66667 * CA - 0.666667 * NF;
  //Note that even with NLO as running,
  //as(K mu)= as(mu)- log(K)*beta0 (as(mu))^2/(2*PI)+N^2LL
  //we note that we then take:
  double K = pow(Exp,-2.0*Gamma1/(beta0));
  //  if(*(params+7)==0){K=1;}
  
  //the factor of 2 is due to the 1/PI we have from the as integral used here.
  //the (2PI) from expanding the as(K mu)
  //and finally the fact that the cusp is defined with the expansion
  //Gamma[as] = (as/PI) Gamma0/4 + (as/(PI))^2 (Gamma1/16) + ....
  //We have gone ahead and lumped the 1/16 into our Gamma1
  //this is different from what is done in the cusp integral functions above, but consistent
  //K will effectively redefine what the hard scale means, our Q
  
  double Pg = (1.0)*(IntToZmaxPgTOgg(1-MinZ, params)+IntToZmaxPgTOqqbar(1-MinZ, params));
  double Pq = (IntToZmaxPqTOqg(1-MinZ, params));

  TotalP = 0;
  //First we must establish the total phase space to radiate
  CanSplit = -1;
  
  for(i=0;i< (emissions->size);i++){    
    if( (PSemissions_get_Death(emissions, i ) == 0)  ){

      //This array is important, it allows us to find for the given ActiveEmission
      //What is the corresponding position in the Total emissions List of all partons in the event
      //That way we do not have to create a list of active emissions
      //We just have a list that tells us how to find the emissions which can split!
      //If it is not on that list, it can't split
      CanSplit = CanSplit + 1;
      TranslateActiveEmissions[CanSplit] = i;
      
      if( PSemissions_get_Flavor(emissions, i) == 0){
	ThePartonPhaseSpace[CanSplit] = Pg; //HAVE TO REMEMBER TO INDEX BY CanSplit here
      }else{
	ThePartonPhaseSpace[CanSplit] = Pq; //HAVE TO REMEMBER TO INDEX BY CanSplit here
      };//IF
      TotalP = TotalP + ThePartonPhaseSpace[CanSplit];
    }//IF
  };//for

  //(CanSplit+1) is the number of partons that can split currently
  //No parton will be able to split
  if( (CanSplit+1) <= 0 ){
    *(EndEvent) = 1;
    return;
  };

  DT = - log( gsl_rng_uniform_pos(r)) /  TotalP ;
  *(t) = *(t) + DT;

  CoinFlipForEmission = gsl_rng_uniform_pos(r);
  
  i=-1;
  while(i < (CanSplit+1) ){
    //Probability for i-th dipole
    i++;
    prob = ThePartonPhaseSpace[i]/TotalP;
    //If probability is great than the coin flip, the parton is chosen
    if(CoinFlipForEmission-prob<0){
      *(TheChosen) = TranslateActiveEmissions[i];//We have chosen the ith active parton, now we find its index where it sits in the master list
      i = (emissions->size)+10;//This will force exiting the while loop
    }else{
      //else, we chew off that bit of probability from coin flip and go to the next parton.
      CoinFlipForEmission = CoinFlipForEmission-prob;
    };//ELSEIF
  };//while

  //These are the quantum numbers of the chosen parton
  AngularPos[0] = PSemissions_get_Polar(emissions, *(TheChosen));
  AngularPos[1] = PSemissions_get_Azimuth(emissions, *(TheChosen));
  CurrentFlavor = PSemissions_get_Flavor(emissions, *(TheChosen));
  Efrac = PSemissions_get_Efrac(emissions, *(TheChosen) );

  //We now compute the splitting directions
  //We note that t is defined to be:
  //t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}]/PI, RF < RI
  //so now we solve for RF to get split angle
  //as is evaluated at NLO
  //Note that using  Q/K rather than Q effectively
  //uses the two loop value for the cusp anomalous dimension
  //in the 1/(1-Z) term of the splitting function
  //This is because:
  //t(RI,RF) = - Int[as(m)/m,{m,Q*RI/K,Q*RF/K}]/PI
  //         = - Int[as(K*m)/m,{m,Q*RI,Q*RF}]/PI
  //         = - Int[as(m)/m,{m,Q*RI,Q*RF}]/PI + beta0*log(K)*Int[as^2(m)/m,{m,Q*RI,Q*RF}]/(2*PI^2)
  // And
  // log(K) = -2*Gamma1/beta0
  //  SplitAngle = ComputeSplitAngleNLO(Q/K, Rmax, *t, params);
  SplitAngle = ComputeSplitAngle(Q/K, Rmax, *t, NF, CA);

  //Now we determine the energy fraction

  if(CurrentFlavor == 0){
    //We chose a gluon! Fun!

      //Now we decide whether to split to gg or qqbar
      //We split to gluons
      //We generate a random number, then invert the integral of the PgTOgg splitting kernel to find the
      //Correct energy fraction
      RND = gsl_rng_uniform_pos(r);
      Z = *(InvertPg + (int)floor(RND*(double)(NumInverts)) );

      Daughter1[2] = Z * Efrac;
      Daughter2[2] = (1-Z) * Efrac;

      if(Z >= 0.5){
	WTA1 = 1;
	WTA2 = 0;
      }else{
	WTA1 = 0;
	WTA2 = 1;     
      }
      
	RND = gsl_rng_uniform_pos(r);
	//Now we decide whether it splits to qqbar or g
	if( RND < ( PgTOgg(Z, params) )/( PgTOgg(Z, params) + PgTOqqbar(Z, params) ) ){
	  AssignFlavor[0] = 0;
	  AssignFlavor[1] = 0;
	}else{
	  AssignFlavor[0] = gsl_rng_uniform_int(r, (int)ceil(NF) ) + 1;
	  AssignFlavor[1] = -AssignFlavor[0];
	};
  } else {
    //We chose a quark! Fun!
    FoundaZ = 0;
    RND = gsl_rng_uniform_pos(r);
    Z = *(InvertPq + (int)floor(RND*(double)(NumInverts)) );

    Daughter1[2] = Z * Efrac;//THIS IS THE GLUON
    Daughter2[2] = (1-Z) * Efrac;//This is the quark or anti-quark

    if(Z >= 0.5){
      WTA1 = 1;
      WTA2 = 0;
    }else{
      WTA1 = 0;
      WTA2 = 1;     
    }

    AssignFlavor[0] = 0;
    AssignFlavor[1] = CurrentFlavor;

  };//IF

  *(Angle) = SplitAngle;
  DetermineSplittingDirections(AngularPos, Daughter1, Daughter2, Z, SplitAngle, r);
  PartonADead = 0;
  PartonBDead = 0;
  if(Daughter1[2]>MinZ){
    PartonADead = 0;
  }else{
    PartonADead = 1;
  }
  if(Daughter2[2]>MinZ){
    PartonBDead = 0;
  }else{
    PartonBDead = 1;
  }

  Daughter1[3] = SplitAngle;
  Daughter2[3] = SplitAngle;
  Daughter1[4] = Efrac;
  Daughter2[4] = Efrac;
  PSemissions_edit(DaughterEmissions, 0, Daughter1, AssignFlavor[0], PartonADead, WTA1); 
  PSemissions_edit(DaughterEmissions, 1, Daughter2, AssignFlavor[1], PartonBDead, WTA2);

  return;

}



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void UpdateEmissionList(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis,
			double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int TheChosen){

  double Daughter1[5], Daughter2[5];
  int Daughter1Flavor, Daughter2Flavor, Daughter1Death, Daughter2Death, Daughter1Label, Daughter2Label;

  Daughter1[0] = PSemissions_get_Polar(DaughterEmissions, 0);
  Daughter1[1] = PSemissions_get_Azimuth(DaughterEmissions, 0);
  Daughter1[2] = PSemissions_get_Efrac(DaughterEmissions, 0);
  Daughter1[3] = PSemissions_get_SplitAngle(DaughterEmissions, 0);
  Daughter1[4] = PSemissions_get_ParentEfrac(DaughterEmissions, 0);
  Daughter1Flavor = PSemissions_get_Flavor(DaughterEmissions, 0);
  Daughter1Death = PSemissions_get_Death(DaughterEmissions, 0);
  Daughter1Label = PSemissions_get_Label(DaughterEmissions, 0);//This is zero or one. If one, this is the WTA axis of the CURRENT splitting

  Daughter2[0] = PSemissions_get_Polar(DaughterEmissions, 1);
  Daughter2[1] = PSemissions_get_Azimuth(DaughterEmissions, 1);
  Daughter2[2] = PSemissions_get_Efrac(DaughterEmissions, 1);
  Daughter2[3] = PSemissions_get_SplitAngle(DaughterEmissions, 1);
  Daughter2[4] = PSemissions_get_ParentEfrac(DaughterEmissions, 1);
  Daughter2Flavor = PSemissions_get_Flavor(DaughterEmissions, 1);
  Daughter2Death = PSemissions_get_Death(DaughterEmissions, 1);
  Daughter2Label = PSemissions_get_Label(DaughterEmissions, 1);//This is zero or one. If one, this is the WTA axis of the CURRENT splitting

  //If the parent is the WTA axis, then one of the daughters becomes the WTA axis.
  //Otherwise, no one is the WTA axis.
  //Daughter 1 replaces the parent, the other daughter is added to the list of emissions.
  PSemissions_edit(emissions, TheChosen, Daughter1, Daughter1Flavor, Daughter1Death, *(ParentLabel) * Daughter1Label); 
  PSemissions_append(emissions, Daughter2, Daughter2Flavor, Daughter2Death, *(ParentLabel) * Daughter2Label);

  if( (*(ParentLabel) == 1)&&(Daughter1[2]>=Daughter2[2]) ){
    CurrentWTAaxis[0] = Daughter1[0];
    CurrentWTAaxis[1] = Daughter1[1];
  }else if(  (*(ParentLabel) == 1)&&(Daughter1[2]<Daughter2[2])  ){
    CurrentWTAaxis[0] = Daughter2[0];
    CurrentWTAaxis[1] = Daughter2[1];
  }

  return;
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void DGLAPDownToAngle(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t, 
		      double *params, double RF, double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int *TheChosen,
		      double *InvertPg, double *InvertPq, int NumInverts, int *EndEvent, int Mode, const gsl_rng *r){
  double FinalT, DT, MinZ,CA,NF,CF,Q,Rmax,MinQ;
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);
  
  //We note that t is defined to be:
  //t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}], RF < RI
  FinalT = AlphasIntegralLO(Q, Rmax, RF, NF, CA);

  SplitSelectedParton(emissions, DaughterEmissions, Angle, t,
		      TheChosen, params, r, EndEvent, Mode,
		      InvertPg, InvertPq, NumInverts);

  *(ParentMomentum) = PSemissions_get_Polar(emissions, *(TheChosen) );
  *(ParentMomentum+1) = PSemissions_get_Azimuth(emissions, *(TheChosen));
  *(ParentMomentum+2) = PSemissions_get_Efrac(emissions, *(TheChosen));
  *(ParentFlavor) = PSemissions_get_Flavor(emissions, *(TheChosen));
  *(ParentLabel) = PSemissions_get_Label(emissions, *(TheChosen));//This is zero or one. If one, this is the current WTA axis

  while(*(t)<FinalT){

    //When we first enter the loop, we update the emissions list
    //replacing the parent with the daughters
    //AFTER this, we generate another splitting
    //We record the parent and structure of the splitting,
    //BUT do not add the daughters to the list of emissions UNLESS we stay in the loop
    //That way the emission list does not contain ANY splittings BELOW the scale FinalT
    UpdateEmissionList(emissions, DaughterEmissions, CurrentWTAaxis,
		       ParentMomentum, ParentFlavor, ParentLabel, *(TheChosen) );
    

    SplitSelectedParton(emissions, DaughterEmissions, Angle, t,
			TheChosen, params, r, EndEvent, Mode,
			InvertPg, InvertPq, NumInverts);
    

    //We do record the parent
    *(ParentMomentum) = PSemissions_get_Polar(emissions, *(TheChosen));
    *(ParentMomentum+1) = PSemissions_get_Azimuth(emissions, *(TheChosen));
    *(ParentMomentum+2) = PSemissions_get_Efrac(emissions, *(TheChosen));
    *(ParentFlavor) = PSemissions_get_Flavor(emissions, *(TheChosen));
    *(ParentLabel) = PSemissions_get_Label(emissions, *(TheChosen));//This is zero or one. If one, this is the current WTA axis


  }//WHILE
  //We exit the While loop once we have a splitting below the angular scale RF
  //Now we have saved all the generated emissions so far, and in particular, the 
  //ParentMomentum and ParentFlavor of the last parton that split, and the structure of the Daughters 
  //the parent split into
  //We DO NOT add the last splitting to the emissions list, though we record its structure
  //in DaughterEmissions, since this splitting went below the resolution scale
  //We also save the MC time when the splitting occured, so that we can reconstruct the angle.

  return;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void DGLAPDownToAngleNoConserve(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t, 
				double *params, double RF, double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int *TheChosen,
				double *InvertPg, double *InvertPq, int NumInverts, int *EndEvent, int Mode, const gsl_rng *r){
  double FinalT, DT, MinZ,CA,NF,CF,Q,Rmax,MinQ;
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);
  
  //We note that t is defined to be:
  //t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}], RF < RI

  SplitSelectedPartonNoConserve(emissions, DaughterEmissions, Angle, t,
				TheChosen, params, r, EndEvent, Mode,
				InvertPg, InvertPq, NumInverts);

  *(ParentMomentum) = PSemissions_get_Polar(emissions, *(TheChosen) );
  *(ParentMomentum+1) = PSemissions_get_Azimuth(emissions, *(TheChosen));
  *(ParentMomentum+2) = PSemissions_get_Efrac(emissions, *(TheChosen));
  *(ParentFlavor) = PSemissions_get_Flavor(emissions, *(TheChosen));
  *(ParentLabel) = PSemissions_get_Label(emissions, *(TheChosen));//This is zero or one. If one, this is the current WTA axis

  while(*(Angle)>RF){

    //When we first enter the loop, we update the emissions list
    //replacing the parent with the daughters
    //AFTER this, we generate another splitting
    //We record the parent and structure of the splitting,
    //BUT do not add the daughters to the list of emissions UNLESS we stay in the loop
    //That way the emission list does not contain ANY splittings BELOW the scale FinalT
    UpdateEmissionList(emissions, DaughterEmissions, CurrentWTAaxis,
		       ParentMomentum, ParentFlavor, ParentLabel, *(TheChosen) );
    

    SplitSelectedPartonNoConserve(emissions, DaughterEmissions, Angle, t,
				  TheChosen, params, r, EndEvent, Mode,
				  InvertPg, InvertPq, NumInverts);
    

    //We do record the parent
    *(ParentMomentum) = PSemissions_get_Polar(emissions, *(TheChosen));
    *(ParentMomentum+1) = PSemissions_get_Azimuth(emissions, *(TheChosen));
    *(ParentMomentum+2) = PSemissions_get_Efrac(emissions, *(TheChosen));
    *(ParentFlavor) = PSemissions_get_Flavor(emissions, *(TheChosen));
    *(ParentLabel) = PSemissions_get_Label(emissions, *(TheChosen));//This is zero or one. If one, this is the current WTA axis


  }//WHILE
  //We exit the While loop once we have a splitting below the angular scale RF
  //Now we have saved all the generated emissions so far, and in particular, the 
  //ParentMomentum and ParentFlavor of the last parton that split, and the structure of the Daughters 
  //the parent split into
  //We DO NOT add the last splitting to the emissions list, though we record its structure
  //in DaughterEmissions, since this splitting went below the resolution scale
  //We also save the MC time when the splitting occured, so that we can reconstruct the angle.

  return;
}


void DGLAPDownToAnglePrimaryOnlyNoConserve(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t, 
					   double *params, double RF, double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int *TheChosen,
					   double *InvertPg, double *InvertPq, int NumInverts, int *EndEvent, int Mode, const gsl_rng *r){
  double FinalT, DT, MinZ,CA,NF,CF,Q,Rmax,MinQ;
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);
  
  //We note that t is defined to be:
  //t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}], RF < RI

  SplitDirectedPartonNoConserve(emissions, DaughterEmissions, Angle, t,
				TheChosen, params, r, EndEvent, Mode,
				InvertPg, InvertPq, NumInverts);

  *(ParentMomentum) = PSemissions_get_Polar(emissions, *(TheChosen) );
  *(ParentMomentum+1) = PSemissions_get_Azimuth(emissions, *(TheChosen));
  *(ParentMomentum+2) = PSemissions_get_Efrac(emissions, *(TheChosen));
  *(ParentFlavor) = PSemissions_get_Flavor(emissions, *(TheChosen));
  *(ParentLabel) = PSemissions_get_Label(emissions, *(TheChosen));//This is zero or one. If one, this is the current WTA axis

  while(*(Angle)>RF){

    //When we first enter the loop, we update the emissions list
    //replacing the parent with the daughters
    //AFTER this, we generate another splitting
    //We record the parent and structure of the splitting,
    //BUT do not add the daughters to the list of emissions UNLESS we stay in the loop
    //That way the emission list does not contain ANY splittings BELOW the scale FinalT
    UpdateEmissionList(emissions, DaughterEmissions, CurrentWTAaxis,
		       ParentMomentum, ParentFlavor, ParentLabel, *(TheChosen) );
    

    SplitDirectedPartonNoConserve(emissions, DaughterEmissions, Angle, t,
				  TheChosen, params, r, EndEvent, Mode,
				  InvertPg, InvertPq, NumInverts);
    

    //We do record the parent
    *(ParentMomentum) = PSemissions_get_Polar(emissions, *(TheChosen));
    *(ParentMomentum+1) = PSemissions_get_Azimuth(emissions, *(TheChosen));
    *(ParentMomentum+2) = PSemissions_get_Efrac(emissions, *(TheChosen));
    *(ParentFlavor) = PSemissions_get_Flavor(emissions, *(TheChosen));
    *(ParentLabel) = PSemissions_get_Label(emissions, *(TheChosen));//This is zero or one. If one, this is the current WTA axis


  }//WHILE
  //We exit the While loop once we have a splitting below the angular scale RF
  //Now we have saved all the generated emissions so far, and in particular, the 
  //ParentMomentum and ParentFlavor of the last parton that split, and the structure of the Daughters 
  //the parent split into
  //We DO NOT add the last splitting to the emissions list, though we record its structure
  //in DaughterEmissions, since this splitting went below the resolution scale
  //We also save the MC time when the splitting occured, so that we can reconstruct the angle.

  return;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void DGLAPDownToAngleNLOThresh(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t, 
			       double *params, double RF, double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int *TheChosen,
			       double *InvertPg, double *InvertPq, int NumInverts, int *EndEvent, int Mode, const gsl_rng *r){
  double FinalT, DT, MinZ,CA,NF,CF,Q,Rmax,MinQ;
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);
  
  SplitSelectedPartonNLOThresh(emissions, DaughterEmissions, Angle, t,
			       TheChosen, params, r, EndEvent, Mode,
			       InvertPg, InvertPq, NumInverts);

  *(ParentMomentum) = PSemissions_get_Polar(emissions, *(TheChosen) );
  *(ParentMomentum+1) = PSemissions_get_Azimuth(emissions, *(TheChosen));
  *(ParentMomentum+2) = PSemissions_get_Efrac(emissions, *(TheChosen));
  *(ParentFlavor) = PSemissions_get_Flavor(emissions, *(TheChosen));
  *(ParentLabel) = PSemissions_get_Label(emissions, *(TheChosen));//This is zero or one. If one, this is the current WTA axis


  while(*(Angle)>RF){

    //When we first enter the loop, we update the emissions list
    //replacing the parent with the daughters
    //AFTER this, we generate another splitting
    //We record the parent and structure of the splitting,
    //BUT do not add the daughters to the list of emissions UNLESS we stay in the loop
    //That way the emission list does not contain ANY splittings BELOW the scale FinalT
    UpdateEmissionList(emissions, DaughterEmissions, CurrentWTAaxis,
		       ParentMomentum, ParentFlavor, ParentLabel, *(TheChosen) );


    SplitSelectedPartonNLOThresh(emissions, DaughterEmissions, Angle, t,
				 TheChosen, params, r, EndEvent, Mode,
				 InvertPg, InvertPq, NumInverts);
  

    //We do record the parent
    *(ParentMomentum) = PSemissions_get_Polar(emissions, *(TheChosen));
    *(ParentMomentum+1) = PSemissions_get_Azimuth(emissions, *(TheChosen));
    *(ParentMomentum+2) = PSemissions_get_Efrac(emissions, *(TheChosen));
    *(ParentFlavor) = PSemissions_get_Flavor(emissions, *(TheChosen));
    *(ParentLabel) = PSemissions_get_Label(emissions, *(TheChosen));//This is zero or one. If one, this is the current WTA axis


  }//WHILE

  //We exit the While loop once we have a splitting below the angular scale RF
  //Now we have saved all the generated emissions so far, and in particular, the 
  //ParentMomentum and ParentFlavor of the last parton that split, and the structure of the Daughters 
  //the parent split into
  //We DO NOT add the last splitting to the emissions list, though we record its structure
  //in DaughterEmissions, since this splitting went below the resolution scale
  //We also save the MC time when the splitting occured, so that we can reconstruct the angle.

  return;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


void WriteToDiskHistgramsFRAG(char filename[1000], int TotalNumEvents, int NumofRadii,
			      int NumBins, double **BINS, double *Radii, double TheMinZBook,
			      int Flavor, double *params, int LogBin){
  int i,j,k,CurrentFlavor;
  double CurrentZ,CurrentBin, MinZ,CA,NF,CF,Q,Rmax,MinQ;
  FILE *OutputFile=NULL;
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);

	OutputFile = fopen(filename,"w");
	if(OutputFile==NULL){
	  printf("NOPE. NO FILE MADE.\n");
	}else{
	  fprintf(OutputFile, "%s", "{{");
	  fprintf(OutputFile, "%d,", TotalNumEvents);
	  fprintf(OutputFile, "%d,", Flavor);
	  fprintf(OutputFile, "%f,", MinZ);
	  fprintf(OutputFile, "%f,", Q);
	  fprintf(OutputFile, "%f},", Rmax);
	  fprintf(OutputFile, "{");
	  
	  /////////////////////////////////////////////////////
	  for(j=0;j<NumofRadii;j++){
	    fprintf(OutputFile, "{%f,", Radii[j]);
      

	    /////////////////////////////////////////////////////
	    fprintf(OutputFile, "{");
	    for(k = 0; k < NumBins;k++){
	      if(LogBin==0){
		CurrentZ = ( (double)(k) )/ ( (double)(NumBins) );
	      }else{
		CurrentZ = pow(TheMinZBook,1.0 -( (double)(k) )/((double)(NumBins)) );
	      }
	      CurrentBin = BINS[j][k]/( (double)(TotalNumEvents) );
	      fprintf(OutputFile, "{%e,%e}", CurrentZ, CurrentBin);
	      
	      if( k<(NumBins-1) ){//more zbins to go!
		fprintf(OutputFile, "," );
	      }//Else we do not need the comma, list is ending
	    }//We have finished looping over energy bins, now time to move to the next radius
	    fprintf(OutputFile, "}");//This ends the Energy/Angle bins
	    /////////////////////////////////////////////////////;

	    fprintf(OutputFile, "}" ); //This ends this radius 
	    if( j<(NumofRadii-1) ){//more radii to go!
	      fprintf(OutputFile, "," );
	    }//Else we have no more radii, so we don't need comma
	  }
	  /////////////////////////////////////////////////////


	  fprintf(OutputFile, "}" );
	  fprintf(OutputFile, "}" );
	  fclose(OutputFile);
	  OutputFile==NULL;
	};

  return;
}

void WriteToDiskHistgramsSIMP(char filename[1000], int TotalNumEvents, int NumofRadii,
			      int NumBins, double *BINS, double *Radii, double TheMinZBook,
			      int Flavor, double *params, int LogBin){
  int i,j,k,CurrentFlavor;
  double CurrentZ,CurrentBin, MinZ,CA,NF,CF,Q,Rmax,MinQ;
  FILE *OutputFile=NULL;
  MinZ = *(params);
  CA = *(params+1);
  NF = *(params+2);
  CF = *(params+3);
  Q = *(params+4);
  Rmax = *(params+5);
  MinQ = *(params+6);

	OutputFile = fopen(filename,"w");
	if(OutputFile==NULL){
	  printf("NOPE. NO FILE MADE.\n");
	}else{
	  fprintf(OutputFile, "%s", "{{");
	  fprintf(OutputFile, "%d,", TotalNumEvents);
	  fprintf(OutputFile, "%d,", Flavor);
	  fprintf(OutputFile, "%f,", MinZ);
	  fprintf(OutputFile, "%f,", Q);
	  fprintf(OutputFile, "%f},", Rmax);
	  fprintf(OutputFile, "{");
	  
	  /////////////////////////////////////////////////////
	  for(j=0;j<NumofRadii;j++){
	    fprintf(OutputFile, "{%f,%f}", Radii[j], BINS[j]/( (double)(TotalNumEvents)));
      	    if( j<(NumofRadii-1) ){//more radii to go!
	      fprintf(OutputFile, "," );
	    }//Else we have no more radii, so we don't need comma
	  }
	  /////////////////////////////////////////////////////


	  fprintf(OutputFile, "}" );
	  fprintf(OutputFile, "}" );
	  fclose(OutputFile);
	  OutputFile==NULL;
	};

  return;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void ReInitialize(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, int Flavor, double *t){
  double NorthPole[4];
  *(t)=0;
  NorthPole[0]=0;
  NorthPole[1]=0;
  NorthPole[2]=1;
  NorthPole[3]=PI/2;//This is the maximum split angle of the seed parton
  CurrentWTAaxis[0] = 0;
  CurrentWTAaxis[1] = 0;

  PSemissions_free(emissions);
  PSemissions_free(DaughterEmissions);

  PSemissions_init(emissions);
  PSemissions_init(DaughterEmissions);

  PSemissions_append(emissions,NorthPole,Flavor,0,1);
  //Daughter emissions is always two partons, this is just to seed
  //the two partons
  //we could in principle put anything we want here as long as they are
  //created
  //they will be immediately overwritten by the first splitting
  PSemissions_append(DaughterEmissions,NorthPole,0,0,0);
  PSemissions_append(DaughterEmissions,NorthPole,0,0,1);

  return;
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


void ReInitializeCustom(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions,
			double *CurrentWTAaxis, int Flavor, double *t, double *NorthPole){
  *(t)=0;
  CurrentWTAaxis[0] = NorthPole[0];
  CurrentWTAaxis[1] = NorthPole[1];


  PSemissions_free(emissions);
  PSemissions_free(DaughterEmissions);

  PSemissions_init(emissions);
  PSemissions_init(DaughterEmissions);

  PSemissions_append(emissions,NorthPole,Flavor,0,1);
  PSemissions_append(DaughterEmissions,NorthPole,0,0,0);
  PSemissions_append(DaughterEmissions,NorthPole,0,0,1);

  return;
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void CleanString(char *str){
  int i;

  for(i = 0;i<33; i++) {
    if(str[i] == '\n') {
      str[i] = '\0';
      return;
    }
  }

}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void LoadShowerParams(char filename[1000], double *params, int *NumRadii, double *JetRadii){

  FILE *InitializationFile=NULL;
  char InitializeFilename[1000];
  ssize_t nread;
  size_t len = 32;
  char *CurrentLine;
  char CleanLine[32];

  CurrentLine = (char *)malloc(len * sizeof(char));
  if( CurrentLine == NULL){
    perror("Unable to allocate buffer");
    exit(1);
  }

  int k;

  InitializationFile = fopen(filename,"r");
  printf("Accessing Initialization File.\n");
  if(InitializationFile==NULL){
    printf("NOPE. NO INITIALIZATION FILE.\n");
    exit(1);
  }else{
    /*
      We note the vales in params:
      *(params) = Minimum Energy Fraction or Zcut;
      *(params+1) = CA;
      *(params+2) = NF;
      *(params+3) = CF;
      *(params+4) = Jet Energy;
      *(params+5) = Jet Initial Opening Angle;
      *(params+6) = Minimum Mass Scale;
      */


    //CUTOFF
    nread = getline(&CurrentLine, &len, InitializationFile);
    if((nread == -1)||(nread == 0)){printf("BAD INITIALIZATION FILE.\n");exit(1);}
    CleanString(CurrentLine);
    *(params) = atof(CurrentLine);
    //CA
    nread = getline(&CurrentLine, &len, InitializationFile);
    if((nread == -1)||(nread == 0)){printf("BAD INITIALIZATION FILE.\n");exit(1);}
    CleanString(CurrentLine);
    *(params+1) = atof(CurrentLine);
    //NF
    nread = getline(&CurrentLine, &len, InitializationFile);
    if((nread == -1)||(nread == 0)){printf("BAD INITIALIZATION FILE.\n");exit(1);}
    CleanString(CurrentLine);
    *(params+2) = atof(CurrentLine);
    //CF
    nread = getline(&CurrentLine, &len, InitializationFile);
    if((nread == -1)||(nread == 0)){printf("BAD INITIALIZATION FILE.\n");exit(1);}
    CleanString(CurrentLine);
    *(params+3) = atof(CurrentLine);
    //MinQ0
    nread = getline(&CurrentLine, &len, InitializationFile);
    if((nread == -1)||(nread == 0)){printf("BAD INITIALIZATION FILE.\n");exit(1);}
    CleanString(CurrentLine);
    *(params+6) = atof(CurrentLine);      
    //Number of Jet Radii
    nread = getline(&CurrentLine, &len, InitializationFile);
    if((nread == -1)||(nread == 0)){printf("BAD INITIALIZATION FILE.\n");exit(1);}
    CleanString(CurrentLine);
    *(NumRadii) = atoi(CurrentLine); 
				       
    //We now read in all the jet radii that the evolution will stop at and book
    k = 0;
    if(*(NumRadii)>99){printf("TOO MANY RADII.\n");exit(1);}
    printf("Fragmentation spectra With %d Radii\n",*(NumRadii));
  
    while(k < *(NumRadii)){
      nread = getline(&CurrentLine, &len, InitializationFile);
      if((nread == -1)||(nread == 0)){printf("BAD INITIALIZATION FILE.\n");exit(1);}
      CleanString(CurrentLine);
      *(JetRadii+k) = atof(CurrentLine);	
      printf("Radius: %f\n",*(JetRadii+k));
      k++;
    };   
	     
  	

	fclose(InitializationFile);
	InitializationFile==NULL;
      
  }
  return;
}
