#ifndef Parton_Shower_Lib
#define Parton_Shower_Lib
#define as 0.1187
#define MZ 91.87
#define PI 3.14159265359
#define DIPOLE_INITIAL_CAPACITY 200000

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

typedef struct {
  int size;      // slots used so far
  int capacity;  // total available slots
  int *emissions_Label;
  int *Flavor;
  int *Dead;
  double *Polar;     
  double *Azimuth;     
  double *Efrac;     
  double *ParentSplitAngle;
  double *ParentEfrac;
} PSEmissionsList;


extern void PSemissions_init(PSEmissionsList *emissions);
extern void PSemissions_append(PSEmissionsList *emissions, double *value, int PartonFlavor, int PartonDead, int Label);
extern void PSemissions_edit(PSEmissionsList *emissions, int index, double *value, int PartonFlavor, int PartonDead, int Label);
extern void RetrieveParton(PSEmissionsList *emissions, int TheChosen, double *Momentum, int *Flavor, int *Death, int *Label);
extern void PSemissions_free(PSEmissionsList *emissions);
extern int  WriteNewPSEmission(PSEmissionsList *emissions, double *Momentum, int PartonFlavor,
				   int PartonDead, int *EmissionsLabel);

////////////////////////////////////////////////////////////////


//We solve for RF in the equation t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}]/PI, RF < RI
extern double ComputeSplitAngle(double Q, double RI, double t, double  NF, double CA);

//We solve for KT in the equation t(RI,RF) = - Int[as(m)/m,{m,Q,KT}]/PI, RF < RI
extern double ComputeKTRunningCoupling(double Q, double RI, double t, double  NF, double CA);

extern double Alphas(double MU, double NF, double TheCA);
extern double AlphasNLO(double MU, double *params);
extern double AlphasIntegral(double MUF, double MUI, double *params);
extern double AlphasSquIntegral(double MUF, double MUI, double *params);
extern double AlphasLogIntegral(double MUF, double MUI, double *params);
extern double AlphasSquLogIntegral(double MUF, double MUI, double *params);
extern double CuspIntegralLL(double Q, double MUF, double MUI, double *params);
extern double CuspIntegralNLL(double Q, double MUF, double MUI, double *params);
extern double NoLogCuspIntegralNLL(double MUF, double MUI, double *params);

extern double LambertW(const double z);

////////////////////////////////////////////////////////////////

extern double ComputeVirtuality(double Efrac, double Q, double Z, double SplitAngle);
extern double ComputeKT(double Efrac, double Q, double Z, double SplitAngle);
extern double MinZFunction(double Efrac, double Q, double PhysSplitAngle, double MinQ, double p, double q );

////////////////////////////////////////////////////////////////

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

extern void SplitSelectedParton(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
				int *TheChosen, double *params, const gsl_rng *r, int *EndEvent, int Mode,
				double *InvertPg, double *InvertPq, int NumInverts);

extern void SplitSelectedPartonNoConserve(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
					  int *TheChosen, double *params, const gsl_rng *r, int *EndEvent, int Mode,
					  double *InvertPg, double *InvertPq, int NumInverts);

extern void SplitSelectedPartonNLOThresh(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
					 int *TheChosen, double *params, const gsl_rng *r, int *EndEvent, int Mode,
					 double *InvertPg, double *InvertPq, int NumInverts);

extern void BuildSplitFunctionInversion(double *InvertPg, double *InvertPq, int NumInverts, double *params, const gsl_rng *r);

extern void DGLAPDownToAngle(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t, 
			     double *params, double RF, double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int *TheChosen,
			     double *InvertPg, double *InvertPq, int NumInverts, int *EndEvent, int Mode, const gsl_rng *r);

extern void DGLAPDownToAngleNoConserve(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t, 
				       double *params, double RF, double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int *TheChosen,
				       double *InvertPg, double *InvertPq, int NumInverts, int *EndEvent, int Mode, const gsl_rng *r);

extern void DGLAPDownToAngleNLOThresh(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t, 
				      double *params, double RF, double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int *TheChosen,
				      double *InvertPg, double *InvertPq, int NumInverts, int *EndEvent, int Mode, const gsl_rng *r);


extern void WriteToDiskHistgramsFRAG(char filename[1000], int TotalNumEvents, int NumofRadii,
				     int NumBins, double **BINS, double *Radii, double TheMinZBook,
				     int Flavor, double *params, int LogBin);

extern void WriteToDiskHistgramsSIMP(char filename[1000], int TotalNumEvents, int NumofRadii,
				     int NumBins, double *BINS, double *Radii, double TheMinZBook,
				     int Flavor, double *params, int LogBin);

extern void ReInitialize(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, int Flavor, double *t);

extern void ReInitializeCustom(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions,
			       double *CurrentWTAaxis, int Flavor, double *t, double *NorthPole);

extern void LoadShowerParams(char filename[1000], double *params, int *NumRadii, double *JetRadii);

extern void UpdateEmissionList(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis,
			       double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int TheChosen);


extern void SplitDirectedPartonNoConserve(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
					  int *TheChosen, double *params, const gsl_rng *r, int *EndEvent, int Mode,
					  double *InvertPg, double *InvertPq, int NumInverts);


void DGLAPDownToAnglePrimaryOnlyNoConserve(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t, 
					   double *params, double RF, double *ParentMomentum, int *ParentFlavor, int *ParentLabel, int *TheChosen,
					   double *InvertPg, double *InvertPq, int NumInverts, int *EndEvent, int Mode, const gsl_rng *r);

////////////////////////////////////////////////////////////////

#endif
