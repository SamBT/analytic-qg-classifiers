#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "Parton_Shower_Lib.h"

#define NUMEVENTSNoteWorking 1000
#define NUMEVENTSSave 50000
#define InvertQgluon -2.0
#define InvertQquark -0.5
#define SetNumInverts 10000000
#define NumTries 10000
#define TheNumEbins 800
#define MinZBook 0.00001


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void JetFragmentation(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *JetRadii, 
		      int NumRadii, double *params, int NumEbins, 
		      double **LeadingJetSpectra, double **EventWideSpectra, double **EventWideSpectraLogBin, double **EEC,
		      double **KTconserve, double **BroadeningSpectra, double **ThrustSpectra, int NumEvent, double *t, double *InvertPg, double *InvertPq, int NumInverts,  
		      int Mode, const gsl_rng *r){

  int i,j,k,Ebin,LogEbin,TheChosen,CurrentFlavor,CurrentLabel,Death,EndEvent, CurrentMaxEbin,MadeNewEmission,CurrentFlavor2,Death2,CurrentLabel2;
  double NumEbinsFloat = (double)(NumEbins);
  double Angle,RI,RF,FindLogEbinStart, MinZ =*(params);
  double ParentMomentum[5],Momentum[5],Momentum2[5];
  int ParentFlavor,ParentLabel;
  double JetMomentum[3], Broadening,Thrust,AngleIJ,Xk,Xj,Yk,Yj,Zk,Zj;
  /*
    Momentum[4]={Polar,Azimuth,Efrac,SplitAngle}
   */
  
  *(t) = 0;
  i = 0;
  j = 0;
  MadeNewEmission = 0;
  //Before we start evolution, we log the histograms of the zero-th
  //order emission
  if((emissions->size)!=1){printf("What, more than one emission to start?\n");}
  Angle = *(params+5);

  while(i<NumRadii){
    RF = *(JetRadii+i);

    if(Angle>RF){
      DGLAPDownToAngle(emissions, DaughterEmissions, CurrentWTAaxis, &Angle, t, 
		       params, RF, ParentMomentum, &ParentFlavor, &ParentLabel, &TheChosen,
		       InvertPg, InvertPq, NumInverts, &EndEvent, Mode, r);
      MadeNewEmission = 1;
    }
    //Now we update the event wide histograms
    CurrentMaxEbin = 0;
    JetMomentum[0]=0;
    JetMomentum[1]=0;
    JetMomentum[2]=0;
    Broadening = 0;
    Thrust = 0;
    for(j=0; j < (emissions->size); j++){
      RetrieveParton(emissions, j, &Momentum, &CurrentFlavor, &Death, &CurrentLabel);
      Xj=sin(Momentum[0])*cos(Momentum[1]);
      Yj=sin(Momentum[0])*sin(Momentum[1]);
      Zj=cos(Momentum[0]);
      JetMomentum[0] = JetMomentum[0] + Momentum[2]*Xj;
      JetMomentum[1] = JetMomentum[1] + Momentum[2]*Yj;
      JetMomentum[2] = JetMomentum[2] + Momentum[2]*Zj;
      Broadening = Broadening + Momentum[2]*sqrt(sin(Momentum[0])*sin(Momentum[0]));
      for(k=j; k < (emissions->size); k++){
	RetrieveParton(emissions, k, &Momentum2, &CurrentFlavor2, &Death2, &CurrentLabel2);
	Xk=sin(Momentum2[0])*cos(Momentum2[1]);
	Yk=sin(Momentum2[0])*sin(Momentum2[1]);
	Zk=cos(Momentum2[0]);
	AngleIJ = 0.5-(Xk*Xj+Yk*Yj+Zk*Zj)/2.0;
	Ebin = (int)floor(NumEbinsFloat * AngleIJ);
	EEC[i][Ebin] = EEC[i][Ebin] + 2.0*Momentum[2]*Momentum2[2];
      }
      
      Ebin = (int)floor(NumEbinsFloat * Momentum[2] );
      if((Ebin <1)||(Ebin >NumEbins)){Ebin=1;};
      if(Ebin>CurrentMaxEbin){
	CurrentMaxEbin = Ebin;
      }; 
      FindLogEbinStart = NumEbinsFloat*(1 - log(Momentum[2])/log(MinZBook) );
      if(FindLogEbinStart > 0){
	LogEbin = (int)floor(FindLogEbinStart);
      }else{
	LogEbin = 0;
      }	
      
      
      EventWideSpectra[i][Ebin] = EventWideSpectra[i][Ebin] + pow(NumEbinsFloat,1.0);
      EventWideSpectraLogBin[i][LogEbin] = EventWideSpectraLogBin[i][LogEbin] + pow(pow(MinZBook,1.0 -( (double)(LogEbin+1) )/NumEbinsFloat)
										    -pow(MinZBook,1.0 -( (double)(LogEbin) )/NumEbinsFloat),-1.0);


    };//for
    
    Ebin = (int)floor(NumEbinsFloat * (1.0-JetMomentum[2]) );
    if((Ebin <1)||(Ebin >NumEbins)){Ebin=1;};
    ThrustSpectra[i][Ebin] = ThrustSpectra[i][Ebin] + pow(NumEbinsFloat,1.0);

    Ebin = (int)floor(NumEbinsFloat * Broadening );
    if((Ebin <1)||(Ebin >NumEbins)){Ebin=1;};
    /*    if((1.0-JetMomentum[2])<pow(NumEbinsFloat,-1.0)){Ebin=0;};
    if(Broadening<pow(NumEbinsFloat,-1.0)){Ebin=0;};*/
    BroadeningSpectra[i][Ebin] = BroadeningSpectra[i][Ebin] + pow(NumEbinsFloat,1.0);

    Ebin = (int)floor(NumEbinsFloat *  sqrt(JetMomentum[0]*JetMomentum[0]+JetMomentum[1]*JetMomentum[1]) );
    if((Ebin <1)||(Ebin >NumEbins)){Ebin=1;};
    KTconserve[i][Ebin] = KTconserve[i][Ebin] + pow(NumEbinsFloat,1.0);

    LeadingJetSpectra[i][CurrentMaxEbin] = LeadingJetSpectra[i][CurrentMaxEbin] + pow(NumEbinsFloat,1.0);
    
    i=i+1;
    //The very last splitting was not folded into the emissions.
    //we do that now, and resume evolution.
    if(MadeNewEmission == 1){
      UpdateEmissionList(emissions, DaughterEmissions, CurrentWTAaxis,
			 &ParentMomentum, &ParentFlavor, &ParentLabel, TheChosen );
      /*      RetrieveParton(DaughterEmissions, 0, &Momentum, &CurrentFlavor, &Death, &CurrentLabel);
      RetrieveParton(DaughterEmissions, 1, &Momentum2, &CurrentFlavor2, &Death2, &CurrentLabel2);
      Xj=sin(Momentum[0])*cos(Momentum[1]);
      Yj=sin(Momentum[0])*sin(Momentum[1]);
      Zj=cos(Momentum[0]);
      Xk=sin(Momentum2[0])*cos(Momentum2[1]);
      Yk=sin(Momentum2[0])*sin(Momentum2[1]);
      Zk=cos(Momentum2[0]);
      AngleIJ = 0.5-(Xk*Xj+Yk*Yj+Zk*Zj)/2.0;
      Ebin = (int)floor(NumEbinsFloat * AngleIJ);
      EEC[i][Ebin] = EEC[i][Ebin] + 2.0*Momentum[2]*Momentum2[2];*/
      MadeNewEmission = 0;
    }

  };
  return;
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[]){

  int RANDOMSEED = atoi(argv[1]);
  int NUMEVENTS = atoi(argv[2]);
  double Q = atof(argv[3]);
  double Rmax = atof(argv[4]);
  int Flavor = atoi(argv[5]);
  int ProgramMode = atoi(argv[6]);

  int i,j;

  char filenameThrustSpectra[1000];
  char filenameEEC[1000];
  char filenameConserveSpectra[1000];
  char filenameBroadeningSpectra[1000];
  char filenameInit[1000];
  char filenameLeadingJetSpectra[1000];
  char filenameEventWideSpectra[1000];
  char filenameEventWideSpectraLogBin[1000];
  clock_t begin, end;
  double time_spent;
  double **LeadingJetSpectra;
  double **EventWideSpectra;
  double **EventWideSpectraLogBin;
  double **KTSpectra,**ThrustSpectra,**BroadeningSpectra,**EEC;
  begin = clock();
  printf("\n");
  printf("Random seed: %d.\n",RANDOMSEED);

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Initialize Random number generator
  gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);//gsl_rng_taus  
  gsl_rng_set(r,RANDOMSEED);

  //Initialize values for shower
  int NumRadii;
  double ActualJetRadii[100], ActualParams[10];
  double *params, *JetRadii;
  JetRadii = &ActualJetRadii;
  params = &ActualParams;
  
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

  
  sprintf(filenameInit,argv[7]); 
  LoadShowerParams(filenameInit, params, &NumRadii, JetRadii);

  *(params+4) = Q;
  *(params+5) = Rmax;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  sprintf(filenameLeadingJetSpectra, "LeadingJetSpectra_rseed%d_Q%f_Rmax%f_ProgramMode%d_Flavor%d_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);
  sprintf(filenameConserveSpectra, "KTconservation_rseed%d_Q%f_Rmax%f_ProgramMode%d_Flavor%d_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);
  sprintf(filenameBroadeningSpectra, "KTBroadening_rseed%d_Q%f_Rmax%f_ProgramMode%d_Flavor%d_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);
  sprintf(filenameThrustSpectra, "ThrustSpectra_rseed%d_Q%f_Rmax%f_ProgramMode%d_Flavor%d_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);
  sprintf(filenameEEC, "EEC_rseed%d_Q%f_Rmax%f_ProgramMode%d_Flavor%d_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);

  sprintf(filenameEventWideSpectra, "InclusiveSpectra_rseed%d_Q%f_Rmax%f_ProgramMode%d_Flavor%d_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);
  sprintf(filenameEventWideSpectraLogBin, "InclusiveSpectraLogBin_rseed%d_Q%f_Rmax%f_ProgramMode%d_Flavor%d_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  PSEmissionsList emissions, DaughterEmissions, EmissionsWithinJet; 
  double t,AveT,RND,Z;
  double CurrentWTAaxis[2];

  PSemissions_init(&emissions);
  PSemissions_init(&DaughterEmissions);
  PSemissions_init(&EmissionsWithinJet);

  //This wipes out all the EmissionsList's, reseeds emissions with a parton pointed at the z axis with energy fraction 1, and flavor specified
  ReInitialize(&emissions, &DaughterEmissions, &CurrentWTAaxis, Flavor, &t);

  AveT = 0;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Build Splitting Function Inversion
  printf("Computing Inversion of Splitting Functions.\n" , 1 );
  double *InvertPg,*InvertPq;
  int NumInverts = SetNumInverts;
  InvertPg = (double *) malloc( (NumInverts+1) * sizeof(double ));
  InvertPq = (double *) malloc( (NumInverts+1) * sizeof(double ));

  BuildSplitFunctionInversion(InvertPg, InvertPq, NumInverts, params, r);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Initialize Histograms

  LeadingJetSpectra = (double **) malloc( (NumRadii+1) * sizeof(double *));
  EEC = (double **) malloc( (NumRadii+1) * sizeof(double *));
  EventWideSpectra = (double **) malloc( (NumRadii+1) * sizeof(double *));
  BroadeningSpectra = (double **) malloc( (NumRadii+1) * sizeof(double *));
  ThrustSpectra = (double **) malloc( (NumRadii+1) * sizeof(double *));
  KTSpectra = (double **) malloc( (NumRadii+1) * sizeof(double *));
  EventWideSpectraLogBin = (double **) malloc( (NumRadii+1) * sizeof(double *));
  for(j=0;j<(NumRadii+1);j++){
    LeadingJetSpectra[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
    EEC[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
    EventWideSpectra[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
    BroadeningSpectra[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
    ThrustSpectra[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
    KTSpectra[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
    EventWideSpectraLogBin[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
  };
  
  for(j=0; j<(NumRadii+1); j++){
    for(i=0; i<(TheNumEbins+1); i++){
      EEC[j][i]=0;
      LeadingJetSpectra[j][i]=0;
      EventWideSpectra[j][i]=0;
      BroadeningSpectra[j][i]=0;
      ThrustSpectra[j][i]=0;
      KTSpectra[j][i]=0;
      EventWideSpectraLogBin[j][i]=0;
    };
  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  
  printf("BEGIN!\n");
  i=0;
  while(i<NUMEVENTS){
      t=0;
      //Generate An Event
      JetFragmentation(&emissions, &DaughterEmissions, &CurrentWTAaxis, JetRadii, 
		       NumRadii, params, TheNumEbins, LeadingJetSpectra, EventWideSpectra, EEC,
		       EventWideSpectraLogBin, KTSpectra, BroadeningSpectra, ThrustSpectra, i+1, &t, InvertPg, InvertPq, NumInverts, ProgramMode, r);

      //Write jet
      FILE *fout = NULL;
      fout = fopen("ps_output.txt","a");
      fprintf(fout,"BEGIN EVENT %i\n",i+1);
      for (int k = 0; k < emissions.size; k++) {
        fprintf(fout,"%i",emissions.Flavor[k]);
        fprintf(fout," %f",emissions.Efrac[k]);
        fprintf(fout," %f",emissions.Polar[k]);
        fprintf(fout," %f\n",emissions.Azimuth[k]);
      }
      fprintf(fout,"END EVENT %i\n",i+1);
      fclose(fout);
      fout==NULL;

      AveT = AveT + t;
      //This wipes out all the EmissionsList's, reseeds emissions with a parton pointed at the z axis with energy fraction 1, and flavor specified
      ReInitialize(&emissions, &DaughterEmissions, &CurrentWTAaxis, Flavor, &t);
      
      if(0==((i+1)%NUMEVENTSNoteWorking)){
	printf("Working on Event %d\n" , (i+1) );
      };
      if(0==((i+1)%NUMEVENTSSave)){
	  end = clock();
	  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
          //printf("Saving Histograms at Event %d\n",(i+1));
	  printf("Average MC time: %f\n",AveT/( (double)(i+1) ));
	  printf("Time spent so far seconds: %f\n",time_spent);
	  printf("Time spent so far minutes: %f\n",time_spent/60);
	  printf("Time spent so far hours: %f\n",time_spent/(60*60));
	  printf("Time spent per event in seconds: %f\n",time_spent/(i+1));
	  //SAVE EVENTS SO FAR
	    

	  /*WriteToDiskHistgramsFRAG(filenameLeadingJetSpectra, (i+1), NumRadii,
				   TheNumEbins, LeadingJetSpectra, JetRadii, MinZBook,
				   Flavor, params, 0);
	  WriteToDiskHistgramsFRAG(filenameEventWideSpectra, (i+1), NumRadii,
				   TheNumEbins, EventWideSpectra, JetRadii, MinZBook,
				   Flavor, params, 0);
	  WriteToDiskHistgramsFRAG(filenameEventWideSpectraLogBin, (i+1), NumRadii,
				   TheNumEbins, EventWideSpectraLogBin, JetRadii, MinZBook,
				   Flavor, params, 1);

	  WriteToDiskHistgramsFRAG(filenameConserveSpectra, (i+1), NumRadii,
				   TheNumEbins, KTSpectra, JetRadii, MinZBook,
				   Flavor, params, 0);	  
	  WriteToDiskHistgramsFRAG(filenameBroadeningSpectra, (i+1), NumRadii,
				   TheNumEbins,  BroadeningSpectra,  JetRadii, MinZBook,
				   Flavor, params, 0);
	  WriteToDiskHistgramsFRAG(filenameThrustSpectra, (i+1), NumRadii,
				   TheNumEbins, ThrustSpectra, JetRadii, MinZBook,
				   Flavor, params, 0);
	  WriteToDiskHistgramsFRAG(filenameEEC, (i+1), NumRadii,
				   TheNumEbins, EEC, JetRadii, MinZBook,
				   Flavor, params, 0);*/

      };
    i = i + 1;
    };//WHILE




  for(j=0;j<(NumRadii+1);j++){
    free(LeadingJetSpectra[j]);
    free(EventWideSpectra[j]);
    free(EventWideSpectraLogBin[j]);
    free(EEC[j]);
    free(BroadeningSpectra[j]);
    free(ThrustSpectra[j]);
    free(KTSpectra[j]);
  };
  free(LeadingJetSpectra);
  free(EventWideSpectra);
  free(EventWideSpectraLogBin);
  free(EEC);
  free(BroadeningSpectra);
  free(ThrustSpectra);
  free(KTSpectra);

  PSemissions_free(&emissions);
  PSemissions_free(&DaughterEmissions);

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("\n");
  printf("Time spent seconds: %f\n",time_spent);
  printf("Time spent minutes: %f\n",time_spent/60);
  printf("Time spent hours: %f\n",time_spent/(60*60));
  printf("Time spent per event in seconds: %f\n",time_spent/NUMEVENTS);

  printf("\n");
  return 0;
}

