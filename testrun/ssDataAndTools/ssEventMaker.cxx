#include <stdio.h>
#include <stdlib.h>
#include "Adu5Pat.h"
#include "TFile.h"
#include "TTree.h"

int main()
{
  
  TFile *fGpsFile=0;
  TTree *fAdu5aPatTree=0; 
  Adu5Pat *fPat=0;
  Adu5Pat *fAdu5APatPtr=0;

  int run=262;
  char filename[150];
  sprintf(filename,"/u/home/agoodhue/ssEventMaker/state_%d.txt",run);
  FILE *fp1 = fopen(filename,"r");
  int ncols=0;
  double latTemp, lonTemp, altTemp, headingTemp, realtimeTemp, runTemp, elevAngleTemp, ssTemp, interpTemp,
    elevAngleRealTemp;
  double latArray[20000];
  double lonArray[20000];
  double altArray[20000];
  double headingArray[20000];
  double realTimeArray[20000];
  for (int i=0;i<20000;i++){
    latArray[i]=0;
    lonArray[i]=0;
    altArray[i]=0;
    headingArray[i]=0;
    realTimeArray[i]=0;
  }

  int nEntries=0;
  //get ss info from file
  for (int i=0;i<20000;i++)
    {
      ncols=fscanf(fp1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
		   &runTemp,&realtimeTemp,&lonTemp,&latTemp,&altTemp,&headingTemp,&elevAngleTemp,
		   &ssTemp, &interpTemp,&elevAngleRealTemp);
      latArray[i]=latTemp;
      lonArray[i]=lonTemp;
      altArray[i]=altTemp;
      headingArray[i]=headingTemp;
      realTimeArray[i]=realtimeTemp;
      if (realtimeTemp==realTimeArray[i-1]) break;
      //printf("i: %d realtime: %lf \n",i,realtimeTemp);
      if (realTimeArray[i]!=0) nEntries++;
    }

  printf("nentries in SS file: %d \n",nEntries);

  //create outputfile
  sprintf(filename,"/u/home/agoodhue/ssEventMaker/gpsSsFile%d.root",run);
  TFile *rootfile = new TFile(filename,"RECREATE");
  TTree *adu5PatTreeNew=new TTree("adu5PatTree","adu5PatTree");
  adu5PatTreeNew->Branch("pat","Adu5Pat",&fPat);

  //get gps file for this run
  char gpsName[FILENAME_MAX];
  sprintf(gpsName,"/u/osgstorage/anita/data/flight0809/uclRoot/run%d/gpsFile%d.root",run, run);
  fGpsFile = TFile::Open(gpsName);
  if(!fGpsFile) {
    printf("Couldn't open gps file\n");
  }
  fAdu5aPatTree = (TTree*) fGpsFile->Get("adu5PatTree");
  if(!fAdu5aPatTree) {
    printf("Couldn't get adu5aPatTree\n");
  }
  else {
    fAdu5aPatTree->SetBranchAddress("pat",&fAdu5APatPtr);
  }
  fAdu5aPatTree->BuildIndex("realTime");
  int nentriesgps=fAdu5aPatTree->GetEntries();
  printf("nentries in gps tree: %d \n",nentriesgps);

  int ctr=0;

  //go through each second
  for (int i=0;i<nEntries;i++){
    if (fPat) delete fPat;
    fPat = new Adu5Pat();
    
    int patEntry = fAdu5aPatTree->GetEntryNumberWithBestIndex(int(realTimeArray[i]));
    fAdu5aPatTree->GetEntry(patEntry);

    if ((fAdu5APatPtr->realTime - realTimeArray[i]) < 30 &&
	(fAdu5APatPtr->realTime - realTimeArray[i]) > -30 ){//use gps data
      fPat->run=fAdu5APatPtr->run;
      fPat->realTime=fAdu5APatPtr->realTime;
      fPat->latitude=fAdu5APatPtr->latitude;
      fPat->longitude=fAdu5APatPtr->longitude;
      fPat->altitude=fAdu5APatPtr->altitude;
      fPat->heading=fAdu5APatPtr->heading;
      fPat->pitch=fAdu5APatPtr->pitch;
      fPat->roll=fAdu5APatPtr->roll;
      fPat->mrms=fAdu5APatPtr->mrms;
      fPat->brms=fAdu5APatPtr->brms;
      fPat->attFlag=fAdu5APatPtr->attFlag;
      fPat->intFlag=fAdu5APatPtr->intFlag;   
      fPat->readTime=fAdu5APatPtr->readTime;
      fPat->payloadTime=fAdu5APatPtr->payloadTime;
      fPat->payloadTimeUs=fAdu5APatPtr->payloadTimeUs;
      fPat->timeOfDay=fAdu5APatPtr->timeOfDay;
    }
    
    else { //use sunsensor data
      //printf("using ss data: ctr %d \n",ctr);
      fPat->run=run;
      fPat->realTime=int(realTimeArray[i]);//arb
      fPat->latitude=latArray[i];
      fPat->longitude=lonArray[i];
      fPat->altitude=altArray[i];
      fPat->heading=headingArray[i];//arb?? can we get a heading?
      fPat->pitch=0;//arb
      fPat->roll=0;//arb
      fPat->mrms=0;//arb
      fPat->brms=0;//arb
      fPat->attFlag=0;//arb
      fPat->intFlag=0;//arb
      fPat->readTime=0;
      fPat->payloadTime=0;
      fPat->payloadTimeUs=0;
      fPat->timeOfDay=0;
      ctr++;
    }

    adu5PatTreeNew->Fill();
    
  }
  
  adu5PatTreeNew->AutoSave();
  rootfile->Close();
}
