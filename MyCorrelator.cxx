//////////////////////////////////////////////////////////////////////////////
/////  MyCorrelator.cxx                                                 /////
/////                                                                    /////
/////  Description:                                                      /////
/////     Class for making pretty event canvases for ANITA-II            /////
/////  Author: ACG,class structure from RJN                              /////
//////////////////////////////////////////////////////////////////////////////
//System includes
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <sys/timeb.h>
#include <time.h>
#include <vector>

//Magic Display Includes
#include "MyCorrelator.h"

//Event Reader Includes
#include "AnitaConventions.h"
#include "UsefulAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "TurfRate.h"
#include "SurfHk.h"
#include "SummedTurfRate.h"
#include "AveragedSurfHk.h"
#include "AnitaGeomTool.h"
#include "AcqdStart.h"
#include "AnitaEventCalibrator.h"
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h"
#include "CalibratedAnitaEvent.h"
#include "Antarctica.h"
#include "Vector.h"
#include "BedmapReader.h"
#include "RampdemReader.h"

//ROOT Includes
#include "TROOT.h"
#include "TRandom3.h" 
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TButton.h"
#include "TGroupButton.h"
#include <TGClient.h>
#include "TStyle.h"
#include "TPostScript.h"
#include "TTree.h"
#include "math.h"
#include "TText.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "FFTtools.h"
#include "TProfile.h"
#include "Math/Interpolator.h"
#include "TImage.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "THStack.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"

// HEALPix C++ module includes.
#include "healpix_map.h"
#include "healpix_base.h"
#include "chealpix.h"
#include "datatypes.h"
#include "healpix_data_io.h"


//#include <Healpix_cxx/cxxsupport/pointing.h>
//#include <Healpix_cxx/cxxsupport/vec3.h>


using namespace std;
class MyCorrelator;
int eventStartedFlag=0;
int eventEntryGottenFlag=0;
int readBaselineFlag=0;
int readRampDemFlag=-1;
TGraph *grEv[NUM_ANTS_WITH_NADIRS];
TGraph *grEv_backup[NUM_ANTS_WITH_NADIRS];
TGraph *grEvHoriz[NUM_ANTS_WITH_NADIRS];
TGraph *grEvUnfiltered[NUM_ANTS_WITH_NADIRS];
TGraph *grEvHorizUnfiltered[NUM_ANTS_WITH_NADIRS];
TGraph *grVertBaseline[NUM_ANTS_WITH_NADIRS];
TGraph *grHorizBaseline[NUM_ANTS_WITH_NADIRS];

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<double> >+;
#endif


TGraph *grVertCoherentBaseline=0;
TGraph *grHorizCoherentBaseline=0;
TGraph *grCor[NUM_ANTS_WITH_NADIRS][NUM_ANTS_WITH_NADIRS];
TGraph *grCoherentWaveformVert=0;
TGraph *grCoherentWaveformHoriz=0;
TGraph *grCoherentWaveformDeconVert=0;
TGraph *grCoherentWaveformDeconHoriz=0;
double sinArray[3600], cosArray[3600], tanArray[3600];
double phiAnt[NUM_ANTS_WITH_NADIRS], rAnt[NUM_ANTS_WITH_NADIRS], zAnt[NUM_ANTS_WITH_NADIRS];


vector< vector<double> > baselinearray (NUM_ANTS_WITH_NADIRS, vector<double>(130));
vector< vector<double> > baselinearray_dB (NUM_ANTS_WITH_NADIRS, vector<double>(130));
vector< vector<double> > baselinearrayHoriz (NUM_ANTS_WITH_NADIRS, vector<double>(130));
vector< vector<double> > baselinearrayHoriz_dB (NUM_ANTS_WITH_NADIRS, vector<double>(130));






vector<double> SNR_ant(NUM_ANTS_WITH_NADIRS);
vector<double> SNR_ant_triggered(NUM_ANTS_WITH_NADIRS);
vector<double> SNR_ant_closest(NUM_ANTS_WITH_NADIRS);

vector<double> PowerCut(40,0);
vector<double> CoherentAnts(9,0);

vector< vector<int> > antenna_group_holder;
vector<int> unique_phis;
int antenna_groups_size;
int groupFlag=0;

double SNR_ant_coherent;
double rmsNoiseCoherent=0;
double SNRpeak;
int passed;
double percentage;

double baselineFreq[2000]={0.};


int ant_interested=0;
int what_antenna=20;
double distance_from_source=0.;
vector< vector<double> > theFFTarray (300, vector<double>(2));


//TRandom3 random; // for generating random numbers

vector< vector<double> > theFFTarrayHoriz (300, vector<double>(2));

int whichAntennasCoherent[NUM_ANTS_WITH_NADIRS];
int thetaArrayIndex, phi1ArrayIndex, phi2ArrayIndex;
int firstEventFromList[262],lastEventFromList[262],numEventsFromList[262];
int saturatedChannels[80];
int thermalFlagUniversal=0;
int allGPTimesFromList[1000000][2];
double baseLatitudeAndLongitude[1000][2];
std::string baseNames[1000];
double baseHeights[1000];
double sigmaTheta[15]={0.365945,0.365945,0.365945,0.357472,0.336054,0.313131,0.275286,0.258882,0.250037,0.208697,0.193625,0.200588,0.213024,0.221101,0.202444};
double sigmaPhi[15]={1.333438,1.333438,1.333438,0.9549,1.0007,0.946765,0.783229,0.746181,0.754061,0.537165,0.476507,0.632205,0.543501,0.551328,0.421605};
//double sigmaTheta[15];
//double sigmaPhi[15];
int sigmaIndex;
unsigned int realTimeForGroundPulser, nsForGroundPulser;
double snrAfterFilter=0.;
double noiseBeforeFilter=0.;
 double CWheight=0.;

int gpsBadFlag=0;
int strongCWFlag=0;
int polToggle=1;//vertical pol=0
double mcLatEvent, mcLonEvent, mcAltEvent, mcWeightEvent, mcExponentEvent;
int printAndresFlag=0;
double step_size_cos_fine;
double max_theta_cos;
double min_theta_cos;
int notchFilterFlag=3;//0 for no-fill,1 for rayleigh, 2 for wiener, 3 for interpolated, 4 no-notch,5 for sine subtraction
int phase_flag=3;//0=old phase, 1 = new (random) phase, 2 = interp phase, 3= geometric, 4=simple shift to zero mean
int thermalSample=0;
int cos_reconstruction_flag=1;
int window_flag=0;//0=256 points,1=512 points,2=512 points, match 1 freq
double phase_slope = -.193;//Radians/MHz
double phase_slope_hypoth = -.193;//Radians/MHz
double sineFreq=0.;
double threshold =0.1;
TStyle* babarStyle();
TStyle *babar=babarStyle();

MyCorrelator*  MyCorrelator::fgInstance = 0;
//Leave these as global variables for now

TStyle* babarStyle(){
  TStyle *babarStyle= new TStyle("BABAR","BaBar approved plots style");
#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                          // save the global style reference
  
  gStyle = babarStyle;
#endif
  // use plain black on white colors
babarStyle->SetFrameBorderMode(0);
babarStyle->SetCanvasBorderMode(0);
babarStyle->SetPadBorderMode(0);
babarStyle->SetPadColor(0);
babarStyle->SetCanvasColor(0);
babarStyle->SetStatColor(0);
//babarStyle->SetFillColor(0);

// set the paper & margin sizes
//babarStyle->SetPaperSize(20,26);
babarStyle->SetPadTopMargin(0.03);
babarStyle->SetPadRightMargin(0.15);
babarStyle->SetPadBottomMargin(0.16);
babarStyle->SetPadLeftMargin(0.12);

// use large Times-Roman fonts
babarStyle->SetTextFont(132);
babarStyle->SetTextSize(0.08);
babarStyle->SetLabelFont(132,"x");
babarStyle->SetLabelFont(132,"y");
babarStyle->SetLabelFont(132,"z");
babarStyle->SetLabelSize(0.05,"x");
babarStyle->SetTitleSize(0.06,"x");
babarStyle->SetLabelSize(0.05,"y");
babarStyle->SetTitleSize(0.06,"y");
babarStyle->SetLabelSize(0.05,"z");
babarStyle->SetTitleSize(0.06,"z");

// use bold lines and markers
babarStyle->SetMarkerStyle(20);
 babarStyle->SetHistLineWidth(2);//1.85
babarStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

// get rid of X error bars and y error bar caps
//babarStyle->SetErrorX(0.001);

// do not display any of the standard histogram decorations
babarStyle->SetOptTitle(0);
babarStyle->SetOptStat(0);
babarStyle->SetOptFit(0000);

// put tick marks on top and RHS of plots
babarStyle->SetPadTickX(1);
babarStyle->SetPadTickY(1);

 return babarStyle;
}




MyCorrelator::MyCorrelator()
{
  //Default constructor
  fRawEventPtr=0;
  fCalEventPtr=0;
  fHeadPtr=0;
  fHkPtr=0;
  fUsefulEventPtr=0;
  fUsefulAdu5Ptr=0;
  fSurfPtr=0;
  fTurfPtr=0;
  fUseCalibratedEventFile=0;
  fUseEventFile=0;
  fAvgSurfPtr=0;
  fSumTurfPtr=0;
  fgInstance=this;
  fCurrentRun=0;
  fEventTree=0;
  tgroundPulser=0;
  myCally=0;
  fUPGeomTool=0;
  fRampdemReader=0;
  fAntarctica=0;
  fFFTTools=0;
  fSurfHkTree=0;
  fTurfRateTree=0;
  fAvgSurfHkTree=0;
  fSumTurfRateTree=0;
  fHeadFile=0;
  fEventFile=0;
  fTurfRateFile=0;
  fSumTurfRateFile=0;
  fSurfHkFile=0;
  fAvgSurfHkFile=0;
  fAcqPtr=0;
  fGpsFile=0;
  for (int i=0;i<NUM_ANTS_WITH_NADIRS;i++){
    for (int j=0;j<NUM_ANTS_WITH_NADIRS;j++){
      grCor[i][j]=0;
    }
    grEv[i]=0;
    grEv_backup[i]=0;
    grEvHoriz[i]=0;
    grEvUnfiltered[i]=0;
    grEvHorizUnfiltered[i]=0;

  }
  for (int i=0;i<1000000;i++){
    allGPTimesFromList[i][0]=0;
    allGPTimesFromList[i][1]=0;
  }
  for (int i=0;i<1000;i++){
    baseLatitudeAndLongitude[i][0]=0;
    baseLatitudeAndLongitude[i][1]=0;
  }
  
}


MyCorrelator::~MyCorrelator()
{
  //Default destructor
}


MyCorrelator::MyCorrelator(char *baseDir, int run, WaveCalType::WaveCalType_t calType)
{
  //Offline constructor
  fCalEventPtr=0;
   fUseCalibratedEventFile=0;
   fHeadFile=0;
   fEventFile=0;
   fTurfRateFile=0;
   fSumTurfRateFile=0;
   fSurfHkFile=0;
   fAvgSurfHkFile=0;
   fRawEventPtr=0;
   fHeadPtr=0;
   fHkPtr=0;
   fAvgSurfPtr=0;
   fSumTurfPtr=0;
   myCally=0;
   fUsefulEventPtr=0;
   fUsefulAdu5Ptr=0;
   fTurfRateTree=0;
   fTurfPtr=0;
   fSurfPtr=0;
   fSurfHkEntry=0;
   fSurfIDEntry=0;
   fSurfHkTree=0;
   fgInstance=this;
   fUPGeomTool=0;
   fRampdemReader=0;
   fAntarctica=0;
   fFFTTools=0;
   //  cout << "MyCorrelator::MyCorrelator(" << baseDir << " , " << run
   //   << ")" << endl;
   fCurrentRun=run;
   strncpy(fCurrentBaseDir,baseDir,179);
   fCalType=calType;
   fEventTree=0;
   tgroundPulser=0;
   fSurfHkTree=0;
   fTurfRateTree=0;
   fAvgSurfHkTree=0;
   fSumTurfRateTree=0;
   fAdu5aPatTree=0;
   fAdu5aPatEntry=0;
   fAdu5APatPtr=0;
   
   
    for (int i=0;i<NUM_ANTS_WITH_NADIRS;i++){
    for (int j=0;j<NUM_ANTS_WITH_NADIRS;j++){
      grCor[i][j]=0;
    }
    grEv[i]=0;
    grEv_backup[i]=0;
    grEvHoriz[i]=0;
    grEvUnfiltered[i]=0;
    grEvHorizUnfiltered[i]=0;

  }
  for (int i=0;i<1000000;i++){
    allGPTimesFromList[i][0]=0;
    allGPTimesFromList[i][1]=0;
  }
  for (int i=0;i<1000;i++){
    baseLatitudeAndLongitude[i][0]=0;
    baseLatitudeAndLongitude[i][1]=0;
  }



}

TStyle* MyCorrelator::RootStyle() {


  TStyle *RootStyle = new TStyle("Root-Style","The Perfect Style for Plots ;-)");

#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                          // save the global style reference

  gStyle = RootStyle;
#endif
// otherwise you need to call TROOT::SetStyle("Root-Style")

  // Paper size

  RootStyle->SetPaperSize(TStyle::kUSLetter);

  // Canvas

  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);
  RootStyle->SetCanvasDefX      (10);
  RootStyle->SetCanvasDefY      (10);

  // Pads

  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderSize  (10);
  RootStyle->SetPadBorderMode  (0);
  //  RootStyle->SetPadBottomMargin(0.13);
  RootStyle->SetPadBottomMargin(0.16);
  RootStyle->SetPadTopMargin   (0.08);
  RootStyle->SetPadLeftMargin  (0.18);
  RootStyle->SetPadRightMargin (0.05);
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // Frames

  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);
  RootStyle->SetFrameBorderMode( 0);


  // Histograms

  RootStyle->SetHistFillColor(0);
  RootStyle->SetHistFillStyle(1);
  RootStyle->SetHistLineColor(1);
  RootStyle->SetHistLineStyle(0);
  RootStyle->SetHistLineWidth(2);

  // Functions

  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(1);

  //Legends 

  RootStyle->SetStatBorderSize(1);
  RootStyle->SetStatFont      (42);
  //RootStyle->SetOptStat       (111111);
  RootStyle->SetOptStat       (0);
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatX         (0.93);
  RootStyle->SetStatY         (0.90);
  RootStyle->SetStatFontSize  (0.07);
  //  RootStyle->SetStatW         (0.2);
  //  RootStyle->SetStatH         (0.15);

  // Labels, Ticks, and Titles

  RootStyle->SetTickLength ( 0.015,"X");
  RootStyle->SetTitleSize  ( 0.05,"X");
  RootStyle->SetTitleOffset( 1.20,"X");
  RootStyle->SetTitleBorderSize(0);
  //  RootStyle->SetTitleFontSize((double)3.);
  RootStyle->SetLabelOffset( 0.015,"X");
  RootStyle->SetLabelSize  ( 0.050,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");

  RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.05,"Y");
  RootStyle->SetTitleOffset( 1.600,"Y");
  RootStyle->SetLabelOffset( 0.015,"Y");
  RootStyle->SetLabelSize  ( 0.050,"Y");
  RootStyle->SetLabelFont  ( 42   ,"Y");

  RootStyle->SetTitleFont  (42,"XY");
  RootStyle->SetTitleColor  (1);

  // Options

  RootStyle->SetOptFit     (0);

  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(0.4);

  return RootStyle;
}
////////////////////////////////////////
Double_t fitSine(Double_t *x, Double_t *par){
  
  Double_t fitval = (par[0])*sin(2.0*TMath::Pi()*par[2]*x[0] + par[1]);//-2*TMath::Pi()*par[2]*10);
    return fitval;
}
Double_t sineWave(const double *xx, double *par)
{
  const Double_t t = xx[0];
  const Double_t freq = par[0];
  const Double_t mag = par[1];
  const Double_t phase = par[2];
  return mag*TMath::Sin( 2*TMath::Pi()  * freq * t + phase);
} 
////////////////////////////////////////
int MyCorrelator::startEachEvent(int myEventNumber) //does everything except get the event waveforms
{
  if (!fEventTree) initialize();
  fEventEntry=fHeadTree->GetEntryNumberWithIndex(myEventNumber);
  if(printFlag==1) cout<<"fEventEntry is "<<fEventEntry<<" for event "<<myEventNumber<<"\n";
  int retValHeader=getHeaderEntry();
  retValHeader=0;
  int patEntry = fAdu5aPatTree->GetEntryNumberWithBestIndex(fHeadPtr->realTime);
  if(printFlag==1){ 
    cout<<"patEntry is "<<patEntry<<"\n";
    cout<<"realTime is "<<fHeadPtr->realTime<<"\n";
    cout<<"eventNumber is "<<myEventNumber<<"\n";
  }
  if (patEntry==-1) patEntry=0;
  gpsBadFlag=0;
  fAdu5aPatTree->GetEntry(patEntry);
  if(printFlag==1) cout<<"lat is "<<fAdu5APatPtr->latitude<<" lon is "<<fAdu5APatPtr->longitude<<"\n";
  if ((fAdu5APatPtr->latitude>-70 || fAdu5APatPtr->latitude<-91 || fAdu5APatPtr->longitude>180 || 
       fAdu5APatPtr->longitude<-180 || fAdu5APatPtr->heading>361 || fAdu5APatPtr->heading<-1) && patEntry!=0){
    if (printFlag==1) cout<<"bad packet: using previous, event: "<<myEventNumber<<endl;
    fAdu5aPatTree->GetEntry(patEntry-1);
  }
  if ((fAdu5APatPtr->latitude>-70 || fAdu5APatPtr->latitude<-91 || fAdu5APatPtr->longitude>180 || 
       fAdu5APatPtr->longitude<-180
       || fAdu5APatPtr->heading>361 || fAdu5APatPtr->heading<-1) && patEntry!=fAdu5aPatTree->GetEntries()-1){
    if (printFlag==1) cout<<"bad packet: using next"<<endl;
    fAdu5aPatTree->GetEntry(patEntry+1);
  }
  if (fAdu5APatPtr->latitude>-70 || fAdu5APatPtr->latitude<-91 || fAdu5APatPtr->longitude>180 || 
      fAdu5APatPtr->longitude<-180 || fAdu5APatPtr->heading>361 || fAdu5APatPtr->heading<-1){
    gpsBadFlag=1;
  }
  if ((fAdu5APatPtr->realTime-fHeadPtr->realTime<=30 && fAdu5APatPtr->realTime-fHeadPtr->realTime>=0) 
      || (fHeadPtr->realTime-fAdu5APatPtr->realTime<=30 && fHeadPtr->realTime-fAdu5APatPtr->realTime>=0)){}
  else gpsBadFlag=1;
  if (printFlag==1 && gpsBadFlag==1) cout<<"bad/old gps event. realtime: "<<fAdu5APatPtr->realTime<<
    ", realtime header: "<<fHeadPtr->realTime<<endl;
  
  if (fUsefulAdu5Ptr) delete fUsefulAdu5Ptr;
  fUsefulAdu5Ptr = new UsefulAdu5Pat(fAdu5APatPtr);
  
  if (printFlag==1){ 
    cout<<"heading: "<<fAdu5APatPtr->heading<<", lat: "
	<<fAdu5APatPtr->latitude<<", lon: "<<fAdu5APatPtr->longitude
	<<", alt: "<<fAdu5APatPtr->altitude<<endl;
  }
  if (printFlag==1) cout<<"headerTime: "<<fHeadPtr->realTime<<", gpsTime: "<<fAdu5APatPtr->realTime<<endl;
  //cout<<"time to nearest gps: "<<fabs(fAdu5APatPtr->realTime-fHeadPtr->realTime)<<endl;
  return myEventNumber;
  
}
///////////////////////////////////
void MyCorrelator::initialize()
{
  int eventTree=loadEventTree();
  eventTree=0;
  int gps=loadGpsTrees();
  gps=0;
  if (!fUPGeomTool) fUPGeomTool=AnitaGeomTool::Instance();
  cout<<"did ANITAGEROMTOOL \n";
  if (!fRampdemReader)fRampdemReader =RampdemReader::Instance();
  if (!myCally) myCally=AnitaEventCalibrator::Instance();
  cout<<"did ANITAEVENTCALIBRATOR \n";
  setupCosSinTanArray();
  getPositionsOfEachAntenna();
  readInFirstLastEventsOfEachRun();
}
//////////////////////////////
void MyCorrelator::initializeAntarctica()
{
  if (!fAntarctica){
    std::cout<<"Please wait while I create the continent...\n";
    fAntarctica= new Antarctica();
    std::cout<<"Finished constructing Antarctica.  Thank you for your patience.\n";
  }
}
//////////////////////////////////////////
void MyCorrelator::initializeBaseList()
{  
  cout<<"initializing base list..."<<endl;
  if (!fAntarctica) initializeAntarctica();
  std::string name;
  double input_lat;
  double input_lon;
  char lat_sign; //A character: this one should always be 'S'.
  char lon_sign; //A character: either 'E' or 'W'
  int misc_field;
  int ctr=0;
  std::string filename="baseLocations/all_base_locations_new.txt";
  std::ifstream base_file(filename.c_str());
  if (base_file.fail())
    {
      std::cerr<<"Error!  Could not open "
	       <<("all_base_locations_new.txt")<<" to read in base locations!\n";
    } //end if
  while (base_file >> name >> input_lat >> lat_sign >> input_lon >> lon_sign >> misc_field )
    {

      if (lat_sign == 'S')
	input_lat *= -1.;
      if (lon_sign == 'W')
	input_lon = (input_lon*-1.);//+360;
      baseLatitudeAndLongitude[ctr][0]=input_lat;
      baseLatitudeAndLongitude[ctr][1]=input_lon;
      if (ctr<328 || ctr>432) baseNames[ctr]=name;
      else baseNames[ctr]="AWS_"+name;
      if (ctr!=261)
	baseHeights[ctr]=fAntarctica->surfaceAboveGeoid(input_lon,input_lat);
      else baseHeights[ctr]=0;
      //cout<<baseNames[ctr]<<"\t"<<baseLatitudeAndLongitude[ctr][0]<<"\t"<<lat_sign<<"\t"<<baseLatitudeAndLongitude[ctr][1]
      //  <<"\t"<<lon_sign<<"\t"<<baseHeights[ctr]<<endl;
      ctr++;
    }
  
  //now do pseudo bases
  //std::string filenamePseudo="${ANITA_ANALYSIS_AGOODHUE}/baseLocations/pseudoBases.txt";
  std::string filenamePseudo="/rh5stuff/64bit/src/anita/analysis/agoodhue/baseLocations/pseudoBases.txt";
  std::ifstream base_filePseudo(filenamePseudo.c_str());
  if (base_filePseudo.fail())
    {
      std::cerr<<"Error!  Could not open "
	       <<("pseudoBases.txt")<<" to read in base locations!\n";
    } //end if
  while (base_filePseudo >> name >> input_lat >> lat_sign >> input_lon >> lon_sign >> misc_field )
    {
      
      if (lat_sign == 'S')
	input_lat *= -1.;
      if (lon_sign == 'W')
	input_lon = (input_lon*-1.);//+360;
      baseLatitudeAndLongitude[ctr][0]=input_lat;
      baseLatitudeAndLongitude[ctr][1]=input_lon;
      baseNames[ctr]=name;
      baseHeights[ctr]=fAntarctica->surfaceAboveGeoid(input_lon,input_lat);
      //cout<<baseNames[ctr]<<"\t"<<baseLatitudeAndLongitude[ctr][0]<<"\t"<<lat_sign<<"\t"<<baseLatitudeAndLongitude[ctr][1]
      //  <<"\t"<<lon_sign<<"\t"<<baseHeights[ctr]<<endl;
      ctr++;
    }

  //now there is room here to add my own hard-coded bases or traverses like stephen did
  cout<<"done initializing base list"<<endl;
}
//////////////////////////
void MyCorrelator::setupCosSinTanArray()
{
  double phi;
  for(int i=0;i<3600;i++){
    phi=(i/10.);
    cosArray[i] = cos(phi*deg2rad);
    sinArray[i] = sin(phi*deg2rad);
    tanArray[i] = tan(phi*deg2rad);
    // cout<<tanArray[i]<<", "<<i<<endl;
  }
}
/////////////////////////////
void MyCorrelator::readInFirstLastEventsOfEachRun()
{
  //FILE *fp2 = fopen("/rh5stuff/64bit/src/anita/analysis/agoodhue/firstLastEventNumbersEachRun.txt","r");
  fstream fin;
  char * pEnv;
  char path[256];
  pEnv = getenv("ANITA_ANALYSIS_AGOODHUE");
  sprintf(path,"%s/firstLastEventNumbersEachRun.txt",pEnv);  
  //fin.open("firstLastEventNumbersEachRun.txt",ios::in);
  fin.open(path,ios::in);
  if( !fin.is_open() )
    {
      cout<<"file firstLastEventNumbersEachRun.txt doesn't exist in ${ANITA_ANALYSIS_AGOODHUE}"<<endl;
      //cout<<"file firstLastEventNumbersEachRun.txt exists in ${ANITA_ANALYSIS_AGOODHUE}"<<endl;
      //exit;      
    }

  FILE *fp2 = fopen(path,"r");//"firstLastEventNumbersEachRun.txt","r");
  int run, firstEvent,lastEvent,numEvents;
  int ncols=0;
  for (int i=0;i<262;i++){

    //cout<<i<<"\t"<<" run, event "<<run<<" "<<firstEvent<<" "<<lastEvent<<" "<<numEvents<<endl;
    ncols=fscanf(fp2,"%d\t%d\t%d\t%d",&run,&firstEvent,&lastEvent,&numEvents);
    firstEventFromList[run]=firstEvent;
    lastEventFromList[run]=lastEvent;
    numEventsFromList[run]=numEvents;
  }
  
  fclose(fp2);
}
//////////////////////////
void MyCorrelator::getPositionsOfEachAntenna()
{
  //get positions of each antenna
  for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    phiAnt[ant]=fUPGeomTool->getAntPhiPositionRelToAftFore(ant)*rad2deg;
    
    rAnt[ant]=fUPGeomTool->getAntR(ant);
    zAnt[ant]=fUPGeomTool->getAntZ(ant);
    //cout<<"rAnt["<<ant<<"] is "<<rAnt[ant]<<" phi is "<<phiAnt[ant]<<" z is "<<zAnt[ant]<<"\n";
    //cout<<"phi sector is "<<fUPGeomTool->getPhiSector(ant)<<"\n";
  }
  
}
////////////////////////////////
void MyCorrelator::getGraphsThisEvent(int windowWaveformFlag, double &snrPeak, double &maxSignalPeak, int &peakAnt)
{
  

  int nantennas=NUM_ANTS_WITH_NADIRS;
  double snrPeak_max=0.;
  double snrPeak_test=0.;
  Double_t deltaTInt=1./(2.6);
  double maxPeakVal=0;
  snrPeak=0;
  double rmsNoise;
  
  int filter_flag=0;
 
  int filtered_switch=0.;
 
  for (int ant=0;ant<nantennas;ant++){
    filtered_switch=0;
   
    TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
    TGraph *gr1Horiz=fUsefulEventPtr->getGraph(ant,AnitaPol::kHorizontal);
   
    TGraph *grInterp=FFTtools::getInterpolatedGraph(gr1,deltaTInt);
    TGraph *grInterpHoriz=FFTtools::getInterpolatedGraph(gr1Horiz, deltaTInt);//interpolate
   
   
    TGraph *gr_changing = new TGraph(grInterp->GetN(),grInterp->GetX(),grInterp->GetY());
    TGraph *gr_changingHoriz = new TGraph(grInterpHoriz->GetN(),grInterpHoriz->GetX(),grInterpHoriz->GetY());

    for(int n=0;n<40;n++){
      PowerCut[n]=integrateTDPower(grInterp);
    }
    if(notchFilterFlag==5 && ant!=1){
    
     
      TGraph *gfit = fitSineWave(gr_changing,0,ant);
      
      Double_t tdpower = integrateTDPower(gfit);
      Double_t oldtdpower = 4*tdpower;
      Double_t oldtdpower2 = 2*tdpower;
      
      TGraph *gpower = new TGraph();
      gpower->SetPoint(0, 0, tdpower);
      Int_t c = 1;
      filter_flag=1;
      cout<<"original power is "<<tdpower<<"\n";
      while( abs((tdpower-oldtdpower2)/tdpower) >= threshold ){
	oldtdpower = tdpower;
	oldtdpower2 = oldtdpower;
	gfit = fitSineWave(gfit, 0,ant);
	tdpower = integrateTDPower( gfit );
	gpower->SetPoint(c, c, tdpower);
	cout <<"Power trial" <<  c++ << "\t new power,old power, ratio: " << tdpower <<"\t" << oldtdpower << "\t" << (tdpower-oldtdpower)/tdpower << "\t" << threshold <<  std::endl;
	
      }
      
      
      delete gr_changing;
      gr_changing =new TGraph(gfit->GetN(),gfit->GetX(),gfit->GetY());
      delete gfit;
      delete gpower;
      
    }

    if(notchFilterFlag==5 && ant!=1){
     
     
      TGraph *gfit = fitSineWave(gr_changingHoriz,0,ant);

      Double_t tdpower = integrateTDPower(gfit);
      Double_t oldtdpower = 4*tdpower;
      Double_t oldtdpower2 = 2*tdpower;
      
      TGraph *gpower = new TGraph();
      gpower->SetPoint(0, 0, tdpower);
      Int_t c = 1;
      filter_flag=1;
      cout<<"original power is "<<tdpower<<"\n";
      while( abs((tdpower-oldtdpower2)/tdpower) >= threshold ){
	oldtdpower = tdpower;
	oldtdpower2 = oldtdpower;
	gfit = fitSineWave(gfit, 0,ant);
	tdpower = integrateTDPower( gfit );
	gpower->SetPoint(c, c, tdpower);
	cout <<"Power trial" <<  c++ << "\t new power,old power, ratio: " << tdpower <<"\t" << oldtdpower << "\t" << (tdpower-oldtdpower)/tdpower << "\t" << threshold <<  std::endl;
	
      }
      
      
      delete gr_changingHoriz;
      gr_changingHoriz =new TGraph(gfit->GetN(),gfit->GetX(),gfit->GetY());
      delete gfit;
      delete gpower;
      
    }
    

    TGraph *grFiltered;
    TGraph *grFilteredHoriz;
   
   
    if(filter_flag==1){
      grFiltered = gr_changing;
      grFilteredHoriz=gr_changingHoriz;
    }
    else{
       grFiltered = Resizeplots(grInterp);//force plot to 256 points
       grFilteredHoriz= Resizeplots(grInterpHoriz);//force plot to 256 points
    }
    //zero mean the graphs
   
    Double_t mean=grFiltered->GetMean(2);
    Double_t *volts=grFiltered->GetY();
    Double_t *times=grFiltered->GetX();
    // cout<<"here before Horiz \n";
    Double_t meanHoriz=grFilteredHoriz->GetMean(2);
    Double_t *voltsHoriz=grFilteredHoriz->GetY();
    Double_t *timesHoriz=grFilteredHoriz->GetX();
    
    for(int i=0;i<grFiltered->GetN();i++){
      volts[i]-=mean;

    }
    for(int i=0;i<grFilteredHoriz->GetN();i++){
      voltsHoriz[i]-=meanHoriz;
    }
   
    //get peak time
    Int_t peakBin=FFTtools::getPeakBin(grFiltered);
    double tVal, dummyPeakVal;
    grFiltered->GetPoint(peakBin,tVal,dummyPeakVal);
    //cout<<"ant is "<<ant<<" tVal is "<<tVal<<" peak is "<<dummyPeakVal<<"\n";
    //get peak to peak
    double peak2peak=getPeak2Peak(grFiltered);
    double corVal=peak2peak/2.;
    
    snrPeak_test = getSNR(grFiltered,rmsNoise);
    
    //get max snr, antenna
    // if (corVal>maxPeakVal){
    if(snrPeak_test > snrPeak_max){
      maxPeakVal=corVal;
      snrPeak_max = snrPeak_test;
      snrPeak=getSNR(grFiltered,rmsNoise);
      maxSignalPeak=maxPeakVal;
      peakAnt=ant;
    }
   
    //window the waveform if need be
    if (windowWaveformFlag==1){
      for (int i=0;i<gr1->GetN();i++){
	if (times[i]>tVal-15 && times[i]<tVal+20) volts[i]=volts[i];
	else volts[i]=0;
      }
    }
   
    //delete the graphs
   
    if (grEv[ant]) delete grEv[ant];
    if (grEvHoriz[ant]) delete grEvHoriz[ant];
    if (grEvUnfiltered[ant]) delete grEvUnfiltered[ant];
    if (grEvHorizUnfiltered[ant]) delete grEvHorizUnfiltered[ant];

    grEv[ant]=new TGraph(grFiltered->GetN(),times,volts);
    grEvHoriz[ant]=new TGraph(grFilteredHoriz->GetN(),timesHoriz,voltsHoriz);
    grEvUnfiltered[ant]=new TGraph(grFiltered->GetN(),times,volts);
    grEvHorizUnfiltered[ant]=new TGraph(grFilteredHoriz->GetN(),timesHoriz,voltsHoriz);
   
    delete gr1;
    delete gr1Horiz;
    delete grFilteredHoriz;
   
    delete grFiltered;
    
     
   
  }//ant  
  noiseBeforeFilter = rmsNoise;
  // myfile.close();
  cout<<"noiseBeforeFilter is "<<noiseBeforeFilter<<"\n";
  cout<<"peakAntenna is "<<peakAnt<<"\n";
}
//////////////////////////
TGraph *MyCorrelator::Resizeplots(TGraph *grWave)
{ 
 
  TGraph *gr1 = new TGraph(grWave->GetN(),grWave->GetX(),grWave->GetY());
 
  double *oldY = gr1->GetY();
  double *oldX = gr1->GetX();
  double deltaT=oldX[1]-oldX[0];
  
  
  
  if(window_flag!=0){
    int length=512;
    int offset = 128;
    double Xarray[length];
    double Yarray[length];
    for(int j=offset;j>=0;j--){
      if(j==offset){
	Xarray[j]=oldX[0];
	Yarray[j]=oldY[0];
      }
      else
	Xarray[j]=Xarray[j+1]-deltaT;
      Yarray[j]=0.;
    }
    
    
    for(int j=0;j<length-offset;j++){
      
      if(j<gr1->GetN()){
	Xarray[j+offset]=oldX[j];
	Yarray[j+offset]=oldY[j];
	//cout<<"index is "<<j-128<<" normal: x is "<<oldX[j-128]<<" and y is "<<oldY[j-128]<<"\n";
      }
      else{
	Xarray[j+offset]=Xarray[j-1+offset]+deltaT;
	Yarray[j+offset]=0;
	//cout<<"outside of normal: x is "<<Xarray[j]<<" and y is "<<Yarray[j]<<"\n";
      }
    }//for j
    grWave=new TGraph(length,Xarray,Yarray);
    delete gr1;
  
  }
  
  if(window_flag==0){
    /* cout<<"HACK! CHANGED RESIZE TO 220, rather than 256! \n";
       int length=220;*/
    int length=256;
    double Xarray[length];
    double Yarray[length];
    for(int j=0;j<length;j++){
      
      if(j<gr1->GetN()){
	Xarray[j]=oldX[j];
	Yarray[j]=oldY[j];
	
      }
      else{
	Xarray[j]=Xarray[j-1]+deltaT;
	Yarray[j]=0;
	
      }
    }//for j
    grWave=new TGraph(length,Xarray,Yarray);
    delete gr1;
  
  }
  
 
  return grWave;
  
}

//////////////////////////
void MyCorrelator::processEventsFromAList(int drawMaps, int rfOnlyFlag, int whichMcMFlag, int whichPolarization)
{
  if (!fEventTree) initialize();

  //some variables we want
  int wantToTraceBackFlag=1;
  int ctr=0;
  struct timeb timebuffer;
  char *timeline;
  double peakThetaFinal, peakPhiFinal;
  int eventPointedFlag, eventTracedFlag;
  double sourceLat, sourceLon, sourceAlt;
  int myEventNumber;
  int quietBaseFlagIndex;
  int eventNumber;
  double anitaLatitude, anitaLongitude, anitaAltitude;
  int realTime;

  int xCorPassFlag=0;
  int ratioOfPeaksPassFlag=0;
  int elevationAnglePassFlag=0;
  int peakCrossCorrFlag=0;
  int polFractionFlag=0;
  int peakHilbertFlag=0;
  int triggerFlag=0;
  double finaltheta;
  
  ftime( &timebuffer );
  timeline = ctime( & ( timebuffer.time ) );
  eventTracedFlag=0;
  eventPointedFlag=0;
  
  ///setup output files
  char filename[150];
  sprintf(filename,"rootOutputs/outputMcM%d.root",
	  fCurrentRun);
  if (printFlag==1) cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");
  
  TNtuple *ndata = new TNtuple("ndata","stuff to plot","deltaTheta:deltaPhi:deltamcmTheta:deltamcmPhi:thetaMap:phiMap:mapSNR:peakVal:ratioFirstToSecondPeak:snrCoherent:snrPeakAnt:maxSignalPeak:distanceTD:peakHilbertCoherent");
  TNtuple *ndata2= new TNtuple("ndata2","stuff to plot 2","deltaTTD:snrPeakAfterFilter:didIFilter:triggerOrPhiMaskFlag:thetaTD:phiTD:thetaWilly:phiWilly:hwTriggerAngle:thisPhiMask:distanceMcM:secondTheta:secondPhi:strongCWFlag");
  TNtuple *ndata3= new TNtuple("ndata3","stuff to plot 3","headingOfThisEvent:nadirFlag:thirdTheta:thirdPhi:varnerFlag:varnerFlag2:pitch:roll:heading:phiMaskFlag:hwTriggerFlag:ratioOfPeaksPassFlag:elevationAnglePassFlag:xCorPassFlag");
  TNtuple *ndata4= new TNtuple("ndata4","stuff to plot 4","payloadBlastFlag:polAngleCoherent:polFractionCoherent:didIFilterAboveSatellite:didIFilterHoriz:didIFilterAboveSatelliteHoriz:meanFreqVert:meanFreqHoriz");
  //TNtuple *ndataTracing=new TNtuple("ndataTracing","stuff to plot for tracing","eventNumber:peakThetaFinal:peakPhiFinal:sourceLon:sourceLat:sourceHeight:eventTracedFlag:anitaLatitude:anitaLongitude:anitaAltitude:quietBaseFlagIndex");
  
  TTree *tdataEvent = new TTree("tdataEvent","eventInfo");
  tdataEvent->Branch("eventNumber", &eventNumber,"eventNumber/I");
  
  TTree *tdataTracing = new TTree("tdataTracing","stuff to plot for tracing");
  tdataTracing->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tdataTracing->Branch("peakThetaFinal",&peakThetaFinal,"peakThetaFinal/D");
  tdataTracing->Branch("peakPhiFinal",&peakPhiFinal,"peakPhiFinal/D");
  tdataTracing->Branch("sourceLon",&sourceLon,"sourceLon/D");
  tdataTracing->Branch("sourceLat",&sourceLat,"sourceLat/D");
  tdataTracing->Branch("sourceHeight",&sourceAlt,"sourceHeight/D");
  tdataTracing->Branch("eventTracedFlag",&eventTracedFlag,"eventTracedFlag/I");
  tdataTracing->Branch("anitaLatitude",&anitaLatitude,"anitaLatitude/D");
  tdataTracing->Branch("anitaLongitude",&anitaLongitude,"anitaLongitude/D");
  tdataTracing->Branch("anitaAltitude",&anitaAltitude,"anitaAltitude/D");
  tdataTracing->Branch("quietBaseFlagIndex",&quietBaseFlagIndex,"quietBaseFlagIndex/I");
  tdataTracing->Branch("realTime",&realTime,"realTime/I");

  //get input file of events
  int ncols=0;
  int ngpEvents=111586;
  FILE *fp2;
  if (whichMcMFlag==0){ 
    fp2=fopen("groundPulser/BHEvents.txt","r");
    ngpEvents=384993;
  }
  if (whichMcMFlag==1){ 
    fp2=fopen("groundPulser/SeaveyVEvents.txt","r");
    ngpEvents=8500;
  }
  if (whichMcMFlag==2){ 
    fp2=fopen("groundPulser/SeaveyHEvents.txt","r");
    ngpEvents=5586;
  }
  if (whichMcMFlag==3){ 
    fp2=fopen("groundPulser/Seavey45Events.txt","r");
    ngpEvents=162070;
  }
  
  for (int i=0;i<ngpEvents;i++){
    ncols=fscanf(fp2,"%d",&myEventNumber);
    eventNumber=myEventNumber;
    if (myEventNumber>=firstEventFromList[fCurrentRun] && myEventNumber<=lastEventFromList[fCurrentRun]){
      //get the event
      eventStartedFlag=startEachEvent(myEventNumber);
      anitaLatitude=fAdu5APatPtr->latitude;
      anitaLongitude=fAdu5APatPtr->longitude;
      anitaAltitude=fAdu5APatPtr->altitude;
      realTime=fHeadPtr->realTime;
      
      cout<<"Is This a Borehole Or Seavey From the Times List? : "<<isMcMBoreholeOrSeaveyFromList(myEventNumber)<<endl;
      
      if (rfOnlyFlag==1){
	if (fHeadPtr->trigType&(1<<0)){ 
	  if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	  if (isChannelSaturated(myEventNumber)==-1 && !isSyncSlip(myEventNumber)){
	    
	    ftime( &timebuffer );
	    timeline = ctime( & ( timebuffer.time ) );
	    printf("Event: %d, the time is %.19s.%hu \n",myEventNumber, timeline, 
		   timebuffer.millitm);
	    
	    eventPointedFlag=pointThisEvent(myEventNumber,drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	    tdataEvent->Fill();
	    if (printFlag==1)cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
							 <<" ...about to trace back!"<<endl;
	    
	    if (wantToTraceBackFlag==1)
	      eventTracedFlag=traceBackToContinent_Brian(myEventNumber, 
						   peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	    
	    if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	      drawEventOnContinent(myEventNumber,sourceLon,sourceLat);
	    
	    ctr++;
	    cout<<"Event Included, counter: "<<ctr<<endl;
	    // ndataTracing->Fill(myEventNumber,peakThetaFinal*rad2deg,peakPhiFinal*rad2deg,
	    //	       sourceLon,sourceLat,sourceAlt,eventTracedFlag,fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,
	    //	       quietBaseFlagIndex);
	    tdataTracing->Fill();
	    
	  }
	}
      }
      if (rfOnlyFlag==0){
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(myEventNumber)==-1 && !isSyncSlip(myEventNumber)){
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  printf("Event: %d, the time is %.19s.%hu \n",myEventNumber, timeline, 
		 timebuffer.millitm);
	  
	  eventPointedFlag=pointThisEvent(myEventNumber,drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;
	  
	  if (wantToTraceBackFlag==1)
	    eventTracedFlag=traceBackToContinent_Brian(myEventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(myEventNumber,sourceLon,sourceLat);
	  
	  ctr++;
	  if (printFlag==1) cout<<"Event Included, counter: "<<ctr<<endl;
	  //ndataTracing->Fill(myEventNumber,peakThetaFinal*rad2deg,peakPhiFinal*rad2deg,
	  //	     sourceLon,sourceLat,sourceAlt,eventTracedFlag,fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude, 
	  //	     quietBaseFlagIndex);
	  tdataTracing->Fill();
	  
	}
      }
      
    }
  }

  fclose(fp2); 
  if (writeData==1){
    rootfile->Write("tdataTracing");
  }
  
}
/////////////////////////////////
void MyCorrelator::loopOverEvents(int eventCtrStart, int eventCtrEnd, int drawMaps, int taylorFlag, int thermalFlag, int whichPolarization, std::string current_dir, int filter_number,int phase_number)
{
  cout<<"filter_num, phase num are "<<filter_number<<","<<phase_number<<"\n";
  //  initializeBaseList();
  int wantToTraceBackFlag=1;
  if (!fEventTree) initialize();
  int ctr=0; 
  //int ctr_brian;
  
  if (eventCtrEnd>fEventTree->GetEntries()) eventCtrEnd=fEventTree->GetEntries();
  if (eventCtrStart>fEventTree->GetEntries()) eventCtrStart=fEventTree->GetEntries();
  fHeadTree->GetEntry(eventCtrStart);
  
  // fEventTree->GetEvent(eventCtrStart);
  int firstEvent=(int)fHeadPtr->eventNumber;
  int lastEvent=fHeadPtr->eventNumber+eventCtrEnd-eventCtrStart;
  if (printFlag==1) cout<<"first Event Number: "<<firstEvent<<", last Event Number: "<<lastEvent
			<<", number of Events: "<<lastEvent-firstEvent<<endl;
  notchFilterFlag=filter_number;
  phase_flag = phase_number;
  if(thermalFlag>0) thermalSample=1;
  thermalSample=0;
  cout<<"notchFilterFlag is "<<notchFilterFlag<<"\n";
  cout<<"thermalSample is "<<thermalSample<<"\n";
  struct timeb timebuffer;
  char *timeline;
  double peakThetaFinal, peakPhiFinal;
  double finaltheta;
  int eventNumber;
  int eventPointedFlag, eventTracedFlag;
  double sourceLat=0;
  double sourceLon=0;
  double sourceAlt=0;
  string baseName;
  int quietBaseFlag[4]={0,0,0,0};
  int quietBaseFlagIndex=-1; 
  double anitaLatitude, anitaLongitude, anitaAltitude, anitaHeading;
  int realTime;
  int mcmflag,tdflag,shorttraceflag,rfonlyflag,calpulserflag,varnerflag,channelsaturatedflag,syncslipflag,mainrfcmflag,nadirrfcmflag,tdreflectionflag, dcoffsetflag, bigenoughpeakflag, varnerflag2, payloadblastflag;
  double thetaMap, phiMap;

  int xCorPassFlag=0;
  int ratioOfPeaksPassFlag=0;
  int elevationAnglePassFlag=0;
  int peakCrossCorrFlag=0;
  int polFractionFlag=0;
  int peakHilbertFlag=0;
  int triggerFlag=0;
 
 
 
  int phase_info_flag = phase_flag;
  
  char filename[150];//simCWdata/largesample/CW2/abby/
  cout<<"current_dir is "<<current_dir<<"\n";
  // sprintf(filename,"/data/anita/btdailey/rootOutputs/notchFilter/output%d_%d_%d.root",fCurrentRun,whichPolarization,eventCtrStart);//change this
  sprintf(filename,"%s/output%d_%d_%d.root",current_dir.c_str(),fCurrentRun,whichPolarization,eventCtrStart);//change this
  cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");
  
  TNtuple *ndata = new TNtuple("ndata","stuff to plot","deltaTheta:deltaPhi:deltamcmTheta:deltamcmPhi:thetaMap:phiMap:mapSNR:peakVal:ratioFirstToSecondPeak:snrCoherent:snrPeakAnt:maxSignalPeak:distanceTD:peakHilbertCoherent");
  TNtuple *ndata2= new TNtuple("ndata2","stuff to plot 2","deltaTTD:snrPeakAfterFilter:didIFilter:triggerOrPhiMaskFlag:thetaTD:phiTD:thetaWilly:phiWilly:hwTriggerAngle:thisPhiMask:distanceMcM:secondTheta:secondPhi:strongCWFlag");
  TNtuple *ndata3= new TNtuple("ndata3","stuff to plot 3","headingOfThisEvent:nadirFlag:thirdTheta:thirdPhi:varnerFlag:varnerFlag2:pitch:roll:heading:phiMaskFlag:hwTriggerFlag:ratioOfPeaksPassFlag:elevationAnglePassFlag:xCorPassFlag");
  TNtuple *ndata4= new TNtuple("ndata4","stuff to plot 4","payloadBlastFlag:polAngleCoherent:polFractionCoherent:didIFilterAboveSatellite:didIFilterHoriz:didIFilterAboveSatelliteHoriz:meanFreqVert:meanFreqHoriz");
  
  TTree *tdataflags=new TTree("ndataflags","flags");
  tdataflags->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tdataflags->Branch("tdflag",&tdflag,"tdflag/I");
  tdataflags->Branch("rfonlyflag",&rfonlyflag,"rfonlyflag/I");
  tdataflags->Branch("shorttraceflag",&shorttraceflag,"shorttraceflag/I");
  tdataflags->Branch("mcmflag",&mcmflag,"mcmflag/I");
  tdataflags->Branch("calpulserflag",&calpulserflag,"calpulserflag/I");
  tdataflags->Branch("varnerflag",&varnerflag,"varnerflag/I");
  tdataflags->Branch("varnerflag2",&varnerflag2,"varnerflag2/I");
  tdataflags->Branch("channelsaturatedflag",&channelsaturatedflag,"channelsaturatedflag/I");
  tdataflags->Branch("syncslipflag",&syncslipflag,"syncslipflag/I");
  tdataflags->Branch("mainrfcmflag",&mainrfcmflag,"mainrfcmflag/I");
  tdataflags->Branch("nadirrfcmflag",&nadirrfcmflag,"nadirrfcmflag/I");
  tdataflags->Branch("tdreflectionflag",&tdreflectionflag,"tdreflectionflag/I");
  tdataflags->Branch("dcoffsetflag",&dcoffsetflag,"dcoffsetflag/I");
  tdataflags->Branch("bigenoughpeakflag",&bigenoughpeakflag,"bigenoughpeakflag/I");
  tdataflags->Branch("gpsBadFlag",&gpsBadFlag,"gpsBadFlag/I");
  tdataflags->Branch("payloadblastflag",&payloadblastflag,"payloadblastflag/I");

  tdataflags->Branch("xCorPassFlag",&xCorPassFlag,"xCorPassFlag/I");
  tdataflags->Branch("ratioOfPeaksPassFlag",&ratioOfPeaksPassFlag,"ratioOfPeaksPassFlag/I");
  tdataflags->Branch("elevationAnglePassFlag",&elevationAnglePassFlag,"elevationAnglePassFlag/I");
  tdataflags->Branch("peakCrossCorrFlag",&peakCrossCorrFlag,"peakCrossCorrFlag/I");
  tdataflags->Branch("polFractionFlag",&polFractionFlag,"polFractionFlag/I");
  tdataflags->Branch("peakHilbertFlag",&peakHilbertFlag,"peakHilbertFlag/I");
  tdataflags->Branch("eventTracedFlag",&eventTracedFlag,"eventTracedFlag/I");
  tdataflags->Branch("triggerFlag",&triggerFlag,"triggerFlag/I");
  
  tdataflags->Branch("eventPointedFlag",&eventPointedFlag,"eventPointedFlag/I");
  tdataflags->Branch("phase_info_flag",&phase_info_flag,"phase_info_flag/I");
  TTree *tdataEvent = new TTree("tdataEvent","eventInfo");
  tdataEvent->Branch("eventNumber", &eventNumber,"eventNumber/I");
 

  TTree *tdataTracing = new TTree("tdataTracing","stuff to plot for tracing");
  tdataTracing->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tdataTracing->Branch("peakThetaFinal",&peakThetaFinal,"peakThetaFinal/D");
  tdataTracing->Branch("peakPhiFinal",&peakPhiFinal,"peakPhiFinal/D");
  tdataTracing->Branch("sourceLon",&sourceLon,"sourceLon/D");
  tdataTracing->Branch("sourceLat",&sourceLat,"sourceLat/D");
  tdataTracing->Branch("sourceHeight",&sourceAlt,"sourceHeight/D");
  tdataTracing->Branch("heading",&anitaHeading,"heading/D");
  tdataTracing->Branch("eventTracedFlag",&eventTracedFlag,"eventTracedFlag/I");
  tdataTracing->Branch("anitaLatitude",&anitaLatitude,"anitaLatitude/D");
  tdataTracing->Branch("anitaLongitude",&anitaLongitude,"anitaLongitude/D");
  tdataTracing->Branch("anitaAltitude",&anitaAltitude,"anitaAltitude/D");
  tdataTracing->Branch("quietBaseFlagIndex",&quietBaseFlagIndex,"quietBaseFlagIndex/I");
  tdataTracing->Branch("realTime",&realTime,"realTime/I");
  tdataTracing->Branch("thetaMap",&thetaMap,"thetaMap/D");
  tdataTracing->Branch("snrPeakAfterFilter",&snrAfterFilter,"snrPeakAfterFilter/D");
  tdataTracing->Branch("rmsNoiseCoherent",&rmsNoiseCoherent,"rmsNoiseCoherent/D");
  tdataTracing->Branch("noiseBeforeFilter",&noiseBeforeFilter,"noiseBeforeFilter/D");
  tdataTracing->Branch("phiMap",&phiMap,"phiMap/D");
  tdataTracing->Branch("eventPointedFlag",&eventPointedFlag,"eventPointedFlag/I");
  tdataTracing->Branch("xCorPassFlag",&xCorPassFlag,"xCorPassFlag/I");

  tdataTracing->Branch("ratioOfPeaksPassFlag",&ratioOfPeaksPassFlag,"ratioOfPeaksPassFlag/I");
  tdataTracing->Branch("elevationAnglePassFlag",&elevationAnglePassFlag,"elevationAnglePassFlag/I");
  tdataTracing->Branch("peakCrossCorrFlag",&peakCrossCorrFlag,"peakCrossCorrFlag/I");
  tdataTracing->Branch("polFractionFlag",&polFractionFlag,"polFractionFlag/I");
  tdataTracing->Branch("peakHilbertFlag",&peakHilbertFlag,"peakHilbertFlag/I");
  tdataTracing->Branch("eventTracedFlag",&eventTracedFlag,"eventTracedFlag/I");
  tdataTracing->Branch("triggerFlag",&triggerFlag,"triggerFlag/I");
  tdataTracing->Branch("finaltheta",&finaltheta,"finaltheta/D");
  tdataTracing->Branch("CWheight",&CWheight,"CWheight/D");
  tdataTracing->Branch("phase_info_flag",&phase_info_flag,"phase_info_flag/I");
  tdataTracing->Branch("SNR_ant",&SNR_ant);
  tdataTracing->Branch("SNR_ant_triggered",&SNR_ant_triggered);
  tdataTracing->Branch("SNR_ant_closest",&SNR_ant_closest);
  tdataTracing->Branch("SNR_ant_coherent",&SNR_ant_coherent,"SNR_ant_coherent/D");
  tdataTracing->Branch("PowerCut",&PowerCut);
  tdataTracing->Branch("CoherentAnts",&CoherentAnts);
  tdataTracing->Branch("distance_from_source",&distance_from_source,"distance_from_source/D");
  if (!fUseCalibratedEventFile && !fUseEventFile){
    tdataTracing->Branch("mcWeightEvent",&mcWeightEvent,"mcWeightEvent/D");
    tdataTracing->Branch("mcExponentEvent",&mcExponentEvent,"mcExponentEvent/D");
    tdataTracing->Branch("mcLatEvent",&mcLatEvent,"mcLatEvent/D");
    tdataTracing->Branch("mcLonEvent",&mcLonEvent,"mcLonEvent/D");
    tdataTracing->Branch("mcAltEvent",&mcAltEvent,"mcAltEvent/D");
  }

  TTree *tdataPointed = new TTree("tdataPointed","only good events that passed cuts!");
  tdataPointed->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tdataPointed->Branch("sourceLon",&sourceLon,"sourceLon/D");
  tdataPointed->Branch("sourceLat",&sourceLat,"sourceLat/D");
  tdataPointed->Branch("sourceHeight",&sourceAlt,"sourceHeight/D");
  tdataPointed->Branch("heading",&anitaHeading,"heading/D");
  tdataPointed->Branch("anitaLatitude",&anitaLatitude,"anitaLatitude/D");
  tdataPointed->Branch("anitaLongitude",&anitaLongitude,"anitaLongitude/D");
  tdataPointed->Branch("anitaAltitude",&anitaAltitude,"anitaAltitude/D");
  tdataPointed->Branch("thetaMap",&thetaMap,"thetaMap/D");
  tdataPointed->Branch("phiMap",&phiMap,"phiMap/D");
  tdataPointed->Branch("snrPeakAfterFilter",&snrAfterFilter,"snrPeakAfterFilter/D");
  tdataPointed->Branch("distance_from_source",&distance_from_source,"distance_from_source/D");
  if (!fUseCalibratedEventFile && !fUseEventFile){
    tdataPointed->Branch("mcWeightEvent",&mcWeightEvent,"mcWeightEvent/D");
    tdataPointed->Branch("mcExponentEvent",&mcExponentEvent,"mcExponentEvent/D");
    tdataPointed->Branch("mcLatEvent",&mcLatEvent,"mcLatEvent/D");
    tdataPointed->Branch("mcLonEvent",&mcLonEvent,"mcLonEvent/D");
    tdataPointed->Branch("mcAltEvent",&mcAltEvent,"mcAltEvent/D");
  }
  
 
  thermalFlagUniversal=thermalFlag;
  
  ftime( &timebuffer );
  timeline = ctime( & ( timebuffer.time ) );
  //cout<<"LOOPING OVER EVERY 100 EVENTS!!!! \n\n\n";
  for (int eventctr=eventCtrStart;eventctr<eventCtrEnd;eventctr+=1){//LOOP START
    if(notchFilterFlag!=filter_number) notchFilterFlag=filter_number;
    passed=0;
    percentage=0.;
    eventTracedFlag=0;
    eventPointedFlag=0;
    //get the event
    fHeadTree->GetEntry(eventctr);
    eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
    eventNumber=fHeadPtr->eventNumber;
    anitaLatitude=fAdu5APatPtr->latitude;
    anitaLongitude=fAdu5APatPtr->longitude;
    anitaAltitude=fAdu5APatPtr->altitude;
    realTime=fHeadPtr->realTime;
    anitaHeading=fAdu5APatPtr->heading;
 
    //set flags to -10
    mcmflag=-10;
    tdflag=-10;
    tdreflectionflag=-10;
    syncslipflag=-10;
    rfonlyflag=-10;
    calpulserflag=-10;
    mainrfcmflag=-10;
    nadirrfcmflag=-10;
    channelsaturatedflag=-10;
    dcoffsetflag=-10;
    shorttraceflag=-10;
    varnerflag=-10;
    payloadblastflag=-10;

    xCorPassFlag=0;
    ratioOfPeaksPassFlag=0;
    elevationAnglePassFlag=0;
    peakCrossCorrFlag=0;
    polFractionFlag=0;
    peakHilbertFlag=0;
    triggerFlag=0;
    
    ///////////taylor dome
    if (taylorFlag==1){
      if (isTaylor(fHeadPtr->eventNumber)==1){
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isMainRFCMOn(fHeadPtr->eventNumber) && isShortTrace(fHeadPtr->eventNumber)==0){
	 
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber,drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;
	  
	  if (wantToTraceBackFlag==1){
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  
	    //eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber,peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  }
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1){ 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  //if (eventTracedFlag==1) checkIfNearAnyBase(fHeadPtr->eventNumber, 
	  //				     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt, baseName);
	  }
	  ctr++;
	 
	  if (printFlag==2) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  //  ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,peakPhiFinal*rad2deg,
	  //	     sourceLon,sourceLat,sourceAlt,eventTracedFlag,fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,
	  //	     quietBaseFlagIndex);
	  tdataTracing->Fill();
	  thetaMap=-1.*rad2deg*peakThetaFinal;
	  phiMap=rad2deg*peakPhiFinal;
	  if ((eventTracedFlag==1 || eventTracedFlag==2)) tdataPointed->Fill();//filling for all events that traced back to continent!
	}
      }
    }
    
    ////////tdreflection
     if (taylorFlag==2){
      if (isTaylorReflection(fHeadPtr->eventNumber)==1){
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isMainRFCMOn(fHeadPtr->eventNumber) && isShortTrace(fHeadPtr->eventNumber)==0){
	  
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber,drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;
	  
	  if (wantToTraceBackFlag==1)
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  
	  ctr++;
	  if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  //  ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,peakPhiFinal*rad2deg,
	  //	     sourceLon,sourceLat,sourceAlt,eventTracedFlag,fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,
	  //	     quietBaseFlagIndex);

	  tdataTracing->Fill();
	}
      }
    }


    ///////////thermals/SW
    if (thermalFlag==1){
      if ((fHeadPtr->trigType&(1<<3)) && isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber)==0
	  && isTaylor(fHeadPtr->eventNumber)==0 && isTaylorReflection(fHeadPtr->eventNumber)==0
	  && isCalPulser(fHeadPtr->eventNumber)==0 && isMainRFCMOn(fHeadPtr->eventNumber)){ //&& fHeadPtr->triggerTimeNs>1000000){ 
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber)
	    && isShortTrace(fHeadPtr->eventNumber)==0){
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;
	  
	  if (wantToTraceBackFlag==1)
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);  
	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  
	  ctr++;
	  if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  // ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,
	  //	     peakPhiFinal*rad2deg,sourceLon,sourceLat,sourceAlt,eventTracedFlag,
	  //	     fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,
	  //	     quietBaseFlagIndex);
	  tdataTracing->Fill();
	  
	}
      }
    }
    
    if (thermalFlag==2){//all non-rf 
      if ((fHeadPtr->trigType&(1<<0))==0 && isCalPulser(fHeadPtr->eventNumber)==0 
	  && isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber)==0 && isTaylor(fHeadPtr->eventNumber)==0
	  && isTaylorReflection(fHeadPtr->eventNumber)==0 && isMainRFCMOn(fHeadPtr->eventNumber)){ 
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isShortTrace(fHeadPtr->eventNumber)==0 ){
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;
	  
	  if (wantToTraceBackFlag==1)
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);  
	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  
	  ctr++;
	  if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  //ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,
	  //	     peakPhiFinal*rad2deg,sourceLon,sourceLat,sourceAlt,eventTracedFlag,
	  //	     fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,quietBaseFlagIndex );
	  tdataTracing->Fill();
	  
	}
      }
    }
    
    if (thermalFlag==3){//upward pointing rf events 
      if ((fHeadPtr->trigType&(1<<0)) && isCalPulser(fHeadPtr->eventNumber)==0 && isTaylor(fHeadPtr->eventNumber)==0
	  && isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber)==0
	  && isTaylorReflection(fHeadPtr->eventNumber)==0 && isMainRFCMOn(fHeadPtr->eventNumber)){
	  //&& isVarnerEvent(fHeadPtr->eventNumber)==0){ 
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isShortTrace(fHeadPtr->eventNumber)==0 ){
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  cout<<"about to point to an event! \n";
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  if (peakThetaFinal<0) tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal<<endl;	  
	  if (peakThetaFinal<0){//upward pointing
	    printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		   timebuffer.millitm);
	    

	    ctr++;
	    if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	    //  ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,
	    //	       peakPhiFinal*rad2deg,0,0,0,0,
	    //	       fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude, quietBaseFlagIndex);
	    tdataTracing->Fill();
	  }
	}
      }
    }
    
    if (thermalFlag==4){// rf events near a small base with no events from ANITA-I 
      if ((fHeadPtr->trigType&(1<<0)) && isCalPulser(fHeadPtr->eventNumber)==0 && isTaylor(fHeadPtr->eventNumber)==0
	  && isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber)==0
	  && isTaylorReflection(fHeadPtr->eventNumber)==0){
	  // && isVarnerEvent(fHeadPtr->eventNumber)==0){ 
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isMainRFCMOn(fHeadPtr->eventNumber) && isShortTrace(fHeadPtr->eventNumber)==0){
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;	  
	  if (wantToTraceBackFlag==1){
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  }
	  if (printFlag==1) cout<<"event traced flag: "<<eventTracedFlag<<endl;
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  for (int i=0;i<4;i++){
	    quietBaseFlag[i]=0;
	  }
	  if (eventTracedFlag==1){
	    quietBaseFlag[0]=checkIfNearSpecificBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt,"AWS_Zoe_(Mega_B)",-1);
	    quietBaseFlag[1]=checkIfNearSpecificBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt,"AWS_AGO_4",-1);
	    quietBaseFlag[2]=checkIfNearSpecificBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt,"Kohnen_Station_(EPICA)[Germany]",-1);
	    quietBaseFlag[3]=checkIfNearSpecificBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt,"A._de_Navajo_Sobral[Argentina]",-1);  
	  }
	  for (int i=0;i<4;i++){
	    if (quietBaseFlag[i]==1) quietBaseFlagIndex=i;
	  }
	  
	  if ((quietBaseFlag[0]==1 || quietBaseFlag[1]==1 ||quietBaseFlag[2]==1 ||quietBaseFlag[3]==1) && eventTracedFlag==1){
	    ctr++;
	    printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		   timebuffer.millitm);
	    if (printFlag==2) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	    // ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,
	    //	       peakPhiFinal*rad2deg,sourceLon,sourceLat,sourceAlt,eventTracedFlag,
	    //	       fAdu5APatPtr->latitude,fAdu5APatPtr->longitude, fAdu5APatPtr->altitude, quietBaseFlagIndex);
	    tdataTracing->Fill();
	  }
	}
      }
    }
    cout<<"thermal and taylor flag is "<<thermalFlag<<" "<<taylorFlag<<"\n";
    if (thermalFlag==0 && taylorFlag==0){//all events that pass basic cuts!!!!
      if (fHeadPtr->trigType&(1<<1) || fHeadPtr->trigType&(1<<2) || fHeadPtr->trigType&(1<<3) ||
	  !fHeadPtr->trigType&(1<<0)){
	rfonlyflag=0;
	
	cout<<"eventnumber for forced trigger is "<<eventNumber<<"\n";
      }
      else rfonlyflag=1;
      
      calpulserflag=isCalPulser(fHeadPtr->eventNumber);
      tdflag=isTaylor(fHeadPtr->eventNumber);
      tdreflectionflag=isTaylorReflection(fHeadPtr->eventNumber);
      mcmflag=isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber);
      syncslipflag=isSyncSlip(fHeadPtr->eventNumber);
      
      if (rfonlyflag!=0 && calpulserflag==0 && tdflag==0
	  && mcmflag==0 && tdreflectionflag==0 && !syncslipflag && gpsBadFlag==0){

	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	
	channelsaturatedflag=isChannelSaturated(fHeadPtr->eventNumber);
	cout<<"here after saturated \n";
	mainrfcmflag=isMainRFCMOn(fHeadPtr->eventNumber);   
	bigenoughpeakflag=isBigEnoughPeakToPeak(fHeadPtr->eventNumber); 
	dcoffsetflag=isDCOffsetLarge(fHeadPtr->eventNumber);
	shorttraceflag=isShortTrace(fHeadPtr->eventNumber);
	varnerflag2=isVarnerEvent2(fHeadPtr->eventNumber);
	varnerflag=isVarnerEvent(fHeadPtr->eventNumber);
	nadirrfcmflag=isNadirRFCMOn(fHeadPtr->eventNumber);
	payloadblastflag=isPayloadBlast(fHeadPtr->eventNumber);
	cout<<"channel, mainrfcm, bigenough, dcoffset, shorttraceflag is "<<channelsaturatedflag<<","<<mainrfcmflag<<","<<bigenoughpeakflag<<","<<dcoffsetflag<<","<<shorttraceflag<<"\n";
	if (channelsaturatedflag<=3 && mainrfcmflag && shorttraceflag==0 && dcoffsetflag==0 && bigenoughpeakflag==1){//ACG used to be channelsaturatedflag==-1
	  //int pass=1;
	  //if(pass==1){
	ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  cout<<"inside the if statement. should point next \n";
	  //cout<<"event number header: "<<fHeadPtr->eventNumber<<", event number event: "<<fHeadPtr->eventNumber<<endl;
	 
	    eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);

	  
	   cout<<"eventnumber is "<<eventNumber<<" "<<passed<<" passed, "<<percentage<<"% \n";
	 
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;	  
	  if (wantToTraceBackFlag==1){
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  }//want to traceBackFlag
	  if (printFlag==1) cout<<"event traced flag: "<<eventTracedFlag<<endl;
	  if (eventNumber==2462929) cout<<"ratio, peakCross,peakHilbert,polFraction,xCor,elevation,eventTraced,trigger are "<<ratioOfPeaksPassFlag<<" "<<peakCrossCorrFlag<<" "<<peakHilbertFlag<<" "<<polFractionFlag<<" "<<xCorPassFlag<<" "<<elevationAnglePassFlag<<" "<<eventTracedFlag<<" "<<triggerFlag<<"\n";
	  if (eventTracedFlag==1 || eventTracedFlag==2) checkIfNearAnyBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt, baseName);
	  thetaMap=-1.*rad2deg*peakThetaFinal;
	  phiMap=rad2deg*peakPhiFinal;

	  if(eventNumber==2462929){
	    cout<<"thetaMap is "<<thetaMap<<" phiMap is "<<phiMap<<"\n";

	  }
	  if ((eventTracedFlag==1 || eventTracedFlag==2)){// && eventPointedFlag==1){
	    //rootfile output for clustering algorithm
	    //events that pointed and traced to continent
	    tdataPointed->Fill();
	    
	  }//eventTracedFlag	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  
	  tdataTracing->Fill();
	}//channelsaturated 
      }//all flags ==0
      
     
      tdataflags->Fill();
      
    }//thermal&taylordome==0
   
    ///////////simulated events loop
    if (thermalFlag==-1 && taylorFlag==-1){	 
      ftime( &timebuffer );
      timeline = ctime( & ( timebuffer.time ) );
      printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
	     timebuffer.millitm);
      
      eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber,drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
      tdataEvent->Fill();
      if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
			    <<" ...about to trace back!"<<endl;
      
      if (wantToTraceBackFlag==1)
	eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
					     peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
      
      if ((eventTracedFlag==1) && drawMaps==1) 
	drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
      
      ctr++;
      if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
      tdataTracing->Fill();
      thetaMap=-1.*rad2deg*peakThetaFinal;
      phiMap=rad2deg*peakPhiFinal;
      if ((eventTracedFlag==1) && eventPointedFlag==1) tdataPointed->Fill();            
    }
    
  }
  cout<<"writing data" <<endl;
  
  if (writeData==1){    
    //tdataTracing->Print();
    // rootfile->cd();
    //tdataTracing->Write();
    rootfile=tdataTracing->GetCurrentFile();
    rootfile->Write("tdataTracing");
    //rootfile->Write("ndata");
    //rootfile->Write("ndata2");
    rootfile->Close();
  }
  cout<<"exiting loopoverevents \n";
}
///////////////////////

/////////////////////////////////
void MyCorrelator::loopOverThirdEvents(int whichThird, int drawMaps, int taylorFlag, int thermalFlag, int whichPolarization, std::string current_dir)
{
 
  //  initializeBaseList();
  int wantToTraceBackFlag=1;
  if (!fEventTree) initialize();
  int ctr=0; 
  //int ctr_brian;
  int eventCtrStart;
  int eventCtrEnd;
  int thirds = (int) floor(fEventTree->GetEntries()/3.);
  
  eventCtrEnd = (int)floor(whichThird*thirds);
  eventCtrStart = eventCtrEnd - thirds;
  cout<<"eventCtrStart is "<<eventCtrStart<<" eventCtrEnd is "<<eventCtrEnd<<"\n";

  if(whichThird >1){
    eventCtrStart++;
  }
  //hack!!  
  //cout<<"HACKING! CHANGED NUMBER OF EVENTS \n";
  // eventCtrEnd = eventCtrStart+50;
  //end hack!!
  if (eventCtrEnd>fEventTree->GetEntries()) eventCtrEnd=fEventTree->GetEntries();
  if (eventCtrStart>fEventTree->GetEntries()) eventCtrStart=fEventTree->GetEntries();
  fHeadTree->GetEntry(eventCtrStart);
  
  // fEventTree->GetEvent(eventCtrStart);
  int firstEvent=(int)fHeadPtr->eventNumber;
  int lastEvent=fHeadPtr->eventNumber+eventCtrEnd-eventCtrStart;
  if (printFlag==1) cout<<"first Event Number: "<<firstEvent<<", last Event Number: "<<lastEvent
			<<", number of Events: "<<lastEvent-firstEvent<<endl;

  struct timeb timebuffer;
  char *timeline;
  double peakThetaFinal, peakPhiFinal;
  double finaltheta;
  int eventNumber;
  int eventPointedFlag, eventTracedFlag;
  double sourceLat=0;
  double sourceLon=0;
  double sourceAlt=0;
  string baseName;
  int quietBaseFlag[4]={0,0,0,0};
  int quietBaseFlagIndex=-1; 
  double anitaLatitude, anitaLongitude, anitaAltitude, anitaHeading;
  int realTime;
  int mcmflag,tdflag,shorttraceflag,rfonlyflag,calpulserflag,varnerflag,channelsaturatedflag,syncslipflag,mainrfcmflag,nadirrfcmflag,tdreflectionflag, dcoffsetflag, bigenoughpeakflag, varnerflag2, payloadblastflag;
  double thetaMap, phiMap;

  int xCorPassFlag=0;
  int ratioOfPeaksPassFlag=0;
  int elevationAnglePassFlag=0;
  int peakCrossCorrFlag=0;
  int polFractionFlag=0;
  int peakHilbertFlag=0;
  int triggerFlag=0;
  
 
  vector< vector<double> > notched_freqs(NUM_ANTS_WITH_NADIRS,vector<double> (10));
  vector< vector<double> > notched_freqsHoriz(NUM_ANTS_WITH_NADIRS,vector<double> (10));
 
 
  char filename[150];//simCWdata/largesample/CW2/abby/
  cout<<"current_dir is "<<current_dir<<"\n";
  // sprintf(filename,"/data/anita/btdailey/rootOutputs/notchFilter/output%d_%d_%d.root",fCurrentRun,whichPolarization,eventCtrStart);//change this
  sprintf(filename,"%s/output%d_%d_%d.root",current_dir.c_str(),fCurrentRun,whichPolarization,eventCtrStart);//change this
  cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");
  
  TNtuple *ndata = new TNtuple("ndata","stuff to plot","deltaTheta:deltaPhi:deltamcmTheta:deltamcmPhi:thetaMap:phiMap:mapSNR:peakVal:ratioFirstToSecondPeak:snrCoherent:snrPeakAnt:maxSignalPeak:distanceTD:peakHilbertCoherent");
  TNtuple *ndata2= new TNtuple("ndata2","stuff to plot 2","deltaTTD:snrPeakAfterFilter:didIFilter:triggerOrPhiMaskFlag:thetaTD:phiTD:thetaWilly:phiWilly:hwTriggerAngle:thisPhiMask:distanceMcM:secondTheta:secondPhi:strongCWFlag");
  TNtuple *ndata3= new TNtuple("ndata3","stuff to plot 3","headingOfThisEvent:nadirFlag:thirdTheta:thirdPhi:varnerFlag:varnerFlag2:pitch:roll:heading:phiMaskFlag:hwTriggerFlag:ratioOfPeaksPassFlag:elevationAnglePassFlag:xCorPassFlag");
  TNtuple *ndata4= new TNtuple("ndata4","stuff to plot 4","payloadBlastFlag:polAngleCoherent:polFractionCoherent:didIFilterAboveSatellite:didIFilterHoriz:didIFilterAboveSatelliteHoriz:meanFreqVert:meanFreqHoriz");
  
  TTree *tdataflags=new TTree("ndataflags","flags");
  tdataflags->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tdataflags->Branch("tdflag",&tdflag,"tdflag/I");
  tdataflags->Branch("rfonlyflag",&rfonlyflag,"rfonlyflag/I");
  tdataflags->Branch("shorttraceflag",&shorttraceflag,"shorttraceflag/I");
  tdataflags->Branch("mcmflag",&mcmflag,"mcmflag/I");
  tdataflags->Branch("calpulserflag",&calpulserflag,"calpulserflag/I");
  tdataflags->Branch("varnerflag",&varnerflag,"varnerflag/I");
  tdataflags->Branch("varnerflag2",&varnerflag2,"varnerflag2/I");
  tdataflags->Branch("channelsaturatedflag",&channelsaturatedflag,"channelsaturatedflag/I");
  tdataflags->Branch("syncslipflag",&syncslipflag,"syncslipflag/I");
  tdataflags->Branch("mainrfcmflag",&mainrfcmflag,"mainrfcmflag/I");
  tdataflags->Branch("nadirrfcmflag",&nadirrfcmflag,"nadirrfcmflag/I");
  tdataflags->Branch("tdreflectionflag",&tdreflectionflag,"tdreflectionflag/I");
  tdataflags->Branch("dcoffsetflag",&dcoffsetflag,"dcoffsetflag/I");
  tdataflags->Branch("bigenoughpeakflag",&bigenoughpeakflag,"bigenoughpeakflag/I");
  tdataflags->Branch("gpsBadFlag",&gpsBadFlag,"gpsBadFlag/I");
  tdataflags->Branch("payloadblastflag",&payloadblastflag,"payloadblastflag/I");

  tdataflags->Branch("xCorPassFlag",&xCorPassFlag,"xCorPassFlag/I");
  tdataflags->Branch("ratioOfPeaksPassFlag",&ratioOfPeaksPassFlag,"ratioOfPeaksPassFlag/I");
  tdataflags->Branch("elevationAnglePassFlag",&elevationAnglePassFlag,"elevationAnglePassFlag/I");
  tdataflags->Branch("peakCrossCorrFlag",&peakCrossCorrFlag,"peakCrossCorrFlag/I");
  tdataflags->Branch("polFractionFlag",&polFractionFlag,"polFractionFlag/I");
  tdataflags->Branch("peakHilbertFlag",&peakHilbertFlag,"peakHilbertFlag/I");
  tdataflags->Branch("eventTracedFlag",&eventTracedFlag,"eventTracedFalg/I");
  tdataflags->Branch("triggerFlag",&triggerFlag,"triggerFlag/I");
  
  tdataflags->Branch("eventPointedFlag",&eventPointedFlag,"eventPointedFlaf/I");

  TTree *tdataEvent = new TTree("tdataEvent","eventInfo");
  tdataEvent->Branch("eventNumber", &eventNumber,"eventNumber/I");
  

  TTree *tdataTracing = new TTree("tdataTracing","stuff to plot for tracing");
  tdataTracing->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tdataTracing->Branch("peakThetaFinal",&peakThetaFinal,"peakThetaFinal/D");
  tdataTracing->Branch("peakPhiFinal",&peakPhiFinal,"peakPhiFinal/D");
  tdataTracing->Branch("sourceLon",&sourceLon,"sourceLon/D");
  tdataTracing->Branch("sourceLat",&sourceLat,"sourceLat/D");
  tdataTracing->Branch("sourceHeight",&sourceAlt,"sourceHeight/D");
  tdataTracing->Branch("heading",&anitaHeading,"heading/D");
  tdataTracing->Branch("eventTracedFlag",&eventTracedFlag,"eventTracedFlag/I");
  tdataTracing->Branch("anitaLatitude",&anitaLatitude,"anitaLatitude/D");
  tdataTracing->Branch("anitaLongitude",&anitaLongitude,"anitaLongitude/D");
  tdataTracing->Branch("anitaAltitude",&anitaAltitude,"anitaAltitude/D");
  tdataTracing->Branch("quietBaseFlagIndex",&quietBaseFlagIndex,"quietBaseFlagIndex/I");
  tdataTracing->Branch("realTime",&realTime,"realTime/I");
  tdataTracing->Branch("thetaMap",&thetaMap,"thetaMap/D");
  tdataTracing->Branch("snrPeakAfterFilter",&snrAfterFilter,"snrPeakAfterFilter/D");
  tdataTracing->Branch("phiMap",&phiMap,"phiMap/D");
  tdataTracing->Branch("eventPointedFlag",&eventPointedFlag,"eventPointedFlag/I");
  tdataTracing->Branch("xCorPassFlag",&xCorPassFlag,"xCorPassFlag/I");

  tdataTracing->Branch("ratioOfPeaksPassFlag",&ratioOfPeaksPassFlag,"ratioOfPeaksPassFlag/I");
  tdataTracing->Branch("elevationAnglePassFlag",&elevationAnglePassFlag,"elevationAnglePassFlag/I");
  tdataTracing->Branch("peakCrossCorrFlag",&peakCrossCorrFlag,"peakCrossCorrFlag/I");
  tdataTracing->Branch("polFractionFlag",&polFractionFlag,"polFractionFlag/I");
  tdataTracing->Branch("peakHilbertFlag",&peakHilbertFlag,"peakHilbertFlag/I");
  tdataTracing->Branch("eventTracedFlag",&eventTracedFlag,"eventTracedFalg/I");
  tdataTracing->Branch("triggerFlag",&triggerFlag,"triggerFlag/I");
  tdataTracing->Branch("finaltheta",&finaltheta,"finaltheta/D");
  if (!fUseCalibratedEventFile && !fUseEventFile){
    tdataTracing->Branch("mcWeightEvent",&mcWeightEvent,"mcWeightEvent/D");
    tdataTracing->Branch("mcExponentEvent",&mcExponentEvent,"mcExponentEvent/D");
    tdataTracing->Branch("mcLatEvent",&mcLatEvent,"mcLatEvent/D");
    tdataTracing->Branch("mcLonEvent",&mcLonEvent,"mcLonEvent/D");
    tdataTracing->Branch("mcAltEvent",&mcAltEvent,"mcAltEvent/D");
  }

  TTree *tdataPointed = new TTree("tdataPointed","only good events that passed cuts!");
  tdataPointed->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tdataPointed->Branch("sourceLon",&sourceLon,"sourceLon/D");
  tdataPointed->Branch("sourceLat",&sourceLat,"sourceLat/D");
  tdataPointed->Branch("sourceHeight",&sourceAlt,"sourceHeight/D");
  tdataPointed->Branch("heading",&anitaHeading,"heading/D");
  tdataPointed->Branch("anitaLatitude",&anitaLatitude,"anitaLatitude/D");
  tdataPointed->Branch("anitaLongitude",&anitaLongitude,"anitaLongitude/D");
  tdataPointed->Branch("anitaAltitude",&anitaAltitude,"anitaAltitude/D");
  tdataPointed->Branch("thetaMap",&thetaMap,"thetaMap/D");
  tdataPointed->Branch("phiMap",&phiMap,"phiMap/D");
  tdataPointed->Branch("snrPeakAfterFilter",&snrAfterFilter,"snrPeakAfterFilter/D");
  if (!fUseCalibratedEventFile && !fUseEventFile){
    tdataPointed->Branch("mcWeightEvent",&mcWeightEvent,"mcWeightEvent/D");
    tdataPointed->Branch("mcExponentEvent",&mcExponentEvent,"mcExponentEvent/D");
    tdataPointed->Branch("mcLatEvent",&mcLatEvent,"mcLatEvent/D");
    tdataPointed->Branch("mcLonEvent",&mcLonEvent,"mcLonEvent/D");
    tdataPointed->Branch("mcAltEvent",&mcAltEvent,"mcAltEvent/D");
  }
  

  
  thermalFlagUniversal=thermalFlag;
  
  ftime( &timebuffer );
  timeline = ctime( & ( timebuffer.time ) );
  for (int eventctr=eventCtrStart;eventctr<eventCtrEnd;eventctr+=1){//LOOP START
    eventTracedFlag=0;
    eventPointedFlag=0;
    //get the event
    fHeadTree->GetEntry(eventctr);
    eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
    eventNumber=fHeadPtr->eventNumber;
    anitaLatitude=fAdu5APatPtr->latitude;
    anitaLongitude=fAdu5APatPtr->longitude;
    anitaAltitude=fAdu5APatPtr->altitude;
    realTime=fHeadPtr->realTime;
    anitaHeading=fAdu5APatPtr->heading;
 
    //set flags to -10
    mcmflag=-10;
    tdflag=-10;
    tdreflectionflag=-10;
    syncslipflag=-10;
    rfonlyflag=-10;
    calpulserflag=-10;
    mainrfcmflag=-10;
    nadirrfcmflag=-10;
    channelsaturatedflag=-10;
    dcoffsetflag=-10;
    shorttraceflag=-10;
    varnerflag=-10;
    payloadblastflag=-10;

    xCorPassFlag=0;
    ratioOfPeaksPassFlag=0;
    elevationAnglePassFlag=0;
    peakCrossCorrFlag=0;
    polFractionFlag=0;
    peakHilbertFlag=0;
    triggerFlag=0;
    
    ///////////taylor dome
    if (taylorFlag==1){
      if (isTaylor(fHeadPtr->eventNumber)==1){
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isMainRFCMOn(fHeadPtr->eventNumber) && isShortTrace(fHeadPtr->eventNumber)==0){
	 
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber,drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;
	  
	  if (wantToTraceBackFlag==1)
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  //if (eventTracedFlag==1) checkIfNearAnyBase(fHeadPtr->eventNumber, 
	  //				     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt, baseName);
	  
	  ctr++;
	 
	  if (printFlag==2) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  //  ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,peakPhiFinal*rad2deg,
	  //	     sourceLon,sourceLat,sourceAlt,eventTracedFlag,fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,
	  //	     quietBaseFlagIndex);
	  tdataTracing->Fill();
	  thetaMap=-1.*rad2deg*peakThetaFinal;
	  phiMap=rad2deg*peakPhiFinal;
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && eventPointedFlag==1) tdataPointed->Fill();
	}
      }
    }
    
    ////////tdreflection
     if (taylorFlag==2){
      if (isTaylorReflection(fHeadPtr->eventNumber)==1){
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isMainRFCMOn(fHeadPtr->eventNumber) && isShortTrace(fHeadPtr->eventNumber)==0){
	  
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber,drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;
	  
	  if (wantToTraceBackFlag==1)
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  
	  ctr++;
	  if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  //  ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,peakPhiFinal*rad2deg,
	  //	     sourceLon,sourceLat,sourceAlt,eventTracedFlag,fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,
	  //	     quietBaseFlagIndex);

	  tdataTracing->Fill();
	}
      }
    }


    ///////////thermals/SW
    if (thermalFlag==1){
      if ((fHeadPtr->trigType&(1<<3)) && isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber)==0
	  && isTaylor(fHeadPtr->eventNumber)==0 && isTaylorReflection(fHeadPtr->eventNumber)==0
	  && isCalPulser(fHeadPtr->eventNumber)==0 && isMainRFCMOn(fHeadPtr->eventNumber)){ //&& fHeadPtr->triggerTimeNs>1000000){ 
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber)
	    && isShortTrace(fHeadPtr->eventNumber)==0){
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;
	  
	  if (wantToTraceBackFlag==1)
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);  
	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  
	  ctr++;
	  if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  // ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,
	  //	     peakPhiFinal*rad2deg,sourceLon,sourceLat,sourceAlt,eventTracedFlag,
	  //	     fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,
	  //	     quietBaseFlagIndex);
	  tdataTracing->Fill();
	  
	}
      }
    }
    
    if (thermalFlag==2){//all non-rf 
      if ((fHeadPtr->trigType&(1<<0))==0 && isCalPulser(fHeadPtr->eventNumber)==0 
	  && isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber)==0 && isTaylor(fHeadPtr->eventNumber)==0
	  && isTaylorReflection(fHeadPtr->eventNumber)==0 && isMainRFCMOn(fHeadPtr->eventNumber)){ 
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isShortTrace(fHeadPtr->eventNumber)==0 ){
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;
	  
	  if (wantToTraceBackFlag==1)
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);  
	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  
	  ctr++;
	  if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  //ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,
	  //	     peakPhiFinal*rad2deg,sourceLon,sourceLat,sourceAlt,eventTracedFlag,
	  //	     fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,quietBaseFlagIndex );
	  tdataTracing->Fill();
	  
	}
      }
    }
    
    if (thermalFlag==3){//upward pointing rf events 
      if ((fHeadPtr->trigType&(1<<0)) && isCalPulser(fHeadPtr->eventNumber)==0 && isTaylor(fHeadPtr->eventNumber)==0
	  && isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber)==0
	  && isTaylorReflection(fHeadPtr->eventNumber)==0 && isMainRFCMOn(fHeadPtr->eventNumber)){
	  //&& isVarnerEvent(fHeadPtr->eventNumber)==0){ 
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isShortTrace(fHeadPtr->eventNumber)==0 ){
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
      
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  if (peakThetaFinal<0) tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal<<endl;	  
	  if (peakThetaFinal<0){//upward pointing
	    printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		   timebuffer.millitm);
	    

	    ctr++;
	    if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	    //  ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,
	    //	       peakPhiFinal*rad2deg,0,0,0,0,
	    //	       fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude, quietBaseFlagIndex);
	    tdataTracing->Fill();
	  }
	}
      }
    }
    
    if (thermalFlag==4){// rf events near a small base with no events from ANITA-I 
      if ((fHeadPtr->trigType&(1<<0)) && isCalPulser(fHeadPtr->eventNumber)==0 && isTaylor(fHeadPtr->eventNumber)==0
	  && isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber)==0
	  && isTaylorReflection(fHeadPtr->eventNumber)==0){
	  // && isVarnerEvent(fHeadPtr->eventNumber)==0){ 
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isMainRFCMOn(fHeadPtr->eventNumber) && isShortTrace(fHeadPtr->eventNumber)==0){
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
      
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;	  
	  if (wantToTraceBackFlag==1){
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  }
	  if (printFlag==1) cout<<"event traced flag: "<<eventTracedFlag<<endl;
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  for (int i=0;i<4;i++){
	    quietBaseFlag[i]=0;
	  }
	  if (eventTracedFlag==1){
	    quietBaseFlag[0]=checkIfNearSpecificBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt,"AWS_Zoe_(Mega_B)",-1);
	    quietBaseFlag[1]=checkIfNearSpecificBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt,"AWS_AGO_4",-1);
	    quietBaseFlag[2]=checkIfNearSpecificBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt,"Kohnen_Station_(EPICA)[Germany]",-1);
	    quietBaseFlag[3]=checkIfNearSpecificBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt,"A._de_Navajo_Sobral[Argentina]",-1);  
	  }
	  for (int i=0;i<4;i++){
	    if (quietBaseFlag[i]==1) quietBaseFlagIndex=i;
	  }
	  
	  if ((quietBaseFlag[0]==1 || quietBaseFlag[1]==1 ||quietBaseFlag[2]==1 ||quietBaseFlag[3]==1) && eventTracedFlag==1){
	    ctr++;
	    printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		   timebuffer.millitm);
	    if (printFlag==2) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	    // ndataTracing->Fill(fHeadPtr->eventNumber,peakThetaFinal*rad2deg,
	    //	       peakPhiFinal*rad2deg,sourceLon,sourceLat,sourceAlt,eventTracedFlag,
	    //	       fAdu5APatPtr->latitude,fAdu5APatPtr->longitude, fAdu5APatPtr->altitude, quietBaseFlagIndex);
	    tdataTracing->Fill();
	  }
	}
      }
    }

    if (thermalFlag==0 && taylorFlag==0){//all events that pass basic cuts!!!!
      if (fHeadPtr->trigType&(1<<1) || fHeadPtr->trigType&(1<<2) || fHeadPtr->trigType&(1<<3) ||
	  !fHeadPtr->trigType&(1<<0)){
	rfonlyflag=0;
	
	cout<<"eventnumber for forced trigger is "<<eventNumber<<"\n";
      }
      else rfonlyflag=1;
      
      calpulserflag=isCalPulser(fHeadPtr->eventNumber);
      tdflag=isTaylor(fHeadPtr->eventNumber);
      tdreflectionflag=isTaylorReflection(fHeadPtr->eventNumber);
      mcmflag=isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber);
      syncslipflag=isSyncSlip(fHeadPtr->eventNumber);
      
      if (rfonlyflag!=0 && calpulserflag==0 && tdflag==0
	  && mcmflag==0 && tdreflectionflag==0 && !syncslipflag && gpsBadFlag==0){

	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	
	channelsaturatedflag=isChannelSaturated(fHeadPtr->eventNumber);
	cout<<"here after saturated \n";
	mainrfcmflag=isMainRFCMOn(fHeadPtr->eventNumber);   
	bigenoughpeakflag=isBigEnoughPeakToPeak(fHeadPtr->eventNumber); 
	dcoffsetflag=isDCOffsetLarge(fHeadPtr->eventNumber);
	shorttraceflag=isShortTrace(fHeadPtr->eventNumber);
	varnerflag2=isVarnerEvent2(fHeadPtr->eventNumber);
	varnerflag=isVarnerEvent(fHeadPtr->eventNumber);
	nadirrfcmflag=isNadirRFCMOn(fHeadPtr->eventNumber);
	payloadblastflag=isPayloadBlast(fHeadPtr->eventNumber);
	cout<<"channel, mainrfcm, bigenough, dcoffset, shorttraceflag is "<<channelsaturatedflag<<","<<mainrfcmflag<<","<<bigenoughpeakflag<<","<<dcoffsetflag<<","<<shorttraceflag<<"\n";
	if (channelsaturatedflag<=3 && mainrfcmflag && shorttraceflag==0 && dcoffsetflag==0 && bigenoughpeakflag==1){//ACG used to be channelsaturatedflag==-1
	  ftime( &timebuffer );
	  timeline = ctime( & ( timebuffer.time ) );
	  cout<<"inside the if statement. should point next \n";
	  //cout<<"event number header: "<<fHeadPtr->eventNumber<<", event number event: "<<fHeadPtr->eventNumber<<endl;
	  
	  eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber, drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);

	  
	  
	  cout<<"eventnumber is "<<eventNumber<<" "<<passed<<" passed, "<<percentage<<"% \n";
	 
	  tdataEvent->Fill();
	  if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
						       <<" ...about to trace back!"<<endl;	  
	  if (wantToTraceBackFlag==1){
	    eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
						 peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
	  }//want to traceBackFlag
	  if (printFlag==1) cout<<"event traced flag: "<<eventTracedFlag<<endl;
	  if (eventNumber==2462929) cout<<"ratio, peakCross,peakHilbert,polFraction,xCor,elevation,eventTraced,trigger are "<<ratioOfPeaksPassFlag<<" "<<peakCrossCorrFlag<<" "<<peakHilbertFlag<<" "<<polFractionFlag<<" "<<xCorPassFlag<<" "<<elevationAnglePassFlag<<" "<<eventTracedFlag<<" "<<triggerFlag<<"\n";
	  if (eventTracedFlag==1 || eventTracedFlag==2) checkIfNearAnyBase(fHeadPtr->eventNumber, 
						     peakThetaFinal, peakPhiFinal, sourceLon, sourceLat, sourceAlt, baseName);
	  thetaMap=-1.*rad2deg*peakThetaFinal;
	  phiMap=rad2deg*peakPhiFinal;

	  if(eventNumber==2462929){
	    cout<<"thetaMap is "<<thetaMap<<" phiMap is "<<phiMap<<"\n";

	  }
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && eventPointedFlag==1){
	    //rootfile output for clustering algorithm
	    //events that pointed and traced to continent
	    tdataPointed->Fill();
	  }//eventTracedFlag	  
	  if ((eventTracedFlag==1 || eventTracedFlag==2) && drawMaps==1) 
	    drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
	  
	  printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
		 timebuffer.millitm);
	  if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
	  
	  tdataTracing->Fill();
	}//channelsaturated 
      }//all flags ==0
      
      tdataflags->Fill();
      
    }//thermal&taylordome==0
   
    ///////////simulated events loop
    if (thermalFlag==-1 && taylorFlag==-1){	 
      ftime( &timebuffer );
      timeline = ctime( & ( timebuffer.time ) );
      printf("Event: %d, the time is %.19s.%hu \n",fHeadPtr->eventNumber, timeline, 
	     timebuffer.millitm);
      
      eventPointedFlag=pointThisEvent(fHeadPtr->eventNumber,drawMaps, ndata, ndata2, ndata3, ndata4, peakThetaFinal,peakPhiFinal,whichPolarization,xCorPassFlag,ratioOfPeaksPassFlag,elevationAnglePassFlag,peakCrossCorrFlag,polFractionFlag,peakHilbertFlag,triggerFlag,finaltheta);
      tdataEvent->Fill();
      if (printFlag==1) cout<<"Event Pointed To "<<peakThetaFinal<<", "<<peakPhiFinal
			    <<" ...about to trace back!"<<endl;
      
      if (wantToTraceBackFlag==1)
	eventTracedFlag=traceBackToContinent_Brian(fHeadPtr->eventNumber, 
					     peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
      
      if ((eventTracedFlag==1) && drawMaps==1) 
	drawEventOnContinent(fHeadPtr->eventNumber,sourceLon,sourceLat);
      
      ctr++;
      if (printFlag==1) cout<<"Event Included, Counter "<<ctr<<", "<<eventctr<<endl;
      tdataTracing->Fill();
      thetaMap=-1.*rad2deg*peakThetaFinal;
      phiMap=rad2deg*peakPhiFinal;
      if ((eventTracedFlag==1) && eventPointedFlag==1) tdataPointed->Fill();            
    }
   
  }
  cout<<"writing data" <<endl;
  
  if (writeData==1){    
    //tdataTracing->Print();
    // rootfile->cd();
    //tdataTracing->Write();
    rootfile=tdataTracing->GetCurrentFile();
    rootfile->Write("tdataTracing");
    //rootfile->Write("ndata");
    //rootfile->Write("ndata2");
    rootfile->Close();
  }
  cout<<"exiting loopoverthirdevents \n";
}
///////////////////////



void MyCorrelator::readBaselineFFTs()
{

  TFile *fread = new TFile("outputcreateBaselineRun191_test.root");
  TTree *tBaseline= (TTree*)fread->Get("tBaseline");
  int ant;
  const int nPoints=129;//124
  int nPointsForTree;
  float frequencyArray[nPoints];
  float vertFFT[nPoints];
  float horizFFT[nPoints];

  tBaseline->SetBranchAddress("ant",&ant);
  tBaseline->SetBranchAddress("nPointsForTree",&nPointsForTree);
  tBaseline->SetBranchAddress("frequencyArray",frequencyArray);
  tBaseline->SetBranchAddress("vertFFT",vertFFT);
  tBaseline->SetBranchAddress("horizFFT",horizFFT);

  float vertCoherentFFT[nPoints];
  float horizCoherentFFT[nPoints];
  
  for (int i=0;i<nPoints;i++){
    vertCoherentFFT[i]=0.;//zero
    horizCoherentFFT[i]=0.;//zero 
  }

  //create an array of FFT -> vertCoherentFFT[i][j].  i for antennas, j for freq. Needs to be able to be used in other functions.
 
  for (int i=0;i<NUM_ANTS_WITH_NADIRS;i++){//for all antennas, get coherent FFT. This is the part I would need to move to later
    tBaseline->GetEntry(i);
   
    for (int j=0;j<nPointsForTree;j++){
     
      baselineFreq[j]=frequencyArray[j];
      if (frequencyArray[j]>=200 && frequencyArray[j]<=1200){
	if (i!=1) vertCoherentFFT[j]+=pow(10,vertFFT[j]/10.);//get rid of 2V
	horizCoherentFFT[j]+=pow(10,horizFFT[j]/10.);
	baselinearray[i][j]=pow(10,vertFFT[j]/10.);//not in dB
	
	baselinearrayHoriz[i][j]=pow(10,horizFFT[j]/10.);
      }
    }
    
    grVertBaseline[i]=new TGraph(nPointsForTree,frequencyArray,vertFFT);
    grHorizBaseline[i]=new TGraph(nPointsForTree,frequencyArray,vertFFT);//change to horiz
     
  }
  //change to dB
  for (int j=0;j<nPointsForTree;j++){
    if (frequencyArray[j]>=200 && frequencyArray[j]<=1200){
      vertCoherentFFT[j]=10*log10(vertCoherentFFT[j]/(double(NUM_ANTS_WITH_NADIRS)-1));
      horizCoherentFFT[j]=10*log10(horizCoherentFFT[j]/double(NUM_ANTS_WITH_NADIRS));
      
    }
    else{
      vertCoherentFFT[j]=-1000;
      horizCoherentFFT[j]=-1000;
      }
  }

  for(int i1=0;i1<NUM_ANTS_WITH_NADIRS;i1++){
    for(int j1=0;j1<nPointsForTree;j1++){
      baselinearray_dB[i1][j1]=10*log10(baselinearray[i1][j1]);
      baselinearrayHoriz_dB[i1][j1]=10*log10(baselinearrayHoriz[i1][j1]);
    }
  }

  



  float vertCoherentFFTAvg[nPoints];
  float horizCoherentFFTAvg[nPoints];
  
  //now do some averaging.  Fo shizzle.
  //use new array for averaging on each antenna
  for (int j=0;j<nPointsForTree;j++){
    vertCoherentFFTAvg[j]=vertCoherentFFT[j];
    horizCoherentFFTAvg[j]=horizCoherentFFT[j];
    if (frequencyArray[j-3]>=260 && frequencyArray[j+3]<=1200){
      vertCoherentFFTAvg[j]=(vertCoherentFFT[j]+vertCoherentFFT[j-1]+vertCoherentFFT[j+1]+vertCoherentFFT[j-2]+vertCoherentFFT[j+2]
			     +vertCoherentFFT[j-3]+vertCoherentFFT[j+3])/7.;
      horizCoherentFFTAvg[j]=(horizCoherentFFT[j]+horizCoherentFFT[j-1]+horizCoherentFFT[j+1]+horizCoherentFFT[j-2]+horizCoherentFFT[j+2]
			      +horizCoherentFFT[j-3]+horizCoherentFFT[j+3])/7.;
    }
    else if (frequencyArray[j-2]>=260 && frequencyArray[j+2]<=1200){
      vertCoherentFFTAvg[j]=(vertCoherentFFT[j]+vertCoherentFFT[j-1]+vertCoherentFFT[j+1]+vertCoherentFFT[j-2]+vertCoherentFFT[j+2])/5.;
      horizCoherentFFTAvg[j]=(horizCoherentFFT[j]+horizCoherentFFT[j-1]+horizCoherentFFT[j+1]+horizCoherentFFT[j-2]+horizCoherentFFT[j+2])/5.;
    }
    else if (frequencyArray[j-1]>=200 && frequencyArray[j+1]<=1200){
      vertCoherentFFTAvg[j]=(vertCoherentFFT[j]+vertCoherentFFT[j-1]+vertCoherentFFT[j+1])/3.;
      horizCoherentFFTAvg[j]=(horizCoherentFFT[j]+horizCoherentFFT[j-1]+horizCoherentFFT[j+1])/3.;
    }

  }//averaging over j

  /*  double baselinearray_avg[NUM_ANTS_WITH_NADIRS][nPoints];

  for(int i=0;i<NUM_ANTS_WITH_NADIRS;i++){
    for(int j=0;j<nPointsForTree;j++){
      baselinearray_avg[i][j]=baselinearray[i][j];
       if (frequencyArray[j-3]>=260 && frequencyArray[j+3]<=1200){
      baselinearray_avg[i][j]=(baselinearray[i][j]+baselinearray[i][j-1]+baselinearray[i][j+1]+baselinearray[i][j-2]
			    +baselinearray[i][j+2] +baselinearray[i][j-3]+baselinearray[i][j+3])/7.;
       }
       else if (frequencyArray[j-2]>=260 && frequencyArray[j+2]<=1200){
       baselinearray_avg[i][j]=(baselinearray[i][j]+baselinearray[i][j-1]+baselinearray[i][j+1]+baselinearray[i][j-2]
			    +baselinearray[i][j+2])/5.;
       }
       else if (frequencyArray[j-1]>=200 && frequencyArray[j+1]<=1200){
	  baselinearray_avg[i][j]=(baselinearray[i][j]+baselinearray[i][j-1]+baselinearray[i][j+1])/3.;
       }
    }//j

  }//average i

  */

  //now do some smoothing over the satellite hump.
   for (int j=0;j<nPointsForTree;j++){
     if (frequencyArray[j]>=230 && frequencyArray[j]<=310){
       vertCoherentFFTAvg[j]=(vertCoherentFFTAvg[30]-vertCoherentFFTAvg[21])*(j-21)/9.+vertCoherentFFTAvg[21];
       horizCoherentFFTAvg[j]=(horizCoherentFFTAvg[30]-horizCoherentFFTAvg[21])*(j-21)/9.+horizCoherentFFTAvg[21];
     }
   }

   /*  for(int i=0;i<NUM_ANTS_WITH_NADIRS;i++){
      for(int j=0;j<nPointsForTree;j++){
	 if (frequencyArray[j]>=230 && frequencyArray[j]<=310){
	    baselinearray_avg[i][j]=(baselinearray_avg[i][30]-baselinearray_avg[i][21])*(j-21)/9.+baselinearray_avg[i][21];
	 }
      }
    }
   */
   //can I add the spectra back together now?
   /*  if (!grVertCoherentBaseline)
    grVertCoherentBaseline=new TGraph(nPointsForTree,frequencyArray,vertCoherentFFTAvg);
  if (!grHorizCoherentBaseline)
    grHorizCoherentBaseline=new TGraph(nPointsForTree,frequencyArray,horizCoherentFFTAvg);
   */
  readBaselineFlag=1;
  fread->Close();
  
}
//////////////////////////////////
void MyCorrelator::GetBaselineperPhi(int pol, double *baseline, int nantennasToUse, vector<int>& whichAntennasToUse ){
 
  int nPoints=129;//124
  int ant;

  for (int i=0;i<nantennasToUse;i++){
    ant = whichAntennasToUse[i];
    for(int j=0;j<nPoints;j++){
      if(pol==1){
	baseline[j]+=baselinearrayHoriz[ant][j];
      }
      else{
	baseline[j]+=baselinearray[ant][j];//not in dB
      }
    }//j
  }//i

  //change to dB
  for(int j1=0;j1<nPoints;j1++){
    if (baselineFreq[j1]>=200 && baselineFreq[j1]<1200){
      baseline[j1]=10*log10(baseline[j1]/(double(nantennasToUse)));
    }
    else
      baseline[j1]=-1000;
  }


  //time for some averaging
  double baselinearray_avg[nPoints];

  for(int j=0;j<nPoints;j++){
    baselinearray_avg[j]=baseline[j];
    if (baselineFreq[j-3]>=260 && baselineFreq[j+3]<=1200){
      baselinearray_avg[j]=(baseline[j]+baseline[j-1]+baseline[j+1]+baseline[j-2]
			    +baseline[j+2] +baseline[j-3]+baseline[j+3])/7.;
    }
    else if (baselineFreq[j-2]>=260 && baselineFreq[j+2]<=1200){
      baselinearray_avg[j]=(baseline[j]+baseline[j-1]+baseline[j+1]+baseline[j-2]+baseline[j+2])/5.;
    }
    else if (baselineFreq[j-1]>=200 && baselineFreq[j+1]<1200){
      baselinearray_avg[j]=(baseline[j]+baseline[j-1]+baseline[j+1])/3.;
    }
  }//j
  
  

  //average over satellite bump
  for(int j=0;j<nPoints;j++){
    if (baselineFreq[j]>=230 && baselineFreq[j]<=310){
      baselinearray_avg[j]=(baselinearray_avg[30]-baselinearray_avg[21])*(j-21)/9.+baselinearray_avg[21];
    }
  }
  
  for(int j=0;j<nPoints;j++){
    baseline[j]=baselinearray_avg[j];
  }

    //can I add the spectra back together now?
  //if (!grVertCoherentBaseline)
  if(pol==1){
    grHorizCoherentBaseline=new TGraph(nPoints,baselineFreq,baseline);
  }
  else{
    grVertCoherentBaseline=new TGraph(nPoints,baselineFreq,baseline);
  }
  // readBaselineFlag=1;
  



}//filter


/////////////////////////////////
double MyCorrelator::getMaximum(int n, double *array, int &index){
  double max;
  max=array[0];
  index=0;
  for (int i=1;i<n;i++){
    if (array[i]>max){ 
      max=array[i];
      index=i;
    }
  }
  return max;
}
////////////////////////////////////
double MyCorrelator::getRMS(TGraph *gr)
{
  double *y=gr->GetY();
  int numPoints=gr->GetN();
  double rms=TMath::RMS(numPoints,y);
  return rms;
}
////////////////////////////
double MyCorrelator::getSNR(TGraph *gr,double &rmsNoise)
{
  //double *y=gr->GetY();
  int numPoints=gr->GetN();
  //rmsNoise=TMath::RMS(numPoints/5,y); //first 1/5th of the graph
  //cout<<"old rms is "<<rmsNoise<<" rms ignoring first bin is "<<getRMSOfRange(gr,1,numPoints/5)<<"\n";
  rmsNoise = getRMSOfRange(gr,1,numPoints/5);
  //cout<<"rmsNoise is "<<rmsNoise<<"\n";
  double peak2peak=getPeak2Peak(gr);
  //cout<<"peak2peak is "<<peak2peak<<"\n";
  double SNR=peak2peak/(rmsNoise*2);
  //cout<<"SNR is "<<SNR<<"\n";
  //cout<<"SNR ignoring 1st bin is "<<peak2peak/(getRMSOfRange(gr,1,numPoints/5)*2)<<"\n";
  return SNR;
}
//////////////////////////////
double MyCorrelator::getRMSOfRange(TGraph *gr, double xLow, double xHigh)
{
  double *y=gr->GetY();
  //double *x=gr->GetX();
  int numPoints=gr->GetN();
  int numNew=0;
  double yNew[2000];
  int j=0;
  
  for (int i=0;i<numPoints;i++){
    // if (x[i]>=xLow && x[i]<=xHigh){
    if (i>=xLow && i<xHigh){
      numNew++;
      yNew[j]=y[i];
      j++;
      //      cout<<"numnew: "<<numNew<<", j: "<<j<<", ynew: "<<yNew[j]<<endl;
    }
  }
  double rms=TMath::RMS(numNew,yNew);
  return rms;
 
}
///////////////////////////////
int MyCorrelator::getPeakAntenna(int myEventNumber, int nantennas)
{
  if (!fEventTree) initialize();
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  
  int antThis=-1;
  double rmsNoiseThis=0;
  double maxPeakValAfterFilter=0;
  for (int ant=0;ant<nantennas;ant++){
    Int_t peakBin=FFTtools::getPeakBin(grEv[ant]);
    double tVal, dummyPeakVal;
    grEv[ant]->GetPoint(peakBin,tVal,dummyPeakVal);
    if (dummyPeakVal>maxPeakValAfterFilter){ 
      maxPeakValAfterFilter=dummyPeakVal;
      snrAfterFilter=getSNR(grEv[ant],rmsNoiseThis);    
      antThis=ant;
    }
  }
  return antThis;
}
////////////////////////////////////
double MyCorrelator::getPeak2Peak(TGraph *gr)
{
  double posPeakVal=0;
  double negPeakVal=0;
  double peak2peak;

  double *y=gr->GetY();
  int numPoints=gr->GetN();
  for (int i=0;i<numPoints;i++){
    if (y[i]>posPeakVal) posPeakVal=y[i];
    if (y[i]<negPeakVal) negPeakVal=y[i];
  }
  peak2peak=posPeakVal+fabs(negPeakVal);
  return peak2peak;

}
////////////////////////////////
double MyCorrelator::getPeakHilbert(TGraph *gr)
{
  Double_t peakTimeHilbert=0;
  Double_t peakValHilbert=0;	
  TGraph *ghilbert=FFTtools::getHilbertEnvelope(gr);
  Int_t peakBinHilbert=FFTtools::getPeakBin(ghilbert);
  ghilbert->GetPoint(peakBinHilbert,peakTimeHilbert,peakValHilbert);
  delete ghilbert;
  return peakValHilbert;
  
}
///////////////////////////////
TGraph *MyCorrelator::simpleNotchFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq)
{

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);

    int newLength=(length/2)+1;

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e3; //MHz
    //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;

    double tempF=0;
    for(int i=0;i<newLength;i++) {
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      if(tempF>minFreq && tempF<maxFreq) {
	theFFT[i].re=0;
	theFFT[i].im=0;
      }      
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      tempF+=deltaF;
    }

    double *filteredVals = FFTtools::doInvFFT(length,theFFT);

    TGraph *grFilteredNotch=new TGraph(length,oldX,filteredVals);  
    // }
    delete [] theFFT;
    delete [] filteredVals;
    return grFilteredNotch;

}
/////////////////////////////////////
TGraph *MyCorrelator::complicatedNotchFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq, int ant, int pol,double *baseY)
{
  
  //TRandom3 random; // for generating random numbers
  //random.SetSeed(fHeadPtr->eventNumber);
  //cout<<"minFreq is "<<minFreq<< " max freq is "<<maxFreq<<"\n";
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();//256
  double frequencyArray[2000];
  double magFFT[2000];

  for(int j0=0;j0<2000;j0++){
    frequencyArray[j0]=0;
    magFFT[j0]=0;
  }
 
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  
  int newLength=(length/2)+1;
  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz

    //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;
  double mag=0.;
  double tempF=0.;
  int temperature=340;//tsys+tice in kelvin
  double gain=75.;  
  double meanamp=sqrt(1.38e-23*double(temperature)*deltaF*1e6*50*pow(10.,gain/10.)/2.)*sqrt(pi/2); //sigma of a gaussian distribution of imaginary and real components
  float u1, u2, v1, v2,z1,z2,rsquared;
  double x,y,phi;
  float u3,u4,v4;
  
  char filename[256];
  string rayleigh_string;
  if(pol==0){
    rayleigh_string ="Vert";
  }
  if(pol==1){
    rayleigh_string ="Horiz";
  }

  sprintf(filename,"/home/dailey.110/analysis/rayleigh%s_run191.root",rayleigh_string.c_str());
  
  //TFile *rayleighfile = new TFile(filename);
  
  
  
  double parameters[20][2];
  double index;
  double sigma_rayleigh;
  double sigma;
  double sigma_squared;
  int antctr;
  const char* dataname;
  
  double phi_old;
  double meantosigma = sqrt(2/ TMath::Pi());
  double phase_intercept=0.;
  double last_phi=0.;
  double phi_interp=0.;
  if(pol==0){
    dataname = "Vertdata"; 
  }
  if(pol==1){
    dataname = "Horizdata";  
  }
  /*  TTree *data = (TTree*)rayleighfile->Get(dataname);
  if(newnotchflag==1){
    data->SetBranchAddress("ant",&antctr);
    data->SetBranchAddress("parameters",&parameters);
    
  }
  */
  for(int i=0;i<newLength;i++) {
    //std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
    //std::cout<<"minFreq is "<<minFreq<< " max freq is "<<maxFreq<<"\n";

    theFFTarray[i][0]=theFFT[i].re;
    theFFTarray[i][1]=theFFT[i].im;
    
    
    if(tempF>minFreq && tempF<maxFreq ) {
     
      while(1){
	
	u1=gRandom->Rndm();
	u2=gRandom->Rndm();
	v1=2*u1-1;
	v2=2*u2-1;
	rsquared=v1*v1+v2*v2;
	if (rsquared<=1) break;
	else continue;
      }

      if(newnotchflag==0){
	cout<<"using Abby filter \n";
	z1=meanamp*v1*sqrt(-2*log(rsquared)/rsquared);
	z2=meanamp*v2*sqrt(-2*log(rsquared)/rsquared);
	//cout<<"meanamp is "<<meanamp<<" rsquared is "<<rsquared<<"\n";
	//cout<<"real part was "<<theFFT[i].re<<" is now "<<z1<<"\n";
	//cout<<"im part was "<<theFFT[i].im<<" is now "<<z2<<"\n";
	
	
	// cout<<"now z1 is "<<z1<<" and z2 is "<<z2<<"\n";
      
	z1*=1E3;
	z2*=1E3;
	
	theFFT[i].re=z1;//random
	theFFT[i].im=z2;//random
	
      }//oldway
      else{
	
	//data->GetEvent(ant);//PUT WHICH ANTENNA INTO FUNCTION CALL

	/*	if(ant!=antctr){
	  cout<<"wrong antenna! ant is "<<ant<<" and event ant is "<<antctr<<"\n";
	}
	*/
	index = tempF/50.;
	index = index-4.;
	index = floor(index);

	u3 = gRandom->Rndm();
	u4 = gRandom->Rndm();
	v4 = 2*u4-1;

	if(u3==1){
	  cout<<"picked 1. BAD! choose again \n";
	  u3 = gRandom->Rndm();
	}

	if(index<20){
	  //sigma_rayleigh = parameters[(int) index][0];
	  
	  sigma = baseY[i]*meantosigma;
	  sigma_squared = pow(sigma,2);
	 
	}
	else{
	  sigma=0;
	}

	

	mag = 1-u3;
	if(mag==0){
	  cout<<"mag was 0! \n";
	  mag=pow(10,-10);
	}
	mag = log(mag);//base e
	
	mag = -2*sigma_squared*mag;

	mag = sqrt(mag);
	
	phi_old = atan2(theFFT[i].im,theFFT[i].re);//old phase
	last_phi = atan2(theFFT[i-1].im,theFFT[i-1].re);
	phi = v4 * TMath::Pi();
	phase_intercept = last_phi-phase_slope*(tempF - deltaF);
	phi_interp = phase_slope*tempF + phase_intercept;
	if(phase_flag==1){
	  x = mag*cos(phi);
	  y = mag*sin(phi);
	}
	else if(phase_flag==2){
	  x = mag*cos(phi_interp);
	  y = mag*sin(phi_interp);
	}
	else{
	  x = mag*cos(phi_old);
	  y = mag*sin(phi_old);
	 
	}

	
	  theFFT[i].re=x;//fill in with rayleigh noise
	  theFFT[i].im=y;//fill in with rayleigh nosie
	
      }//else oldway
     
    }//in inside freq range      
    //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
    tempF+=deltaF;
  }
  
  //rayleighfile->Close();
  
  double *filteredVals = FFTtools::doInvFFT(length,theFFT);
  TGraph *grFilteredNotch=new TGraph(length,oldX,filteredVals);  //oldX

  // }
 
  delete [] theFFT;
 
  delete [] filteredVals;
 
  //delete rayleighfile;
  //rayleighfile->Close();
  return grFilteredNotch; 

}
//////////////

TGraph *MyCorrelator::nofillNotchFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq)
{
  
  //TRandom3 random; // for generating random numbers
  //random.SetSeed(fHeadPtr->eventNumber);
  //cout<<"minFreq is "<<minFreq<< " max freq is "<<maxFreq<<"\n";
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();//256
  double frequencyArray[2000];
  double magFFT[2000];
  float u4,v4;
  for(int j0=0;j0<2000;j0++){
    frequencyArray[j0]=0;
    magFFT[j0]=0;
  }
 
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  
  int newLength=(length/2)+1;
  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  double tempF=0.;
  double phi_old=0.;
  double phi=0.;
  double mag=3E-10;
  double x=0.;
  double y=0.;

  double phase_intercept=0.;
  double last_phi=0.;
  double phi_interp=0.;
  for(int i=0;i<newLength;i++) {
   
    if(tempF>minFreq && tempF<maxFreq ) {
      phi_old = atan2(theFFT[i].im,theFFT[i].re);//old phase
      last_phi = atan2(theFFT[i-1].im,theFFT[i-1].re);
      u4 = gRandom->Rndm();
      v4 = 2*u4-1;
      phi = v4 * TMath::Pi();
      phase_intercept = last_phi-phase_slope*(tempF - deltaF);
      phi_interp = phase_slope*tempF + phase_intercept;
       if(phase_flag==1){
	 x = mag*cos(phi);
	 y = mag*sin(phi);
       }
       else if(phase_flag==2){
	 x = mag*cos(phi_interp);
	 y = mag*sin(phi_interp);
       }
       else{
	 x = mag*cos(phi_old);
	 y = mag*sin(phi_old);
       }
      theFFT[i].re=x;
      theFFT[i].im=y;
     
    }
    tempF+=deltaF;
  }

  double *filteredVals = FFTtools::doInvFFT(length,theFFT);
  TGraph *grFilteredNotch=new TGraph(length,oldX,filteredVals);  //oldX

  // }
 
  delete [] theFFT;
 
  delete [] filteredVals;
  return grFilteredNotch; 

}

/////////////////////////////////////
TGraph *MyCorrelator::wienerFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq, int ant, int pol,double *baseY)
{
  
  //TRandom3 random; // for generating random numbers
  //random.SetSeed(fHeadPtr->eventNumsber);
  //cout<<"minFreq is "<<minFreq<< " max freq is "<<maxFreq<<"\n";
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();//256
  double frequencyArray[2000];
  double magFFT[2000];
  
  for(int j0=0;j0<2000;j0++){
    frequencyArray[j0]=0;
    magFFT[j0]=0;
  }
  
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  
  int newLength=(length/2)+1;
  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz

    //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;
  double mag=0.;
  double mag_old=0.;
  double tempF=0.;
  int temperature=340;//tsys+tice in kelvin
  double gain=75.;  
  double meanamp=sqrt(1.38e-23*double(temperature)*deltaF*1e6*50*pow(10.,gain/10.)/2.)*sqrt(pi/2); //sigma of a gaussian distribution of imaginary and real components
  float u1, u2, v1, v2,z1,z2,rsquared;
  double x,y,phi,phi_old;
  float u3,u4,v3,v4;
  
  double wiener_filter;
  double noise_value;
  double min_diff=1000;
  int j_baseline=0;
  double phase_intercept=0.;
  double last_phi=0.;
  double phi_interp=0.;
  for(int i=0;i<newLength;i++) {
    min_diff=1000;
    if(tempF>minFreq && tempF<maxFreq) {
      for(int j=0;j<124;j++){
	if(abs(baselineFreq[j]-tempF)<min_diff){
	  min_diff = abs(baselineFreq[j]-tempF);
	  j_baseline = j;
	}
      }
     
	mag_old = sqrt(pow(theFFT[i].re,2)+pow(theFFT[i].im,2));
	if(mag_old > baseY[j_baseline]){
	  noise_value =sqrt( pow(mag_old,2) - pow(baseY[j_baseline],2));
	}
	else{
	  noise_value=0;
	}

	if(mag_old >0){
	  wiener_filter = pow(baseY[i],2)/pow(mag_old,2);
	}
	else
	  wiener_filter=0;
	phi_old = atan2(theFFT[i].im,theFFT[i].re);//old phase
	last_phi = atan2(theFFT[i-1].im,theFFT[i-1].re);
	u4 = gRandom->Rndm();
	v4 = 2*u4-1;
	phi = v4 * TMath::Pi();
	phase_intercept = last_phi-phase_slope*(tempF - deltaF);
	phi_interp = phase_slope*tempF + phase_intercept;
	if(phase_flag==1){
	  x = mag_old*wiener_filter*cos(phi);
	  y = mag_old*wiener_filter*sin(phi);
	}
	else if(phase_flag==2){
	  x = mag_old*wiener_filter*cos(phi_interp);
	  y = mag_old*wiener_filter*sin(phi_interp);
	}
	else{
	  x = mag_old*wiener_filter*cos(phi_old);
	  y = mag_old*wiener_filter*sin(phi_old);
	 
	}
     
	theFFT[i].re=x;//
	theFFT[i].im=y;//
	
      }//in inside freq range      
    //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
    tempF+=deltaF;
  }
 
  double *filteredVals = FFTtools::doInvFFT(length,theFFT);
  TGraph *grFilteredNotch=new TGraph(length,oldX,filteredVals);  //oldX
 
  delete [] theFFT;
 
  delete [] filteredVals;
 
  return grFilteredNotch; 

}

/////////////////////////////////////

TGraph *MyCorrelator::interpolatedFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq)
{
  if(maxFreq>1200){
    maxFreq=1200;
  }
  //TRandom3 random; // for generating random numbers
  //random.SetSeed(fHeadPtr->eventNumber);
  //cout<<"minFreq is "<<minFreq<< " max freq is "<<maxFreq<<"\n";
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();//256
  double frequencyArray[2000];
  double magFFT[2000];

  for(int j0=0;j0<2000;j0++){
    frequencyArray[j0]=0;
    magFFT[j0]=0;
  }
 
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  
  int newLength=(length/2)+1;
  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  double tempF=0.;
  double startFreq=0.;
  double endFreq=0.;
  int i_start=0;
  int i_end=0;
  for(int i=0;i<newLength;i++){
    //cout<<"minFreq is "<<minFreq<<" maxFreq is "<<maxFreq<<" tempF is "<<tempF<<"\n";
    if(tempF >= (minFreq-deltaF) && tempF <=(minFreq) ){
      //  cout<<"setting starts \n";
      startFreq = tempF;
      i_start = i;
    }
    if(tempF >= (maxFreq) && tempF <= (maxFreq+deltaF)){
      //  cout<<"Setting ends \n";
      endFreq = tempF;
      i_end=i;
    }
    tempF+=deltaF;  
  }
  tempF=0;
  double slope;
  double intercept;
  //cout<<"i_start is "<<i_start<<" i_end is "<<i_end<<" pow(theFFT[i_start].re,2) is "<<pow(theFFT[i_start].re,2)<<" pow(theFFT[i_start].im,2) "<<pow(theFFT[i_start].im,2)<<"\n";
  double mag_start = sqrt(pow(theFFT[i_start].re,2) + pow(theFFT[i_start].im,2));
  double mag_end = sqrt(pow(theFFT[i_end].re,2) + pow(theFFT[i_end].im,2));
  slope = (mag_end - mag_start)/(endFreq - startFreq);
  intercept = mag_start - (slope*startFreq);
  // cout<<"minFreq,maxFreq are "<<minFreq<<" "<<maxFreq<<"\n";
  //cout<<"mag_start, end, freqstart, end, slope, intercept are "<<mag_start<<" "<<mag_end<<" "<<startFreq<<" "<<endFreq<<" "<<slope<<" "<<intercept<<"\n";
  //cout<<"\n\n";
  float u4,v4;
  double mag=0.;
  double phi_old=0.;
  double phi=0.;
  double x=0.;
  double y=0.;
 
  double phase_intercept=0.;
  double last_phi=0.;
  double phi_interp=0.;
  for(int i=0;i<newLength;i++) {
   
    if(tempF>startFreq && tempF<endFreq ) {
      phi_old = atan2(theFFT[i].im,theFFT[i].re);//old phase
      last_phi = atan2(theFFT[i-1].im,theFFT[i-1].re);
      mag = slope*tempF + intercept;
      // cout<<"mag is "<<mag<<"\n";
      u4 = gRandom->Rndm();
      v4 = 2*u4-1;
      phase_intercept = last_phi-phase_slope*(tempF - deltaF);
      phi_interp = phase_slope*tempF + phase_intercept;
      //  cout<<"freq is "<<tempF<<" phi_interp is "<<phi_interp<<"\n";
       phi = v4 * TMath::Pi();
      // phi=ImpulsePhase[i];
      // cout<<"phi_old is "<<phi_old<<" phi is "<<phi<<"\n";
      if(phase_flag==1){
	x = mag*cos(phi);
	y = mag*sin(phi);
      }
      else if(phase_flag==2){
 
	x = mag*cos(phi_interp);
	y = mag*sin(phi_interp);
      }
      else{ 
	x = mag*cos(phi_old);
	y = mag*sin(phi_old);
      }
      theFFT[i].re = x;
      theFFT[i].im = y;
    }//in inside freq range 
   
    tempF+=deltaF;
  }
  
  double *filteredVals = FFTtools::doInvFFT(length,theFFT);
  TGraph *grFilteredNotch=new TGraph(length,oldX,filteredVals);  //oldX

  // }
 
  delete [] theFFT;
 
  delete [] filteredVals;
 
  return grFilteredNotch; 

}

/////////////////////////////////////

double MyCorrelator::getCoherentPolarization(TGraph *grV, TGraph *grH, double &polarizationFraction)//returns polarization angle in degrees
{
  double polarizationAngle;
  double Q,U,I;
  double V0sq=0;
  double V90sq=0;
  double V0V90=0;
  double tValV, tValH, vValV, vValH;
  
  int nentries=grV->GetN();
  if (grH->GetN() < nentries) nentries=grH->GetN();

  for (int i=0;i<nentries;i++){
    grV->GetPoint(i,tValV,vValV);
    grH->GetPoint(i,tValH,vValH);
    vValH*=pow(10,-1.*sqrt(5)/10)*(0.863);//what is this?
    V0sq+=vValH*vValH;
    V90sq+=vValV*vValV;
    V0V90+=vValV*vValH;
  }
  V0sq/=nentries;
  V90sq/=nentries;
  V0V90/=nentries;

  Q=V0sq-V90sq;
  I=V0sq+V90sq;
  U=2*V0V90;
  
  polarizationFraction=sqrt(Q*Q+U*U)/I;
  polarizationAngle=atan(U/Q)/2*rad2deg;
  if (Q<0 && polarizationAngle>=0) polarizationAngle=90-polarizationAngle;
  if (Q<0 && polarizationAngle<0) polarizationAngle=-90-polarizationAngle;

  //get the peaktopeaks of each graph.
  //double peak2peakVert=getPeak2Peak(grV);
  //double peak2peakHoriz=getPeak2Peak(grH); 
  
  //peak2peakHoriz*=pow(10,-1.*sqrt(5)/10)*(0.863);//adjust down by 5dB to get relative gains out (rfcm and splitter) and extra from noise study
  //double polarization=peak2peakVert/peak2peakHoriz;
  //double polarizationAngle=atan(polarization)*rad2deg;
  
  // cout<<"Q: "<<Q<<", I: "<<I<<", U: "<<U<<endl;
  if (printFlag==1)  cout<<"Polarization angle using Coherent waveforms after filter: "<<polarizationAngle<<", fraction: "<<polarizationFraction<<endl;
  return polarizationAngle;

}
////////////////////////////
int MyCorrelator::isChannelSaturated(int myEventNumber){
  if (!fEventTree) initialize();
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  
  int nantennas=NUM_ANTS_WITH_NADIRS;
  int channelSaturated=-1;
  int saturationCut=1500;//mV
  int nsaturatedV=0;
  int nsaturatedH=0;
  int channelSaturatedFlag=0;
  
  for (int i=0;i<nantennas*2;i++){
    saturatedChannels[i]=0;
  }
  //cout<<"nantennas is "<<nantennas<<"\n";

  for (int ant=0;ant<nantennas;ant++){
    TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
    Double_t *volts=gr1->GetY();
   
    for(int i=0;i<gr1->GetN();i++){ 
        
      
      if (volts[i]>saturationCut || volts[i]<-1.*saturationCut){
	channelSaturated=fUPGeomTool->getChanIndexFromAntPol(ant, AnitaPol::kVertical);    
	channelSaturatedFlag=1;
	saturatedChannels[ant]=1;
	
	break;
      }
    }
    delete gr1;
    if (channelSaturatedFlag==1){ 
      nsaturatedV++;
      channelSaturatedFlag=0;
    }
  }

  for (int ant=0;ant<nantennas;ant++){
    TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kHorizontal);
    Double_t *volts=gr1->GetY();
    for(int i=0;i<gr1->GetN();i++){ 
      if (volts[i]>saturationCut || volts[i]<-1.*saturationCut){
	channelSaturated=fUPGeomTool->getChanIndexFromAntPol(ant, AnitaPol::kHorizontal);    
	channelSaturatedFlag=1;
	saturatedChannels[ant+nantennas]=1;
	break;
      }
    }
    delete gr1;
    if (channelSaturatedFlag==1){ 
      nsaturatedH++;
      channelSaturatedFlag=0;
    }
  }
  //cout<<"about to return \n";
  if (nsaturatedV==0 && nsaturatedH==0) return -1;
  //else if (polToggle==0) return nsaturatedV;
  //else return nsaturatedH;
  if (nsaturatedV>nsaturatedH) return nsaturatedV;
  else return nsaturatedH;
  
}
//////////////////////////////////
int MyCorrelator::isPayloadBlast(int myEventNumber)
{
  if (!fEventTree) initialize();
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  
  int nantennas=NUM_ANTS_WITH_NADIRS;
  int nchannels=0;
  int phisectors[16];
  int nphisectors=0;
  int blastCut=400;//mV
  double peak2peakV, peak2peakH;
  int thisphisector;
  
  for(int i=0;i<16;i++){
    phisectors[i]=0;
  }
  
  for (int ant=0;ant<nantennas;ant++){
    TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
    peak2peakV=getPeak2Peak(gr1);
    delete gr1; 
    TGraph *gr2 = fUsefulEventPtr->getGraph(ant,AnitaPol::kHorizontal);
    peak2peakH=getPeak2Peak(gr2);
    if (peak2peakH>blastCut || peak2peakV>blastCut){
      nchannels++;
      thisphisector=fUPGeomTool->getPhiSector(ant);
      phisectors[thisphisector]=1;
    }
    delete gr2;
  }
  
  for(int i=0;i<16;i++){
    if (phisectors[i]==1) nphisectors++;
  }


  if (printFlag==1) cout<<"nchannels for payload blast: "<<nchannels<<", nphi: "<<nphisectors<<endl;
  if (nchannels>15 || nphisectors>9) return 1;
  else return 0;
  
}
//////////////////////////
int MyCorrelator::isVarnerEvent2(int myEventNumber)
{
  if (!fEventTree) initialize();
  int varnerFlag=0;
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  double peak2peak;
  double peak2peakNadirMax=0;
  double peak2peakTopMax=0;
  double peak2peakMax=0;

  for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
    peak2peak=getPeak2Peak(gr1);
    if (ant<16){ //top ring
      if (peak2peak>peak2peakTopMax) peak2peakTopMax=peak2peak;
    }
    if (ant==34 || ant==33){//ones that get blasted by varner events
      if (peak2peak>peak2peakNadirMax) peak2peakNadirMax=peak2peak;
    }
    if (peak2peak>peak2peakMax) peak2peakMax=peak2peak;
    delete gr1;
  }
  for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kHorizontal);
    peak2peak=getPeak2Peak(gr1);
    if (peak2peak>peak2peakMax) peak2peakMax=peak2peak;
    if (ant==34 || ant==33) if (peak2peak>peak2peakNadirMax) peak2peakNadirMax=peak2peak;
    if (ant<16){ //top ring
      if (peak2peak>peak2peakTopMax) peak2peakTopMax=peak2peak;
    }
    delete gr1;
  }

  if (printFlag==1) cout<<"Inside IsVarnerEvent2, peak2peakNadir/Peak2PeakTop: "<<peak2peakNadirMax/peak2peakTopMax<<endl;
  if (peak2peakTopMax<(peak2peakNadirMax/2) && fHeadPtr->trigType&(1<<0) && peak2peakMax>500) varnerFlag=1;
  
  return varnerFlag;
}
///////////////////////////
int MyCorrelator::isVarnerEvent(int myEventNumber)
{
  if (!fEventTree) initialize();
  int varnerFlag=0;
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  double peak2peak;
  double peak2peakNadirMax=0;
  double peak2peakTopMax=0;
  
  for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
    peak2peak=getPeak2Peak(gr1);
    if (ant<16){ //top ring
      if (peak2peak>peak2peakTopMax) peak2peakTopMax=peak2peak;
    }
    if (ant>31){//nadir ring
      if (peak2peak>peak2peakNadirMax) peak2peakNadirMax=peak2peak;
    }
    delete gr1;
  }
  
  if (printFlag==1) cout<<"Inside IsVarnerEvent, peak2peakNadir/Peak2PeakTop: "<<peak2peakNadirMax/peak2peakTopMax<<endl;
  if (peak2peakTopMax<(peak2peakNadirMax/2) && fHeadPtr->trigType&(1<<0)) varnerFlag=1;

  return varnerFlag;
}
//////////////////////////
int MyCorrelator::isSyncSlip(int myEventNumber)
{
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (fHeadPtr->surfSlipFlag) return fHeadPtr->surfSlipFlag;
  if (fHeadPtr->errorFlag) return fHeadPtr->errorFlag;
  else return 0;
}
////////////////////
int MyCorrelator::isNadirRFCMOn(int eventNumber)
{
  int nantennas=NUM_ANTS_WITH_NADIRS;
  int myEventNumber=eventNumber;
  double peak2peak;
  
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (!(fHeadPtr->calibStatus & (1<<4))) return 0;
  
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  
  for (int ant=32;ant<nantennas;ant++){
    TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
    peak2peak=getPeak2Peak(gr1);
    delete gr1;
    if (peak2peak<20) return 0;
  }
  return 1;
}
////////////////////
int MyCorrelator::isDCOffsetLarge(int eventNumber)
{
  double mean;
  int myEventNumber=eventNumber;
  int nantennas=NUM_ANTS_WITH_NADIRS;
  
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  for (int ant=0;ant<nantennas;ant++){
    if (ant!=1){
      TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
      mean=gr1->GetMean(2);
      delete gr1;
      if (fabs(mean)>100) return 1;
    }
  }
  return 0;
}
//////////////////////////
int MyCorrelator::isMainRFCMOn(int eventNumber)
{
  //  int nantennas=NUM_ANTS_NO_NADIRS;
  int myEventNumber=eventNumber;
  // double peak2peak;
  
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (!(fHeadPtr->calibStatus & (1<<0)) || !(fHeadPtr->calibStatus & (1<<1)) || 
      !(fHeadPtr->calibStatus & (1<<2)) || !(fHeadPtr->calibStatus & (1<<3))){
    return 0;
  }  
  return 1;
}
///////////////////////
int MyCorrelator::isBigEnoughPeakToPeak(int eventNumber)
{
  int nantennas=NUM_ANTS_NO_NADIRS;
  int myEventNumber=eventNumber;
  double peak2peak;
  double maxpeak2peak=-1;
  
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  for (int ant=0;ant<nantennas;ant++){
    if (ant!=1){
      TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
      peak2peak=getPeak2Peak(gr1);
      if(peak2peak > maxpeak2peak) maxpeak2peak=peak2peak;
      delete gr1;
      //if (peak2peak<20) return 0;
    }
  }
  if(maxpeak2peak < 20) return 0;
  else return 1;
  
}
///////////////////
int MyCorrelator::isShortTrace(int myEventNumber)
{
  int nantennas=NUM_ANTS_WITH_NADIRS;
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  for (int ant=0;ant<nantennas;ant++){
    TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
   
    if (gr1->GetN()<240) return 1;
    delete gr1;
    gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kHorizontal);
    if (gr1->GetN()<240) return 1;
    delete gr1;
  } 
  return 0;
}
////////////////////////////////////////////////////////
int MyCorrelator::isTaylor(int eventNumber)
//gives 1 if event is a TD event (by a time cut)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  //get taylor dome location
  double x2,y2,z2; 
  getTDxyz(x2,y2,z2);
  
  double x1,y1,z1;
  LatLonAlt2xyz(fAdu5APatPtr->latitude*(-1.),fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,x1,y1,z1);
  double distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  double timeOfFlight=distance/(C_light);//in nanoseconds
  double delay=timeOfFlight-40e3;     
  int trig_time=fHeadPtr->triggerTimeNs;
  
  if (((trig_time-delay)<100 && (trig_time-delay)>-500) && distance<800e3){  
    if (printFlag==1) cout<<"Event "<<eventNumber<<" is a Taylor Dome Event."<<endl;
    return 1;
  }
  else return 0;

}
//////////////////////////
int MyCorrelator::isTaylorReflection(int eventNumber)
//gives 1 if event is a TD event (by a time cut)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  //get taylor dome location
  double x2,y2,z2; 
  getTDxyz(x2,y2,z2);
  
  double x1,y1,z1;
  LatLonAlt2xyz(fAdu5APatPtr->latitude*(-1.),fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,x1,y1,z1);
  double distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  double timeOfFlight=distance/(C_light);//in nanoseconds
  double delay=timeOfFlight-40e3;     
  int trig_time=fHeadPtr->triggerTimeNs;
  
  if (((trig_time-delay)<1000 && (trig_time-delay)>100) && distance<800e3){  
    if (printFlag==1) cout<<"Event "<<eventNumber<<" is a Taylor Dome Reflection Event."<<endl;
    return 1;
  }
  else return 0;

}
//////////////////////////////////////
int MyCorrelator::isMcMSeaveyByTime(int eventNumber)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  double x2,y2,z2; 
  getWillyxyz(x2,y2,z2);  
  double x1,y1,z1;
  LatLonAlt2xyz(fAdu5APatPtr->latitude*(-1.),fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,x1,y1,z1);
  double distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  
  int trig_time=fHeadPtr->triggerTimeNs;
  int delay=1200; //payload to mcm offset
  
  if ((fabs(trig_time+delay-50e6)<1e3 || fabs(trig_time+delay-250e6)<1e3 
       || fabs(trig_time+delay-450e6)<1e3 || fabs(trig_time+delay-650e6)<1e3) && distance<800e3){  
    if (printFlag==1) cout<<"Event "<<eventNumber<<" is a Willy Seavey Event."<<endl;
    return 1;
  }
  else return 0;

}
//////////////////////////////////////
int MyCorrelator::isMcMBoreholeByTime(int eventNumber)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  double x2,y2,z2; 
  getWillyxyz(x2,y2,z2);  
  double x1,y1,z1;
  LatLonAlt2xyz(fAdu5APatPtr->latitude*(-1.),fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,x1,y1,z1);
  double distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  
  int trig_time=fHeadPtr->triggerTimeNs;
  int delay=1200; //payload to mcm offset
  
  if ((fabs(trig_time+delay-150e6)<1e3 || fabs(trig_time+delay-350e6)<1e3 
       || fabs(trig_time+delay-550e6)<1e3 || fabs(trig_time+delay-750e6)<1e3) && distance<800e3){  
    if (printFlag==1) cout<<"Event "<<eventNumber<<" is a Willy Borehole Event."<<endl;
    return 1;
  }
  else return 0;
  
}
///////////////////////
int MyCorrelator::isMcMBoreholeOrSeaveyFromList(int eventNumber)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  unsigned int ns=fHeadPtr->triggerTimeNs;
  unsigned int realTime=fHeadPtr->realTime;
  int nentries;
  
  if (!tgroundPulser){
    char filename[150];
    sprintf(filename,"groundPulser/ground_pulser_arrival_times.root"); 
    TFile *fGroundPulserFile;
    fGroundPulserFile = TFile::Open(filename);
    if (fGroundPulserFile){
      TFile *rootfile=new TFile(filename,"READ");
      tgroundPulser=(TTree*) rootfile->Get("ground_pulser");
      nentries=tgroundPulser->GetEntries();
      cout<<"Number of entries in the ground pulser tree: "<<nentries<<endl;
      tgroundPulser->BuildIndex("arrivalTime","arrivalTimeNs");
      tgroundPulser->SetBranchAddress("arrivalTime",&realTimeForGroundPulser);
      tgroundPulser->SetBranchAddress("arrivalTimeNs",&nsForGroundPulser);
    }
  }
  
  nentries=tgroundPulser->GetEntries();
  int gpEntry = tgroundPulser->GetEntryNumberWithBestIndex(realTime,ns);
  tgroundPulser->GetEntry(gpEntry);

  if (realTimeForGroundPulser==realTime && fabs(int(nsForGroundPulser)-int(ns))<1000) return 1;
  else if (realTimeForGroundPulser==realTime+1 && 1e9-int(ns)+int(nsForGroundPulser)<1000
	   && 1e9-int(ns)+int(nsForGroundPulser)>0) return 1;
  else if (realTimeForGroundPulser==realTime-1 && 1e9-int(nsForGroundPulser)+int(ns)<1000
	   && 1e9-int(nsForGroundPulser)+int(ns)>0) return 1;
  else if (gpEntry<nentries-1){
    tgroundPulser->GetEntry(gpEntry+1);
    // cout<<"realTime: "<<realTime<<"realTimeForGP: "<<realTimeForGroundPulser<<", ns: "
    //<<ns<<", nsForGP: "<<nsForGroundPulser<<endl;
    
    if (realTimeForGroundPulser==realTime && fabs(int(nsForGroundPulser)-int(ns))<1000) return 1;
    else if (realTimeForGroundPulser==realTime+1 && 1e9-int(ns)+int(nsForGroundPulser)<1000
	     && 1e9-int(ns)+int(nsForGroundPulser)>0) return 1;
    else if (realTimeForGroundPulser==realTime-1 && 1e9-int(nsForGroundPulser)+int(ns)<1000
	     && 1e9-int(nsForGroundPulser)+int(ns)>0) return 1;
    else return 0;
  }
  else return 0;
  
  /*
  int ngpEvents=256968;
      if (allGPTimesFromList[0][0]==0){
      
      FILE *fp2 = fopen("/home/agoodhue/anita/clusterCode/groundPulser/all_ground_pulser_times.txt","r");
      int ncols=0;
      for (int i=0;i<ngpEvents;i++){
      ncols=fscanf(fp2,"%d.%d",&realTimeFromFile,&nsFromFile);
      allGPTimesFromList[i][0]=realTimeFromFile;
      allGPTimesFromList[i][1]=nsFromFile*10;
      }
      fclose(fp2);
      }
  
  for (int i=0;i<ngpEvents;i++){
    if (allGPTimesFromList[i][0]==realTime && fabs(allGPTimesFromList[i][1]-ns)<1000) return 1;
    if (allGPTimesFromList[i][0]==realTime+1 && allGPTimesFromList[i][1]+(1e9-ns)<1000 && allGPTimesFromList[i][1]+(1e9-ns)>0) return 1;
    if (allGPTimesFromList[i][0]==realTime-1 && (1e9-allGPTimesFromList[i][1])+ns<1000 && (1e9-allGPTimesFromList[i][1])+ns>0) return 1;
  }
  
  return 0;
  */
}
////////////////////////////////////////
int MyCorrelator::isCalPulser(int eventNumber)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  
  int trig_time=fHeadPtr->triggerTimeNs;  
  if ((fHeadPtr->calibStatus & (1<<6)) && trig_time<300){
    
    if (printFlag==1) cout<<"Event "<<eventNumber<<" is a Cal Pulser Event."<<endl;
    return 1;
  }
  
  else return 0;
  
}
/////////////////////////////////////

////////////////////////////////////////
int MyCorrelator::isTaylor(int eventNumber, double &distance)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  //get taylor dome location
  double x2,y2,z2; 
  getTDxyz(x2,y2,z2);
  
  double x1,y1,z1;
  LatLonAlt2xyz(fAdu5APatPtr->latitude*(-1.),fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,x1,y1,z1);
  distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  double timeOfFlight=distance/(C_light);//in nanoseconds
  double delay=timeOfFlight-40e3;     
  int trig_time=fHeadPtr->triggerTimeNs;
  
  if (((trig_time-delay)<100 && (trig_time-delay)>-500) && distance<800e3){  
    // cout<<"Event "<<eventNumber<<" is a Taylor Dome Event."<<endl;
    return 1;
  }
  else return 0;
}
////////////////////////////////////////
int MyCorrelator::isTaylor(int eventNumber, double &distance, double &deltaT)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  //get taylor dome location
  double x2,y2,z2; 
  getTDxyz(x2,y2,z2);
  
  double x1,y1,z1;
  LatLonAlt2xyz(fAdu5APatPtr->latitude*(-1.),fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,x1,y1,z1);
  distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  double timeOfFlight=distance/(C_light);//in nanoseconds
  double delay=timeOfFlight-40e3;     
  int trig_time=fHeadPtr->triggerTimeNs;
  deltaT=trig_time-delay;

  if (((trig_time-delay)<100 && (trig_time-delay)>-500) && distance<800e3){  
    // cout<<"Event "<<eventNumber<<" is a Taylor Dome Event."<<endl;
    return 1;
  }
  else return 0;
}
///////////////////////////////
void MyCorrelator::xyz2LatLonAlt(double &lat, double &lon, double &alt, double x, double y, double z)
{
  double a2=6378137.0;//radius of earth in meters
  double f=1./298.257223563;//flattening factor
  double b2=a2-a2*f;
  double p2=sqrt(x*x+y*y);
  double th=atan(z*a2/(p2*b2));
  double eprime=(a2*a2-b2*b2)/(b2*b2);
  double epsilon=2*f-f*f;

  lon=atan2(y,x);
  lat=atan((z+eprime*b2*sin(th)*sin(th)*sin(th))/(p2-epsilon*a2*cos(th)*cos(th)*cos(th)));
  double C2=pow(cos(lat)*cos(lat)+(1-f)*(1-f)*sin(lat)*sin(lat),-0.5);
  alt=p2/cos(lat)-a2*C2;

  lat*=rad2deg;
  lon*=rad2deg;

}
/////////////////////////////
void MyCorrelator::getTDxyz(double &x, double &y, double &z)
//gets the x,y,z coordinates of Taylor Dome
{
  double lat2=latTaylor;
  double lon2=lonTaylor;
  double height2=heightTaylor;
  LatLonAlt2xyz(lat2,lon2,height2,x,y,z); 

}
///////////////////////////////////////////////////////
void MyCorrelator::getWillyxyz(double &x, double &y, double &z)
//gets the x,y,z coordinates of Taylor Dome
{
  double lat2=latWilly;
  double lon2=lonWilly;
  double height2=heightWilly;
  LatLonAlt2xyz(lat2,lon2,height2,x,y,z); 

}

/////////////////////
double MyCorrelator::getMcMDistance(int eventNumber)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  //get taylor dome location
  double x2,y2,z2, distance; 
  getWillyxyz(x2,y2,z2);
  
  double x1,y1,z1;
  LatLonAlt2xyz(fAdu5APatPtr->latitude*(-1.),fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,x1,y1,z1);
  distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  
  return distance;
  
}
/////////////////////////////////////////////////////////////////////BEGIN STUFF FOR POINTING TO TD
void MyCorrelator::LatLonAlt2xyz(double lat, double lon, double alt, double &x, double &y, double &z)
//Calculates x,y,z given a latitude and longitude.
{
  lat*=deg2rad;
  lon*=deg2rad;
  double r=6378137.0;//radius of earth in meters
  double f=1./298.257223563;//flattening factor
  //calculate x,y,z coordinates
  double C2=pow(cos(lat)*cos(lat)+(1-f)*(1-f)*sin(lat)*sin(lat),-0.5);
  double Q2=(1-f)*(1-f)*C2;
  x=(r*C2+alt)*cos(lat)*cos(lon);
  y=(r*C2+alt)*cos(lat)*sin(lon);
  z=(r*Q2+alt)*sin(lat);
}
////////////////////////////////
double MyCorrelator::findDistanceBetweenTwoThings(double lat1, double lon1, 
						  double alt1, double lat2, double lon2, double alt2)
{
  double x1,y1,z1;
  double x2,y2,z2, distance; 
  LatLonAlt2xyz(lat1,lon1,alt1,x1,y1,z1);
  LatLonAlt2xyz(lat2,lon2,alt2,x2,y2,z2);
  distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  
  return distance;

}
///////////////////////////
void MyCorrelator::getRelXYFromLatLong(double latitude, double longitude, double &x, double &y)
{
    //Negative longitude is west
  //    //All latitudes assumed south
  double altitude=0;
  double z=0;
  LatLonAlt2xyz(latitude, longitude, altitude, x,y,z);

  double temp=x;
  x=y;
  y=temp;
  
  y*=scale;
  y+=yOffset;
  y/=ySize;
  x*=scale;
  x+=xOffset;
  x/=xSize;

    
}
////////////////////////////////////
void MyCorrelator::getClosestNAntennas(int nantennasToUse, double peakPhi, vector<int>& whichAntennasToUse, 
				       int nadirFlag)
{

  int nantennasTotal=NUM_ANTS_NO_NADIRS;
  if (nadirFlag==1) nantennasTotal=NUM_ANTS_WITH_NADIRS;
  if (nadirFlag==0) nantennasTotal=NUM_ANTS_NO_NADIRS;
  
  double deltaPhiArray[nantennasToUse];
  double deltaPhiThisOne;  
  int index;
  double phi_ant[nantennasTotal];

  
  for (int i=0;i<nantennasToUse;i++){
    deltaPhiArray[i]=360;//set array full of large values
  }
  
  for (int ant=0;ant<nantennasTotal;ant++){
    //cout<<"ant is "<<ant<<"\n";
    if (ant!=1 && saturatedChannels[ant+polToggle*NUM_ANTS_WITH_NADIRS]==0){
      phi_ant[ant]=fUPGeomTool->getAntPhiPositionRelToAftFore(ant);
      deltaPhiThisOne=phi_ant[ant]*rad2deg-peakPhi;
     
      if (deltaPhiThisOne>180.) deltaPhiThisOne-=360.;
      if (deltaPhiThisOne<-180.) deltaPhiThisOne+=360.;
      if (deltaPhiThisOne<0) deltaPhiThisOne*=-1.;
      
      if (deltaPhiThisOne<getMaximum(nantennasToUse, deltaPhiArray, index)){//find maximum of array/replace if new value is lower
	deltaPhiArray[index]=deltaPhiThisOne;
	whichAntennasToUse[index]=ant;
      }
      }//ant!=1
  }
  
}
/////////////////////////
void MyCorrelator::getGroupsofAntennas(int nAntennasToUse, int nadirFlag){
  vector<int> antenna_group (nAntennasToUse,0);

  int old_group_switch=0;
 
  for(int runoverphi=0;runoverphi<360;runoverphi++){//get all possible antenna groups
   
    old_group_switch=0;
    getClosestNAntennas(nAntennasToUse, runoverphi, antenna_group, nadirFlag);   
    sort(antenna_group.begin(),antenna_group.end());
    
    for(int i0=0;i0<(int)antenna_group_holder.size();i0++){
      if(antenna_group==antenna_group_holder[i0]){
	old_group_switch=1;
	break;
      }//check same  
    }//i0 
    
    
    if(old_group_switch==0){
      antenna_group_holder.push_back(antenna_group);
      unique_phis.push_back(runoverphi);
     
    }
   
  }
  
  antenna_groups_size=(int) antenna_group_holder.size();

}

////////////////////////////////////

double MyCorrelator::getHeadingOfEvent(int eventNumber, double peakPhiMapDegrees)
{
  if(!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  double headingOfEvent=fAdu5APatPtr->heading-peakPhiMapDegrees;
  if (headingOfEvent>=360) headingOfEvent-=360;
  if (headingOfEvent<0) headingOfEvent+=360;
  return headingOfEvent;
}

///////////////////////////////////////////////////////////////////////////////////
void MyCorrelator::getTriggeredPhi(RawAnitaHeader *hdPtr,int triggeredPhi[NUM_PHI_SECTORS]){
  for(int phi=0;phi<NUM_PHI_SECTORS;phi++){
    if(hdPtr->l3TrigPattern & (1 << phi)) triggeredPhi[phi]=1;
    else triggeredPhi[phi]=0;
  }
}
///////////////////////////
void MyCorrelator::getTriggeredL2Phi(RawAnitaHeader *hdPtr,int triggeredPhi[NUM_PHI_SECTORS]){
  for(int phi=0;phi<NUM_PHI_SECTORS;phi++){
    if(hdPtr->upperL2TrigPattern & (1 << phi)) triggeredPhi[phi]=1;
    else triggeredPhi[phi]=0;
    if (hdPtr->lowerL2TrigPattern & (1 << phi)) triggeredPhi[phi]=1;
    if (hdPtr->nadirL2TrigPattern & (1 << phi)) triggeredPhi[phi]=1;
    if (hdPtr->l3TrigPattern & (1 << phi)) triggeredPhi[phi]=1;
 
  }
}
////////////////////////////////////////////////////////////////////////////////////
void MyCorrelator::getTriggeredAnt(int triggeredPhi[NUM_PHI_SECTORS],int triggeredAnt[NUM_ANTS_WITH_NADIRS]){
  int phi;
  for(int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    phi=fUPGeomTool->getPhiFromAnt(ant);
    if(triggeredPhi[phi]) triggeredAnt[ant]=1;
    else triggeredAnt[ant]=0;
  }
}
///////////////////////////////////////////////////////////////
void MyCorrelator::getTriggeredAntpm1PhiSector(int triggeredPhi[NUM_PHI_SECTORS],int triggeredAnt[NUM_ANTS_WITH_NADIRS])
{
  int phi;
  int triggeredPhiDummy[NUM_PHI_SECTORS];
  for(int i=0;i<NUM_PHI_SECTORS;i++){
    triggeredPhiDummy[i]=0;
  }
  for(int i=0;i<NUM_PHI_SECTORS;i++){
   
    int left=i-1;
    int right=i+1;
    if (left<0) left+=NUM_PHI_SECTORS;
    if (right>15) right-=NUM_PHI_SECTORS;
    if (triggeredPhi[i]){ 
      triggeredPhiDummy[i]=1;
      triggeredPhiDummy[left]=1;
      triggeredPhiDummy[right]=1;
    } 
  }
  for(int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    phi=fUPGeomTool->getPhiFromAnt(ant);
    if(triggeredPhiDummy[phi])  triggeredAnt[ant]=1;
    else triggeredAnt[ant]=0;
  }
  for(int i=0;i<NUM_PHI_SECTORS;i++){
    triggeredPhi[i]=triggeredPhiDummy[i];
    
  }

}
/////////////////////////////////////////////////////////////
void MyCorrelator::getTriggeredAntOf3PhiSectors(int triggeredPhi[NUM_PHI_SECTORS],int triggeredAnt[NUM_ANTS_WITH_NADIRS]){
  int phi;
  for (int phiSector=0;phiSector<NUM_PHI_SECTORS;phiSector++){
    if (phiSector!=0 && phiSector!=15 &&
	triggeredPhi[phiSector]==1 && triggeredPhi[phiSector-1]==0 && triggeredPhi[phiSector+1]==0){
      triggeredPhi[phiSector+1]=1;
      triggeredPhi[phiSector-1]=1;
    }
    if (phiSector==0 && triggeredPhi[phiSector]==1 && triggeredPhi[15]==0 && triggeredPhi[1]==0){
      triggeredPhi[1]=1;
      triggeredPhi[15]=1;
    }
    if (phiSector==15 && triggeredPhi[phiSector]==1 && triggeredPhi[14]==0 && triggeredPhi[0]==0){
      triggeredPhi[0]=1;
      triggeredPhi[14]=1;
    }
  }
  
  for(int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    phi=fUPGeomTool->getPhiFromAnt(ant);
    if(triggeredPhi[phi]) triggeredAnt[ant]=1;
    else triggeredAnt[ant]=0;
  }
}
/////////////////////////////////////
int MyCorrelator::allowedPhisPairOfAntennas(double &lowerAngle, double &higherAngle, double &centerTheta1, double &centerTheta2, double &centerPhi1, double &centerPhi2, int ant1, int ant2)
{

  int phi1=fUPGeomTool->getPhiFromAnt(ant1);
  int phi2=fUPGeomTool->getPhiFromAnt(ant2);
  int allowedFlag=0;
  
  int upperlimit=phi2+2;//2 phi sectors on either side
  int lowerlimit=phi2-2;

  if(upperlimit>NUM_PHI_SECTORS-1)upperlimit-=NUM_PHI_SECTORS;
  if(lowerlimit<0)lowerlimit+=NUM_PHI_SECTORS;

  if (upperlimit>lowerlimit){
    if (phi1<=upperlimit && phi1>=lowerlimit){//within 2 phi sectors of eachother
      allowedFlag=1;
    }
  }
  if (upperlimit<lowerlimit){
    if (phi1<=upperlimit || phi1>=lowerlimit){
      allowedFlag=1;

    }
  }
  
  double centerAngle1, centerAngle2;
  if (allowedFlag==1){
    centerAngle1=phi1*PHI_SECTOR_ANGLE-ADU5_FORE_PHI;
    centerAngle2=phi2*PHI_SECTOR_ANGLE-ADU5_FORE_PHI;

    // cout<<"centerAngle1 for ant "<<ant1<<" is "<<centerAngle1<<"\n";
   

    if (centerAngle2>centerAngle1){//NUM_DEGREES =75
      lowerAngle=centerAngle2-NUM_DEGREES_OFF_CENTER;//lowest angle both antennas can see
      higherAngle=centerAngle1+NUM_DEGREES_OFF_CENTER;//highest angle both antennas can see
    }
    if (centerAngle1>centerAngle2){
      lowerAngle=centerAngle1-NUM_DEGREES_OFF_CENTER;
      higherAngle=centerAngle2+NUM_DEGREES_OFF_CENTER; 
    }
    if (lowerAngle<0) lowerAngle+=360;
    if (higherAngle>360) higherAngle-=360;
    
  }
  centerTheta1=10;//degrees down
  centerTheta2=10;//degrees down
  centerPhi1=centerAngle1;
  centerPhi2=centerAngle2;
  
  return allowedFlag;

}
////////////////////////////////
int MyCorrelator::isPhiMaskingOn(int eventNumber, int phiMaskedArray[NUM_PHI_SECTORS])
{  
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  int isPhiMask=0;

  for (int i=0;i<NUM_PHI_SECTORS;i++){
    int thisPhiMask=(fHeadPtr->phiTrigMask & (1<<i));
    if (thisPhiMask>0 && isPhiMask==0) isPhiMask=1;
    phiMaskedArray[i]=thisPhiMask;
  }
  return isPhiMask;
}
////////////////////////////////////////////////////////////////////////////////////
void MyCorrelator::getPhiMaskedAnts(int phiMaskArray[NUM_PHI_SECTORS],int phiMaskedAnts[NUM_ANTS_WITH_NADIRS]){
  int phi;
  for(int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    phi=fUPGeomTool->getPhiFromAnt(ant);
    if(phiMaskArray[phi]) phiMaskedAnts[ant]=1;
    else phiMaskedAnts[ant]=0;
  }
}
/////////////////////////////
void MyCorrelator::getPhiMaskedAntspm1(int phiMaskArray[NUM_PHI_SECTORS],int phiMaskedAnts[NUM_ANTS_WITH_NADIRS])
{
  
  int phi;
  int phiMaskedDummy[NUM_PHI_SECTORS];
  for(int i=0;i<NUM_PHI_SECTORS;i++){
    phiMaskedDummy[i]=0;
  }
  for(int i=0;i<NUM_PHI_SECTORS;i++){
    int left=i-1;
    int right=i+1;
    if (left<0) left+=NUM_PHI_SECTORS;
    if (right>15) right-=NUM_PHI_SECTORS;
    if (phiMaskArray[i]){ 
      phiMaskedDummy[i]=1;
      phiMaskedDummy[left]=1;
      phiMaskedDummy[right]=1;
    } 
  }
  for(int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    phi=fUPGeomTool->getPhiFromAnt(ant);
    if(phiMaskedDummy[phi])  phiMaskedAnts[ant]=1;
    else phiMaskedAnts[ant]=0;
  }
  for(int i=0;i<NUM_PHI_SECTORS;i++){
    phiMaskArray[i]=phiMaskedDummy[i]; 
  }
}
//////////////////////////////////
int MyCorrelator::getBestPhiSector(double phiWaveRadians)
{
  double phiWaveDeg=phiWaveRadians*rad2deg;
  int phiSector=int((phiWaveDeg+22.5+11.25)/22.5);
  if (phiSector>15) phiSector-=16;
  if (phiSector<0) phiSector+=16;
  
  return phiSector;
}
///////////////////////
float MyCorrelator::getHWTriggerAngle(int eventNumber)
{
  int triggeredPhi[NUM_PHI_SECTORS];
  float currentAngle;
  float lastIterationAngle;
  float iterationAngle=0;
  int avgAngleCtr=0;
  
  if (!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  
  if((fHeadPtr->trigType&(1<<0))){
    //if(!(fHeadPtr->trigType&(1<<3)) && !(fHeadPtr->trigType&(1<<2)) && !(fHeadPtr->trigType&(1<<1))){
    getTriggeredPhi(fHeadPtr,triggeredPhi);  
    for (int i=0;i<NUM_PHI_SECTORS;i++){
      if (triggeredPhi[i]){
	currentAngle=i*22.5-22.5;
	if (currentAngle<0) currentAngle+=360;
	if (avgAngleCtr==0){
	  lastIterationAngle=currentAngle;
	  iterationAngle=currentAngle;
	}
	else{
	    if (fabs(currentAngle-lastIterationAngle)<=180) 
	      iterationAngle=(currentAngle+avgAngleCtr*lastIterationAngle)/(1+avgAngleCtr);
	    if ((currentAngle-lastIterationAngle)>180) 
	      iterationAngle=(currentAngle-360+avgAngleCtr*lastIterationAngle)/(1+avgAngleCtr);
	    if ((currentAngle-lastIterationAngle)<-180) 
	      iterationAngle=(currentAngle+avgAngleCtr*(lastIterationAngle-360))/(1+avgAngleCtr);
	    
	    lastIterationAngle=iterationAngle;
	}
	avgAngleCtr++;
      }
  
    }
    if (avgAngleCtr!=0){
      if (iterationAngle>=0) return iterationAngle;
      else return iterationAngle+360;
    }
    else{
      getTriggeredL2Phi(fHeadPtr,triggeredPhi);
      for (int i=0;i<NUM_PHI_SECTORS;i++){
	if (triggeredPhi[i]){
	  currentAngle=i*22.5-22.5;
	  if (currentAngle<0) currentAngle+=360;
	  if (avgAngleCtr==0) lastIterationAngle=currentAngle;
	  else{
	    if (fabs(currentAngle-lastIterationAngle)<=180) 
	      iterationAngle=(currentAngle+avgAngleCtr*lastIterationAngle)/(1+avgAngleCtr);
	    if ((currentAngle-lastIterationAngle)>180) 
	      iterationAngle=(currentAngle-360+avgAngleCtr*lastIterationAngle)/(1+avgAngleCtr);
	    if ((currentAngle-lastIterationAngle)<-180) 
	      iterationAngle=(currentAngle+avgAngleCtr*(lastIterationAngle-360))/(1+avgAngleCtr);
	    
	    lastIterationAngle=iterationAngle;
	  }
	  avgAngleCtr++;
	}
	
      }
      if (avgAngleCtr!=0){
	if (iterationAngle>=0) return iterationAngle;
	else return iterationAngle+360;
      }
      else{
	return -1;
      }
    }
  }
  else return -1;
}
//////////////////////
int MyCorrelator::isPeakPhiTriggeredOrMasked(RawAnitaHeader *hdPtr, double peakPhiInterp)
{  
  int triggeredPhi[NUM_PHI_SECTORS];
  int triggeredAnt[NUM_ANTS_WITH_NADIRS];
  int phiMaskArray[NUM_PHI_SECTORS];
  int phiMaskedAnts[NUM_ANTS_WITH_NADIRS]; 
  int bestPhiSector;
  int phiOK;

  getTriggeredL2Phi(hdPtr,triggeredPhi);
  getTriggeredAntpm1PhiSector(triggeredPhi,triggeredAnt);
  int thisPhiMask=isPhiMaskingOn(hdPtr->eventNumber,phiMaskArray); 
  thisPhiMask=0;
  getPhiMaskedAntspm1(phiMaskArray,phiMaskedAnts);
  
  for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    if (phiMaskedAnts[ant] || triggeredAnt[ant]) triggeredAnt[ant]=1;
  }
  bestPhiSector=getBestPhiSector(peakPhiInterp*deg2rad);
  if (triggeredPhi[bestPhiSector]|| phiMaskArray[bestPhiSector]) phiOK=1;
  else phiOK=0;

  return phiOK;
}
////////////////////////
int MyCorrelator::isPeakPhiMasked(RawAnitaHeader *hdPtr, double peakPhiInterp)
{
  int phiMaskArray[NUM_PHI_SECTORS];
  int bestPhiSector;
  int phiOK;
  
  int thisPhiMask=isPhiMaskingOn(hdPtr->eventNumber,phiMaskArray); 
  thisPhiMask=0;
  bestPhiSector=getBestPhiSector(peakPhiInterp*deg2rad);
  if (phiMaskArray[bestPhiSector]) phiOK=1;
  else phiOK=0;

  return phiOK;
}
//////////////////////
int MyCorrelator::isPeakTriggeredpm1(RawAnitaHeader *hdPtr, double peakPhiInterp)
{
  int triggeredPhi[NUM_PHI_SECTORS];
  int triggeredAnt[NUM_ANTS_WITH_NADIRS];
  int bestPhiSector;
  int phiOK;

  getTriggeredL2Phi(hdPtr,triggeredPhi);
  getTriggeredAntpm1PhiSector(triggeredPhi,triggeredAnt);
  
  bestPhiSector=getBestPhiSector(peakPhiInterp*deg2rad);
  if (triggeredPhi[bestPhiSector]) phiOK=1;
  else phiOK=0;
  
  return phiOK;
}
//////////////////////////////
void MyCorrelator::adaptiveFilterCoherent(TGraph *grCoherent,int pol, double dBCut, int nfreq, 
					  double *frequencies, int drawFlag, double bandWidth)// double &mean)
{
  if (readBaselineFlag!=1){ 
    readBaselineFFTs();
    if (printFlag==1) cout<<"Read in Baseline For Adaptive Filter"<<endl;
  }
  
  double magFFT[2000];
  double frequencyArray[2000];
  double deltaT, deltaF;
  int length;
  int newLength;
 
  for (int i=0;i<2000;i++){
    magFFT[i]=0;
    frequencyArray[i]=-1;
  }
  
  // for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){      
  // if (pol!=0 || ant!=1){//get rid of 2V
  double *Y = grCoherent->GetY();
  double *X = grCoherent->GetX();
  deltaT=X[1]-X[0];
  length=grCoherent->GetN();
 
  newLength=(length/2)+1;
      deltaF=1/(deltaT*length); //Hz
      deltaF*=1e3; //MHz   
   
      FFTWComplex *theFFT=FFTtools::doFFT(length,Y);
      
      for(int i=0;i<newLength;i++) {
	//if (ant==0){
	if (i==0) frequencyArray[i]=0;
	if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;
	//	}
	if (frequencyArray[i]>=200 && frequencyArray[i]<=1200)
	  magFFT[i]=theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im;
	else magFFT[i]=-1000;    

      }    
      delete [] theFFT;
      //  }
      // }

  for(int i=0;i<newLength;i++) {
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      magFFT[i]=10*log10(sqrt(magFFT[i])/10.);///double(NUM_ANTS_WITH_NADIRS-1)/10.);
     
      //if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i])/double(NUM_ANTS_WITH_NADIRS)/10.);
    }    
    else magFFT[i]=-1000;
    // cout<<"aaaa  "<<magFFT[i]<<endl;
  }
  
  //get baseline average so we can bump FFT around
  double meanBaseline=0;
  double mean=0;
  double firstHalfMeanBaseline=0;
  double secondHalfMeanBaseline=0;
  double firstHalfMean=0;
  double secondHalfMean=0;
  int navg=0;
  int nfirsthalfavg=0;
  int nsecondhalfavg=0;

  double *bX, *bY;
  int n;

  if (pol==0){//vertical
    bY=grVertCoherentBaseline->GetY();
    bX=grVertCoherentBaseline->GetX();
    n=grVertCoherentBaseline->GetN();
  }
  if (pol==1){ //horizontal
    bY=grHorizCoherentBaseline->GetY();
    bX=grHorizCoherentBaseline->GetX();
    n=grHorizCoherentBaseline->GetN();
  }
  
  for (int i=0;i<n;i++){
    if (bX[i]>=200 && bX[i]<1200){ 
      meanBaseline+=bY[i];
      navg++;
    }
    if (bX[i]>=200 && bX[i]<700){ 
      firstHalfMeanBaseline+=bY[i];
      nfirsthalfavg++;
    }
    if (bX[i]>=700 && bX[i]<1200){ 
      secondHalfMeanBaseline+=bY[i];
      nsecondhalfavg++;
    }	
  }
  meanBaseline=meanBaseline/double(navg);
  firstHalfMeanBaseline=firstHalfMeanBaseline/double(nfirsthalfavg);
  secondHalfMeanBaseline=secondHalfMeanBaseline/double(nsecondhalfavg);
  
  navg=0;
  nfirsthalfavg=0;
  nsecondhalfavg=0;
  
  //get average of graph in question
  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>=200 && frequencyArray[i]<1200){ 
      mean+=magFFT[i];
      navg++;
    }
    if (frequencyArray[i]>=200 && frequencyArray[i]<700){
      firstHalfMean+=magFFT[i];
      nfirsthalfavg++;
    }
    if (frequencyArray[i]>=700 && frequencyArray[i]<1200){
      secondHalfMean+=magFFT[i];
      nsecondhalfavg++;
    }

  } 
  mean=mean/double(navg);
  firstHalfMean=firstHalfMean/double(nfirsthalfavg);
  secondHalfMean=secondHalfMean/double(nsecondhalfavg);
  
  //now bump the average to the baseline average and apply a tilt correction to baseline
  double deltaMean=mean-meanBaseline;
  double deltaMeanFirst=firstHalfMean-firstHalfMeanBaseline-deltaMean;
  double deltaMeanSecond=secondHalfMean-secondHalfMeanBaseline-deltaMean;
  double slope=(deltaMeanFirst-deltaMeanSecond)/500.;
  //cout<<deltaMean<<", "<<deltaMeanFirst<<", "<<deltaMeanSecond<<" , "<<slope<<", "<<newLength<<endl;
  
  for (int i=0;i<newLength;i++){
    magFFT[i]=magFFT[i]-deltaMean;
  }
  for (int ctr=0;ctr<n;ctr++){
    if (bX[ctr]>=200 && bX[ctr]<1200){
      bY[ctr]=bY[ctr]-slope*(bX[ctr]-700.);
      //magFFT[i]=magFFT[i]+slope*(frequencyArray[i]-700.);
    }
    //  cout<<magFFT[i]<<endl;
  }
  
  //now make the graph
  /*
  if (drawFlag==1){
    TGraph *grCoherentMeanAdjusted=new TGraph(newLength,frequencyArray,magFFT);
    TGraph *grbaselinehere=new TGraph(n,bX,bY);
    
    TH2F *haxes=new TH2F("haxes","haxes",10,0,1400,10,-20,40);
    gStyle->SetOptStat(kFALSE);
    TCanvas *cBaseline3=new TCanvas("cBaseline3","cBaseline3",800,800);
    cBaseline3->cd(0);  
    haxes->Draw();
    grbaselinehere->Draw("l");
    grCoherentMeanAdjusted->Draw("l+");
  }
  */
  //now see if any peaks are ndB above the baseline.
  double deltaMag[newLength];

  int j;
  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>210 && frequencyArray[i]<1190){
      for (j=0;j<n;j++){
	if (bX[j]>frequencyArray[i]) break;
      }
      deltaMag[i]=magFFT[i]-bY[j];
    }
    else deltaMag[i]=-1000;
  }
  
  for (int i=0;i<nfreq;i++){
    int index;
    double maxDelta=getMaximum(newLength, deltaMag, index);
    if (maxDelta>dBCut){
      frequencies[i]=frequencyArray[index];
    }
    else frequencies[i]=-1;
    for (int i=0;i<newLength;i++){
      if (fabs(frequencyArray[i]-frequencyArray[index])<=bandWidth)
	deltaMag[i]=-1000;
    }
  }
  
}
  /* double *Y = grCoherent->GetY();
  double *X = grCoherent->GetX();
  double deltaT=X[1]-X[0];
  int length=grCoherent->GetN();
  FFTWComplex *theFFT=FFTtools::doFFT(length,Y);
  int newLength=(length/2)+1;
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  double magFFT[newLength];
  double frequencyArray[newLength];

  for(int i=0;i<newLength;i++) {
    if (i==0) frequencyArray[i]=0;
    if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200)
      magFFT[i]=10*log10(sqrt(theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im)/10.); 
    else magFFT[i]=-1000;    
  }
  
  //get baseline average so we can bump FFT around
  double meanBaseline=0;
  double mean=0;
  double firstHalfMeanBaseline=0;
  double secondHalfMeanBaseline=0;
  double firstHalfMean=0;
  double secondHalfMean=0;
  int navg=0;
  int nfirsthalfavg=0;
  int nsecondhalfavg=0;

  double *bX, *bY;
  int n;

  if (pol==0){//vertical
    bY=grVertCoherentBaseline->GetY();
    bX=grVertCoherentBaseline->GetX();
    n=grVertCoherentBaseline->GetN();
  }
  if (pol==1){ //horizontal
    bY=grHorizCoherentBaseline->GetY();
    bX=grHorizCoherentBaseline->GetX();
    n=grHorizCoherentBaseline->GetN();
  }
  
  for (int i=0;i<n;i++){
    if (bX[i]>=200 && bX[i]<1200){ 
      meanBaseline+=bY[i];
      navg++;
    }
    if (bX[i]>=200 && bX[i]<700){ 
      firstHalfMeanBaseline+=bY[i];
      nfirsthalfavg++;
    }
    if (bX[i]>=700 && bX[i]<1200){ 
      secondHalfMeanBaseline+=bY[i];
      nsecondhalfavg++;
    }	
  }
  meanBaseline=meanBaseline/double(navg);
  firstHalfMeanBaseline=firstHalfMeanBaseline/double(nfirsthalfavg);
  secondHalfMeanBaseline=secondHalfMeanBaseline/double(nsecondhalfavg);
  
  navg=0;
  nfirsthalfavg=0;
  nsecondhalfavg=0;
  
  //get average of graph in question

  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>=200 && frequencyArray[i]<1200){ 
      mean+=magFFT[i];
      navg++;
    }
    if (frequencyArray[i]>=200 && frequencyArray[i]<700){
      firstHalfMean+=magFFT[i];
      nfirsthalfavg++;
    }
    if (frequencyArray[i]>=700 && frequencyArray[i]<1200){
      secondHalfMean+=magFFT[i];
      nsecondhalfavg++;
    }

  } 
  mean=mean/double(navg);
  firstHalfMean=firstHalfMean/double(nfirsthalfavg);
  secondHalfMean=secondHalfMean/double(nsecondhalfavg);
  
  //now bump the average to the baseline average and apply a tilt correction to baseline
  double deltaMean=mean-meanBaseline;
  double deltaMeanFirst=firstHalfMean-firstHalfMeanBaseline-deltaMean;
  double deltaMeanSecond=secondHalfMean-secondHalfMeanBaseline-deltaMean;
  double slope=(deltaMeanFirst-deltaMeanSecond)/500.;
  
  for (int i=0;i<newLength;i++){
    magFFT[i]=magFFT[i]-deltaMean;
    if (frequencyArray[i]>=200 && frequencyArray[i]<1200){
      magFFT[i]=magFFT[i]+slope*(frequencyArray[i]-700.);
    }
  }

  //now make the graph
  TGraph *grCoherentMeanAdjusted=new TGraph(newLength,frequencyArray,magFFT);
  
  if (drawFlag==1){
    TH2F *haxes=new TH2F("haxes","haxes",10,0,1400,10,-30,20);
    gStyle->SetOptStat(kFALSE);
    TCanvas *cBaseline3=new TCanvas("cBaseline3","cBaseline3",800,800);
    cBaseline3->cd(0);  
    haxes->Draw();
    grVertCoherentBaseline->Draw("l");
    grCoherentMeanAdjusted->Draw("l+");
  }
  //now see if any peaks are ndB above the baseline.

  double deltaMag[newLength];

  int j;
  for (int i=0;i<newLength;i++){
    for (j=0;j<n;j++){
      if (bX[j]>frequencyArray[i]) break;
    }
    deltaMag[i]=magFFT[i]-bY[j];
  }
  
  int index;
  double maxDelta=getMaximum(newLength, deltaMag, index);
  delete [] theFFT;
  
  if (maxDelta>dBCut){
    return frequencyArray[index];
  }
  else return -1;
  */

///////////////////////////////////////

void MyCorrelator::adaptiveFilter(int pol, double dBCut, int nfreq, double *frequencies, int drawFlag, double bandWidth, float &mean_freq, double *freqArray, double *FFTarray, double *baseX, double *baseY, int baseN, int myEventNumber)// double &mean)
{
  if (readBaselineFlag!=1){ 
    readBaselineFFTs();
    if (printFlag==1) cout<<"Read in Baseline For Adaptive Filter"<<endl;
  }  
  
  double magFFT[2000];
  double frequencyArray[2000];
  double deltaT, deltaF;
  int length;
  int newLength;
  double *X, *Y;
  
  for (int i=0;i<2000;i++){
    magFFT[i]=0;
    frequencyArray[i]=-1;
  }
 
  for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){      
    if ((pol!=0 || ant!=1) && saturatedChannels[ant+pol*NUM_ANTS_WITH_NADIRS]==0){//get rid of 2V
      if (pol==0) Y = grEv[ant]->GetY();//set X and Y values
      else{ Y = grEvHoriz[ant]->GetY(); }
      if (pol==0) X = grEv[ant]->GetX();
      else{ X = grEvHoriz[ant]->GetX(); }
      
      deltaT=X[1]-X[0];
      if (pol==0) length=grEv[ant]->GetN();
      else{length=grEvHoriz[ant]->GetN();}
      //length=256;
      //cout<<"length is "<<length<<"\n";
      newLength=(length/2)+1;
      deltaF=1/(deltaT*length); //GHz
      deltaF*=1e3; //MHz
      
      //cout<<"deltaF is "<<deltaF<<"\n";
     
      FFTWComplex *theFFT=FFTtools::doFFT(length,Y);
      
     
      
      for(int i=0;i<newLength;i++) {
	if (ant==0){
	  if (i==0) frequencyArray[i]=0;
	  if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;//make a freq Array. DOES THIS CHANGE FOR EVERY ANT GRAPH??
	 
	}
	if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
	  magFFT[i]+=theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im;//adding squares of all antennas together
	   //magFFT[i]+=sqrt(theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im);
	}
	else magFFT[i]=-1000;    
	
      }   
      delete [] theFFT;
      
    }//pol!=0
  }//ant==0
 
  for(int i=0;i<newLength;i++) {
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      if (pol==0) magFFT[i]=10*log10(sqrt(magFFT[i]/double(NUM_ANTS_WITH_NADIRS-1))/10.);//correct RMS
      if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i]/double(NUM_ANTS_WITH_NADIRS))/10.);//correct RMS
      // if (pol==0) magFFT[i]=10*log10(sqrt(magFFT[i])/double(NUM_ANTS_WITH_NADIRS-1)/10.);//old way
      //if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i])/double(NUM_ANTS_WITH_NADIRS)/10.);//old way
    }    
    else magFFT[i]=-1000;
    // cout<<"aaaa  "<<magFFT[i]<<endl;
  }
  for(int k1=0;k1<newLength;k1++){
    //cout<<"frequencyArray[k1] is "<<frequencyArray[k1]<<"\n";
    freqArray[k1]=frequencyArray[k1];
    FFTarray[k1]=magFFT[k1];
    //cout<<"freqArray is now "<<freqArray[k1]<<" and frequencyaRray is "<<frequencyArray[k1]<<"\n";
    //cout<<"fftarray["<<k1<<"] is now "<<FFTarray[k1]<<" and magFFT is "<<magFFT[k1]<<"\n";
  }
  //get baseline average so we can bump FFT around
  double meanBaseline=0;
  double mean=0;
  double firstHalfMeanBaseline=0;
  double secondHalfMeanBaseline=0;
  double firstHalfMean=0;
  double secondHalfMean=0;
  int navg=0;
  int nfirsthalfavg=0;
  int nsecondhalfavg=0;

  double *bX, *bY;
  vector<double> bX2;
  vector<double> bY2;
  int n;

  if (pol==0){//vertical
    bY=grVertCoherentBaseline->GetY();
    bX=grVertCoherentBaseline->GetX();
    n=grVertCoherentBaseline->GetN();
  }
  if (pol==1){ //horizontal
    bY=grHorizCoherentBaseline->GetY();
    bX=grHorizCoherentBaseline->GetX();
    n=grHorizCoherentBaseline->GetN();
  }
  
  for (int i=0;i<n;i++){
    if (bX[i]>=200 && bX[i]<1200){ 
      meanBaseline+=bY[i];
      navg++;
    }
    if (bX[i]>=200 && bX[i]<700){ 
      firstHalfMeanBaseline+=bY[i];
      nfirsthalfavg++;
    }
    if (bX[i]>=700 && bX[i]<1200){ 
      secondHalfMeanBaseline+=bY[i];
      nsecondhalfavg++;
    }	
  }
  meanBaseline=meanBaseline/double(navg);
  firstHalfMeanBaseline=firstHalfMeanBaseline/double(nfirsthalfavg);
  secondHalfMeanBaseline=secondHalfMeanBaseline/double(nsecondhalfavg);
  
  navg=0;
  nfirsthalfavg=0;
  nsecondhalfavg=0;
  
  //get average of graph in question
  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>=200 && frequencyArray[i]<1200){ 
      mean+=magFFT[i];
      navg++;
    }
    if (frequencyArray[i]>=200 && frequencyArray[i]<700){
      firstHalfMean+=magFFT[i];
      nfirsthalfavg++;
    }
    if (frequencyArray[i]>=700 && frequencyArray[i]<1200){
      secondHalfMean+=magFFT[i];
      nsecondhalfavg++;
    }

  } 
  mean=mean/double(navg);
  firstHalfMean=firstHalfMean/double(nfirsthalfavg);
  secondHalfMean=secondHalfMean/double(nsecondhalfavg);
  
  //now bump the average to the baseline average and apply a tilt correction to baseline
  double deltaMean=mean-meanBaseline;
  double deltaMeanFirst=firstHalfMean-firstHalfMeanBaseline-deltaMean;
  double deltaMeanSecond=secondHalfMean-secondHalfMeanBaseline-deltaMean;
  double slope=(deltaMeanFirst-deltaMeanSecond)/500.;
  //cout<<deltaMean<<", "<<deltaMeanFirst<<", "<<deltaMeanSecond<<" , "<<slope<<", "<<newLength<<endl;
  
  for (int i=0;i<newLength;i++){
    //cout<<"magFFT before is "<<magFFT[i]<<"\n";
    magFFT[i]=magFFT[i]-deltaMean;
    //cout<<"magFFT after is "<<magFFT[i]<<"\n";
    
  }
  for (int ctr=0;ctr<n;ctr++){
    if (bX[ctr]>=200 && bX[ctr]<1200){
      bY[ctr]=bY[ctr]-slope*(bX[ctr]-700.);
      //magFFT[i]=magFFT[i]+slope*(frequencyArray[i]-700.);
    }
    //  cout<<magFFT[i]<<endl;
  }

  for(int ctr1=0;ctr1<n;ctr1++){
    //  cout << "before bY, bY2 are " << bY[ctr1] << " " << bY2[ctr1] << "\n";
    // bY2[ctr1] = bY[ctr1]+2.;
    bX2.push_back(bX[ctr1]);
    bY2.push_back(bY[ctr1]+2.);
    baseX[ctr1]=bX[ctr1];
    baseY[ctr1]=bY[ctr1];
    //cout<<"baseline mean is "<<meanBaseline<<" and mean of event is "<<mean<<"\n";
    //cout<<"baseY before was "<<baseY[ctr1]<<"\n";
    baseY[ctr1]=baseY[ctr1]-meanBaseline+mean;
    //cout<<"baseY is now "<<baseY[ctr1]<<"\n";
    //cout<<" /////////////////////// \n";
	
    //cout << "after bY, bY2 are " << bY[ctr1] << " " << bY2[ctr1] << "\n";
  }
  baseN=n;
  //now make the graph
  /*
  if (drawFlag==1){
    if (pol==0){
      TGraph *grCoherentMeanAdjusted=new TGraph(newLength,frequencyArray,magFFT);
      TGraph *grbaselinehere=new TGraph(n,bX,bY);
      TGraph *grbaselinehere2=new TGraph(n,&bX2[0],&bY2[0]);
     
      
      TH2F *haxes1=new TH2F("haxes1",";Freq(MHz);dB",10,0,1400,10,-20,40);
      gStyle->SetOptStat(kFALSE);
      TCanvas *cBaseline3=new TCanvas("cBaseline3","cBaseline3",800,800);
      cBaseline3->cd(0);  
      haxes1->Draw();
      grbaselinehere->Draw("l");
      grbaselinehere2->SetLineColor(kRed);
      grbaselinehere2->Draw("l");
      grCoherentMeanAdjusted->Draw("l+");
      cBaseline3->Print("fft.eps");
      char printer[256];
  
      sprintf(printer,"FFT_%i.png",myEventNumber);
      cBaseline3->Print(printer);
      //cBaseline3->Print("cBaseline3.png");
    }
    else {
      TGraph *grCoherentMeanAdjustedH=new TGraph(newLength,frequencyArray,magFFT);
      TGraph *grbaselinehereH=new TGraph(n,bX,bY);
      TGraph *grbaselinehereH2=new TGraph(n,&bX2[0],&bY2[0]);

      TH2F *haxesH=new TH2F("haxesH",";Freq (MHz);dB",10,0,1400,10,-20,40);
      gStyle->SetOptStat(kFALSE);
      TCanvas *cBaseline3H=new TCanvas("cBaseline3H","cBaseline3H",800,800);
      cBaseline3H->cd(0);  
      haxesH->Draw();
      grbaselinehereH->Draw("l");
      grbaselinehereH2->SetLineColor(kRed);
      grbaselinehereH2->Draw("l+");
     
      grCoherentMeanAdjustedH->Draw("l+");
      cBaseline3H->Print("fft.eps");
      char printer[256];
  
      sprintf(printer,"FFTHoriz_%i.png",myEventNumber);
      cBaseline3H->Print(printer);
      // cBaseline3H->Print("cBaseline3H.png");
    }
    
  }
  */
  //find mean frequency
  mean_freq=0;
  double cum_power=0;
  double avg_power=0;
  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>200 && frequencyArray[i]<1200){
      avg_power+=pow(10,magFFT[i]/10);
    }
  }

  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>200 && frequencyArray[i]<1200){
      if (cum_power<avg_power/2.) cum_power+=pow(10,magFFT[i]/10);
      else{
	mean_freq=frequencyArray[i];
	break;
      }
    }
  }
  if (printFlag==1) cout<<"Mean frequency before filtering: "<<mean_freq<<endl;
  
  //now see if any peaks are ndB above the baseline.
  double deltaMag[newLength];

  int j;
  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>210 && frequencyArray[i]<1190){
      for (j=0;j<n;j++){
	if (bX[j]>frequencyArray[i]) break;
      }
      deltaMag[i]=magFFT[i]-bY[j];
    }
    else deltaMag[i]=-1000;
  }
  
  for (int i=0;i<nfreq;i++){
    int index;
    double maxDelta=getMaximum(newLength, deltaMag, index);
    if (maxDelta>dBCut){
      frequencies[i]=frequencyArray[index];
    }
    else frequencies[i]=-1;
    if (i==0) strongCWFlag=0;
    if (i==0 && maxDelta>5) strongCWFlag=int(maxDelta);//set a flag if CW is 5dB above baseline
    for (int i=0;i<newLength;i++){
      if (fabs(frequencyArray[i]-frequencyArray[index])<=bandWidth)
	deltaMag[i]=-1000;
    }
  }
  
}
///////////////////////////////////////
void MyCorrelator::adaptiveFilterPartialPayload(int pol, double dBCut, int nfreq, 
						double *frequencies,double *bandwidth,double *magPeak, int drawFlag, double bandWidth, double peakPhi, 
						int nantennasToUse,vector<int>& whichAntennasToUse, float &mean_freq,
						int myEventNumber, int nadirFlag, int antenna_groups)
{
  if (readBaselineFlag!=1){ 
    readBaselineFFTs();
    if (printFlag==1) cout<<"Read in Baseline For Adaptive Filter"<<endl;
  }  
  
  

  const int nPoints = 129;//124
  double magFFT[2000];
  double frequencyArray[2000];
  double deltaT, deltaF;
  int length;
  int newLength;
  double *X, *Y;
  double bY2dB[2000]={0};
  //int phi_sector =(int) round(peakPhi/7.5);
  int phi_sector = antenna_groups;
  //int whichAntennasToUse[nantennasToUse];
  int ant;
  for (int i=0;i<2000;i++){
    magFFT[i]=0;
   
    frequencyArray[i]=-1;
  }
 

  double baseline[nPoints]={0};
  //getClosestNAntennas(nantennasToUse, peakPhi, whichAntennasToUse, nadirFlag);
 
  GetBaselineperPhi(pol,baseline,nantennasToUse, whichAntennasToUse);//take n antennas and get the avg baseline for those
  if(printFlag==1) cout<<"peakPhi is "<<peakPhi<<"\n";
  for(int i0=0;i0<nantennasToUse;i0++){
    //cout<<"antennas to use are "<<whichAntennasToUse[i0]<<"\n";

  }
  // for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){  
  for (int ctr=0;ctr<nantennasToUse;ctr++){
    ant=whichAntennasToUse[ctr];
    
    if ((pol!=0 || ant!=1) && saturatedChannels[ant+pol*NUM_ANTS_WITH_NADIRS]==0){//get rid of 2V
      if (pol==0) Y = grEv[ant]->GetY();//set X and Y values
      else{ Y = grEvHoriz[ant]->GetY(); }
      if (pol==0) X = grEv[ant]->GetX();
      else{ X = grEvHoriz[ant]->GetX(); }
      
      deltaT=X[1]-X[0];
      if (pol==0) length=grEv[ant]->GetN();
      else{length=grEvHoriz[ant]->GetN();}
      //length=256;
      //cout<<"length is "<<length<<"\n";
      newLength=(length/2)+1;
      deltaF=1/(deltaT*length); //GHz
      deltaF*=1e3; //MHz
      
      // cout<<"length, deltaT, deltaF are "<<length<<" "<<deltaT<<" "<<deltaF<<" "<<newLength<<"\n";
      //cout<<"deltaF is "<<deltaF<<"\n";
     
      FFTWComplex *theFFT=FFTtools::doFFT(length,Y);
      
     
      
      for(int i=0;i<newLength;i++) {
	if (ctr==0){
	  if (i==0) frequencyArray[i]=0;
	  if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;//make a freq Array. DOES THIS CHANGE FOR EVERY ANT GRAPH??
	 
	}
	if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
	  magFFT[i]+=theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im;//adding squares of all antennas together
	 
	  //magFFT[i]+=sqrt(theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im);
	}
	else {
	  magFFT[i]=-1000;
	 
	}
	
      }   
      delete [] theFFT;
      
    }//pol!=0
  }//ctr==0
 
  for(int i=0;i<newLength;i++) {
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      if (pol==0) magFFT[i]=10*log10(sqrt(magFFT[i]/double(nantennasToUse))/10.);//correct RMS
      if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i]/double(nantennasToUse))/10.);//correct RMS
     
      //if (pol==0) magFFT[i]=10*log10(sqrt(magFFT[i])/double(NUM_ANTS_WITH_NADIRS-1)/10.);//old way
      //if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i])/double(NUM_ANTS_WITH_NADIRS)/10.);//old way
    }    
    else {
      magFFT[i]=-1000;
     
    }
   
  }

  
 
  //get baseline average so we can bump FFT around
  double meanBaseline=0;
  double mean=0;
  double firstHalfMeanBaseline=0;
  double secondHalfMeanBaseline=0;
  double firstHalfMean=0;
  double secondHalfMean=0;
  int navg=0;
  int nfirsthalfavg=0;
  int nsecondhalfavg=0;

  // double *bX, *bY;
  double *bX_set,*bY_set;
  vector<double> bX2;
  vector<double> bY2;
  int n;
  int n_set;
  double old_mean=-1000.;

  vector<int> index_skip (newLength,1);
  
  if (pol==0){//vertical
    bY_set=grVertCoherentBaseline->GetY();//should be set to partial payload
    bX_set=grVertCoherentBaseline->GetX();//should be set to partial payload
    n_set=grVertCoherentBaseline->GetN();//should be set to partial payload
  }
  if (pol==1){ //horizontal
    bY_set=grHorizCoherentBaseline->GetY();
    bX_set=grHorizCoherentBaseline->GetX();
    n_set=grHorizCoherentBaseline->GetN();
  }
  double bX[n_set];
  double bY[n_set];

  vector< double> baseline_normal(n_set,0);
   vector< double> baseline_tilt(n_set,0);
   vector <double> baseline_shift(n_set,0);
   vector <double> baseline_shift1(n_set,0);
 int n_max=10;
 
 for(int n0=0;n0<n_max;n0++){//do iterative process
  
   if( fabs (mean-old_mean) <= 0.01){
     n0 = n_max-1;
   }
   meanBaseline=0.;
   mean=0.;
   firstHalfMeanBaseline=0.;
   secondHalfMeanBaseline=0.;
   firstHalfMean=0.;
   secondHalfMean=0.;
   navg=0;
   nfirsthalfavg=0;
   nsecondhalfavg=0;
   
   for(int i0=0;i0<n_set;i0++){
     baseline_normal[i0]=0.;
     baseline_tilt[i0]=0.;
     baseline_shift[i0]=0.;
     bY[i0] = bY_set[i0];
     bX[i0] = bX_set[i0];
     n = n_set;
   }
  
   for (int i=0;i<n;i++){
     baseline_normal[i]=bY[i];
     if (bX[i]>=200 && bX[i]<1200){
       meanBaseline+=bY[i];
       navg++;
     }
     if (bX[i]>=200 && bX[i]<700){ 
       firstHalfMeanBaseline+=bY[i];
       nfirsthalfavg++;
     }
     if (bX[i]>=700 && bX[i]<1200){ 
       secondHalfMeanBaseline+=bY[i];
       nsecondhalfavg++;
     }	
   }
   meanBaseline=meanBaseline/double(navg);
  firstHalfMeanBaseline=firstHalfMeanBaseline/double(nfirsthalfavg);
  secondHalfMeanBaseline=secondHalfMeanBaseline/double(nsecondhalfavg);
  
  navg=0;
  nfirsthalfavg=0;
  nsecondhalfavg=0;
  
 

  //get average of graph in question

  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>=200 && frequencyArray[i]<1200 && index_skip[i]==1){ 
       mean+=magFFT[i];
      
      navg++;
    }
    if (frequencyArray[i]>=200 && frequencyArray[i]<700 && index_skip[i]==1){
      firstHalfMean+=magFFT[i];
     
      nfirsthalfavg++;
    }
    if (frequencyArray[i]>=700 && frequencyArray[i]<1200 && index_skip[i]==1){
      secondHalfMean+=magFFT[i];
     
      nsecondhalfavg++;
    }

  } 
  mean=mean/double(navg);
  firstHalfMean=firstHalfMean/double(nfirsthalfavg);
  secondHalfMean=secondHalfMean/double(nsecondhalfavg);
  // cout<<"mean is "<<mean<<" old mean is "<<old_mean<<"\n";
  old_mean = mean;

  //now bump the average to the baseline average and apply a tilt correction to baseline
  double deltaMean=mean-meanBaseline;
  double deltaMeanFirst=firstHalfMean-firstHalfMeanBaseline-deltaMean;
  double deltaMeanSecond=secondHalfMean-secondHalfMeanBaseline-deltaMean;
  double slope=(deltaMeanFirst-deltaMeanSecond)/500.;
 
  for(int ctr0=0;ctr0<n;ctr0++){
    baseline_shift[ctr0]=baseline_normal[ctr0]+deltaMean;
    //baseline_shift1[ctr0]=0.;
  }


  for (int ctr=0;ctr<n;ctr++){
    if (bX[ctr]>=200 && bX[ctr]<1200){
     
      bY[ctr]=bY[ctr]-slope*(bX[ctr]-700.);
      baseline_tilt[ctr]=bY[ctr];
     
      //magFFT[i]=magFFT[i]+slope*(frequencyArray[i]-700.);
    }
    else{
      baseline_tilt[ctr]=-1000;
    }
    //  cout<<magFFT[i]<<endl;
  }
  //cout<<"ant is "<<ant<<"\n";
 
  for(int i=0;i<n;i++){
    bY[i]= bY[i]+deltaMean;
    //bY[i]=bY[i]-5;
  }
 for(int ctr1=0;ctr1<n;ctr1++){
   
    bX2.push_back(bX[ctr1]);
    bY2.push_back(bY[ctr1]+2.);
   
    
    bY2dB[ctr1] = bY[ctr1]+2.;
   
  }//ctr1
 
  //now make the graph
  if(printFlag==1) cout<<"phi_sector is "<<phi_sector<<"\n";

  //find mean frequency
  mean_freq=0;
  double cum_power=0;
  double avg_power=0;
  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>200 && frequencyArray[i]<1200){
      avg_power+=pow(10,magFFT[i]/10);
    }
  }

  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>200 && frequencyArray[i]<1200){
      if (cum_power<avg_power/2.) cum_power+=pow(10,magFFT[i]/10);
      else{
      	mean_freq=frequencyArray[i];
	break;
      }
    }
  }
  if (printFlag==1) cout<<"Mean frequency before filtering: "<<mean_freq<<endl;
  
  //now see if any peaks are ndB above the baseline.
  double deltaMag[newLength];

  int j;
  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>210 && frequencyArray[i]<1190){
      for (j=0;j<n;j++){
	if (bX[j]>frequencyArray[i]) break;
      }
      deltaMag[i]=magFFT[i]-bY[j];
    }
    else deltaMag[i]=-1000;
  }

  
  int j_index;
  int k_index;
   for(int i=0;i<nfreq;i++){
    int index;
    double maxDelta=getMaximum(newLength,deltaMag,index);
    if(maxDelta>CWheight){
      CWheight=maxDelta;
    }
    if(maxDelta>dBCut){
     
      for(int j=index;j>0;j--){

	if(deltaMag[j]<dBCut-1){
	  j_index=j;
	 
	  break;
	}
      }
     
      for(int k=index;k<newLength;k++){
	
	if(deltaMag[k]<dBCut-1){
	  k_index=k;
	 
	  break;
	}
      }
      
      index = (int)ceil((k_index+j_index)/2);
      frequencies[i]=frequencyArray[index];
      bandwidth[i]=frequencyArray[k_index] - frequencyArray[j_index];
      
      magPeak[i]=maxDelta;
      if(bandwidth[i]<deltaF){
	bandwidth[i]=deltaF;
      }

      for(int reset=j_index;reset<=k_index;reset++){
	index_skip[reset]=0;
      }

      // cout<<"index is "<<index<<" freq to cut is "<<frequencyArray[index]<<" bandwidth is "<<bandwidth[i]<<"\n";
      for(int clear=j_index;clear<=k_index;clear++){
	deltaMag[clear]=-1000;
      }
    }
    else{
      magPeak[i]=-1;
      frequencies[i]=-1;
    }
   }//nfreq
  
   /*
  
   

  //static baseline method
  
    for (int i=0;i<nfreq;i++){
    int index;
    double maxDelta=getMaximum(newLength, deltaMag, index);
    cout<<"maxDelta is "<<maxDelta<<"\n";
    if (maxDelta>dBCut){
      frequencies[i]=frequencyArray[index];
      bandwidth[i]=bandWidth;
      
      for(int reset=0;reset<=newLength;reset++){
	if(frequencyArray[reset] >= (frequencyArray[index] - bandWidth) && frequencyArray[reset] <= (frequencyArray[index] + bandWidth)){
	 
	  index_skip[reset]=0;
	}
      }
      magPeak[i]=maxDelta;
      if(printFlag==0) cout<<"freq need to be cut is "<<frequencyArray[index]<<"\n";
    }
    else{
      frequencies[i]=-1;
      magPeak[i]=-1;
    }
    if (i==0) strongCWFlag=0;
    if (i==0 && maxDelta>5) strongCWFlag=int(maxDelta);//set a flag if CW is 5dB above baseline
    for (int i=0;i<newLength;i++){
      if (fabs(frequencyArray[i]-frequencyArray[index])<=bandWidth)
	deltaMag[i]=-1000;
    }
  }
   */
 
 }//n=0
  if(printFlag==1) cout<<"about to delete grVert/Horiz Coherent \n";
  delete grVertCoherentBaseline;
  delete grHorizCoherentBaseline;

  grVertCoherentBaseline=0;
  grHorizCoherentBaseline=0;

   if(printFlag==1) cout<<" deleted grVert/Horiz Coherent \n";
}
///////////////////////////////
void MyCorrelator::GetFrequenciestoCut(int antenna,vector< vector<double> > &antennaFreq,vector< vector<double> > &bandwidth, vector< vector<double> > &PeakMag, vector<double> &uniquefreqs,vector<double> &uniquebandwidth, int nfreq, vector<double> &uniquePhase, vector<double> &uniquePhase_bandwidth){
  double max=0.;
  int maxbin=0;
  double bandWidth=0.;
  vector< vector<double> > freq_mag(50, vector<double>(3));
  int new_freq_flag=0;
  int non_zero_flag=0;
  uniquefreqs.clear();
  uniquebandwidth.clear();
  //get freqs to cut

  for(int n=0;n<50;n++){
    max=0.;
    maxbin=0;
    for(int j=0;j<50;j++){
      if(PeakMag[antenna][j]>max){
	max=PeakMag[antenna][j];
	maxbin =j;
      }
    }
    if(max>0){
     
      freq_mag[n][0]=antennaFreq[antenna][maxbin];
      freq_mag[n][1]=max;
      freq_mag[n][2]=bandwidth[antenna][maxbin];
    }
    PeakMag[antenna][maxbin]=-1000;
  }//n
 
  for(int k0=0;k0<50;k0++){
    new_freq_flag=0;
   
    if(non_zero_flag==0){
      if(freq_mag[k0][0]> 1 && freq_mag[k0][0]<1210){
	uniquefreqs.push_back(freq_mag[k0][0]); 
	non_zero_flag=1;
      }
    }
    for(int k1=0;k1<(int)uniquefreqs.size();k1++){
      if((int)freq_mag[k0][0]==(int)uniquefreqs[k1]){
	new_freq_flag=0;
	break;
      }
      else{
	new_freq_flag=1;
      }
    }//k1
    
    if(new_freq_flag==1 && (int)uniquefreqs.size()<nfreq  && freq_mag[k0][0]> 190 && freq_mag[k0][0]<1210){
      uniquefreqs.push_back(freq_mag[k0][0]); 
     
    }
    
  }//k0
  vector<double> uniqueMag;
  double mag;
  //vector<double> uniquePhase;
  //vector<double> uniquePhase_bandwidth;
  
  if((int)uniquefreqs.size() >0){
    sort(uniquefreqs.begin(),uniquefreqs.end());//sort all the freqs for this antenna

    for(int k0=0;k0<(int)uniquefreqs.size();k0++){
      bandWidth=0.;
      for(int k1=0;k1<50;k1++){
	if(uniquefreqs[k0]==freq_mag[k1][0]){
	  if(freq_mag[k1][2]>bandWidth){
	    bandWidth=freq_mag[k1][2];
	    mag = freq_mag[k1][1];
	  }
	}
      }
      if(bandWidth >1) {
	
	uniquebandwidth.push_back(bandWidth);
	uniqueMag.push_back(mag);
      }
    }
    
    uniquePhase.push_back(uniquefreqs[0]);
    uniquePhase_bandwidth.push_back(uniquebandwidth[0]);
    int kappa=0;
    for(int trial0=0;trial0<(int)uniquefreqs.size();trial0++){
      
      if(uniquefreqs[trial0] < uniquePhase[kappa]-15 || uniquefreqs[trial0]>uniquePhase[kappa]+15){
	if(uniquefreqs[trial0] >200 && uniquefreqs[trial0]<1200){
	  uniquePhase.push_back(uniquefreqs[trial0]);
	  uniquePhase_bandwidth.push_back(uniquebandwidth[trial0]);
	  kappa++;
	}
      }
    }

    
    double middle_freq_value;
    
    double new_bandwidth;
    
    double lower1,lower2;
    double higher1,higher2;
  
    for(int tester=0;tester<(int)uniquefreqs.size()-1;tester++){
      lower1 = uniquefreqs[tester]-uniquebandwidth[tester];
      higher1 = uniquefreqs[tester]+uniquebandwidth[tester];
      
      lower2 = uniquefreqs[tester+1]-uniquebandwidth[tester+1];
      higher2 = uniquefreqs[tester+1]+uniquebandwidth[tester+1];
     
      if(lower1<=lower2-1 && higher1>=higher2+1){//2nd freq to cut lies inside first's bandwidth
	//	cout<<"first scenario! \n\n";
	uniquefreqs.erase (uniquefreqs.begin()+tester+1);	    
	uniquebandwidth.erase (uniquebandwidth.begin()+tester+1);	    
	tester=tester-1;
      }
      else if(lower2<=lower1-1 && higher2>=higher1+1){//1st freq to cut lies inside 2nds bandwidth
	//	cout<<"second scenario! \n\n";
	uniquefreqs.erase (uniquefreqs.begin()+tester);	    
	uniquebandwidth.erase (uniquebandwidth.begin()+tester);	    
	tester=tester-1;
      }
      
      
      else if(higher1 >= lower2-1){//overlap. combine.
	//cout<<"third scenario! \n\n";
	middle_freq_value = (lower1 + higher2)/2;
	new_bandwidth = middle_freq_value - lower1;
	//cout<<"middle_freq_value is "<<middle_freq_value<<" new_bandwidth is "<<new_bandwidth<<"\n";


	uniquefreqs[tester]=middle_freq_value;
	uniquebandwidth[tester]=new_bandwidth;

	uniquefreqs.erase (uniquefreqs.begin()+tester+1);
	
	uniquebandwidth.erase (uniquebandwidth.begin()+tester+1);

	tester=tester-1;

      }
    }
  }//uniquefreqs.size()>0




 
}



///////////////////////////////
void MyCorrelator::GetSatelliteFrequenciestoCut(int antenna,vector<double> &uniquefreqs,vector<double> &uniquebandwidth, double &satellite_freq,double &satellite_bandwidth,
						double &satellite_freq2, double &satellite_bandwidth2, int &satellite_flag){
  double lowest_freq;
  double highest_freq;
  double lowest_satellite;
  double highest_satellite;

  if(uniquefreqs.size() >0){
    lowest_freq=uniquefreqs[0]-uniquebandwidth[0];
    highest_freq=uniquefreqs[0]+uniquebandwidth[0];
    lowest_satellite =  satellite_freq-satellite_bandwidth;
    highest_satellite = satellite_freq + satellite_bandwidth;
    if(lowest_freq >= lowest_satellite && lowest_freq <=highest_satellite && highest_freq >= highest_satellite){
      
      satellite_freq = (lowest_freq+lowest_satellite)/2;
      satellite_bandwidth = satellite_freq -lowest_satellite;
    }
    else if(lowest_freq <= lowest_satellite && highest_freq >= lowest_satellite && highest_freq <=highest_satellite){
      
      satellite_freq = (highest_freq+highest_satellite)/2;
      satellite_bandwidth  =highest_satellite - satellite_freq;
    }
    else if(lowest_freq > lowest_satellite && highest_freq < highest_satellite){
      
      satellite_flag=2;
      satellite_freq = (lowest_freq +lowest_satellite)/2;
      satellite_bandwidth = (satellite_freq - lowest_satellite);
      
      satellite_freq2 = (highest_satellite+highest_freq)/2;
      satellite_bandwidth2 = highest_satellite - satellite_freq2;
    }
    else if(lowest_freq < lowest_satellite && highest_freq > highest_satellite){
      
      satellite_flag=1;
    }
  } 

}
///////////////////////////////
void MyCorrelator::GetFFTandBaseline(int nantennasToUse,vector<int>& whichAntennasToUse, int pol, double *frequencyArray,
				     double *baseX, double *baseY2, double *magFFT2, int nfreq, int dBCut){

 
  double deltaT, deltaF;
  int length;
  int newLength;
  double *X, *Y;
  int ant;
  double magFFT[2000]={0.};
  double phase[2000]={0.};
  for (int ctr=0;ctr<nantennasToUse;ctr++){
    ant = whichAntennasToUse[ctr];
    
   if ((pol!=0 || ant!=1) && saturatedChannels[ant+pol*NUM_ANTS_WITH_NADIRS]==0){//get rid of 2V
    
     if (pol==0) Y = grEv[ant]->GetY();//set X and Y values
     else{ Y = grEvHoriz[ant]->GetY(); }
     if (pol==0) X = grEv[ant]->GetX();
     else{ X = grEvHoriz[ant]->GetX(); }
     
     deltaT=X[1]-X[0];
     if (pol==0) length=grEv[ant]->GetN();
     else{length=grEvHoriz[ant]->GetN();}
    
     //length=256;
     // cout<<"length is "<<length<<"\n";
     newLength=(length/2)+1;
     deltaF=1/(deltaT*length); //GHz
     deltaF*=1e3; //MHz
     
     //cout<<"deltaF is "<<deltaF<<"\n";
     
     FFTWComplex *theFFT=FFTtools::doFFT(length,Y);
     
    
     for(int i=0;i<newLength;i++) {
       if (ctr==0){
	 if (i==0) frequencyArray[i]=0.;
	 if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;//make a freq Array. DOES THIS CHANGE FOR EVERY ANT GRAPH??
	 
       }
       if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
	
	 magFFT[i]+=(theFFT[i].re*theFFT[i].re) + (theFFT[i].im*theFFT[i].im);//adding squares of all antennas together
	 magFFT2[i]+=(theFFT[i].re*theFFT[i].re) + (theFFT[i].im*theFFT[i].im);//adding squares of all antennas together
	 //magFFT[i]+=sqrt(theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im);
	 phase[i]= atan2(theFFT[i].im,theFFT[i].re);
       }
       else{
	 magFFT[i]=-40;
	 magFFT2[i]=0;
       }    
       
     }   
     delete [] theFFT;
     
   }//pol!=0
 }//ctr
 
  for(int i=0;i<newLength;i++) {
   
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      if (pol==0) magFFT[i]=10*log10(sqrt(magFFT[i]/double(nantennasToUse))/10.);//correct RMS
      if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i]/double(nantennasToUse))/10.);//correct RMS
      
      magFFT2[i]=sqrt(magFFT2[i]/double(nantennasToUse));
      // if (pol==0) magFFT[i]=10*log10(sqrt(magFFT[i])/double(NUM_ANTS_WITH_NADIRS-1)/10.);//old way
      //if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i])/double(NUM_ANTS_WITH_NADIRS)/10.);//old way
    }    
    else {
      magFFT[i]=-20;
      magFFT2[i]=0;
    }
    // cout<<"aaaa  "<<magFFT[i]<<endl;
  }
 
  //get baseline average so we can bump FFT around
 double meanBaseline=0;
  double mean=0.;
  double firstHalfMeanBaseline=0;
  double secondHalfMeanBaseline=0;
  double firstHalfMean=0;
  double secondHalfMean=0;
  int navg=0;
  int nfirsthalfavg=0;
  int nsecondhalfavg=0;

  double *bX_set,*bY_set;
  double bY2dB[2000]={0};
  int n;
  int n_set;
 
  vector<int> index_skip (newLength,1);
  int n_max=10;
  
  double old_mean=-1000;
  if (pol==0){//vertical
    bY_set=grVertCoherentBaseline->GetY();//should be set to partial payload
    bX_set=grVertCoherentBaseline->GetX();//should be set to partial payload
    n_set=grVertCoherentBaseline->GetN();//should be set to partial payload
  }
  if (pol==1){ //horizontal
    bY_set=grHorizCoherentBaseline->GetY();
    bX_set=grHorizCoherentBaseline->GetX();
    n_set=grHorizCoherentBaseline->GetN();
  }
  double bX[n_set];
  double bY[n_set];
  n = n_set;
  for(int n0=0;n0<n_max;n0++){
    if(fabs(mean-old_mean)<0.01){
      n0=n_max-1;
    }
    meanBaseline=0.;
    mean=0.;
    firstHalfMeanBaseline=0.;
    secondHalfMeanBaseline=0.;
    firstHalfMean=0.;
    secondHalfMean=0.;
    navg=0;
    nfirsthalfavg=0;
    nsecondhalfavg=0;
  for(int i0=0;i0<n_set;i0++){
    bY[i0] = bY_set[i0];
    bX[i0] = bX_set[i0];
    n = n_set;
  }

  for (int i=0;i<n;i++){
   
    if (bX[i]>=200 && bX[i]<1200){
      meanBaseline+=bY[i];
      navg++;
    }
    if (bX[i]>=200 && bX[i]<700){ 
      firstHalfMeanBaseline+=bY[i];
      nfirsthalfavg++;
    }
    if (bX[i]>=700 && bX[i]<1200){ 
      secondHalfMeanBaseline+=bY[i];
      nsecondhalfavg++;
    }	
  }
  meanBaseline=meanBaseline/double(navg);
  firstHalfMeanBaseline=firstHalfMeanBaseline/double(nfirsthalfavg);
  secondHalfMeanBaseline=secondHalfMeanBaseline/double(nsecondhalfavg);
  
  navg=0;
  nfirsthalfavg=0;
  nsecondhalfavg=0;
  
 

  //get average of graph in question

  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>=200 && frequencyArray[i]<1200 && index_skip[i]==1){ 
       mean+=magFFT[i];
       
      navg++;
    }
    if (frequencyArray[i]>=200 && frequencyArray[i]<700 && index_skip[i]==1){
      firstHalfMean+=magFFT[i];
      
      nfirsthalfavg++;
    }
    if (frequencyArray[i]>=700 && frequencyArray[i]<1200 && index_skip[i]==1){
      secondHalfMean+=magFFT[i];
     
      nsecondhalfavg++;
    }

  }
  if(navg <1) navg=1;
  if(nfirsthalfavg <1) nfirsthalfavg=1;
  if(nsecondhalfavg <1) nsecondhalfavg=1;

  mean=mean/double(navg);
  firstHalfMean=firstHalfMean/double(nfirsthalfavg);
  secondHalfMean=secondHalfMean/double(nsecondhalfavg);
 
  old_mean = mean;

  //now bump the average to the baseline average and apply a tilt correction to baseline
  double deltaMean=mean-meanBaseline;
  double deltaMeanFirst=firstHalfMean-firstHalfMeanBaseline-deltaMean;
  double deltaMeanSecond=secondHalfMean-secondHalfMeanBaseline-deltaMean;
  double slope=(deltaMeanFirst-deltaMeanSecond)/500.;

  for (int ctr=0;ctr<n;ctr++){
    if (bX[ctr]>=200 && bX[ctr]<1200){
     
      bY[ctr]=bY[ctr]-slope*(bX[ctr]-700.);
     
    } 
    else bY[ctr]=-15;
  }
 
  for(int i=0;i<n;i++){
    bY[i]= bY[i]+deltaMean;
   
  }

 for(int i=0;i<n;i++){
   if(bX[i]>=200 && bX[i]<=1200){
     baseY2[i]=10.*pow(10.,bY[i]/10.);//used for sigma
   }
   else
     baseY2[i]=0;
 }

 for(int ctr1=0;ctr1<n;ctr1++){
   baseX[ctr1]=bX[ctr1]; 
   bY2dB[ctr1] = bY[ctr1]+2.;
   
  }//ctr1
 
  //now see if any peaks are ndB above the baseline.
  double deltaMag[newLength];

  int j;
  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>210 && frequencyArray[i]<1190){
      for (j=0;j<n;j++){
	if (bX[j]>frequencyArray[i]) break;
      }
      deltaMag[i]=magFFT[i]-bY[j];
    }
    else deltaMag[i]=-1000;
  }

  
  int j_index;
  int k_index;
  int nfreq=5;
   for(int i=0;i<nfreq;i++){
    int index;
    double maxDelta=getMaximum(newLength,deltaMag,index);
    if(maxDelta>dBCut){
     
      for(int j=index;j>0;j--){
	
	if(deltaMag[j]<1){
	  j_index=j-1;
	  break;
	}
      }
     
      for(int k=index;k<newLength;k++){
	if(deltaMag[k]<1){
	  k_index=k+1;
	  break;
	}
      }

      for(int reset=j_index;reset<=k_index;reset++){
	index_skip[reset]=0;
      }

      for(int clear=j_index;clear<=k_index;clear++){
	deltaMag[clear]=-1000;
      }
    }
   
   }//nfreq
  
  }//n=0
 

}
///////////////////////////////
void MyCorrelator::GetFFTandBaseline_old(int nantennasToUse,vector<int>& whichAntennasToUse, int pol, double *frequencyArray,
				     double *baseX, double *baseY2, double *magFFT2){

 
  double deltaT, deltaF;
  int length;
  int newLength;
  double *X, *Y;
  int ant;
  double magFFT[2000];
  for (int ctr=0;ctr<nantennasToUse;ctr++){
    ant = whichAntennasToUse[ctr];
    
   if ((pol!=0 || ant!=1) && saturatedChannels[ant+pol*NUM_ANTS_WITH_NADIRS]==0){//get rid of 2V
    
     if (pol==0) Y = grEv[ant]->GetY();//set X and Y values
     else{ Y = grEvHoriz[ant]->GetY(); }
     if (pol==0) X = grEv[ant]->GetX();
     else{ X = grEvHoriz[ant]->GetX(); }
     
     deltaT=X[1]-X[0];
     if (pol==0) length=grEv[ant]->GetN();
     else{length=grEvHoriz[ant]->GetN();}
    
     //length=256;
     //cout<<"length is "<<length<<"\n";
     newLength=(length/2)+1;
     deltaF=1/(deltaT*length); //GHz
     deltaF*=1e3; //MHz
     
     //cout<<"deltaF is "<<deltaF<<"\n";
     
     FFTWComplex *theFFT=FFTtools::doFFT(length,Y);
     
    
     for(int i=0;i<newLength;i++) {
       if (ctr==0){
	 if (i==0) frequencyArray[i]=0.;
	 if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;//make a freq Array. DOES THIS CHANGE FOR EVERY ANT GRAPH??
	 
       }
       if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
	 magFFT[i]+=theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im;//adding squares of all antennas together
	 magFFT2[i]+=theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im;//adding squares of all antennas together
	 
	 //magFFT[i]+=sqrt(theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im);
       }
       else{
	 magFFT[i]=-1000;
	 magFFT2[i]=0;
       }    
       
     }   
     delete [] theFFT;
     
   }//pol!=0
 }//ctr
 
  for(int i=0;i<newLength;i++) {
   
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      if (pol==0) magFFT[i]=10*log10(sqrt(magFFT[i]/double(nantennasToUse))/10.);//correct RMS
      if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i]/double(nantennasToUse))/10.);//correct RMS
     
      magFFT2[i]=sqrt(magFFT2[i]/double(nantennasToUse));
      // if (pol==0) magFFT[i]=10*log10(sqrt(magFFT[i])/double(NUM_ANTS_WITH_NADIRS-1)/10.);//old way
      //if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i])/double(NUM_ANTS_WITH_NADIRS)/10.);//old way
    }    
    else {
      magFFT[i]=-1000;
      magFFT2[i]=0;
    }
    // cout<<"aaaa  "<<magFFT[i]<<endl;
  }
 
  //get baseline average so we can bump FFT around
  double meanBaseline=0;
  double mean=0;
  double firstHalfMeanBaseline=0;
  double secondHalfMeanBaseline=0;
  double firstHalfMean=0;
  double secondHalfMean=0;
  int navg=0;
  int nfirsthalfavg=0;
  int nsecondhalfavg=0;

  double *bX, *bY;
 
  double bY2dB[2000]={0};
  int n;
 
 
  if (pol==0){//vertical
    bY=grVertCoherentBaseline->GetY();//should be set to partial payload
    bX=grVertCoherentBaseline->GetX();//should be set to partial payload
    n=grVertCoherentBaseline->GetN();//should be set to partial payload
  }
  if (pol==1){ //horizontal
    bY=grHorizCoherentBaseline->GetY();
    bX=grHorizCoherentBaseline->GetX();
    n=grHorizCoherentBaseline->GetN();
  }
 
  double baseline_normal[n];
  for (int i=0;i<n;i++){
    if (bX[i]>=200 && bX[i]<1200){ 
      meanBaseline+=bY[i];
      //cout<<"meanBaseline is "<<meanBaseline<<"\n";
      navg++;
    }
    if (bX[i]>=200 && bX[i]<700){ 
      firstHalfMeanBaseline+=bY[i];
      nfirsthalfavg++;
    }
    if (bX[i]>=700 && bX[i]<1200){ 
      secondHalfMeanBaseline+=bY[i];
      nsecondhalfavg++;
    }	
  }
  meanBaseline=meanBaseline/double(navg);
  firstHalfMeanBaseline=firstHalfMeanBaseline/double(nfirsthalfavg);
  secondHalfMeanBaseline=secondHalfMeanBaseline/double(nsecondhalfavg);
  
  navg=0;
  nfirsthalfavg=0;
  nsecondhalfavg=0;
 
  //get average of graph in question
  for (int i=0;i<newLength;i++){
    if (frequencyArray[i]>=200 && frequencyArray[i]<1200){ 
      mean+=magFFT[i];
      navg++;
    }
    if (frequencyArray[i]>=200 && frequencyArray[i]<700){
      firstHalfMean+=magFFT[i];
      nfirsthalfavg++;
    }
    if (frequencyArray[i]>=700 && frequencyArray[i]<1200){
      secondHalfMean+=magFFT[i];
      nsecondhalfavg++;
    }

  } 
  mean=mean/double(navg);
  firstHalfMean=firstHalfMean/double(nfirsthalfavg);
  secondHalfMean=secondHalfMean/double(nsecondhalfavg);
 
  //now bump the average to the baseline average and apply a tilt correction to baseline
  double deltaMean=mean-meanBaseline;
  double deltaMeanFirst=firstHalfMean-firstHalfMeanBaseline-deltaMean;
  double deltaMeanSecond=secondHalfMean-secondHalfMeanBaseline-deltaMean;
  double slope=(deltaMeanFirst-deltaMeanSecond)/500.;
  
  /* for (int i=0;i<newLength;i++){
    magFFT[i]=magFFT[i]-deltaMean;
    }*/
  //cout<<" ** ant is "<<ant<<" deltaMean is "<<deltaMean<<"\n";
  
  for (int ctr=0;ctr<n;ctr++){
    if (bX[ctr]>=200 && bX[ctr]<1200){
      bY[ctr]=bY[ctr]-slope*(bX[ctr]-700.);
      //magFFT[i]=magFFT[i]+slope*(frequencyArray[i]-700.);
    }
    //  cout<<magFFT[i]<<endl;
  }
  
  for(int i=0;i<n;i++){
    bY[i]= bY[i]+deltaMean;
  }
  

  for(int i=0;i<n;i++){
    baseline_normal[i]=bY[i];
    if(bX[i]>=200 && bX[i]<=1200){
      baseY2[i]=10.*pow(10.,bY[i]/10.);//used for sigma
    }
    else
      baseY2[i]=0;
  }

   


  for(int ctr1=0;ctr1<n;ctr1++){
    //  cout << "before bY, bY2 are " << bY[ctr1] << " " << bY2[ctr1] << "\n";
    // bY2[ctr1] = bY[ctr1]+2.;
   
    baseX[ctr1]=bX[ctr1];
    
    bY2dB[ctr1]=bY[ctr1]+2;
   
    //baseY[ctr1]=baseY[ctr1]-meanBaseline+mean;//puts baseY mean = mean of FFT
    
  }
 
  navg=0;
   for (int i=0;i<n;i++){
    if (bX[i]>=200 && bX[i]<1200){ 
      meanBaseline+=bY[i];
      navg++;
    }
    	
  }
  meanBaseline=meanBaseline/double(navg);
 
   //cout<<"abouttoexit \n";

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MyCorrelator::applyAdaptiveFilter(double centerFrequency, double bandWidth, int polFlag)
{
  int nantennas=NUM_ANTS_WITH_NADIRS;
  double baseY[100];
  if (polFlag==0){
    if (centerFrequency!=-1){
      for (int ant=0;ant<nantennas;ant++){
	TGraph *gr1 = new TGraph(grEv[ant]->GetN(),grEv[ant]->GetX(),grEv[ant]->GetY());
	delete grEv[ant];
	//TGraph *grNotch=simpleNotchFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth);
	//cout<<"center freq is "<<centerFrequency<<" and bandwidth is "<<bandWidth<<"\n";
	TGraph *grNotch=complicatedNotchFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth,ant,polFlag,baseY);
	grEv[ant]=new TGraph(grNotch->GetN(),grNotch->GetX(),grNotch->GetY());
	delete gr1;
	delete grNotch;
      }//ant
    }//centerfreq
  }//polFlag
  else{
      if (centerFrequency!=-1){
      for (int ant=0;ant<nantennas;ant++){
	TGraph *gr1Horiz=new TGraph(grEvHoriz[ant]->GetN(),grEvHoriz[ant]->GetX(),grEvHoriz[ant]->GetY());
	delete grEvHoriz[ant];
	//TGraph *grNotchHoriz=simpleNotchFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth);
	TGraph *grNotchHoriz=complicatedNotchFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth,ant,polFlag,baseY);
	grEvHoriz[ant]=new TGraph(grNotchHoriz->GetN(),grNotchHoriz->GetX(),grNotchHoriz->GetY());
	delete gr1Horiz;
	delete grNotchHoriz;
    
      }//ant 
    
      }//centerfreq
  }//else
 
}//main

///////////////////////////////
void MyCorrelator::applyAdaptiveFilter_singleAnt(double centerFrequency, double bandWidth, int polFlag,int ant,double *baseY)
{
  
  if (polFlag==0){
    if (centerFrequency!=-1){
      TGraph *gr1 = new TGraph(grEv[ant]->GetN(),grEv[ant]->GetX(),grEv[ant]->GetY());
      delete grEv[ant];
      //TGraph *grNotch=simpleNotchFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth);
      //cout<<"center freq is "<<centerFrequency<<" and bandwidth is "<<bandWidth<<"\n";
      if(notchFilterFlag==0){
	cout<<"using noFill Notch Filter! \n";
	TGraph *grNotch = nofillNotchFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth);
	grEv[ant]=new TGraph(grNotch->GetN(),grNotch->GetX(),grNotch->GetY());

	delete gr1;
	delete grNotch;
      }
      if(notchFilterFlag==1){
	cout<<"using rayleigh notch filter \n";
	TGraph *grNotch=complicatedNotchFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth,ant, polFlag,baseY);	
	grEv[ant]=new TGraph(grNotch->GetN(),grNotch->GetX(),grNotch->GetY());
	delete gr1;
	delete grNotch;
      }
       if(notchFilterFlag==2){
	cout<<"using wiener notch filter \n";
	TGraph *grNotch=wienerFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth,ant, polFlag,baseY);	
	grEv[ant]=new TGraph(grNotch->GetN(),grNotch->GetX(),grNotch->GetY());
	delete gr1;
	delete grNotch;
      }
      if(notchFilterFlag==3){
	//cout<<"using interpolated notch filter \n";
	TGraph *grNotch=interpolatedFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth);
	grEv[ant]=new TGraph(grNotch->GetN(),grNotch->GetX(),grNotch->GetY());
	delete gr1;
	delete grNotch;
      }
     
    }//centerfreq
  }//polFlag
  else{
    if (centerFrequency!=-1){
      TGraph *gr1Horiz=new TGraph(grEvHoriz[ant]->GetN(),grEvHoriz[ant]->GetX(),grEvHoriz[ant]->GetY());
      delete grEvHoriz[ant];
      //TGraph *grNotchHoriz=simpleNotchFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth);
      if(notchFilterFlag==0){
	TGraph *grNotchHoriz = nofillNotchFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth);
	 grEvHoriz[ant]=new TGraph(grNotchHoriz->GetN(),grNotchHoriz->GetX(),grNotchHoriz->GetY());
	 delete gr1Horiz;
	 delete grNotchHoriz;
      }
      if(notchFilterFlag==1){
	TGraph *grNotchHoriz=complicatedNotchFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth,ant, polFlag,baseY);
	grEvHoriz[ant]=new TGraph(grNotchHoriz->GetN(),grNotchHoriz->GetX(),grNotchHoriz->GetY());
	delete gr1Horiz;
	delete grNotchHoriz;
      }
      if(notchFilterFlag==3){
	TGraph *grNotchHoriz=interpolatedFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth);
	grEvHoriz[ant]=new TGraph(grNotchHoriz->GetN(),grNotchHoriz->GetX(),grNotchHoriz->GetY());
	delete gr1Horiz;
	delete grNotchHoriz;
      }
      
      }//centerfreq
  }//else
  
}//main
////////////////////////////////////////////

void MyCorrelator::applySatelliteFilter(double centerFrequency, double bandWidth,int ant,int pol,double *baseY)
{

  double headingAnt;
  
  if (centerFrequency!=-1){
    // for (int ant=0;ant<nantennas;ant++){
      headingAnt=fAdu5APatPtr->heading-phiAnt[ant]; 
      if (headingAnt>=360) headingAnt-=360;
      if (headingAnt<0) headingAnt+=360;
      if (headingAnt<=120 || headingAnt>=240){//northward facing antennas (plus 30 degrees on either side, so filtering 180+30+30 degrees
	//cout<<"ant being filtered by satellite is "<<ant<<"\n";
	
	if(pol==0){
	  TGraph *gr1 = new TGraph(grEv[ant]->GetN(),grEv[ant]->GetX(),grEv[ant]->GetY());
	  delete grEv[ant];
	 
	  if(notchFilterFlag==0){
	    TGraph *grNotch = nofillNotchFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth);
	    grEv[ant]=new TGraph(grNotch->GetN(),grNotch->GetX(),grNotch->GetY());
	    delete gr1;
	    delete grNotch;
	  }
	  if(notchFilterFlag==1){
	    TGraph *grNotch=complicatedNotchFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth,ant,0,baseY);	   
	    grEv[ant]=new TGraph(grNotch->GetN(),grNotch->GetX(),grNotch->GetY());	   
	    delete gr1;
	    delete grNotch;	   
	  }
	  
	  if(notchFilterFlag==3){
	   
	    TGraph *grNotch=interpolatedFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth);
	    grEv[ant]=new TGraph(grNotch->GetN(),grNotch->GetX(),grNotch->GetY());
	    delete gr1;
	    delete grNotch;
	   
	  }
	 
	}//pol==0
	if(pol==1){
	  TGraph *gr1Horiz=new TGraph(grEvHoriz[ant]->GetN(),grEvHoriz[ant]->GetX(),grEvHoriz[ant]->GetY());
	
	  delete grEvHoriz[ant];
	  if(notchFilterFlag==0){	   
	    TGraph *grNotchHoriz = nofillNotchFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth);
	    grEvHoriz[ant]=new TGraph(grNotchHoriz->GetN(),grNotchHoriz->GetX(),grNotchHoriz->GetY());
	    delete gr1Horiz;
	    delete grNotchHoriz;
	    
	    
	  }
	  if(notchFilterFlag==1){
	    TGraph *grNotchHoriz=complicatedNotchFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth,ant,1,baseY);
	    grEvHoriz[ant]=new TGraph(grNotchHoriz->GetN(),grNotchHoriz->GetX(),grNotchHoriz->GetY());
	    delete gr1Horiz;
	    delete grNotchHoriz;
	    
	  }
	  if(notchFilterFlag==3){
	    TGraph *grNotchHoriz=interpolatedFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth);
	    grEvHoriz[ant]=new TGraph(grNotchHoriz->GetN(),grNotchHoriz->GetX(),grNotchHoriz->GetY());
	    delete gr1Horiz;
	    delete grNotchHoriz;
	  }
	 
	}//pol==1
      }
      // }  
  }
  
}
////////////////////////////////
inline Double_t MyCorrelator::getDeltaTExpected(Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave)
{
  if (thetaWave<0) thetaWave+=360;
  thetaArrayIndex=int(thetaWave*10);  

  phi1ArrayIndex=int((phiWave-phiAnt[ant1])*10);
  if (phi1ArrayIndex<0) phi1ArrayIndex+=3600;
  if (phi1ArrayIndex>=3600) phi1ArrayIndex-=3600;
  phi2ArrayIndex=int((phiWave-phiAnt[ant2])*10);
  if (phi2ArrayIndex<0) phi2ArrayIndex+=3600;
  if (phi2ArrayIndex>=3600) phi2ArrayIndex-=3600;
  
  Double_t part1=zAnt[ant1]*tanArray[thetaArrayIndex] - rAnt[ant1] * cosArray[phi1ArrayIndex];
  Double_t part2=zAnt[ant2]*tanArray[thetaArrayIndex] - rAnt[ant2] * cosArray[phi2ArrayIndex];
  
  double geomDelay=1e9*((cosArray[thetaArrayIndex] * (part1 - part2))/C_LIGHT);    //returns time in ns
  
  //add in off-axis antenna delay here
  if (groupDelayFlag==1){
    Double_t phi1Diff=fUPGeomTool->getPhiDiff(phiWave*deg2rad,phiAnt[ant1]*deg2rad);
    Double_t delay1=getGroupDelay(phi1Diff,thetaWave*deg2rad);
    Double_t phi2Diff=fUPGeomTool->getPhiDiff(phiWave*deg2rad,phiAnt[ant2]*deg2rad);
    Double_t delay2=getGroupDelay(phi2Diff,thetaWave*deg2rad);
    geomDelay+=(delay1-delay2);
  }
  return geomDelay;
}
////////////////////////////////////////
inline Double_t MyCorrelator::getGroupDelay(Double_t phiToAntBoresight, Double_t thetaWave)//in radians, with positive being down
{
  Double_t thetaDeg=thetaWave*TMath::RadToDeg()-10;
  Double_t phiDeg=phiToAntBoresight*TMath::RadToDeg();
  Double_t totalAngleDeg=thetaDeg*thetaDeg+phiDeg*phiDeg;
  if (totalAngleDeg>2500) totalAngleDeg=2500;
  //Double_t delayTime=(totalAngleDeg*totalAngleDeg*totalAngleDeg*totalAngleDeg)*1.303e-8;
  //delayTime+=(totalAngleDeg*totalAngleDeg*totalAngleDeg)*6.544e-8;
  //delayTime-=(totalAngleDeg*totalAngleDeg)*4.770e-6;
  //delayTime+=totalAngleDeg*8.097e-4;

  Double_t delayTime=(totalAngleDeg*totalAngleDeg)*1.45676e-8;
  delayTime-=(totalAngleDeg)*5.01452e-6;

  return delayTime;
}
/////////////////////////////////
TGraph *MyCorrelator::makeCoherentlySummedDeconvolvedWaveform(int myEventNumber, double peakTheta, 
							      double peakPhi, int nadirFlag, int polFlag,
							      int nantennas, int drawFlag)
{
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  TGraph *grEvCoherent[nantennas];
  int upSampleFactor=1;
  Double_t deltaTInt=1./(2.6*upSampleFactor);
  
  vector<int> whichAntennasToUse (nantennas,0);
  getClosestNAntennas(nantennas, peakPhi, whichAntennasToUse,nadirFlag);
  
  
  
  for (int i=0;i<nantennas;i++){
    int ant=whichAntennasToUse[i];
    // grEvCoherent[i]=FFTtools::getInterpolatedGraph(grEv[ant], deltaTInt);
    //cout<<"ant: "<<ant<<endl;
    if (polFlag==0) grEvCoherent[i]=deconvolveWaveformUsingStephens(grEvUnfiltered[ant],0,0);//grEvUnfiltered
    else grEvCoherent[i]=deconvolveWaveformUsingStephens(grEvHorizUnfiltered[ant],1,0);
  }
  
  int npoints=260*upSampleFactor*3;
  double voltsCoherent[npoints];
  double timesCoherent[npoints];
  for (int j=0;j<npoints;j++){
    timesCoherent[j]=0;
    voltsCoherent[j]=0;
  }
  
  int nsToLeaveOut=20;
  int factorToLeaveOut=int(nsToLeaveOut*2.6);

  for (int i=0;i<nantennas;i++){
    int ant=whichAntennasToUse[i];
    //cout<<"antenna "<<whichAntennasToUse[i]+1<<" used for coherently summed waveform"<<endl;
    double timeExpected=getDeltaTExpected(ant, whichAntennasToUse[0], 
							  peakPhi, peakTheta);
    //cout<<"time expected between this and 0: "<<timeExpected<<endl;
    //get the value of graph at this delay.
    double tVal1, tVal2, val1, val2, voltVal;			
    double tVal0,val0, tValend, valend;
    grEvCoherent[i]->GetPoint(0,tVal0,val0);
    grEvCoherent[i]->GetPoint(grEvCoherent[i]->GetN()-1,tValend,valend);
    double lengthOfTrace=tValend-tVal0;
    for (int j=factorToLeaveOut*upSampleFactor;j<npoints-factorToLeaveOut*upSampleFactor;j++){
      int bin1=int((timeExpected+j*deltaTInt-tVal0)/lengthOfTrace*grEvCoherent[i]->GetN());
      int bin2=bin1+1;
      //if (j==factorToLeaveOut*upSampleFactor) cout<<"bin1: "<<bin1<<", bin2: "<<bin2<<endl;
      if (bin1>=0 && bin2<=grEvCoherent[i]->GetN()){
	grEvCoherent[i]->GetPoint(bin1,tVal1,val1);
	grEvCoherent[i]->GetPoint(bin2,tVal2,val2);
	voltVal=(val2-val1)*(timeExpected+j*deltaTInt-tVal1)/(tVal2-tVal1)+val1;//interpolate 2 bins 
	voltsCoherent[j]+=voltVal;
      }
    }
  }
  for (int i=0;i<nantennas;i++){
    delete grEvCoherent[i];
  }
  
  for (int j=factorToLeaveOut*upSampleFactor;j<npoints-factorToLeaveOut*upSampleFactor;j++){
    timesCoherent[j]=deltaTInt*j;
    voltsCoherent[j]=voltsCoherent[j]/nantennas;
  }
  
  //for plotting, get rid of 0's
  double timesCoherentPlot[npoints], voltsCoherentPlot[npoints];
  for (int j=0;j<npoints-2*factorToLeaveOut*upSampleFactor;j++){
    timesCoherentPlot[j]=timesCoherent[j+factorToLeaveOut*upSampleFactor];
    voltsCoherentPlot[j]=voltsCoherent[j+factorToLeaveOut*upSampleFactor];
  }
  
  TGraph *grCoherentWaveformDeconVertUnfilt;
  TGraph *grCoherentWaveformDeconHorizUnfilt;
  //now fill the waveform
  if (polFlag==0){
    if (grCoherentWaveformDeconVert) delete grCoherentWaveformDeconVert;
    grCoherentWaveformDeconVertUnfilt=new TGraph(npoints-(2*upSampleFactor*factorToLeaveOut)-1,
						 timesCoherentPlot,voltsCoherentPlot);
    grCoherentWaveformDeconVert=simpleNotchFilter(grCoherentWaveformDeconVertUnfilt,1250,1270);
  }
  else{
    if (grCoherentWaveformDeconHoriz) delete grCoherentWaveformDeconHoriz;
    grCoherentWaveformDeconHorizUnfilt=new TGraph(npoints-(2*upSampleFactor*factorToLeaveOut)-1,
						  timesCoherentPlot,voltsCoherentPlot);
    grCoherentWaveformDeconHoriz=simpleNotchFilter(grCoherentWaveformDeconHorizUnfilt,1250,1270);
  }
  
  /////////////////////////
  TGraph *grFFTVert;
  TGraph *grFFTHoriz;
  double *oldX, *oldY;
  int length;

  if (polFlag==0){ 
    oldY = grCoherentWaveformDeconVert->GetY();
    oldX = grCoherentWaveformDeconVert->GetX();
    length=grCoherentWaveformDeconVert->GetN();
  }
  else {
    oldY = grCoherentWaveformDeconHoriz->GetY();
    oldX = grCoherentWaveformDeconHoriz->GetX();
    length=grCoherentWaveformDeconHoriz->GetN();
  }
  
  double deltaT=oldX[1]-oldX[0];
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  int newLength=(length/2)+1;
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  
  //define variables for FFT of coherent waveform
  double frequencyCoherent[newLength];//MHz
  double magnitudeCoherent[newLength];
  frequencyCoherent[0]=0;
  
  for(int i=0;i<newLength;i++) {
    if (i>0) frequencyCoherent[i]=frequencyCoherent[i-1]+deltaF;
    magnitudeCoherent[i]=sqrt(theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im);
    magnitudeCoherent[i]=log10(magnitudeCoherent[i])*10;
    magnitudeCoherent[i]=magnitudeCoherent[i]*2; // go to power
    if (frequencyCoherent[i]<200 || frequencyCoherent[i]>1200) magnitudeCoherent[i]=-100;
  }
  
  if (polFlag==0) grFFTVert=new TGraph(newLength,frequencyCoherent,magnitudeCoherent);
  else  grFFTHoriz=new TGraph(newLength,frequencyCoherent,magnitudeCoherent);
  ////////////

  TGraph *grCoherentWaveformDecon;
  if (polFlag==0) grCoherentWaveformDecon=FFTtools::getInterpolatedGraph(grCoherentWaveformDeconVert, deltaTInt/4.);
  else grCoherentWaveformDecon=FFTtools::getInterpolatedGraph(grCoherentWaveformDeconHoriz, deltaTInt/4.);
    
  if (drawFlag==1){//draw coherent waveform
    if (polFlag==0){
      TCanvas *cCoherentDecon=new TCanvas("cCoherentDecon","cCoherentDecon",800,800);
      cCoherentDecon->Divide(1,2);
      cCoherentDecon->cd(1);
      grCoherentWaveformDecon->GetXaxis()->SetRangeUser(100,180);
      grCoherentWaveformDecon->GetXaxis()->SetTitle("Time (ns)");
      grCoherentWaveformDecon->GetYaxis()->SetTitle("Electric Field (mV/m)");
      grCoherentWaveformDecon->Draw("al");
      cCoherentDecon->cd(2);
      grFFTVert->Draw("al");
    }
    else{
      TCanvas *cCoherentDeconHoriz=new TCanvas("cCoherentHorizDecon","cCoherentHorizDecon",800,800);
      cCoherentDeconHoriz->Divide(1,2);
      cCoherentDeconHoriz->cd(1);
      grCoherentWaveformDecon->GetXaxis()->SetRangeUser(100,180);
      grCoherentWaveformDecon->GetXaxis()->SetTitle("Time (ns)");
      grCoherentWaveformDecon->GetYaxis()->SetTitle("Electric Field (mV/m)");
      grCoherentWaveformDecon->Draw("al");
      cCoherentDeconHoriz->cd(2);
      grFFTHoriz->Draw("al");
    }
  }

  //print stuff for Andres
  if (printAndresFlag==1){
    char filenameAndres[150];
    if (polFlag==0) sprintf(filenameAndres,"coherentDeconvolved_V_%d.txt",fHeadPtr->eventNumber);
    else sprintf(filenameAndres,"coherentDeconvolved_H_%d.txt",fHeadPtr->eventNumber);
    ofstream foutAndres(filenameAndres);
    if (polFlag==0){
      double *xAndres=grCoherentWaveformDeconVert->GetX();
      double *yAndres=grCoherentWaveformDeconVert->GetY();
      for (int i=0;i<grCoherentWaveformDeconVert->GetN();i++)
	if (xAndres[i]>=100 && xAndres[i]<=180) 
	  foutAndres<<xAndres[i]<<"\t"<<yAndres[i]<<endl;
    }
    else{
      double *xAndres=grCoherentWaveformDeconHoriz->GetX();
      double *yAndres=grCoherentWaveformDeconHoriz->GetY();
      for (int i=0;i<grCoherentWaveformDeconHoriz->GetN();i++)
	if (xAndres[i]>=100 && xAndres[i]<=180) 
	  foutAndres<<xAndres[i]<<"\t"<<yAndres[i]<<endl;
    }
    foutAndres.close();
  }

  if (polFlag==0)
    return grCoherentWaveformDeconVert;
  else return grCoherentWaveformDeconHoriz;
  
}
////////////////////////////////////
TGraph *MyCorrelator::makeCoherentlySummedWaveform(int myEventNumber, double peakTheta, 
						   double peakPhi, int nadirFlag, int polFlag,
						   int nantennas, int drawFlag)
{
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  TGraph *grEvCoherent[nantennas];
  int upSampleFactor=8;
  Double_t deltaTInt=1./(2.6*upSampleFactor);

  vector<int> whichAntennasToUse (nantennas,0);
  if(printFlag==1) cout<<"peakPhi is "<<peakPhi<<"\n";
  getClosestNAntennas(nantennas, peakPhi, whichAntennasToUse,nadirFlag);
  for(int k=0;k<NUM_ANTS_WITH_NADIRS;k++){
    if(k<nantennas){
      whichAntennasCoherent[k]=whichAntennasToUse[k];
      cout<<"Antennas to use are "<<whichAntennasToUse[k]<<"\n";
      if(printFlag==1) cout<<"Antennas to use are "<<whichAntennasToUse[k]<<"\n";
    }
    else{
      whichAntennasCoherent[k]=-1;
    }
  }//k

  /*
  for (int i=0;i<nantennas;i++){
    int ant=whichAntennasToUse[i];
    if (polFlag==0) grEvCoherent[i]=FFTtools::getInterpolatedGraph(grEvUnfiltered[ant], deltaTInt); //vertical pol grEvUnfiltered
    else grEvCoherent[i]=FFTtools::getInterpolatedGraph(grEvHorizUnfiltered[ant], deltaTInt); //horizontal pol
    }*/
  //Check to see if should use filtered!
  char namer[256];
  for (int i=0;i<nantennas;i++){
    int ant=whichAntennasToUse[i];
    // sprintf(namer,"pre_grEvCoherent_%i",ant);
    //DrawFreqDomain(grEv[ant],myEventNumber,namer);
    if (polFlag==0) grEvCoherent[i]=FFTtools::getInterpolatedGraph(grEv[ant], deltaTInt); //vertical pol grEvUnfiltered
    else grEvCoherent[i]=FFTtools::getInterpolatedGraph(grEvHoriz[ant], deltaTInt); //horizontal pol
    //sprintf(namer,"grEvCoherent_%i",ant);
    //DrawFreqDomain(grEvCoherent[i],myEventNumber,namer);
  }
 

  int npoints=260*upSampleFactor;
  double voltsCoherent[npoints];
  double timesCoherent[npoints];
  for (int j=0;j<npoints;j++){
    timesCoherent[j]=0;
    voltsCoherent[j]=0;
  }
  
  int nsToLeaveOut=20;
  int factorToLeaveOut=int(nsToLeaveOut*2.6);
  //interpolate!!
  for (int i=0;i<nantennas;i++){
    int ant=whichAntennasToUse[i];
    //cout<<"antenna "<<whichAntennasToUse[i]+1<<" used for coherently summed waveform"<<endl;
    double timeExpected=getDeltaTExpected(ant, whichAntennasToUse[0], 
							  peakPhi, peakTheta);
    
    //cout<<"time expected between this and 0: "<<timeExpected<<endl;
    //get the value of graph at this delay.
    double tVal1, tVal2, val1, val2, voltVal;			
    double tVal0,val0, tValend, valend;
    grEvCoherent[i]->GetPoint(0,tVal0,val0);//get start of time and volt(?)
    grEvCoherent[i]->GetPoint(grEvCoherent[i]->GetN()-1,tValend,valend);//get end of time and volt(?)
    double lengthOfTrace=tValend-tVal0;
    for (int j=factorToLeaveOut*upSampleFactor;j<npoints-factorToLeaveOut*upSampleFactor;j++){//start upsample and iterpolate
      int bin1=int((timeExpected+j*deltaTInt-tVal0)/lengthOfTrace*grEvCoherent[i]->GetN());
      int bin2=bin1+1;
      //if (j==factorToLeaveOut*upSampleFactor) cout<<"bin1: "<<bin1<<", bin2: "<<bin2<<endl;
      if (bin1>=0 && bin2<=grEvCoherent[i]->GetN()){
	grEvCoherent[i]->GetPoint(bin1,tVal1,val1);
	grEvCoherent[i]->GetPoint(bin2,tVal2,val2);
	voltVal=(val2-val1)*(timeExpected+j*deltaTInt-tVal1)/(tVal2-tVal1)+val1;//interpolate 2 bins 
	voltsCoherent[j]+=voltVal;//upsample adds bins at 0?
      }
    }
  }
  for (int i=0;i<nantennas;i++){
    delete grEvCoherent[i];
  }
  ofstream myfile;
  myfile.open("CWCoherentWaveform.txt");
  myfile<<"Time (ns) \t Volts (mV) \n";
  for (int j=factorToLeaveOut*upSampleFactor;j<npoints-factorToLeaveOut*upSampleFactor;j++){
    timesCoherent[j]=deltaTInt*j;
    voltsCoherent[j]=voltsCoherent[j]/nantennas;
    myfile<<timesCoherent[j]<<"\t"<<voltsCoherent[j]<<"\n";
  }
  myfile.close();
  //for plotting, get rid of 0's
  double timesCoherentPlot[npoints], voltsCoherentPlot[npoints];
  for (int j=0;j<npoints-2*factorToLeaveOut*upSampleFactor;j++){
    timesCoherentPlot[j]=timesCoherent[j+factorToLeaveOut*upSampleFactor];
    voltsCoherentPlot[j]=voltsCoherent[j+factorToLeaveOut*upSampleFactor];
  }
  
  if (polFlag==0){
    if (grCoherentWaveformVert)delete grCoherentWaveformVert;

    grCoherentWaveformVert=new TGraph(npoints-(2*upSampleFactor*factorToLeaveOut)-1,
				      timesCoherentPlot,voltsCoherentPlot);
  }
  else{
    if (grCoherentWaveformHoriz) delete grCoherentWaveformHoriz;
    grCoherentWaveformHoriz=new TGraph(npoints-(2*upSampleFactor*factorToLeaveOut)-1,
				       timesCoherentPlot,voltsCoherentPlot);
  }
  
  
  
  if (drawFlag==1){//draw coherent waveform
    char printer[256];
    if (polFlag==0){
      TCanvas *cCoherent=new TCanvas("cCoherent","cCoherent",800,300);
      cCoherent->cd(0);
      grCoherentWaveformVert->GetXaxis()->SetTitle("Time (ns)");
      grCoherentWaveformVert->GetYaxis()->SetTitle("Voltage (mV)");
      grCoherentWaveformVert->Draw("al");
      
      
      sprintf(printer,"Coherent_%i.png",myEventNumber);
      cCoherent->Print(printer);
      
    }
    else{
      TCanvas *cCoherentHoriz=new TCanvas("cCoherentHoriz","cCoherentHoriz",800,300);
      cCoherentHoriz->cd(0);
      grCoherentWaveformHoriz->GetXaxis()->SetTitle("Time (ns)");
      grCoherentWaveformHoriz->GetYaxis()->SetTitle("Voltage (mV)");
      grCoherentWaveformHoriz->Draw("al");
      // cCoherentHoriz->Print("cCoherentHoriz.png");
      sprintf(printer,"CoherentHoriz_%i.png",myEventNumber);
      cCoherentHoriz->Print(printer);
    }
  }

  
  //print stuff for Andres
  if (printAndresFlag==1){
    char filenameAndres[150];
    if (polFlag==0) sprintf(filenameAndres,"coherent_V_%d.txt",fHeadPtr->eventNumber);
    else sprintf(filenameAndres,"coherent_H_%d.txt",fHeadPtr->eventNumber);
    ofstream foutAndres(filenameAndres);
    if (polFlag==0){
      double *xAndres=grCoherentWaveformVert->GetX();
      double *yAndres=grCoherentWaveformVert->GetY();
      for (int i=0;i<grCoherentWaveformVert->GetN();i++)
	if (xAndres[i]>=20 && xAndres[i]<=80) 
	  foutAndres<<xAndres[i]<<"\t"<<yAndres[i]<<endl;
    }
    else{
      double *xAndres=grCoherentWaveformHoriz->GetX();
      double *yAndres=grCoherentWaveformHoriz->GetY();
      for (int i=0;i<grCoherentWaveformHoriz->GetN();i++)
	if (xAndres[i]>=20 && xAndres[i]<=80) 
	  foutAndres<<xAndres[i]<<"\t"<<yAndres[i]<<endl;
    }
    foutAndres.close();
  }
   if (polFlag==0)
    return grCoherentWaveformVert;
  else return grCoherentWaveformHoriz;
}
///////////////////////////////////////
int MyCorrelator::pointThisEvent(int eventNumber, int drawMaps, TNtuple *ndata, TNtuple *ndata2, TNtuple *ndata3, TNtuple *ndata4,
				 double &peakThetaFinal, double &peakPhiFinal, int whichPolarization,int &xCorPassFlag,
				 int &ratioOfPeaksPassFlag, int &elevationAnglePassFlag, int &peakCrossCorrFlag,
				 int &polFractionFlag, int &peakHilbertFlag,int &triggerFlag, double &finaltheta)
{ 
  
  gStyle=babar;
  if (!fEventTree) initialize();
  
  //GetHealPixMap();
  //initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  if (eventEntryGottenFlag!=(int)eventNumber) eventEntryGottenFlag=getEventEntry();
  cout<<"evententrygotten flag is "<<eventEntryGottenFlag<<"\n";
  int deconFlag=0;
  polToggle=whichPolarization;
  
  //random.SetSeed(fHeadPtr->eventNumber);
  //define taylor dome theta and phi
  double deltaTheta=0;
  double deltaPhi=0;
  double thetaMap=0;
  double phiMap=0;
  double phi=0;
  double theta=0;
  double mcmphi=0;
  double mcmtheta=0;
  cout<<"finished setting up \n";
  fUsefulAdu5Ptr->getThetaAndPhiWaveWillyBorehole(mcmtheta,mcmphi);
 
  fUsefulAdu5Ptr->getThetaAndPhiWaveTaylorDome(theta,phi);
  if (printFlag==1){ 
    cout<<"WZ Event Number: "<<eventNumber<<endl;
    cout<<"Taylor dome theta: "<<theta*rad2deg<<", phi: "<<phi*rad2deg<<endl;  
    cout<<"Willy theta: "<<mcmtheta*rad2deg<<", phi: "<<mcmphi*rad2deg<<endl;
    cout<<"WZ realtime: "<<fHeadPtr->realTime<<endl;
    cout<<"WZ heading: "<<fAdu5APatPtr->heading<<"\n";
    cout<<"altitude: "<<fAdu5APatPtr->altitude<<"\n";
  }
  
  //get the graphs for this event
  int nadirFlag=0;
  int nantennas=0;
  if (isNadirRFCMOn(eventNumber)){ 
    nadirFlag=1;
    nantennas=NUM_ANTS_WITH_NADIRS;
  }
  else{
    nadirFlag=0;
     nantennas=NUM_ANTS_NO_NADIRS;
  }
 
  int windowWaveformFlag=0;
  double snrPeakAnt=0;
  double maxSignalPeak=0;
  int peakAnt=0;
  double distanceTD=0;
  double deltaTTD=0;
  double thetaTD=0; 
  double phiTD=0; 
  double distanceMcM=0;
  double secondThetaMap=0;
  double secondPhiMap=0;
  double thetaWilly=0;
  double phiWilly=0;
  int phiMaskArray[NUM_PHI_SECTORS];
  thetaTD=theta*rad2deg*-1;
  phiTD=phi*rad2deg;
  thetaWilly=mcmtheta*rad2deg*-1;
  phiWilly=mcmphi*rad2deg;
 
  SNRpeak =-1;
  getGraphsThisEvent(windowWaveformFlag,snrPeakAnt, maxSignalPeak, peakAnt);
  SNRpeak = snrPeakAnt;
  /*
  for(int n=0;n<40;n++){
    PowerCut[n]=integrateTDPower(grEv[n]);
  }
  */

  if (printFlag==1) cout<<"snr peak antenna: "<<snrPeakAnt<<endl;
  cout<<"maxSignalPeak is "<<maxSignalPeak<<"\n";
  if (printFlag==1) cout<<"peakAnt is "<<peakAnt<<"\n";
  int dummyTaylorFlag=isTaylor(eventNumber, distanceTD, deltaTTD);
  int thisPhiMask=isPhiMaskingOn(eventNumber,phiMaskArray); 
  distanceMcM=getMcMDistance(eventNumber);
  if (printFlag==1) cout<<"distance to mcmurdo: "<<distanceMcM<<endl;
  dummyTaylorFlag=0;
  
  //define variables needed for filter
  double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI];
  int myEventNumber=eventNumber;
  double peakVal=0;
  double peakTheta=0;
  double peakPhi=0;
  double mapCorValRefined[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI];
  double peakValRefined=0;
  int peakThetaBin=0;
  int peakPhiBin=0;
  double refinedPeakTheta=0;
  double refinedPeakPhi=0;
  double FWHMTheta=0;
  double FWHMPhi=0;
  double peakThetaInterp=0;
  double peakPhiInterp=0;
  double peakValInterp=0;
  TGraph *grCoherent;
  TGraph *grCoherentHoriz;
  TGraph *grDeconvolvedCoherent;
  TGraph *grDeconvolvedCoherentHoriz;
  double bandWidth=26;//MHz band (really half bandwidth) for adaptive filtering
  double dBCut=4.0;
  cout<<"CHANGED dBCut to 4! \n";
  int nfreq=5;//change this
  double frequencies[nfreq];
  double frequenciesHoriz[nfreq];
  double bandwidth[nfreq];
  double bandwidthHoriz[nfreq];
  double magPeak[nfreq];
  double magPeakHoriz[nfreq];
  
  double baseX[2000];
  double baseY[2000];
  double baseXHoriz[2000];
  double baseYHoriz[2000];
 
  int didIFilter=0;
  int didIFilterHoriz=0;
  int didIFilterAboveSatellite=0;
  int didIFilterAboveSatelliteHoriz=0;
  int triggerOnlyFlag=0;
  // cout<<"trying triggered phis! \n";
  // cout<<"trying triggered phis + -  1! \n";
  //  int nantennasToUse=20;
  int triggerOrPhiMaskFlag=1;
  int phiMaskFlag=0;
  int hwTriggerFlag=1;
 
  /*
  int xCorPassFlag=0;
  int ratioOfPeaksPassFlag=0;
  int elevationAnglePassFlag=0;
  int peakCrossCorrFlag=0;
  int polFractionFlag=0;
  int peakHilbertFlag=0;
  */
  float meanFreqHoriz, meanFreqVert;
  strongCWFlag=0;
  CWheight=0.;
  //check for saturated channels (and set up saturated array of channel indeces)
  int nSaturated=isChannelSaturated(myEventNumber);
  nSaturated*=1;
  cout<<"zeroing some stuff \n";
  vector< vector<double> >antennaFreq(NUM_ANTS_WITH_NADIRS, vector<double>(50));
  vector< vector<double> >antennaFreqHoriz(NUM_ANTS_WITH_NADIRS, vector<double>(50));
  vector< vector<double> >antennaBandwidth(NUM_ANTS_WITH_NADIRS, vector<double>(50));
  vector< vector<double> >antennaBandwidthHoriz(NUM_ANTS_WITH_NADIRS, vector<double>(50));
  vector< vector<double> >PeakMag(NUM_ANTS_WITH_NADIRS, vector<double>(50));
  vector< vector<double> >PeakMag_backup(NUM_ANTS_WITH_NADIRS, vector<double>(50));
  vector< vector<double> >PeakMagHoriz(NUM_ANTS_WITH_NADIRS, vector<double>(50));
 
  //SWITCH TO USE ABBY'S NOTCH FILTER:
  //newnotchflag=0; abby's
  //newnotchflag=1; brian's

 
  int nAntennasToUse=9;//9;
  if(newnotchflag==0){
    nAntennasToUse=39;
  }
  if(printFlag==1) cout<<"nAntennasToUse is "<<nAntennasToUse<<"\n";
  vector<int> whichAntennasToUse (nAntennasToUse,0);
  int antused;
  int n;

  
  if(groupFlag==0) {
    getGroupsofAntennas(nAntennasToUse,nadirFlag);
    groupFlag=1;
  }

  
  if(notchFilterFlag==5){
    notchFilterFlag=4;
  }


  for(int antenna_groups=0;antenna_groups<antenna_groups_size;antenna_groups++){//go through all phi sectors
      
    peakPhi = (double) unique_phis[antenna_groups];
    
    if (printFlag==1) cout<<"antenna_groups is "<<antenna_groups<<" and peak phi is "<<peakPhi<<"\n";
      for(int k=0;k<nAntennasToUse;k++){
	whichAntennasToUse[k]=antenna_group_holder[antenna_groups][k];
	//cout<<"using antenna "<<whichAntennasToUse[k]<<"\n";
	
      }
     
      //find frequencies to cut, put them into frequencies
      adaptiveFilterPartialPayload(0, dBCut, nfreq,frequencies,bandwidth,magPeak,drawMaps, bandWidth, peakPhi, nAntennasToUse,whichAntennasToUse,
				   meanFreqVert,myEventNumber,nadirFlag, antenna_groups);
      adaptiveFilterPartialPayload(1, dBCut, nfreq,frequenciesHoriz,bandwidthHoriz,magPeakHoriz,drawMaps, bandWidth, peakPhi, nAntennasToUse,
				   whichAntennasToUse, meanFreqHoriz, myEventNumber,nadirFlag, antenna_groups);
     
      //put freqs cut into each antenna used

      
      for(int j=0;j<nAntennasToUse;j++){

	antused = whichAntennasToUse[j];
	n=0;
	

	for(int j1=0;j1<50;j1++){
	  if(antennaFreq[antused][j1]<=0 && n<nfreq){ 
	    antennaFreq[antused][j1]=frequencies[n];
	    antennaBandwidth[antused][j1]=bandwidth[n];
	    PeakMag[antused][j1]=magPeak[n];
	    antennaFreqHoriz[antused][j1]=frequenciesHoriz[n];
	    antennaBandwidthHoriz[antused][j1]=bandwidthHoriz[n];
	    PeakMagHoriz[antused][j1]=magPeakHoriz[n];
	  
	    n++;
	  }//antennaFreq==0
	 
	}//j1
	
      }//j
      
      
      
  }//antenna_groups
  vector<double> uniquefreqs;
  vector<double> uniquefreqs1;
  vector<double> uniquefreqsHoriz;
  vector<double> uniquefreqsHoriz1;

  vector<double> uniquebandwidth;
  vector<double> uniquebandwidthHoriz;
   
  vector<double> uniquePhase;
  vector<double> uniquePhaseHoriz;

  vector<double> uniquePhase_bandwidth;
  vector<double> uniquePhaseHoriz_bandwidth;
  int triggeredPhi[NUM_PHI_SECTORS];
  int triggeredAnt[NUM_ANTS_WITH_NADIRS];
  for(int i=0;i<NUM_ANTS_WITH_NADIRS;i++){
    triggeredAnt[i]=-1;
  }
  getTriggeredL2Phi(fHeadPtr,triggeredPhi);
  getTriggeredAnt(triggeredPhi,triggeredAnt);
 
  //do adaptive filter
    
  vector<int> whichAntennasToUse_1 (1,0); //was 40;
  double baseline[2000];
  
  double frequencyArray[2000];
  double magFFT[2000];
  double magFFT2[2000];
  double baselineHoriz[2000];
  double magFFT2Horiz[2000];
  double frequencyArrayHoriz[2000];
  double magFFTHoriz[2000];
  double magFFTHoriz2[2000];
 
 

  ///////

  
  double ratioFirstToSecondPeak;
  double peakHilbertCoherent;
  double polFractionCoherent;
  double deltamcmTheta;
  double deltamcmPhi;
  double mapSNR;
  double snrCoherent;
  float hwTriggerAngle;
  double headingOfThisEvent;
  double tertiaryPeakTheta;
  double tertiaryPeakPhi;
  double polAngleCoherent;


 
  int satellite_flag=0;
  int satellite_flag_Horiz=0;
  double satellite_freq=0.;
  double satellite_bandwidth=0.;
  double satellite_freq2=0.;
  double satellite_bandwidth2=0.;
  double satellite_freq_Horiz=0.;
  double satellite_bandwidth_Horiz=0.;
  double satellite_freq_Horiz2=0.;
  double satellite_bandwidth_Horiz2=0.;
  
  double SNR_test;
  double rmsNoise_test;

  ratioFirstToSecondPeak=0.;
  peakHilbertCoherent=0.;
  polFractionCoherent=0.;
  deltamcmTheta=0.;
  deltamcmPhi=0.;
  mapSNR=0.;
  snrCoherent=0.;
  hwTriggerAngle=0.;
  headingOfThisEvent=0.;
  tertiaryPeakTheta=0.;
  tertiaryPeakPhi=0.;
  polAngleCoherent=0.;
  satellite_flag=0;
     
  for(int i=0;i<40;i++){
    //cout<<" ant is "<<i<<" phi is "<<fUPGeomTool->getPhiFromAnt(i)<<"\n";
    //grEv_backup[i]=new TGraph(grEv[i]->GetN(),grEv[i]->GetX(),grEv[i]->GetY());
    for(int j=0;j<50;j++){
      PeakMag_backup[i][j]=PeakMag[i][j];
    }
  }
    
   for(int i=0;i<40;i++){//loop over each antenna
    if(i!=1 && saturatedChannels[i+whichPolarization*NUM_ANTS_WITH_NADIRS]==0 ){//dont use bad antennas
      //grEv[i] = new TGraph(grEv_backup[i]->GetN(),grEv_backup[i]->GetX(),grEv_backup[i]->GetY());
      for(int j=0;j<50;j++){
	PeakMag[i][j] = PeakMag_backup[i][j];
      }
     
      if(printFlag==1) cout<<"i is "<<i<<"\n";
      uniquefreqs.clear();//zero stuff
      uniquefreqs1.clear();
      uniquefreqsHoriz.clear();
      uniquefreqsHoriz1.clear();
      uniquebandwidth.clear();
      uniquebandwidthHoriz.clear();
      uniquePhase.clear();
      uniquePhaseHoriz.clear();
      uniquePhase_bandwidth.clear();
      uniquePhaseHoriz_bandwidth.clear();
 
      satellite_freq=260;
      satellite_bandwidth=26;
      satellite_freq2=0;
      satellite_bandwidth2=0;
      satellite_freq_Horiz=260;
      satellite_bandwidth_Horiz=26;
      satellite_freq_Horiz2=0;
      satellite_bandwidth_Horiz2=0;
      satellite_flag=0;
      satellite_flag_Horiz=0;
      for(int firstindex=0;firstindex<300;firstindex++){//zero some stuff
	for(int secondindex=0; secondindex<2;secondindex++){
	  theFFTarray[firstindex][secondindex]=0.;
	  theFFTarrayHoriz[firstindex][secondindex]=0.;
	}
      }
      
       for(int arrayctr=0;arrayctr<2000;arrayctr++){
	baseX[arrayctr]=0.;
	baseXHoriz[arrayctr]=0.;
	baseY[arrayctr]=0.;
	baseYHoriz[arrayctr]=0.;
	baseline[arrayctr]=0.;
  
	frequencyArray[arrayctr]=0.;
	magFFT[arrayctr]=0.;
	magFFT2[arrayctr]=0.;
	magFFT2Horiz[arrayctr]=0.;
	baselineHoriz[arrayctr]=0.;
	
	frequencyArrayHoriz[arrayctr]=0.;
	magFFTHoriz[arrayctr]=0.;
	magFFTHoriz2[arrayctr]=0.;
	

 
      }  
      
      
      //Get rid of duplicated freqs.
       if(printFlag==1) cout<<"about to get freqs to cut \n";
       GetFrequenciestoCut(i, antennaFreq,antennaBandwidth,PeakMag,uniquefreqs1,uniquebandwidth,nfreq,uniquePhase,uniquePhase_bandwidth);//fill uniquefreq with freq to cut!
      
       GetFrequenciestoCut(i, antennaFreqHoriz,antennaBandwidthHoriz,PeakMagHoriz,uniquefreqsHoriz1,uniquebandwidthHoriz,nfreq,uniquePhaseHoriz,uniquePhaseHoriz_bandwidth);
       if(printFlag==1) cout<<"got freqs to cut \n";
      
       //GetSatelliteFrequenciestoCut(i,uniquefreqs1,uniquebandwidth,satellite_freq,satellite_bandwidth,satellite_freq2,satellite_bandwidth2,satellite_flag);
       //GetSatelliteFrequenciestoCut(i,uniquefreqsHoriz1,uniquebandwidthHoriz,satellite_freq_Horiz,satellite_bandwidth_Horiz,satellite_freq_Horiz2,satellite_bandwidth_Horiz2,satellite_flag_Horiz);

       
     
      for(int f=0;f<2000;f++){//zero stuff

	frequencyArray[f]=0.;
	frequencyArrayHoriz[f]=0.;					
	
	magFFT2[f]=0.;
	baseline[f]=0.;
	
	magFFTHoriz2[f]=0.;
	baselineHoriz[f]=0.;
      }
      whichAntennasToUse_1[0]=0;


      ///DO some stuff

      whichAntennasToUse_1[0]=i;//need vector to for function (function used in other places)
      
      //set grCoherent to correct antenna
      GetBaselineperPhi(0, baseline,1, whichAntennasToUse_1);
      GetBaselineperPhi(1, baseline,1, whichAntennasToUse_1);
      //
      
      GetFFTandBaseline(1,whichAntennasToUse_1, 0,frequencyArray, baseX,baseY,magFFT2,nfreq,dBCut);
      GetFFTandBaseline(1,whichAntennasToUse_1, 1,frequencyArrayHoriz, baseXHoriz,baseYHoriz,magFFT2Horiz,nfreq,dBCut);

      //UniquePhase hold mostly the same info as uniquefreqs, but the geom method blows up at center of notch, so we need to keep track of those. 
      cout<<"ant is "<<i<<" ";
      if(notchFilterFlag!=4){
	for(int k=0;k<(int)uniquefreqs1.size();k++){
	  cout<<"cutting "<<uniquefreqs1[k]-uniquebandwidth[k]<<" to "<<uniquefreqs1[k]+uniquebandwidth[k]<<"\n";
	if (uniquefreqs1[k]!=-1 && printFlag==1) cout<<"WZ going to cut: "<<uniquefreqs1[k]<<" MHz in Vert for i == "<<i<<endl;
	if (uniquefreqs1[k]!=-1) didIFilter++;
	if (uniquefreqs1[k]!=-1 && uniquefreqs1[k]>286) didIFilterAboveSatellite++;
	bandWidth = uniquebandwidth[k];	

	if(printFlag==1)	cout<<"bandwidth is "<<bandWidth<<"\n";

	if(uniquefreqs1[k]>0) applyAdaptiveFilter_singleAnt(uniquefreqs1[k], bandWidth,0,i,baseY);//apply filter in Vert
	
      }//k
      
      
	for(int k=0;k<(int)uniquefreqsHoriz1.size();k++){
	  //cout<<"cutting "<<uniquefreqsHoriz1[k]-uniquebandwidthHoriz[k]<<" to "<<uniquefreqsHoriz1[k]+uniquebandwidthHoriz[k]<<"\n";
	if (uniquefreqsHoriz1[k]!=-1 && printFlag==1) cout<<"WZ going to cut: "<<uniquefreqsHoriz1[k]<<" MHz in Horiz for i == "<<i<<endl;
	//	if (uniquefreqsHoriz1[k]!=-1) didIFilter++;
	//if (uniquefreqsHoriz1[k]!=-1 && uniquefreqsHoriz1[k]>286) didIFilterAboveSatellite++;
	bandWidth = uniquebandwidthHoriz[k];

	if(uniquefreqsHoriz1[k]>0) applyAdaptiveFilter_singleAnt(uniquefreqsHoriz1[k],bandWidth,1,i,baseYHoriz);//apply filter in Horiz
	
      }//k
      
	if(phase_flag==3 || phase_flag==4){
	  if(uniquefreqs1.size()>0) GeomMethod(i,0,uniquefreqs1, uniquebandwidth,uniquePhase);
	  if(uniquefreqsHoriz1.size()>0) GeomMethod(i,1,uniquefreqsHoriz1, uniquebandwidthHoriz,uniquePhaseHoriz);
	}

	/*

      vector<double> satellite_freq_vec(1,0);
      vector<double> satellite_freq_vec2(1,0);
      
      vector<double> satellite_bandwidth_vec(1,0);
      vector<double> satellite_bandwidth_vec2(1,0);
     
      if(satellite_flag==0){
	applySatelliteFilter(satellite_freq,satellite_bandwidth,i,0,baseY); //apply Satellite Filter
	satellite_freq_vec[0]=satellite_freq;
	satellite_bandwidth_vec[0]=satellite_bandwidth;
	//ShiftPhase1(i,0,satellite_freq_vec, satellite_bandwidth_vec,satellite_freq_vec);
      }
      else if(satellite_flag==2){
	applySatelliteFilter(satellite_freq,satellite_bandwidth,i,0,baseY); //apply Satellite Filter
	satellite_freq_vec[0]=satellite_freq;
	satellite_bandwidth_vec[0]=satellite_bandwidth;
	//ShiftPhase1(i,0,satellite_freq_vec, satellite_bandwidth_vec,satellite_freq_vec);


	applySatelliteFilter(satellite_freq2,satellite_bandwidth2,i,0,baseY); //apply Satellite Filter
	satellite_freq_vec2[0]=satellite_freq2;
	satellite_bandwidth_vec2[0]=satellite_bandwidth2;
	//ShiftPhase1(i,0,satellite_freq_vec2, satellite_bandwidth_vec2,satellite_freq_vec2);
      }
      
      if(satellite_flag_Horiz==0){
	applySatelliteFilter(satellite_freq_Horiz,satellite_bandwidth_Horiz,i,1,baseYHoriz); //apply Satellite Filter
	satellite_freq_vec[0]=satellite_freq_Horiz;
	satellite_bandwidth_vec[0]=satellite_bandwidth_Horiz;
	//ShiftPhase1(i,1,satellite_freq_vec, satellite_bandwidth_vec,satellite_freq_vec);
      }
      else if(satellite_flag==2){
	applySatelliteFilter(satellite_freq_Horiz,satellite_bandwidth_Horiz,i,1,baseYHoriz); //apply Satellite Filter
	satellite_freq_vec[0]=satellite_freq_Horiz;
	satellite_bandwidth_vec[0]=satellite_bandwidth_Horiz;
	//ShiftPhase1(i,1,satellite_freq_vec, satellite_bandwidth_vec,satellite_freq_vec);

	applySatelliteFilter(satellite_freq_Horiz2,satellite_bandwidth_Horiz2,i,1,baseYHoriz); //apply Satellite Filter
	satellite_freq_vec2[0]=satellite_freq_Horiz;
	satellite_bandwidth_vec2[0]=satellite_bandwidth_Horiz;
	//ShiftPhase1(i,1,satellite_freq_vec2, satellite_bandwidth_vec2,satellite_freq_vec2);
      }
	*/
      
      }//notchfilterFlag!=4
      
     
    delete grVertCoherentBaseline;
    delete grHorizCoherentBaseline;
      
    grVertCoherentBaseline=0;
    grHorizCoherentBaseline=0;
      

     
     //Getting SNR of individual Antenna
    SNR_ant[i]=0.;
    SNR_ant_triggered[i]=0.;
    SNR_test = getSNR(grEv[i],rmsNoise_test);
    SNR_ant[i]=SNR_test;
    
    if(triggeredAnt[i]>0){
      SNR_ant_triggered[i]=SNR_test;
    }
    else{
      SNR_ant_triggered[i]=-1;
    }
   
    }//i!=1//dont use antennas that are bad: This is always antenna 1 and any saturated antennas
   }//i
 
  //do correlation map
 
   for(int n=0;n<40;n++){
    PowerCut[n]=PowerCut[n] - integrateTDPower(grEv[n]);
  }
   char namer[256];
   if(drawMaps==1){
     sprintf(namer,"5");
     DrawFreqDomain(grEv[5], myEventNumber,namer);
  }
  
  doCorrelationMap(myEventNumber, mapCorVal,triggerOnlyFlag,nadirFlag,whichPolarization);
  if (drawMaps==1) drawCorrelationMap(mapCorVal,myEventNumber,0);
  
  //get peak of main map
  findPeakOfMap(mapCorVal,peakVal,peakTheta,peakPhi);
  if (printFlag==1) cout<<"peak val: "<<peakVal<<", peak theta: "<<peakTheta<<", peak phi: "<<peakPhi<<endl;
  cout<<"peakTheta before refined is "<<peakTheta<<"\n";
  cout<<"peakPhi is "<<peakPhi<<"\n";
  //now do a refined map around the peak.
  doRefinedMap(myEventNumber, mapCorValRefined,peakTheta,peakPhi, triggerOnlyFlag, nadirFlag,whichPolarization);
  if (drawMaps==1) drawRefinedMap(mapCorValRefined,peakTheta,peakPhi,myEventNumber);

  //get peak of refined map
  findPeakOfRefinedMap(mapCorValRefined,peakValRefined,peakThetaBin,peakPhiBin);
  
  if(cos_reconstruction_flag==1){
    refinedPeakTheta = peakThetaBin*step_size_cos_fine + max_theta_cos;
  }
  else{
    refinedPeakTheta = peakThetaBin*FINE_BIN_SIZE+(peakTheta-thetaStartOffset);
  }
  refinedPeakPhi=peakPhiBin*FINE_BIN_SIZE+(peakPhi-phiStartOffset);
   cout<<"Refined peak theta: "<<refinedPeakTheta<<", peak Phi: "
			<<refinedPeakPhi<<", peakVal: "<<peakValRefined<<endl;
  if (printFlag==1) cout<<"Refined peak theta: "<<refinedPeakTheta<<", peak Phi: "
			<<refinedPeakPhi<<", peakVal: "<<peakValRefined<<endl;
  
  //get the peak using interpolation method
  cout<<"peakTheta is "<<peakTheta<<" refinedPeakTheta is "<<refinedPeakTheta<<"\n";
  cout<<"refinedpeakPhi and peakPhi are "<<refinedPeakPhi<<" "<<peakPhi<<"\n";
  //doInterpolationPeakFinding(mapCorValRefined, peakValInterp, peakThetaInterp, peakPhiInterp, FWHMTheta, FWHMPhi, peakTheta, peakPhi,0); //drawMaps);
  doInterpolationPeakFinding(mapCorValRefined, peakValInterp, peakThetaInterp, peakPhiInterp, FWHMTheta, FWHMPhi, peakTheta, peakPhi,0); //drawMaps);
  if (printFlag==1) cout<<"peak val interp: "<<peakValInterp<<endl;
  if(cos_reconstruction_flag==1){
    cout<<"peakThetaInterp is "<<peakThetaInterp<<"\n";
    peakThetaInterp = acos(peakThetaInterp);
    peakThetaInterp = pi/2 - peakThetaInterp;
    peakThetaInterp = peakThetaInterp*rad2deg;
  }
  cout<<"peakThetaInterp is "<<peakThetaInterp<<"\n";
  cout<<"peakPhiInterp is "<<peakPhiInterp<<"\n";
  headingOfThisEvent=getHeadingOfEvent(myEventNumber, peakPhiInterp);
  if (printFlag==1){
    cout<<"heading of event peak: "<<headingOfThisEvent<<endl;
    cout<<"unixtime: "<<fHeadPtr->realTime<<endl;
    cout<<"theta: "<<-1.*peakThetaInterp<<endl;
    cout<<"anitalat: "<<fAdu5APatPtr->latitude<<endl;
    cout<<"anitalon: "<<fAdu5APatPtr->longitude<<endl;
    cout<<"anitaheight: "<<fAdu5APatPtr->altitude<<endl;
    
    if (mcWeightEvent!=0){
      cout<<"mcWeightEvent: "<<mcWeightEvent<<endl;
      cout<<"mcLatEvent: "<<mcLatEvent<<endl;
      cout<<"mcLonEvent: "<<mcLonEvent<<endl;
      cout<<"mcAltEvent: "<<mcAltEvent<<endl;
    }
  }  
  //get coherent waveforms
  //cout<<"myEventNumber is "<<myEventNumber<<" nadir flag is "<<nadirFlag<<"\n";
  
  grCoherent=makeCoherentlySummedWaveform(myEventNumber, peakThetaInterp, peakPhiInterp, 
					  nadirFlag, 0, 9,drawMaps);//10
  
  grCoherentHoriz=makeCoherentlySummedWaveform(myEventNumber, peakThetaInterp, peakPhiInterp, 
					       nadirFlag, 1, 9,drawMaps);//10
  
  //HACK!

  for(int i=0;i<NUM_ANTS_WITH_NADIRS;i++){
    SNR_ant_closest[i]=-1;
  }
  /* for(int i=0;i<40;i++){
    
    sprintf(namer,"%i",i);
    DrawFreqDomain(grEv[i], myEventNumber,namer);
    }*/
   /* for(int i=0;i<9;i++){
    SNR_ant_closest[whichAntennasCoherent[i]]=SNR_test;
    CoherentAnts[i] = whichAntennasCoherent[i];
    sprintf(namer,"%i",whichAntennasCoherent[i]);
    DrawFreqDomain(grEv[whichAntennasCoherent[i]], myEventNumber,namer);
    }*/
  
 

  if (deconFlag==1) {
    grDeconvolvedCoherent=makeCoherentlySummedDeconvolvedWaveform(myEventNumber, peakThetaInterp, peakPhiInterp, 
								  nadirFlag, 0, 9,drawMaps); //10
    grDeconvolvedCoherentHoriz=makeCoherentlySummedDeconvolvedWaveform(myEventNumber, peakThetaInterp, peakPhiInterp, 
								       nadirFlag, 1, 9,drawMaps);//10
    double peakValDeconvolvedV = FFTtools::getPeakVal(grDeconvolvedCoherent);
    double peakValDeconvolvedH = FFTtools::getPeakVal(grDeconvolvedCoherentHoriz);
    double peak2peakDeconvolvedV = getPeak2Peak(grDeconvolvedCoherent);
    double peak2peakDeconvolvedH = getPeak2Peak(grDeconvolvedCoherentHoriz);
    
    double rmsNoiseDeconV = getRMSOfRange(grDeconvolvedCoherent, 100, 115);
    double rmsNoiseDeconH = getRMSOfRange(grDeconvolvedCoherentHoriz, 100, 115);
    
    if (printFlag==1) cout<<"peak V: "<<peakValDeconvolvedV<<", peak H: "<<peakValDeconvolvedH<<endl;
    if (printFlag==1) cout<<"peak2peak V: "<<peak2peakDeconvolvedV<<", peak2peak H: "<<peak2peakDeconvolvedH<<endl;
    if (printFlag==1) cout<<"rms Noise V: "<<rmsNoiseDeconV<<", rms H: "<<rmsNoiseDeconH<<endl;
  }
  
 
  int peakAntDummy=getPeakAntenna(myEventNumber, nantennas);
  
  peakAntDummy=0;
  sigmaIndex=int(snrAfterFilter);
 
  if (sigmaIndex>14) sigmaIndex=14;
 
  
  snrCoherent=getSNR(grCoherent,rmsNoiseCoherent);
  SNR_ant_coherent=0.;
  SNR_ant_coherent=snrCoherent;
  
  peakHilbertCoherent=0;
  peakHilbertCoherent=getPeakHilbert(grCoherent);
  if (printFlag==1) cout<<"peak hilbert coherent: "<<peakHilbertCoherent<<endl;

  if(drawMaps==1){
    sprintf(namer,"coherent");
    DrawFreqDomain(grCoherent, myEventNumber,namer);
  }


  
  //get stuff to characterize map.
  mapSNR=getSNROfMap(mapCorVal);
  double secondaryPeakVal=0;
  double secondaryPeakTheta=0;
  double secondaryPeakPhi=0;
  findSecondaryPeakOfMap(mapCorVal, secondaryPeakVal, secondaryPeakTheta, secondaryPeakPhi);
  //cout<<"secondary peak val: "<<secondaryPeakVal<<", peakval: "<<peakVal<<endl;
  ratioFirstToSecondPeak=peakVal/secondaryPeakVal;
  secondThetaMap=secondaryPeakTheta*-1;
  secondPhiMap=secondaryPeakPhi;
  if (printFlag==1) cout<<"Ratio of 1st to 2nd Peak in Map: "<<ratioFirstToSecondPeak<<endl;
  if (printFlag==1) cout<<"WZ peakVal is "<<peakVal<<"\n";
  double tertiaryPeakVal=0;
  tertiaryPeakTheta=0;
  tertiaryPeakPhi=0;
  findTertiaryPeakOfMap(mapCorVal, tertiaryPeakVal, tertiaryPeakTheta, tertiaryPeakPhi);
  
  //now calculate how far from TD it is using interpolated peak
  deltaPhi=peakPhiInterp;
  deltaTheta=peakThetaInterp;				
  deltamcmTheta=peakThetaInterp;
  deltamcmPhi=peakPhiInterp;
  double theta_deg=theta*rad2deg;//taylor dome theta
  double phi_deg=phi*rad2deg; //taylor dome phi
  double mcmtheta_deg=mcmtheta*rad2deg;
  double mcmphi_deg=mcmphi*rad2deg;
  //cout<<"theta_deg is "<<theta_deg<<" phi_deg is "<<phi_deg<<"\n";
  cout<<"Theta_interp is "<<peakThetaInterp<<" phi_interp is "<<peakPhiInterp<<"\n";
  deltaPhi-=phi_deg;
  deltaTheta-=theta_deg;
  deltamcmPhi-=mcmphi_deg;
  deltamcmTheta-=mcmtheta_deg;
  if (deltaPhi>180) deltaPhi-=360; 
  if (deltaPhi<-180) deltaPhi+=360;
  if (deltamcmPhi>180) deltamcmPhi-=360; 
  if (deltamcmPhi<-180) deltamcmPhi+=360;
  cout<<"deltaTheta is "<<deltaTheta<<" deltaPhi is "<<deltaPhi<<"\n";
  if (printFlag==1) cout<<"Interp Peak Theta - TD theta: "<<deltaTheta<<", peak phi - TD phi: "
      <<deltaPhi<<endl;
  if (printFlag==1) cout<<"Interp Peak Theta - McM theta: "<<deltamcmTheta<<", peak phi - McM phi: "
      <<deltamcmPhi<<endl;
  if (printFlag==1) if (fabs(deltaPhi)>5) cout<<"eventNumber: "<<eventNumber<<", deltaPhiTD: "<<deltaPhi<<", deltaThetaTD: "
			    <<deltaTheta<<", deltaPhiMcM: "<<deltamcmPhi<<", deltaThetaMcM: "<<deltamcmTheta<<endl;
  
  phiMap=peakPhiInterp;
  thetaMap=-1.*peakThetaInterp;
  
  //cout<<"thetaMap is "<<thetaMap<<"\n";
  peakPhiFinal=phiMap*deg2rad;
  peakThetaFinal=peakThetaInterp*deg2rad;
 
  if (printFlag==1){ 
    cout<<"WZ Phi and Theta final is "<<peakPhiFinal<<" "<<peakThetaFinal<<"\n";
    cout<<"thetaMap is "<<thetaMap*deg2rad<<"\n";
  }


  //get hw trigger angle
  hwTriggerAngle=getHWTriggerAngle(eventNumber);
  if (printFlag==1)    cout<<"hwTriggerAngle: "<<hwTriggerAngle<<endl;
  //see if triggerbits or phimaskbits are in direction of peak
  triggerOrPhiMaskFlag=isPeakPhiTriggeredOrMasked(fHeadPtr,peakPhiInterp);
  phiMaskFlag=isPeakPhiMasked(fHeadPtr,peakPhiInterp);
  hwTriggerFlag=isPeakTriggeredpm1(fHeadPtr,peakPhiInterp);
  triggerFlag=hwTriggerFlag;
  //get polarization
  polFractionCoherent=0;
  polAngleCoherent=getCoherentPolarization(grCoherent, grCoherentHoriz, polFractionCoherent);
  if (printFlag==1) cout<<"polanglecoherent: "<<fabs(polAngleCoherent)<<endl;
  
  if(drawMaps==1){
    double sourceLon,sourceLat,sourceAlt;
   int  eventTracer=traceBackToContinent_Brian(myEventNumber, peakThetaFinal,peakPhiFinal, sourceLon, sourceLat, sourceAlt);
   cout<<"distance from source is "<<distance_from_source<<"\n";
   cout<<"sourceLon, Lat, Alt are "<<sourceLon<<" "<<sourceLat<<" "<<sourceAlt<<"\n";
  }
 
  //delete graphs
  delete grVertCoherentBaseline;
  delete grHorizCoherentBaseline;


  grVertCoherentBaseline=0;
  grHorizCoherentBaseline=0;
 
  for (int ant1=0;ant1<NUM_ANTS_WITH_NADIRS;ant1++){
    for (int ant2=0;ant2<NUM_ANTS_WITH_NADIRS;ant2++){
      delete grCor[ant1][ant2];
      grCor[ant1][ant2]=0;
    }
    delete grEv[ant1];
    delete grEv_backup[ant1];
    delete grEvHoriz[ant1];
    delete grEvUnfiltered[ant1];
    delete grEvHorizUnfiltered[ant1];
    
    delete grVertBaseline[ant1];
    delete grHorizBaseline[ant1];
      
    grEvHoriz[ant1]=0;
    grEv[ant1]=0;
    grEv_backup[ant1]=0;
    grEvUnfiltered[ant1]=0;
    grEvHorizUnfiltered[ant1]=0;
    grVertBaseline[ant1]=0;
    grHorizBaseline[ant1]=0;

  }

  ratioOfPeaksPassFlag=0;
  elevationAnglePassFlag=0;
  xCorPassFlag=0;
  peakCrossCorrFlag=0;
  peakHilbertFlag=0;
  polFractionFlag=0;

  cout<<"ratioFirstToSecondPeak is "<<ratioFirstToSecondPeak<<"\n";
  //set pointing flags
  if (ratioFirstToSecondPeak>(1/(0.9))) ratioOfPeaksPassFlag=1;
  if (thetaMap>-35 && thetaMap<0) elevationAnglePassFlag=1;
  if (peakHilbertCoherent>=-350*peakVal+57.14) xCorPassFlag=1;
  if (peakVal>0.075) peakCrossCorrFlag=1;
  if (peakHilbertCoherent>15) peakHilbertFlag=1;
  if (polFractionCoherent>0.3) polFractionFlag=1;
  

  if (printFlag==1) cout<<" WX event number is: "<<eventNumber<<", -350*peakVal+57.14 is: "<<-350*peakVal+57.14<<", peakHilbert is: "<<peakHilbertCoherent<<", xCorPassFlag: "<<xCorPassFlag<<", ratioOfPeaks: "<<ratioOfPeaksPassFlag<<", elevation: "<<elevationAnglePassFlag<<", peakVal: "<<peakVal<<", peakHilbertCoherent: "<<peakHilbertCoherent<<", polFractionFlag: "<<polFractionFlag<<endl;
  
  if (printFlag==1) cout<<"didIFilter: "<<didIFilter<<", didIFilterHoriz: "<<didIFilterHoriz<<endl;

  //write the data
  if (writeData==1){
    if ((thermalFlagUniversal!=3) || (thermalFlagUniversal==3 && thetaMap>0)){
      ndata->Fill(deltaTheta,deltaPhi,deltamcmTheta, deltamcmPhi, thetaMap,phiMap,mapSNR,
		  peakVal,ratioFirstToSecondPeak,snrCoherent, snrPeakAnt, maxSignalPeak, distanceTD, peakHilbertCoherent);
      ndata2->Fill(deltaTTD, snrAfterFilter,didIFilter,triggerOrPhiMaskFlag,thetaTD,
		   phiTD,thetaWilly,phiWilly, hwTriggerAngle,thisPhiMask, distanceMcM, secondThetaMap, secondPhiMap, strongCWFlag);
      ndata3->Fill(headingOfThisEvent, nadirFlag,tertiaryPeakTheta*-1,tertiaryPeakPhi, isVarnerEvent(eventNumber),
		   isVarnerEvent2(eventNumber),
		   fAdu5APatPtr->pitch, fAdu5APatPtr->roll, fAdu5APatPtr->heading, phiMaskFlag, 
		   hwTriggerFlag, ratioOfPeaksPassFlag,
		   elevationAnglePassFlag,xCorPassFlag);
      ndata4->Fill(isPayloadBlast(eventNumber),polAngleCoherent, polFractionCoherent, didIFilterAboveSatellite,
		   didIFilterHoriz,didIFilterAboveSatelliteHoriz, meanFreqVert, meanFreqHoriz);
    }
  }
  
  if (peakHilbertCoherent>=-350*peakVal+57.14 && ratioFirstToSecondPeak>(1/(0.9))
      && thetaMap>-35 && thetaMap<0 && peakVal>0.075 && peakHilbertCoherent>15 && polFractionCoherent>0.3){
    if (printFlag==1) cout<<"WX WZ returning 1"<<endl;
    return 1; //a first cut for whether it pointed or not.
  }
  else{
    if (printFlag==1) cout<<"WX WZ returning 0"<<endl;
    return 0;
  }
  //return 1;
  
}
///////////////////////////
void MyCorrelator::drawRefinedMap(double mapCorValRefined[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI], 
				  double peakTheta, double peakPhi,int myEventNumber)
{
  // cout<<"Drawing Refined Map"<<endl;
  TH2F *hmap2;
  // TH2F *hmap2=new TH2F("hmap2","Refined MAP!",NUM_BINS_FINE_PHI,peakPhi-phiStartOffset,
  //		       peakPhi+phiStartOffset,
  //		       NUM_BINS_FINE_THETA,-1.*peakTheta-thetaStartOffset,
  //		       -1.*peakTheta+thetaStartOffset);
  if(cos_reconstruction_flag==1){
    hmap2=new TH2F("hmap2","Refined MAP!",NUM_BINS_FINE_PHI,peakPhi-phiStartOffset, peakPhi+phiStartOffset,NUM_BINS_FINE_THETA,-1*min_theta_cos,-1*max_theta_cos);
  }
  else{
    hmap2=new TH2F("hmap2","Refined MAP!",NUM_BINS_FINE_PHI,peakPhi-phiStartOffset,peakPhi+phiStartOffset,NUM_BINS_FINE_THETA,-1.*peakTheta-thetaStartOffset,-1.*peakTheta+thetaStartOffset);
  }
  
  //TH2F *haxes=new TH2F("haxes","haxes",100,peakPhi-5,peakPhi+5,100,peakTheta-5,peakTheta+5);
  
  for (double i=0;i<NUM_BINS_FINE_PHI;i++){
    for (double j=0;j<NUM_BINS_FINE_THETA;j++){
      double phi=(peakPhi-phiStartOffset)+(i+0.5)*FINE_BIN_SIZE;
     
      if(cos_reconstruction_flag==1){
	double abbyTheta_cos =(-1)* ((j+0.5)*step_size_cos_fine+max_theta_cos);
	//cout<<"abbyTheta_cos is "<<abbyTheta_cos<<"\n";
	hmap2->Fill(phi,abbyTheta_cos,mapCorValRefined[int(j)][int(i)]);//hack for now
      }
      else{
	double theta=-1.*((peakTheta-thetaStartOffset)+(j+0.5)*FINE_BIN_SIZE);
	//cout<<"theta is "<<theta<<"\n";
	hmap2->Fill(phi,theta,mapCorValRefined[int(j)][int(i)]);//hack for now
      }
    }
  }
  gStyle->SetPalette(1);
  TCanvas *c2=new TCanvas("c2","c2",800,400);
  c2->cd(0);
  gStyle->SetOptStat(kFALSE);
  //  cout<<"Min of Histogram in phi: "<<peakPhi-5<<", max: "<<peakPhi+5<<endl;
  // cout<<"Min of Histogram in theta: "<<peakTheta-5<<", max: "<<peakTheta+5<<endl;
  hmap2->Draw("colz");
  c2->Print("refinedMap.eps");
  //c2->Print("c2.png");
   char printer[256];
  
  sprintf(printer,"refinedMap_%i.png",myEventNumber);
  c2->Print(printer);
  //  haxes->Draw("");
}
///////////////////////////
void MyCorrelator::doRefinedMap(int myEventNumber, 
				double mapCorValRefined[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI], 
				double peakTheta, double peakPhi, int onlyTriggered, 
				int nadirFlag, int whichPolarization)
{
  
  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  myCally->setClockUpSampleFactor(4);
  Double_t abbyTheta=0;
  Double_t abbyPhi=0;
  int triggeredPhi[NUM_PHI_SECTORS];
  int triggeredAnt[NUM_ANTS_WITH_NADIRS];
  int nantennas;
  if (nadirFlag==0) nantennas=NUM_ANTS_NO_NADIRS;
  if (nadirFlag==1) nantennas=NUM_ANTS_WITH_NADIRS;

  double rms[nantennas];
  double mapCtr[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI];

  for (int i=0;i<NUM_BINS_FINE_THETA;i++){
    for (int j=0;j<NUM_BINS_FINE_PHI;j++){
      mapCorValRefined[i][j]=0;
      mapCtr[i][j]=0;
    }
  }
  if(cos_reconstruction_flag==1){
    cout<<"peakTheta is "<<peakTheta;
    peakTheta = acos(peakTheta);
    cout<<"acos is "<<peakTheta<<"\n";
    peakTheta = pi/2 - peakTheta;
    cout<<"peakTheta is now "<<peakTheta<<"\n";
    peakTheta = peakTheta*rad2deg;
   
    cout<<" now it is "<<peakTheta<<"\n";
  }
  //get the event and triggered antennas
  getTriggeredL2Phi(fHeadPtr,triggeredPhi);
  getTriggeredAnt(triggeredPhi,triggeredAnt);
  if (onlyTriggered==0){
    for (int ant=0;ant<nantennas;ant++){
      triggeredAnt[ant]=1;
    }
  }

  if (onlyTriggered==1){
    for (int ant=0;ant<nantennas;ant++){
      // if (phiMaskedAnts[ant] || triggeredAnt[ant]) triggeredAnt[ant]=1;
      if (triggeredAnt[ant]) triggeredAnt[ant]=1;
    }
  }
  for (int ant=0;ant<nantennas;ant++){
    if (whichPolarization==0) rms[ant]=getRMS(grEv[ant]); 
    else rms[ant]=getRMS(grEvHoriz[ant]); 
  }

  double lowerAngleThis, higherAngleThis, centerTheta1, centerTheta2, centerPhi1, centerPhi2;
  double tVal1, tVal2, corVal1, corVal2, corVal;			
  double tVal0,corVal0, tValend, corValend;
  int allowedFlag;
  double deltaPhi1_deg, deltaPhi2_deg;
  double abbyPhi_deg, abbyTheta_deg;
  double timeExpected;
  double lengthOfTrace;
  int bin1, bin2;
  int npoints;
  Double_t deltaTInt=1./(2.6*4);
  int loop_flag=0;
  //  cout<<"Beginning refined correlation map"<<endl;
  for (int ant1=0;ant1<nantennas;ant1++){
    for (int ant2=0;ant2<nantennas;ant2++){
      loop_flag=0;
       if(onlyTriggered==0){
	if (ant1<ant2 && ant1!=1 && ant2!=1 && saturatedChannels[ant1+whichPolarization*NUM_ANTS_WITH_NADIRS]==0 
	    && saturatedChannels[ant2+whichPolarization*NUM_ANTS_WITH_NADIRS]==0){//take out 2V
	  loop_flag=1;
       }
      }//onlyTrigger==0
      if(onlyTriggered==1){
	if (ant1<ant2 && ant1!=1 && ant2!=1 &&  saturatedChannels[ant1+whichPolarization*NUM_ANTS_WITH_NADIRS]==0 
	    && saturatedChannels[ant2+whichPolarization*NUM_ANTS_WITH_NADIRS]==0){//take out 2V
	  loop_flag=1;
	}
      }//onlyTrigger==1
      if(loop_flag==1){
	allowedFlag=allowedPhisPairOfAntennas(lowerAngleThis,higherAngleThis, 
					      centerTheta1, centerTheta2, centerPhi1, centerPhi2, 
					      ant1,ant2);
	
	if (higherAngleThis<lowerAngleThis) higherAngleThis+=360;
	if (allowedFlag==1 && ((peakPhi>lowerAngleThis && peakPhi<higherAngleThis) 
			       || (peakPhi+360>lowerAngleThis && peakPhi+360<higherAngleThis))){//only use antennas in this direction
	  
	  deltaPhi1_deg=peakPhi-centerPhi1;
	  deltaPhi2_deg=peakPhi-centerPhi2;  
	  if (deltaPhi1_deg>180 && deltaPhi1_deg<540) deltaPhi1_deg-=360;
	  if (deltaPhi2_deg>180 && deltaPhi2_deg<540) deltaPhi2_deg-=360;
	  if (deltaPhi1_deg>=540) deltaPhi1_deg-=720;
	  if (deltaPhi2_deg>=540) deltaPhi2_deg-=720;
	  if (deltaPhi1_deg<-180) deltaPhi1_deg+=360;
	  if (deltaPhi2_deg<-180) deltaPhi2_deg+=360;
  
	  if(cos_reconstruction_flag==1){
	    double cos_angle1;
	    double cos_angle2;
	    double angle1;
	    double angle2;
	   
	    cos_angle1 = cos(peakTheta*deg2rad)*cos(centerTheta1*deg2rad)+sin(peakTheta*deg2rad)*sin(centerTheta1*deg2rad)*cos(deltaPhi1_deg*deg2rad);
	    cos_angle2 = cos(peakTheta*deg2rad)*cos(centerTheta2*deg2rad)+sin(peakTheta*deg2rad)*sin(centerTheta2*deg2rad)*cos(deltaPhi2_deg*deg2rad);
	    
	    // cos_angle1 = cos(peakTheta*deg2rad)*cos(centerTheta1*deg2rad)+sin(peakTheta)*sin(centerTheta1*deg2rad)*cos(deltaPhi1_deg*deg2rad);
	    //cos_angle2 = cos(peakTheta)*cos(centerTheta2*deg2rad)+sin(peakTheta)*sin(centerTheta2*deg2rad)*cos(deltaPhi2_deg*deg2rad);
	    
	    angle1 = acos(cos_angle1);
	    angle2 = acos(cos_angle2);
	    
	    angle1=angle1*rad2deg;
	    angle2=angle2*rad2deg;
	    

	    double peakTheta_shift = 90 - peakTheta;//change from Abby coord to Brian. (0 is straight up)
	   
	    double max_theta = peakTheta_shift + thetaStartOffset; //Theta closest to up
	    double min_theta = peakTheta_shift - thetaStartOffset;// theta closest to down
	    
	    max_theta = max_theta*deg2rad;
	    min_theta = min_theta*deg2rad;
	    max_theta_cos = cos(max_theta);//smaller
	    min_theta_cos = cos(min_theta);//larger
	    
	    step_size_cos_fine = min_theta_cos - max_theta_cos; //positive number
	     step_size_cos_fine = step_size_cos_fine/NUM_BINS_FINE_THETA;
	    
	    
	    /* if ((fabs(deltaPhi1_deg*deltaPhi1_deg+
	       (peakTheta-centerTheta1)*(peakTheta-centerTheta1))<NUM_DEGREES_OFF_CENTER*NUM_DEGREES_OFF_CENTER) &&
	       (fabs(deltaPhi2_deg*deltaPhi2_deg+
	       (peakTheta-centerTheta2)*(peakTheta-centerTheta2))<NUM_DEGREES_OFF_CENTER*NUM_DEGREES_OFF_CENTER)){*/
	     if((angle1 <NUM_DEGREES_OFF_CENTER && angle2 <NUM_DEGREES_OFF_CENTER)||thermalSample==1 ){
	      
	      if (!grCor[ant1][ant2]){ 
		if (whichPolarization==0) grCor[ant1][ant2] = FFTtools::getInterpolatedCorrelationGraph(grEv[ant1],grEv[ant2],deltaTInt);
		else grCor[ant1][ant2] = FFTtools::getInterpolatedCorrelationGraph(grEvHoriz[ant1],grEvHoriz[ant2],deltaTInt);
	      }
	      grCor[ant1][ant2]->GetPoint(0,tVal0,corVal0);
	      grCor[ant1][ant2]->GetPoint(grCor[ant1][ant2]->GetN()-1,tValend,corValend);
	      lengthOfTrace=tValend-tVal0;
	      npoints=grCor[ant1][ant2]->GetN();
	      
	     
	      double abbyTheta_cos;
	      for (int ctrPhi=0;ctrPhi<NUM_BINS_FINE_PHI;ctrPhi++){
		abbyPhi_deg=ctrPhi*FINE_BIN_SIZE+peakPhi-phiStartOffset;
		abbyPhi=abbyPhi_deg*deg2rad;
		
		//do the scan of the sky once we've picked two antennas we are interested in
		for (int ctrTheta=0;ctrTheta<NUM_BINS_FINE_THETA;ctrTheta++){
		  abbyTheta_cos = ctrTheta*step_size_cos_fine+max_theta_cos;
		  //cout<<"abbyTheta_cos is "<<abbyTheta_cos<<"\n";
		 
		  //if(ctrPhi==0) cout<<"ctrTheta is "<<ctrTheta<<" abbyTheta_Cos is "<<abbyTheta_cos<<"\n";
		  //	abbyTheta_deg=ctrTheta*FINE_BIN_SIZE+peakTheta-thetaStartOffset;	      
		  //abbyTheta=abbyTheta_deg*deg2rad;
		  
		  abbyTheta = acos(abbyTheta_cos);
		  abbyTheta =  pi/2 - abbyTheta;
		  abbyTheta_deg = abbyTheta*rad2deg;
		  
		  //now see if the peak bin is within the beamwidth (generous) of both antennas in 3D angles
		  
		  timeExpected=getDeltaTExpected(ant1, ant2, abbyPhi_deg, abbyTheta_deg);
		  
		  //get the value of the cross correlation graph at this delay.
		  bin1=int((timeExpected-tVal0)/lengthOfTrace*npoints);
		  bin2=bin1+1;
		  grCor[ant1][ant2]->GetPoint(bin1,tVal1,corVal1);
		  grCor[ant1][ant2]->GetPoint(bin2,tVal2,corVal2);
		  //cout<<"grCor is "<<grCor[ant1][ant2]<<"\n";
		  corVal=(corVal2-corVal1)*(timeExpected-tVal1)/(tVal2-tVal1)+corVal1;//interpolate between 2 bins	   
		  mapCorValRefined[ctrTheta][ctrPhi]+=corVal/(rms[ant1]*rms[ant2]);
		  mapCtr[ctrTheta][ctrPhi]+=1;
		}
	      }//end phi loop
	    }//end theta loop
	  }//con_reconstruction
	  else{

	    if ((fabs(deltaPhi1_deg*deltaPhi1_deg+
		      (peakTheta-centerTheta1)*(peakTheta-centerTheta1))<NUM_DEGREES_OFF_CENTER*NUM_DEGREES_OFF_CENTER) &&
		(fabs(deltaPhi2_deg*deltaPhi2_deg+
		      (peakTheta-centerTheta2)*(peakTheta-centerTheta2))<NUM_DEGREES_OFF_CENTER*NUM_DEGREES_OFF_CENTER)){
	      
	      if (!grCor[ant1][ant2]){ 
		if (whichPolarization==0) grCor[ant1][ant2] = FFTtools::getInterpolatedCorrelationGraph(grEv[ant1],grEv[ant2],deltaTInt);
		else grCor[ant1][ant2] = FFTtools::getInterpolatedCorrelationGraph(grEvHoriz[ant1],grEvHoriz[ant2],deltaTInt);
	      }
	      grCor[ant1][ant2]->GetPoint(0,tVal0,corVal0);
	      grCor[ant1][ant2]->GetPoint(grCor[ant1][ant2]->GetN()-1,tValend,corValend);
	      lengthOfTrace=tValend-tVal0;
	      npoints=grCor[ant1][ant2]->GetN();
	      
	      for (int ctrPhi=0;ctrPhi<NUM_BINS_FINE_PHI;ctrPhi++){
		abbyPhi_deg=ctrPhi*FINE_BIN_SIZE+peakPhi-phiStartOffset;
		abbyPhi=abbyPhi_deg*deg2rad;
		
		//do the scan of the sky once we've picked two antennas we are interested in
		for (int ctrTheta=0;ctrTheta<NUM_BINS_FINE_THETA;ctrTheta++){
		  abbyTheta_deg=ctrTheta*FINE_BIN_SIZE+peakTheta-thetaStartOffset;	      
		  abbyTheta=abbyTheta_deg*deg2rad;
		  
		  //now see if the peak bin is within the beamwidth (generous) of both antennas in 3D angles
		  
		  timeExpected=getDeltaTExpected(ant1, ant2, abbyPhi_deg, abbyTheta_deg);
		  
		  //get the value of the cross correlation graph at this delay.
		  bin1=int((timeExpected-tVal0)/lengthOfTrace*npoints);
		  bin2=bin1+1;
		  grCor[ant1][ant2]->GetPoint(bin1,tVal1,corVal1);
		  grCor[ant1][ant2]->GetPoint(bin2,tVal2,corVal2);
		  corVal=(corVal2-corVal1)*(timeExpected-tVal1)/(tVal2-tVal1)+corVal1;//interpolate between 2 bins	   
		  mapCorValRefined[ctrTheta][ctrPhi]+=corVal/(rms[ant1]*rms[ant2]);
		  mapCtr[ctrTheta][ctrPhi]+=1;
		}
	      }//end phi loop
	    }//end theta loop
	  }//cos_reconstruction
	}
      }//end trig==1 cut loop
    }//end ant2loop
  }//end ant1loop
  
  double peakRefinedVal;
  int peakThetaRefinedBin,peakPhiRefinedBin;
  findPeakOfRefinedMap(mapCorValRefined,peakRefinedVal,peakThetaRefinedBin,peakPhiRefinedBin);
  
  for (int i=0;i<NUM_BINS_FINE_THETA;i++){
    for (int j=0;j<NUM_BINS_FINE_PHI;j++){
      
      if (mapCtr[i][j]!=0) mapCorValRefined[i][j]=mapCorValRefined[i][j]/mapCtr[i][j];

      if(mapCorValRefined[i][j] != mapCorValRefined[i][j]) mapCorValRefined[i][j]=0.;
    }
  }
 
}
//////////////////////////////////
void MyCorrelator::drawCorrelationMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], int myEventNumber,int trial)
{
  //  cout<<"Drawing Map"<<endl;
  TH2F *hmap1;
  // TH2F *hmap1=new TH2F("hmap1","MAP!",NUM_BINS_ROUGH_PHI,0,360,NUM_BINS_ROUGH_THETA,-90,90);
  if(cos_reconstruction_flag==1){
    hmap1=new TH2F("hmap1","MAP!",NUM_BINS_ROUGH_PHI,0,360,NUM_BINS_PLOT_COS,-1,1);
  }
  else{
    hmap1=new TH2F("hmap1","MAP!",NUM_BINS_ROUGH_PHI,0,360,NUM_BINS_ROUGH_THETA,-90,90);
  }
  //cout<<"step_size_cos_rough is "<<step_size_cos_rough<<" cosHighestTheta is "<<cosHighestTheta<<"\n";
 
  
  for (int i=0;i<NUM_BINS_ROUGH_PHI;i++){
    for (int j=0;j<NUM_BINS_ROUGH_THETA;j++){
      double phi=ROUGH_BIN_SIZE*(i+0.5);
      //double theta=(ROUGH_BIN_SIZE*-1.*(j+0.5))+HIGHEST_THETA;
      if(cos_reconstruction_flag==1){
	double costheta = (-1)*(step_size_cos_rough*(j+0.5)+cosHighestTheta);//cosHighestTheta is negative (costheta =1 is straight down to match abby's)
	
	// if (theta<=HIGHEST_THETA)
	if (costheta < (-1)*cosHighestTheta) hmap1->Fill(phi,costheta,mapCorVal[j][i]);

      }//cos
      else{
	double theta=(ROUGH_BIN_SIZE*-1.*(j+0.5))+HIGHEST_THETA;
	if (theta<=HIGHEST_THETA) hmap1->Fill(phi,theta,mapCorVal[j][i]);
      }
    }//j
  }//i
  gStyle->SetPalette(1);
  TCanvas *c1=new TCanvas("c1","c1",400,400);
  c1->cd(0);
  gStyle->SetOptStat(kFALSE);
  // hmap1->SetMaximum(.15);
  hmap1->Draw("colz");
  hmap1->GetXaxis()->SetTitle("Payload Azimuth (Degrees)");
  hmap1->GetYaxis()->SetTitle("Elevation (cos(#theta))");
  c1->Print("mainMap.eps");
  char printer[256];
  
  sprintf(printer,"mainMap_%i.png",myEventNumber);
  c1->Print(printer);
   // c1->Print("c1.png");
}
////////////////////////////////////////////////////
void MyCorrelator::doCorrelationMap(int myEventNumber, 
				    double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], 
				    int onlyTriggered, int nadirFlag, int whichPolarization)
{
  
  //struct timeb timebuffer;
  //char *timeline;
  
  //ftime( &timebuffer );
  //timeline = ctime( & ( timebuffer.time ) );
  // printf( "The time is %.19s.%hu %s", timeline, timebuffer.millitm, &timeline[20] );

  if (eventStartedFlag!=myEventNumber) eventStartedFlag=startEachEvent(myEventNumber);
  if (eventEntryGottenFlag!=(int)myEventNumber) eventEntryGottenFlag=getEventEntry();
  myCally->setClockUpSampleFactor(4);
  Double_t abbyTheta=0;
  Double_t abbyPhi=0;
  int triggeredPhi[NUM_PHI_SECTORS];
  int triggeredAnt[NUM_ANTS_WITH_NADIRS];
  int phiMaskArray[NUM_PHI_SECTORS];
  int phiMaskedAnts[NUM_ANTS_WITH_NADIRS]; 
  int bestPhiSector;
  int phiOKToUse[360];
 
  int nantennas=0;
  if (nadirFlag==0) nantennas=NUM_ANTS_NO_NADIRS;
  if (nadirFlag==1) nantennas=NUM_ANTS_WITH_NADIRS;

  double rms[nantennas];
  double mapCtr[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI];

  for (int i=0;i<NUM_BINS_ROUGH_THETA;i++){
    for (int j=0;j<NUM_BINS_ROUGH_PHI;j++){
      mapCorVal[i][j]=0;
      mapCtr[i][j]=0;
    }
  }
 
  // if ((fHeadPtr->trigType(1<<0))){//means to look for RF triggers, 1<<1 would be PPSs (I think)

  //get the event and triggered antennas
 
  if (onlyTriggered==1){
    getTriggeredL2Phi(fHeadPtr,triggeredPhi);
    getTriggeredAntpm1PhiSector(triggeredPhi,triggeredAnt);
    int thisPhiMask=isPhiMaskingOn(fHeadPtr->eventNumber,phiMaskArray); 
    thisPhiMask=0;
    getPhiMaskedAntspm1(phiMaskArray,phiMaskedAnts);
  }
  
  if (onlyTriggered==0){
    for (int ant=0;ant<nantennas;ant++) triggeredAnt[ant]=1;   
    for (int i=0;i<360;i++) phiOKToUse[i]=1;
  }
  int nextPhiSector;
  int lastPhiSector;
  if (onlyTriggered==1){
    for (int ant=0;ant<nantennas;ant++){
      // if (phiMaskedAnts[ant] || triggeredAnt[ant]) triggeredAnt[ant]=1;
      if (triggeredAnt[ant]){
	triggeredAnt[ant]=1;
	cout<<"triggered ant is "<<ant<<"\n";
      }
    }
    /* for (int i=0;i<360;i++){
      bestPhiSector=getBestPhiSector(i*deg2rad);
      // cout<<"bestPhiSector is "<<bestPhiSector<<"\n";
      // if (triggeredPhi[bestPhiSector]|| phiMaskArray[bestPhiSector]) phiOKToUse[i]=1;
      nextPhiSector = bestPhiSector++;
      lastPhiSector = bestPhiSector--;
      if(nextPhiSector > 15) nextPhiSector-=16;
      if(lastPhiSector < 0) lastPhiSector+=16;
      
      if (triggeredPhi[bestPhiSector] || triggeredPhi[lastPhiSector] || triggeredPhi[nextPhiSector]){
       	phiOKToUse[i]=1;
	//cout<<"using phi sector "<<bestPhiSector<<"\n";
	
	//i++;
      }
      else phiOKToUse[i]=0;
     
    }
    */
   
    
    for(int phisector=0;phisector<16;phisector++){
      nextPhiSector = phisector+1;
      lastPhiSector = phisector-1;
     
      if(nextPhiSector > 15) nextPhiSector-=16;
      if(lastPhiSector < 0) lastPhiSector+=16;

      if (triggeredPhi[phisector]){// || triggeredPhi[lastPhiSector] || triggeredPhi[nextPhiSector]){
	cout<<"triggered phi sectors are "<<phisector<<"\n";
	 for (int i=0;i<360;i++){
	   bestPhiSector=getBestPhiSector(i*deg2rad);
	   if(bestPhiSector==phisector){
	     phiOKToUse[i]=1;
	   }//bestPhi==phi
	 }//i
	

      }//triggered
    }//phisector

    
  }
 
  for (int ant=0;ant<nantennas;ant++){
    if (whichPolarization==0) rms[ant]=getRMS(grEv[ant]); 
   else  rms[ant]=getRMS(grEvHoriz[ant]); 
  }

  int abbyPhi_deg_this;
  int deltaPhi1_deg, deltaPhi2_deg;
  int centerPhi1Int, centerPhi2Int;
  int centerTheta1Int, centerTheta2Int;
  double tVal1, tVal2, corVal1, corVal2, corVal;
  double tVal0,corVal0, tValend, corValend;
  double timeExpected;
  double lengthOfTrace;
  int bin1, bin2;
  double lowerAngleThis, higherAngleThis, centerTheta1, centerTheta2, centerPhi1, centerPhi2;
  int allowedFlag;
  int npoints;
  Double_t deltaTInt=1./(2.6*4);
  //  TGraph *grCor[ant1][ant2];

  //double max_corval;
  //double ant1_max;
  //double ant2_max;
  // cout<<"Beginning correlation map"<<endl;
  int loop_flag=0;
  for (int ant1=0;ant1<nantennas;ant1++){
    // cout<<"phi sector is for ant "<<ant1<<" is "<<fUPGeomTool->getPhiFromAnt(ant1)<<"\n";
    //cout<<" phi is "<<fUPGeomTool->getPhiFromAnt(ant1)*PHI_SECTOR_ANGLE-ADU5_FORE_PHI<<"\n";
    
    for (int ant2=0;ant2<nantennas;ant2++){
      loop_flag=0;
      if(onlyTriggered==0){
	if (ant1<ant2 && ant1!=1 && ant2!=1 && saturatedChannels[ant1+whichPolarization*NUM_ANTS_WITH_NADIRS]==0 
	    && saturatedChannels[ant2+whichPolarization*NUM_ANTS_WITH_NADIRS]==0){//take out 2V
	  loop_flag=1;
	}
      }//onlyTrigger==0
      if(onlyTriggered==1){
	//	if ((triggeredAnt[ant1]==1 || triggeredAnt[ant2]==1) && ant1<ant2 && ant1!=1 && ant2!=1 &&  saturatedChannels[ant1+whichPolarization*NUM_ANTS_WITH_NADIRS]==0 
	//    && saturatedChannels[ant2+whichPolarization*NUM_ANTS_WITH_NADIRS]==0){//take out 2V
	  if (ant1<ant2 && ant1!=1 && ant2!=1 &&  saturatedChannels[ant1+whichPolarization*NUM_ANTS_WITH_NADIRS]==0 
	      && saturatedChannels[ant2+whichPolarization*NUM_ANTS_WITH_NADIRS]==0){//take out 2V
	  loop_flag=1;
	}
      }//onlyTrigger==1
      if(loop_flag==1){
	allowedFlag=allowedPhisPairOfAntennas(lowerAngleThis,higherAngleThis,
					      centerTheta1, centerTheta2, centerPhi1, 
					      centerPhi2, ant1,ant2);//checks ants are within 2 phi sectors
	centerPhi2Int=int(centerPhi2);//degrees
	centerPhi1Int=int(centerPhi1);//degrees
	centerTheta1Int=int(centerTheta1);
	centerTheta2Int=int(centerTheta2);

	if (higherAngleThis<lowerAngleThis)higherAngleThis+=360;
	if (allowedFlag==1){
	 
	  //now get cross correlation graph
	  //cout<<"grEv has "<<grEv[ant1]->GetN()<<" number of points \n";
	 
	  if (!grCor[ant1][ant2]){
	    if (whichPolarization==0) grCor[ant1][ant2] = FFTtools::getInterpolatedCorrelationGraph(grEv[ant1],grEv[ant2],deltaTInt);
	    else grCor[ant1][ant2] = FFTtools::getInterpolatedCorrelationGraph(grEvHoriz[ant1],grEvHoriz[ant2],deltaTInt);
	  }
	  grCor[ant1][ant2]->GetPoint(0,tVal0,corVal0);//get first point
	  grCor[ant1][ant2]->GetPoint(grCor[ant1][ant2]->GetN()-1,tValend,corValend);//get last point
	  lengthOfTrace=tValend-tVal0;
	  npoints=grCor[ant1][ant2]->GetN(); 


	  double abbyTheta_deg;
	  double abbyTheta_cos;
	  int inside_cone_flag;
	  double angle1;
	  double cos_angle1;
	  double angle2;
	  double cos_angle2;

	  //do the scan of the sky once we've picked two antennas we are interested in
	  for (int abbyPhi_deg=int(lowerAngleThis);abbyPhi_deg<higherAngleThis;abbyPhi_deg+=ROUGH_BIN_SIZE){//go through all phi angles between the two antennas
	   
	    abbyPhi_deg_this=abbyPhi_deg;
	    if (abbyPhi_deg_this>=360) abbyPhi_deg_this-=360;
	    abbyPhi=abbyPhi_deg_this*deg2rad;
	    deltaPhi1_deg=abbyPhi_deg-centerPhi1Int;
	    deltaPhi2_deg=abbyPhi_deg-centerPhi2Int;
	    if (deltaPhi1_deg>180 && deltaPhi1_deg<540) deltaPhi1_deg-=360;
	    if (deltaPhi2_deg>180 && deltaPhi2_deg<540) deltaPhi2_deg-=360;
	    if (deltaPhi1_deg>=540) deltaPhi1_deg-=720;
	    if (deltaPhi2_deg>=540) deltaPhi2_deg-=720;
	    if (deltaPhi1_deg<-180) deltaPhi1_deg+=360;
	    if (deltaPhi2_deg<-180) deltaPhi2_deg+=360;
	    
	    if(cos_reconstruction_flag==1){
	      for (int ctrTheta=0;ctrTheta<=NUM_BINS_ROUGH_COS;ctrTheta++){
		abbyTheta_cos = ctrTheta*step_size_cos_rough + cosHighestTheta;
		
		//cout<<"ctr theta is "<<ctrTheta<<" and abbyTheta_cos is "<<abbyTheta_cos;
		// abbyTheta=abbyTheta_deg*deg2rad;
		abbyTheta = acos(abbyTheta_cos);
		abbyTheta = pi/2 - abbyTheta;
		abbyTheta_deg = abbyTheta*rad2deg;
		
		
		cos_angle1 = cos(abbyTheta)*cos(centerTheta1Int*deg2rad)+sin(abbyTheta)*sin(centerTheta1Int*deg2rad)*cos(deltaPhi1_deg*deg2rad);
		cos_angle2 = cos(abbyTheta)*cos(centerTheta2Int*deg2rad)+sin(abbyTheta)*sin(centerTheta2Int*deg2rad)*cos(deltaPhi2_deg*deg2rad);
		
		angle1 = acos(cos_angle1);
		angle2 = acos(cos_angle2);
		
		angle1=angle1*rad2deg;
		angle2=angle2*rad2deg;
		
		if(angle1 < NUM_DEGREES_OFF_CENTER && angle2 <NUM_DEGREES_OFF_CENTER && phiOKToUse[int(abbyPhi_deg_this)]==1){
		  
		  timeExpected=getDeltaTExpected(ant1, ant2, abbyPhi_deg, abbyTheta_deg);
		 
		  //double timeExpected=0;
		  //get the value of the cross correlation graph at this delay.
		  
		  bin1=int((timeExpected-tVal0)/lengthOfTrace*npoints);
		  bin2=bin1+1;
		  grCor[ant1][ant2]->GetPoint(bin1,tVal1,corVal1);
		  grCor[ant1][ant2]->GetPoint(bin2,tVal2,corVal2);
		  corVal=(corVal2-corVal1)*(timeExpected-tVal1)/(tVal2-tVal1)+corVal1;//interpolate between 2 bins
		  
		  mapCorVal[ctrTheta][int((abbyPhi_deg_this)/ROUGH_BIN_SIZE)]+=corVal/(rms[ant1]*rms[ant2]);//add this correlation to the bin
		  mapCtr[ctrTheta][int((abbyPhi_deg_this)/ROUGH_BIN_SIZE)]+=1;//keep track of the number of points in that bin
		}
		
	      }//end theta loop
	    }//cos_reconstruction
	    else{
	      for (int abbyTheta_deg=-1*HIGHEST_THETA;abbyTheta_deg<LOWEST_THETA; abbyTheta_deg+=ROUGH_BIN_SIZE){ //Highest theta==25, lowest theta=60
	
		abbyTheta=abbyTheta_deg*deg2rad;
		//now see if the point is within the beamwidth (generous) of both antennas in 3D angles
		if ((fabs(deltaPhi1_deg*deltaPhi1_deg+
			  (abbyTheta_deg-centerTheta1Int)*
			  (abbyTheta_deg-centerTheta1Int))<NUM_DEGREES_OFF_CENTER*NUM_DEGREES_OFF_CENTER) &&
		    (fabs(deltaPhi2_deg*deltaPhi2_deg+
			  (abbyTheta_deg-centerTheta2Int)*
			  (abbyTheta_deg-centerTheta2Int))<NUM_DEGREES_OFF_CENTER*NUM_DEGREES_OFF_CENTER)
		    && phiOKToUse[int(abbyPhi_deg_this)]==1){
		  
		  timeExpected=getDeltaTExpected(ant1, ant2, abbyPhi_deg, abbyTheta_deg);
		  
		  if(ant1==0 && (ant2 ==1 || ant2 == 7 || ant2== 9 || ant2 == 15)){
		    // cout<<"timeExpected for these  antennas 0,"<<ant2<<" at phi = "<<abbyPhi_deg<<" theta = "<<abbyTheta_deg<<" is "<<timeExpected<<"\n";
		    
		  }
		  
		  //double timeExpected=0;
		  //get the value of the cross correlation graph at this delay.
		  
		  
		  bin1=int((timeExpected-tVal0)/lengthOfTrace*npoints);
		  bin2=bin1+1;
		  grCor[ant1][ant2]->GetPoint(bin1,tVal1,corVal1);
		  grCor[ant1][ant2]->GetPoint(bin2,tVal2,corVal2);
		  corVal=(corVal2-corVal1)*(timeExpected-tVal1)/(tVal2-tVal1)+corVal1;//interpolate between 2 bins
		  
		  mapCorVal[int((abbyTheta_deg+HIGHEST_THETA)/ROUGH_BIN_SIZE)][int((abbyPhi_deg_this)/ROUGH_BIN_SIZE)]
		    +=corVal/(rms[ant1]*rms[ant2]);//add this correlation to the bin
		  mapCtr[int((abbyTheta_deg+HIGHEST_THETA)/ROUGH_BIN_SIZE)][int((abbyPhi_deg_this)/ROUGH_BIN_SIZE)]+=1;//keep track of the number of points in that bin
		}
		
	      }//end theta loop
	    }//cos_reconstruction
	  }//end phi loop
	  
	}//allowed flag
      }//end trig==1 cut loop
    }//end ant2loop
  }//end ant1loop


  for (int i=0;i<NUM_BINS_ROUGH_THETA;i++){
    for (int j=0;j<NUM_BINS_ROUGH_PHI;j++){
      
      if (mapCtr[i][j]>2) mapCorVal[i][j]=mapCorVal[i][j]/mapCtr[i][j]; //require 3 antennas to contribute to location
      else mapCorVal[i][j]=0;
      
      if(mapCorVal[i][j] != mapCorVal[i][j]){
	mapCorVal[i][j]=0;
      }
      
    }
  }

  

  //  ftime( &timebuffer );
  //timeline = ctime( & ( timebuffer.time ) );
  
  //printf( "The time is %.19s.%hu %s", timeline, timebuffer.millitm, &timeline[20] );

  /*if (drawMaps==1){
    TCanvas *cantennaMap=new TCanvas("cantennaMap","cantennaMap",800,400);
    TH2F *hantennaMap=new TH2F("hantennaMap","hantennaMap",NUM_BINS_ROUGH_PHI,0,360,NUM_BINS_ROUGH_THETA,-90,90);
    
    for (int i=0;i<NUM_BINS_ROUGH_PHI;i++){
      for (int j=0;j<NUM_BINS_ROUGH_THETA;j++){
	double phi=ROUGH_BIN_SIZE*(i+0.5);
	double theta=(ROUGH_BIN_SIZE*-1.*(j+0.5))+HIGHEST_THETA;
	if (theta<=HIGHEST_THETA)
	  if (mapCtr[j][i]>2) hantennaMap->Fill(phi,theta,mapCtr[j][i]);
      }
    }
    gStyle->SetPalette(1);
    cantennaMap->cd(0);
    gStyle->SetOptStat(kFALSE);
    hantennaMap->Draw("colz");
    cantennaMap->Print("antennaNumberMap.eps");
    }*/
  
  
  //}//end if RF trigger
}
///////////////////////////
void MyCorrelator::findPeakOfMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], double &peakVal, 
				 double &peakTheta, double &peakPhi)
{
  peakVal=0;
  for (int i=0;i<NUM_BINS_ROUGH_THETA;i++){
    for (int j=0;j<NUM_BINS_ROUGH_PHI;j++){
      
      if (mapCorVal[i][j]>peakVal){
	peakVal=mapCorVal[i][j];

	peakPhi=ROUGH_BIN_SIZE*j;
	if(cos_reconstruction_flag==1){
	  peakTheta=i*step_size_cos_rough+cosHighestTheta;//ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	}
	else{
	  peakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	}
	/*peakTheta= acos(peakTheta);
	peakTheta = pi/2 - peakTheta ; //acos -> (0,180), but zero for payload is horizontal, so need to shift angle by 90 deg. Negative is toward sky.
	peakTheta = peakTheta * rad2deg;
	*/
	
      }      
    }
  }
}
/////////////////////////////////
void MyCorrelator::findSecondaryPeakOfMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], 
					  double &secondaryPeakVal, 
					  double &secondaryPeakTheta, double &secondaryPeakPhi)
{
  double peakVal;			       
  secondaryPeakVal=0;
  int peakThetaBin;
  int peakPhiBin;
  int numBinsToCut=5;//10 degrees on each side is out

  findPeakOfMapBin(mapCorVal, peakVal, peakThetaBin, peakPhiBin);
  // cout<<"peakphibin: "<<peakPhiBin<<", peakThetaBin: "<<peakThetaBin<<", peakVal: "<<peakVal<<endl;
  int lowThetaBin=peakThetaBin-numBinsToCut;
  int highThetaBin=peakThetaBin+numBinsToCut;
  int lowPhiBin=peakPhiBin-numBinsToCut;
  int highPhiBin=peakPhiBin+numBinsToCut;
  if (lowPhiBin<0) lowPhiBin+=NUM_BINS_ROUGH_PHI;
  if (highPhiBin>NUM_BINS_ROUGH_PHI) highPhiBin-=NUM_BINS_ROUGH_PHI;
  
  for (int i=1;i<NUM_BINS_ROUGH_THETA-1;i++){
    for (int j=1;j<NUM_BINS_ROUGH_PHI-1;j++){
      if (lowPhiBin<highPhiBin){
	if ((i>=lowThetaBin && i<=highThetaBin) && (j>=lowPhiBin && j<=highPhiBin)){}
	else{
	  if (mapCorVal[i][j]>secondaryPeakVal){
	    secondaryPeakVal=mapCorVal[i][j];
	    //secondaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    secondaryPeakPhi=ROUGH_BIN_SIZE*j;
	    
	    if(cos_reconstruction_flag==1){
	      secondaryPeakTheta =i*step_size_cos_rough+cosHighestTheta;
	    }
	    else{
	      secondaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    }
	  }
	}
      }
      if (lowPhiBin>highPhiBin){
	if ((i>=lowThetaBin && i<=highThetaBin) && (j<=highPhiBin || j>=lowPhiBin)){}
	else{
	  //	if ((i<lowThetaBin || i>highThetaBin) && (j<lowPhiBin && j>highPhiBin)){
	  if (mapCorVal[i][j]>secondaryPeakVal){
	    secondaryPeakVal=mapCorVal[i][j];
	    secondaryPeakPhi=ROUGH_BIN_SIZE*j;
	    
	    if(cos_reconstruction_flag==1){
	      secondaryPeakTheta =i*step_size_cos_rough+cosHighestTheta;
	    }
	    else{
	      secondaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    }
	    
	  }
	}
      }
    }
  }
  if (printFlag==1) cout<<"Secondary peak: "<<secondaryPeakVal<<", theta: "<<secondaryPeakTheta
			<<", phi: "<<secondaryPeakPhi<<endl;
  
}
///////////////////////////////////////////
void MyCorrelator::findTertiaryPeakOfMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], 
					  double &tertiaryPeakVal, 
					  double &tertiaryPeakTheta, double &tertiaryPeakPhi)
{
  double peakVal;			       
  tertiaryPeakVal=0;
  int peakThetaBin;
  int peakPhiBin;
  int numBinsToCut=5;//10 degrees on each side is out
  double secondaryPeakTheta,secondaryPeakPhi,secondaryPeakVal;
  
  findPeakOfMapBin(mapCorVal, peakVal, peakThetaBin, peakPhiBin);
  int lowThetaBin=peakThetaBin-numBinsToCut;
  int highThetaBin=peakThetaBin+numBinsToCut;
  int lowPhiBin=peakPhiBin-numBinsToCut;
  int highPhiBin=peakPhiBin+numBinsToCut;
  if (lowPhiBin<0) lowPhiBin+=NUM_BINS_ROUGH_PHI;
  if (highPhiBin>NUM_BINS_ROUGH_PHI) highPhiBin-=NUM_BINS_ROUGH_PHI;
  findSecondaryPeakOfMap(mapCorVal,secondaryPeakVal,secondaryPeakTheta,secondaryPeakPhi);
  int secondaryThetaBin=int((secondaryPeakTheta+HIGHEST_THETA)/ROUGH_BIN_SIZE);
  int secondaryPhiBin=int((secondaryPeakPhi)/ROUGH_BIN_SIZE);
  int lowThetaBinSecond=secondaryThetaBin-numBinsToCut;
  int highThetaBinSecond=secondaryThetaBin+numBinsToCut;
  int lowPhiBinSecond=secondaryPhiBin-numBinsToCut;
  int highPhiBinSecond=secondaryPhiBin+numBinsToCut;
  if (lowPhiBinSecond<0) lowPhiBinSecond+=NUM_BINS_ROUGH_PHI;
  if (highPhiBinSecond>NUM_BINS_ROUGH_PHI) highPhiBinSecond-=NUM_BINS_ROUGH_PHI;
 
  for (int i=1;i<NUM_BINS_ROUGH_THETA-1;i++){
    for (int j=1;j<NUM_BINS_ROUGH_PHI-1;j++){
      if (lowPhiBin<highPhiBin && lowPhiBinSecond<highPhiBinSecond){
	if (((i>=lowThetaBin && i<=highThetaBin) && (j>=lowPhiBin && j<=highPhiBin)) || 
	    ((i>=lowThetaBinSecond && i<=highThetaBinSecond) && (j>=lowPhiBinSecond && j<=highPhiBinSecond))){}
	else{
	  if (mapCorVal[i][j]>tertiaryPeakVal){
	    tertiaryPeakVal=mapCorVal[i][j];
	    
	    tertiaryPeakPhi=ROUGH_BIN_SIZE*j;


	    if(cos_reconstruction_flag==1){
	      tertiaryPeakTheta=i*step_size_cos_rough+cosHighestTheta;//ROUGH_BIN_SIZE*i-HIGHEST_THETA; 
	    }
	    else{
	      tertiaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    }
	   
	  }
	}
      }
      if (lowPhiBin>highPhiBin && lowPhiBinSecond>highPhiBinSecond){
	if (((i>=lowThetaBin && i<=highThetaBin) && (j<=highPhiBin || j>=lowPhiBin)) ||
	    ((i>=lowThetaBinSecond && i<=highThetaBinSecond) && (j<=highPhiBinSecond || j>=lowPhiBinSecond))){}
	else{
	  //	if ((i<lowThetaBin || i>highThetaBin) && (j<lowPhiBin && j>highPhiBin)){
	  if (mapCorVal[i][j]>tertiaryPeakVal){
	    tertiaryPeakVal=mapCorVal[i][j];
	    //tertiaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    tertiaryPeakPhi=ROUGH_BIN_SIZE*j;

	    if(cos_reconstruction_flag==1){
	      tertiaryPeakTheta=i*step_size_cos_rough+cosHighestTheta;//ROUGH_BIN_SIZE*i-HIGHEST_THETA; 
	    }
	    else{
	      tertiaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    }

	  }
	}
      }
      if (lowPhiBin<highPhiBin && lowPhiBinSecond>highPhiBinSecond){
	if (((i>=lowThetaBin && i<=highThetaBin) && (j>=lowPhiBin && j<=highPhiBin)) || 
	    ((i>=lowThetaBinSecond && i<=highThetaBinSecond) && (j>=lowPhiBinSecond || j<=highPhiBinSecond))){}
	else{
	  if (mapCorVal[i][j]>tertiaryPeakVal){
	    tertiaryPeakVal=mapCorVal[i][j];
	    //tertiaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    tertiaryPeakPhi=ROUGH_BIN_SIZE*j;

	    if(cos_reconstruction_flag==1){
	      tertiaryPeakTheta=i*step_size_cos_rough+cosHighestTheta;//ROUGH_BIN_SIZE*i-HIGHEST_THETA; 
	    }
	    else{
	      tertiaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    }

	  }
	}
      }
      if (lowPhiBin>highPhiBin && lowPhiBinSecond<highPhiBinSecond){
	if (((i>=lowThetaBin && i<=highThetaBin) && (j<=highPhiBin || j>=lowPhiBin)) ||
	    ((i>=lowThetaBinSecond && i<=highThetaBinSecond) && (j<=highPhiBinSecond && j>=lowPhiBinSecond))){}
	else{
	  //	if ((i<lowThetaBin || i>highThetaBin) && (j<lowPhiBin && j>highPhiBin)){
	  if (mapCorVal[i][j]>tertiaryPeakVal){
	    tertiaryPeakVal=mapCorVal[i][j];
	    //tertiaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    tertiaryPeakPhi=ROUGH_BIN_SIZE*j;

	    if(cos_reconstruction_flag==1){
	      tertiaryPeakTheta=i*step_size_cos_rough+cosHighestTheta;//ROUGH_BIN_SIZE*i-HIGHEST_THETA; 
	    }
	    else{
	      tertiaryPeakTheta=ROUGH_BIN_SIZE*i-HIGHEST_THETA;
	    }

	  }
	}
      }
    }
  }

  if (printFlag==1) cout<<"Tertiary peak: "<<tertiaryPeakVal<<", theta: "<<tertiaryPeakTheta
			<<", phi: "<<tertiaryPeakPhi<<endl;
 
}
///////////////////////////
void MyCorrelator::findPeakOfMapBin(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], double &peakVal, 
				 int &peakThetaBin, int &peakPhiBin)
{
  peakVal=0;
  for (int i=0;i<NUM_BINS_ROUGH_THETA;i++){
    for (int j=0;j<NUM_BINS_ROUGH_PHI;j++){
      if (mapCorVal[i][j]>peakVal){
	peakVal=mapCorVal[i][j];
	peakThetaBin=i;
	peakPhiBin=j;
      }      
    }
  }
  if(printFlag==1) cout<<"peakVal of map is "<<peakVal<<"\n";
}
//////////////////////////////////////////
void MyCorrelator::findPeakOfRefinedMap(double mapCorVal[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI], double &peakVal, 
					int &peakThetaBin, int &peakPhiBin)
{
  peakVal=0;
  for (int i=0;i<NUM_BINS_FINE_THETA;i++){
    for (int j=0;j<NUM_BINS_FINE_PHI;j++){
      if (mapCorVal[i][j]>peakVal){
	peakVal=mapCorVal[i][j];
	peakThetaBin=i;
	peakPhiBin=j;
      }      
    }
  }
  if(printFlag==1) cout<<"peakVal refined is "<<peakVal<<"\n";
}
///////////////////////////////////////////////////////
void MyCorrelator::doInterpolationPeakFinding(double mapCorVal[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI], 
					      double &peakVal, 
					      double &peakTheta, double &peakPhi, double &FWHMTheta, 
					      double &FWHMPhi, double peakThetaRough, 
					      double peakPhiRough, int drawFlag)
{
 
  
  const int npointsTheta=NUM_BINS_FINE_THETA-19;
  const int npointsPhi=NUM_BINS_FINE_PHI-19;
  
  int peakThetaBin, peakPhiBin;
  findPeakOfRefinedMap(mapCorVal, peakVal, peakThetaBin, peakPhiBin);
 
  //assigns peakTheta,peakPhi,peakVal to peak of refined map bin.  Then we can do the interpolation
  
  //double refinedPeakTheta=peakThetaBin*FINE_BIN_SIZE+(peakThetaRough-thetaStartOffset);
  //double refinedPeakPhi=peakPhiBin*FINE_BIN_SIZE+(peakPhiRough-phiStartOffset);
 
  double thetaArray[npointsTheta];
  double phiArray[npointsPhi];
  double thetaValArray[npointsTheta], phiValArray[npointsPhi];
  cout<<"max_theta_cos is "<<max_theta_cos<<" step_size is "<<step_size_cos_fine<<"\n";
  for (int i=0;i<npointsTheta;i++){//construct an array of 10 values in phi and theta around the peak
    if (int(peakThetaBin-npointsTheta/2.)>=0 && int(peakThetaBin+npointsTheta/2.)<NUM_BINS_FINE_THETA){
      thetaValArray[i]=mapCorVal[int(peakThetaBin-npointsTheta/2+i)][peakPhiBin];
      
      if(cos_reconstruction_flag==1){
	//cout<<"Theta is "<<int(peakThetaBin-npointsTheta/2+i)*step_size_cos_fine+max_theta_cos<<"\n";
	thetaArray[i] =int(peakThetaBin-npointsTheta/2+i)*step_size_cos_fine+max_theta_cos;
	
      }
      else{
	thetaArray[i]=int(peakThetaBin-npointsTheta/2+i)*FINE_BIN_SIZE+(peakThetaRough-thetaStartOffset);
      }
      
    }
    
    else if(int(peakThetaBin-npointsTheta/2.)<0){
      thetaValArray[i]=mapCorVal[i][peakPhiBin];
      if(cos_reconstruction_flag==1){	
	//cout<<"Theta is "<< i*step_size_cos_fine+max_theta_cos<<"\n";
	thetaArray[i] = i*step_size_cos_fine+max_theta_cos;
      }
      else{
	thetaArray[i]=i*FINE_BIN_SIZE+(peakThetaRough-thetaStartOffset);
      }
      
    }
    else{
      thetaValArray[i]=mapCorVal[NUM_BINS_FINE_THETA-(npointsTheta-i)][peakPhiBin];
      
      if(cos_reconstruction_flag==1){
	//cout<<"Theta is "<<(NUM_BINS_FINE_THETA-(npointsTheta-i))*step_size_cos_fine+max_theta_cos<<"\n";
	thetaArray[i]=(NUM_BINS_FINE_THETA-(npointsTheta-i))*step_size_cos_fine+max_theta_cos;
      }
      else{
	thetaArray[i]=(NUM_BINS_FINE_THETA-(npointsTheta-i))*FINE_BIN_SIZE+(peakThetaRough-thetaStartOffset);
      }
     
    }
  }    
  for (int i=0;i<npointsPhi;i++){
    if (int(peakPhiBin-npointsPhi/2.)>=0 && int(peakPhiBin+npointsPhi/2.)<NUM_BINS_FINE_PHI){
      phiValArray[i]=mapCorVal[peakThetaBin][int(peakPhiBin-npointsPhi/2+i)];
      phiArray[i]=int(peakPhiBin-npointsPhi/2+i)*FINE_BIN_SIZE+(peakPhiRough-phiStartOffset);
    }
    
    else if (int(peakPhiBin-npointsPhi/2.)<0){ 
      phiValArray[i]=mapCorVal[peakThetaBin][i];
      phiArray[i]=i*FINE_BIN_SIZE+(peakPhiRough-phiStartOffset);
    }
    else{
      phiValArray[i]=mapCorVal[peakThetaBin][NUM_BINS_FINE_PHI-(npointsPhi-i)];
      phiArray[i]=(NUM_BINS_FINE_PHI-(npointsPhi-i))*FINE_BIN_SIZE+(peakPhiRough-phiStartOffset);
    }
    //cout<<"i: "<<i<<", phi: "<<phiArray[i]<<", phiVal: "<<phiValArray[i]<<endl;
    //cout<<"i: "<<i<<", theta: "<<thetaArray[i]<<", thetaVal: "<<thetaValArray[i]<<endl;
  }

  //do interpolation with akima to get peak and find FWHM or something like that.
  std::vector<double> thetaVector(npointsTheta);
  std::vector<double> thetaValVector(npointsTheta);
  std::vector<double> phiVector(npointsPhi);
  std::vector<double> phiValVector(npointsPhi);
  for (int i=0;i<npointsTheta;i++){
    thetaVector[i]=thetaArray[i];
    thetaValVector[i]=thetaValArray[i];
  }
 for (int i=0;i<npointsPhi;i++){
   phiVector[i]=phiArray[i];
   phiValVector[i]=phiValArray[i];
 }

  ROOT::Math::Interpolator interpolatorTheta(thetaVector.size(), ROOT::Math::Interpolation::kAKIMA );
  interpolatorTheta.SetData(thetaVector,thetaValVector);
  ROOT::Math::Interpolator interpolatorPhi(phiVector.size(), ROOT::Math::Interpolation::kAKIMA );
  interpolatorPhi.SetData(phiVector,phiValVector);
  
  int upsampleFactor=100;
  int interpCtr=0;
  double interpThetaArray[npointsTheta*upsampleFactor], interpThetaValArray[npointsTheta*upsampleFactor];
  double interpPhiArray[npointsPhi*upsampleFactor], interpPhiValArray[npointsPhi*upsampleFactor];
  if(cos_reconstruction_flag==1){
    for (double interp_theta=thetaArray[0];interp_theta<thetaArray[npointsTheta-1];interp_theta+=step_size_cos_fine/double(upsampleFactor)){
      interpThetaArray[interpCtr]=interp_theta;
      interpThetaValArray[interpCtr]=interpolatorTheta.Eval(interp_theta);
      //cout<<"theta, val is "<<interpThetaArray[interpCtr]<<" "<<interpThetaValArray[interpCtr]<<"\n";
      interpCtr++;
    }
  }
  else{
    for (double interp_theta=thetaArray[0];interp_theta<thetaArray[npointsTheta-1];interp_theta+=FINE_BIN_SIZE/double(upsampleFactor)){
      interpThetaArray[interpCtr]=interp_theta;
      interpThetaValArray[interpCtr]=interpolatorTheta.Eval(interp_theta);
      //cout<<"theta, val is "<<interpThetaArray[interpCtr]<<" "<<interpThetaValArray[interpCtr]<<"\n";
      interpCtr++;
    }
  }
  interpCtr=0;
  for (double interp_phi=phiArray[0];interp_phi<phiArray[npointsPhi-1];interp_phi+=FINE_BIN_SIZE/double(upsampleFactor)){
    interpPhiArray[interpCtr]=interp_phi;
    interpPhiValArray[interpCtr]=interpolatorPhi.Eval(interp_phi);
    interpCtr++;
  }


  //make a graph of
  TGraph *gInterpTheta=new TGraph((npointsTheta-1)*upsampleFactor,interpThetaArray,interpThetaValArray);
  TGraph *gInterpPhi=new TGraph((npointsPhi-1)*upsampleFactor,interpPhiArray,interpPhiValArray);

  if (drawFlag==1){ 
    TCanvas *cInterpolation=new TCanvas("cInterpolation","cInterpolation",800,800);
    cInterpolation->Divide(1,2);
    cInterpolation->cd(1);
    gInterpTheta->Draw("aP");
    gInterpTheta->SetMarkerStyle(22);
    gInterpTheta->SetMarkerSize(0.5);
    cInterpolation->cd(2);
    gInterpPhi->SetMarkerStyle(22);
    gInterpPhi->SetMarkerSize(0.5);
    gInterpPhi->Draw("aP");
    cInterpolation->Print("interpolatedPeak.eps");

  }

  //get the peak of the interpolated graph.
  Int_t peakBinTheta=FFTtools::getPeakBin(gInterpTheta);
  Double_t peakValTheta;
  gInterpTheta->GetPoint(peakBinTheta,peakTheta,peakValTheta);
  Int_t peakBinPhi=FFTtools::getPeakBin(gInterpPhi);
  Double_t peakValPhi;
  gInterpPhi->GetPoint(peakBinPhi,peakPhi,peakValPhi);

  delete gInterpTheta;
  delete gInterpPhi;
  
  if (printFlag==1) cout<<"Peak of interpolation in theta: "<<peakTheta<<", and phi: "<<peakPhi<<endl;
  if (peakValTheta<peakValPhi) peakVal=peakValTheta;
  else peakVal=peakValPhi;

  //get the FWHM of the peak
  double halfPeakThetaLeft=0; 
  double halfPeakThetaRight=0;  
  
  for (int i=0;i<(npointsTheta-1)*upsampleFactor;i++){
    if (peakBinTheta-i>0){
      if (interpThetaValArray[peakBinTheta-i]<peakValTheta/2. 
	  && interpThetaValArray[peakBinTheta-(i-1)]>peakValTheta/2.){
	halfPeakThetaLeft=(interpThetaArray[peakBinTheta-i]+interpThetaArray[peakBinTheta-(i-1)])/2.;
      }
    }
    if (halfPeakThetaLeft!=0) break;
  }
  //if (halfPeakThetaLeft==0) halfPeakThetaLeft=interpThetaArray[0];
  
  for (int i=0;i<(npointsTheta-1)*upsampleFactor;i++){
    if (peakBinTheta+i<(npointsTheta-1)*upsampleFactor){
      if (interpThetaValArray[peakBinTheta+i]<peakValTheta/2. 
	  && interpThetaValArray[peakBinTheta+(i-1)]>peakValTheta/2.){
	halfPeakThetaRight=(interpThetaArray[peakBinTheta+i]+interpThetaArray[peakBinTheta+(i-1)])/2.;
      }    
    }
    if (halfPeakThetaRight!=0) break;
  }
  // if (halfPeakThetaRight==0) halfPeakThetaRight=interpThetaArray[(npointsTheta-1)*upsampleFactor];

  double halfPeakPhiLeft=0; 
  double halfPeakPhiRight=0;  
  for (int i=0;i<(npointsPhi-1)*upsampleFactor;i++){
    if (peakBinPhi-i>0){
      if (interpPhiValArray[peakBinPhi-i]<peakValPhi/2. && interpPhiValArray[peakBinPhi-(i-1)]>peakValPhi/2.){
	halfPeakPhiLeft=(interpPhiArray[peakBinPhi-i]+interpPhiArray[peakBinPhi-(i-1)])/2.;
      }
    }
    if (halfPeakPhiLeft!=0) break;
  }
  //if (halfPeakPhiLeft==0) halfPeakPhiLeft=interpPhiArray[0];
  
  for (int i=0;i<(npointsPhi-1)*upsampleFactor;i++){
    if (peakBinPhi+i<(npointsPhi-1)*upsampleFactor){
      if (interpPhiValArray[peakBinPhi+i]<peakValPhi/2. && interpPhiValArray[peakBinPhi+(i-1)]>peakValPhi/2.){
	halfPeakPhiRight=(interpPhiArray[peakBinPhi+i]+interpPhiArray[peakBinPhi+(i-1)])/2.;
      }    
    }
    if (halfPeakPhiRight!=0) break;
  }
  // if (halfPeakPhiRight==0) halfPeakPhiRight=interpPhiArray[(npointsPhi-1)*upsampleFactor];
  
  if (halfPeakThetaRight==0 || halfPeakThetaLeft==0) FWHMTheta=-1;
  else FWHMTheta=halfPeakThetaRight-halfPeakThetaLeft;
  if (halfPeakPhiRight==0 || halfPeakPhiLeft==0) FWHMPhi=-1;
  else FWHMPhi=halfPeakPhiRight-halfPeakPhiLeft;
  if (printFlag==1) cout<<"FWHM Theta: "<<FWHMTheta<<", FWHM Phi: "<<FWHMPhi<<endl;
  
}
//////////////////////////////////////
double MyCorrelator::getSNROfMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI]){
  double peakVal;
  int peakThetaBin, peakPhiBin;
  int numBinsToCut=5; //around the peak on each side in theta and phi=10 degrees
  double RMS=0;
  int numPointsRMS=0;
  
  findPeakOfMapBin(mapCorVal, peakVal, peakThetaBin, peakPhiBin);
  //now cut out a few points around the peak
  int lowThetaBin=peakThetaBin-numBinsToCut;
  int highThetaBin=peakThetaBin+numBinsToCut;
  int lowPhiBin=peakPhiBin-numBinsToCut;
  int highPhiBin=peakPhiBin+numBinsToCut;
  if (lowPhiBin<0) lowPhiBin+=NUM_BINS_ROUGH_PHI;
  if (highPhiBin>NUM_BINS_ROUGH_PHI) highPhiBin-=NUM_BINS_ROUGH_PHI;

  for (int i=0;i<NUM_BINS_ROUGH_THETA;i++){
    for (int j=0;j<NUM_BINS_ROUGH_PHI;j++){
      if (lowPhiBin<highPhiBin) 
	if ((i<lowThetaBin || i>highThetaBin) && (j<lowPhiBin || j<highPhiBin)){
	  RMS+=mapCorVal[i][j]*mapCorVal[i][j];
	  if (mapCorVal[i][j]!=0) numPointsRMS++;
	}
      if (lowPhiBin>highPhiBin)
	if ((i<lowThetaBin || i>highThetaBin) && (j<lowPhiBin && j>highPhiBin)){
	  RMS+=mapCorVal[i][j]*mapCorVal[i][j];      
	  if (mapCorVal[i][j]!=0) numPointsRMS++;
	}
    }
  }
  RMS=sqrt(RMS/numPointsRMS);
  double SNR=peakVal/RMS;
  if (printFlag==1 ) cout<<"For Main Map SNR: peakVal: "<<peakVal<<", RMS: "<<RMS<<", SNR: "<<SNR<<endl;
  return SNR;
  

}
//////////////////////////////////////////
int MyCorrelator::traceBackTo0Altitude(int eventNumber, Double_t thetaWave, Double_t phiWave, Double_t &sourceLon, Double_t &sourceLat){
  if (!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);

  int groundPointingFlag=fUsefulAdu5Ptr->getSourceLonAndLatAltZero(phiWave, thetaWave, sourceLon, sourceLat);

  return groundPointingFlag;

}
/////////////////////////////////////
int MyCorrelator::traceBackToContinent_Brian(int eventNumber, double thetaWave, double phiWave, double &sourceLon, double &sourceLat, double &sourceAlt)
{
  if (!fAntarctica) initializeAntarctica();
  
   double anitaLatitude=fAdu5APatPtr->latitude;
   double anitaLongitude=fAdu5APatPtr->longitude;
   double anitaAltitude=fAdu5APatPtr->altitude;
   double anitaHeading=fAdu5APatPtr->heading;
   

   // cout<<"HACKED TRACED BACK TO CONTINENT!!!! \n\n\n\n";
   //cout<<"lat, lon, alt, heading are "<<anitaLatitude<<" "<<anitaLongitude<<" "<<anitaAltitude<<" "<<anitaHeading<<"\n";
   
   /*anitaLatitude = -77.7154;
   anitaLongitude = 167.905;
   anitaAltitude = 37069.6;
   anitaHeading = 234.971;
   */
   double pitch=-0.29;
   double roll=0.89;
   //cout<<"pitch and roll are "<<pitch<<" "<<roll<<"\n";
   //cout<<"theta is "<<thetaWave<<"\n";
   double heading = anitaHeading;
   //cout<<"heading is "<<heading<<"\n";
   double tempPhiWave=phiWave;
   double tempThetaWave=TMath::PiOver2()+thetaWave;//- thetaWave
   cout<<"lat, lon, alt, heading, thetawave,phiwave are "<<anitaLatitude<<" "<<anitaLongitude<<" "<<anitaAltitude<<" "<<anitaHeading<<" "<<tempThetaWave<<" "<<tempPhiWave<<"\n";
  
   //tempThetaWave = .558195;
   //tempPhiWave = 4.39482;
    //cout<<"tempThetaWave is "<<tempThetaWave<<"\n";
   
   

   //tempPhiWave = 0;
   //tempThetaWave = 110*TMath::DegToRad();
   // cout<<"tempThetaWave is "<<tempThetaWave<<"\n";
   //tempThetaWave=TMath::Pi();
   //Now need to take account of balloon heading


   TVector3 rollAxis,pitchAxis;
   rollAxis=fUPGeomTool->fRollRotationAxis;
   pitchAxis=fUPGeomTool->fPitchRotationAxis;
  
   TVector3 headingAxis;
   headingAxis =fUPGeomTool->fHeadingRotationAxis;

   pitchAxis = -1*pitchAxis; //PitchAxis for positive rotation to the clockwise of Roll (aviation)
   headingAxis = -1*headingAxis;//HeadingAxis for positive rotation is vertically downward

   
   double Anita_Cart[3];
   TVector3 AnitaPos;
   double source[3];
  
 
   TVector3 rollAxis2,pitchAxis2,headingAxis2;
   
   //cout<<"rollAxis is "<<rollAxis[0]<<" "<<rollAxis[1]<<" "<<rollAxis[2]<<"\n";
   fUPGeomTool->getCartesianCoords(anitaLatitude,anitaLongitude,anitaAltitude,Anita_Cart);
   AnitaPos.SetXYZ(Anita_Cart[1],Anita_Cart[0],Anita_Cart[2]);
   
   if(heading>=0 && heading<=360) {
     
     TVector3 arbDir;
     arbDir.SetMagThetaPhi(1,tempThetaWave,tempPhiWave);//-1*tempPhiWave
     //cout<<"ArbDir is "<<arbDir[0]<<" "<<arbDir[1]<<" "<<arbDir[2]<<"\n";
     TVector3 BalloonPos;
     BalloonPos.SetXYZ(Anita_Cart[1],Anita_Cart[0],Anita_Cart[2]);
     double BalloonTheta = BalloonPos.Theta();
     double BalloonPhi =BalloonPos.Phi();
      //AVIATION: to put into correct position, Correct Heading, then Pitch, then Roll
     arbDir.Rotate(heading*TMath::DegToRad(),headingAxis);
     rollAxis.Rotate((heading)*TMath::DegToRad(),headingAxis);
     pitchAxis.Rotate((heading)*TMath::DegToRad(),headingAxis);
     // cout<<"Theta,Phi are "<<arbDir.Theta()<<" "<<arbDir.Phi()<<"\n";
     //cout<<"rollAxis is "<<rollAxis[0]<<" "<<rollAxis[1]<<" "<<rollAxis[2]<<"\n";
     arbDir.Rotate(pitch*TMath::DegToRad(),pitchAxis);
     rollAxis.Rotate(pitch*TMath::DegToRad(),pitchAxis);
     headingAxis.Rotate(pitch*TMath::DegToRad(),pitchAxis);
     //cout<<"Theta,Phi are "<<arbDir.Theta()<<" "<<arbDir.Phi()<<"\n";
     //cout<<"rollAxis is "<<rollAxis[0]<<" "<<rollAxis[1]<<" "<<rollAxis[2]<<"\n";
     arbDir.Rotate(roll*TMath::DegToRad(),rollAxis);//roll and pitch
     headingAxis.Rotate(roll*TMath::DegToRad(),rollAxis);
     pitchAxis.Rotate(roll*TMath::DegToRad(),rollAxis);
     //cout<<"Theta,Phi are "<<arbDir.Theta()<<" "<<arbDir.Phi()<<"\n";
     // cout<<"rollAxis is "<<rollAxis[0]<<" "<<rollAxis[1]<<" "<<rollAxis[2]<<"\n";

      arbDir.RotateY(BalloonTheta);
      headingAxis.RotateY(BalloonTheta);
      rollAxis.RotateY(BalloonTheta);
      pitchAxis.RotateY(BalloonTheta);
     
      arbDir.RotateZ(-1*BalloonPhi);
      headingAxis.RotateZ(-1*BalloonPhi);
      rollAxis.RotateZ(-1*BalloonPhi);
      pitchAxis.RotateZ(-1*BalloonPhi);
     
     
      tempPhiWave=arbDir.Phi();
     
     //std::cout << "tempPhiWave2: " << tempPhiWave << "\n";
     if(tempPhiWave>TMath::TwoPi()) {
       tempPhiWave-=TMath::TwoPi();
     }
     if(tempPhiWave<0) {
       tempPhiWave+=TMath::TwoPi();
       }
     tempThetaWave=arbDir.Theta();
   }
   
   tempPhiWave = -1*tempPhiWave;//Flip back to correct coordinate system
   
   double true_source[3];
   //2: -11.9744 -87.114 2846.97
   //3: -87.114,168.026,2846.97
   //4: 14.3373,-80.6344,3267.95
   double trueSource_lat= -87.114;
   double trueSource_lon= 168.026;
   double trueSource_alt= fRampdemReader->SurfaceAboveGeoid(trueSource_lon,trueSource_lat);
    
   fUPGeomTool->getCartesianCoords(trueSource_lat,trueSource_lon,trueSource_alt,true_source);
   cout<<"true_source is "<<true_source[1]<<" "<<true_source[0]<<" "<<true_source[2]<<"\n";
   //Try to use vectors to find source Location:
   //cout<<" swapped x and y \n\n";
  
   double x_dir = sin(tempThetaWave)*cos(tempPhiWave);
   double y_dir = sin(tempThetaWave)*sin(tempPhiWave);
   double z_dir = cos(tempThetaWave);
  
   double step_size=1000;
   double surface_height;
   
   int step_max = 10000;
   
   source[1]=Anita_Cart[1];
   source[0]=Anita_Cart[0];
   source[2]=Anita_Cart[2];
  
   cout<<"Anita is at "<<Anita_Cart[1]<<" "<<Anita_Cart[0]<<" "<<Anita_Cart[2]<<"\n";
   cout<<"x_dir, y_dir, z_dir is "<<x_dir<<" "<<y_dir<<" "<<z_dir<<"\n";
   //cout<<"tempThetaWave, Phi are "<<tempThetaWave<<" "<<tempPhiWave<<"\n";
   for(int step=0;step<step_max;step++){
    
     source[1] = source[1]+step_size*x_dir;
     source[0] = source[0]+step_size*y_dir;
     source[2] = source[2]+step_size*z_dir;
       
     
     fUPGeomTool->getLatLonAltFromCartesian(source,sourceLat,sourceLon,sourceAlt);
     
    
     if(sourceLat > -60){
       cout<<"MISSED ANTARCTICA! \n";
       return 0;
     }
    
     surface_height =  fRampdemReader->SurfaceAboveGeoid(sourceLon,sourceLat);
    
      if(surface_height <-9990){
       cout<<"MISSED ANTARCTICA! \n";
       return 0;
     }
     if(surface_height > sourceAlt){
       break;
     }
      if(step==step_max-1){
       cout<<"missed the entire time? surface_height is "<<surface_height<<" altitude is "<<sourceAlt<<"\n";
       return 0;
     }
   }

   source[1] = source[1]-step_size*x_dir;
   source[0] = source[0]-step_size*y_dir;
   source[2] = source[2]-step_size*z_dir;
   
   step_size = 100;
   //cout<<"STEP BACK! \n";
   for(int step=0;step<step_max;step++){
     
     source[1] = source[1]+step_size*x_dir;
     source[0] = source[0]+step_size*y_dir;
     source[2] = source[2]+step_size*z_dir;
     
     fUPGeomTool->getLatLonAltFromCartesian(source,sourceLat,sourceLon,sourceAlt);
     // cout<<" Source Lat Lon, alt is "<<sourceLat<<" "<<sourceLon<<" "<<sourceAlt<<"\n";
     surface_height =  fRampdemReader->SurfaceAboveGeoid(sourceLon,sourceLat);
     
     if(surface_height > sourceAlt){
       break;
     }
   }
   
   source[1] = source[1]-step_size*x_dir;
   source[0] = source[0]-step_size*y_dir;
   source[2] = source[2]-step_size*z_dir;
   step_size=1;
   // cout<<"STEP BACK! \n";
   for(int step=0;step<step_max;step++){
     
     source[1] = source[1]+step_size*x_dir;
     source[0] = source[0]+step_size*y_dir;
     source[2] = source[2]+step_size*z_dir;
     
     fUPGeomTool->getLatLonAltFromCartesian(source,sourceLat,sourceLon,sourceAlt);
     //cout<<" Source Lat Lon, alt is "<<sourceLat<<" "<<sourceLon<<" "<<sourceAlt<<"\n";
     surface_height =  fRampdemReader->SurfaceAboveGeoid(sourceLon,sourceLat);
     
     
     if(surface_height > sourceAlt){ 
       distance_from_source = sqrt(pow(source[1]-Anita_Cart[1],2)+ pow(source[0]-Anita_Cart[0],2)+pow(source[2]-Anita_Cart[2],2));
       cout<<"source is "<<source[1]<<" "<<source[0]<<" "<<source[2]<<"\n";
       cout<<"sourceLat,Lon are "<<sourceLat<<" "<<sourceLon<<"\n";
       return 1;  
     }
   }
  
   return 2;//IF 2 returned, something awkward happened!
   
  
} 
//////////////////////////////////
int MyCorrelator::checkIfNearAnyBase(int eventNumber,double peakThetaFinal, double peakPhiFinal, double sourceLon, //positive theta is down, angle in radians
				  double sourceLat, double sourceAlt, std::string &baseName)
{ 
  if (baseLatitudeAndLongitude[0][0]==0) initializeBaseList();
  if (!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  double thetaWave, phiWave;
  double distance;
  double distanceBase;
  int nearBaseFlag=0;
  double loglikelihood;
  
  for (int ctr=0;ctr<1000;ctr++){
    if(eventNumber==3182406) cout<<"ctr is "<<ctr<<"\n";
    if (baseLatitudeAndLongitude[ctr][0]==0) break;
    fUsefulAdu5Ptr->getThetaAndPhiWave(baseLatitudeAndLongitude[ctr][1], baseLatitudeAndLongitude[ctr][0], 
				       baseHeights[ctr], thetaWave, phiWave);//pos theta is down
    distance=findDistanceBetweenTwoThings(sourceLat,sourceLon,sourceAlt, baseLatitudeAndLongitude[ctr][0], 
					  baseLatitudeAndLongitude[ctr][1],baseHeights[ctr]);
    distanceBase=findDistanceBetweenTwoThings(fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude, baseLatitudeAndLongitude[ctr][0], 
					    baseLatitudeAndLongitude[ctr][1],baseHeights[ctr]);
    
    loglikelihood=((thetaWave*rad2deg-peakThetaFinal*rad2deg)*
		   (thetaWave*rad2deg-peakThetaFinal*rad2deg)/(sigmaTheta[sigmaIndex]*sigmaTheta[sigmaIndex])
		   +(phiWave*rad2deg-peakPhiFinal*rad2deg)*(phiWave*rad2deg-peakPhiFinal*rad2deg)
		   /(sigmaPhi[sigmaIndex]*sigmaPhi[sigmaIndex]));
    if ((loglikelihood<40 || distance<40e3) && distanceBase<800e3){
      //   if (((fabs((thetaWave*rad2deg)-peakThetaFinal*rad2deg)<5*sigmaTheta[sigmaIndex] 
      //  && fabs(phiWave*rad2deg-peakPhiFinal*rad2deg)<5*sigmaPhi[sigmaIndex]) || distance<50e3) && distanceBase<800e3){ 
      baseName=baseNames[ctr];
      if (printFlag==1) cout<<"Event is near base: "<<baseName<<endl;
      nearBaseFlag=ctr+1;
      break;
      // return 1;//distance has to be less than 800 km
    }
  }
  
  return nearBaseFlag;
}
///////////////////////////////////////
int MyCorrelator::checkIfNearSpecificBase(int eventNumber,double peakThetaFinal, double peakPhiFinal, double sourceLon, //positive theta is down, angle in radians
					  double sourceLat, double sourceAlt, std::string baseName, int baseIndex)
{ 
  if (baseLatitudeAndLongitude[0][0]==0) initializeBaseList();
  if (!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  double thetaWave=0;
  double phiWave=0;
  double distance=0;
  int nearBaseFlag=0;
  double distanceBase=0;
  
  for (int ctr=0;ctr<1000;ctr++){
    if (baseLatitudeAndLongitude[ctr][0]==0) break;
    if (baseNames[ctr]==baseName || baseIndex==ctr){
      fUsefulAdu5Ptr->getThetaAndPhiWave(baseLatitudeAndLongitude[ctr][1], baseLatitudeAndLongitude[ctr][0], 
					 baseHeights[ctr], thetaWave, phiWave);//pos theta is down
      distance=findDistanceBetweenTwoThings(sourceLat,sourceLon,sourceAlt, baseLatitudeAndLongitude[ctr][0], 
					    baseLatitudeAndLongitude[ctr][1],baseHeights[ctr]);
      distanceBase=findDistanceBetweenTwoThings(fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,fAdu5APatPtr->altitude, baseLatitudeAndLongitude[ctr][0], 
					    baseLatitudeAndLongitude[ctr][1],baseHeights[ctr]);
      
      if (((fabs((thetaWave*rad2deg)-peakThetaFinal*rad2deg)<5*sigmaTheta[sigmaIndex] 
	    && fabs(phiWave*rad2deg-peakPhiFinal*rad2deg)<5*sigmaPhi[sigmaIndex]) || distance<30e3) && distanceBase<800e3){ 
	if (printFlag==1) cout<<"Event is near base: "<<baseName<<", peakThetaFinal: "<<peakThetaFinal*rad2deg<<", peakPhiFinal: "<<peakPhiFinal*rad2deg
			      <<", thetaWave "<<thetaWave*rad2deg<<", phiwave: "<<phiWave*rad2deg<<endl;
	nearBaseFlag=1;
      }
      // return 1;//distance has to be less than 800 km
    }
  }
  
  return nearBaseFlag;
}
/////////////////////////////////
void MyCorrelator::drawEventOnContinent(int eventNumber,double sourceLon,double sourceLat)
{

  if (!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  
  double xAnita,yAnita, xEvent, yEvent;
  // cout<<"Maps anita is at: "<<fAdu5APatPtr->latitude<<", "<<fAdu5APatPtr->longitude<<endl;

  getRelXYFromLatLong(fAdu5APatPtr->latitude,fAdu5APatPtr->longitude,xAnita,yAnita);
  getRelXYFromLatLong(static_cast<float>(sourceLat),static_cast<float>(sourceLon),xEvent,yEvent);

  gStyle->SetMarkerColor(kBlack);
  gStyle->SetTextSize(0.02);

  //this is located in  /rh5stuff/64bit/src/anita/eventCorrelator/.  Root *should* search for this properly.
  TImage *map = TImage::Open("macros/antarcticaIceMapBW.png");
  TCanvas *canMap=(TCanvas*)gROOT->FindObject("canMap");
  if(!canMap){
    canMap = new TCanvas("canMap","canMap",(int)xSize,(int)ySize);
    canMap->Clear();
    canMap->SetLogz();
    canMap->SetTopMargin(0);
    canMap->SetBottomMargin(0);
    canMap->SetLeftMargin(0);
    canMap->SetRightMargin(0);
    map->Draw("");
  }
 
  TMarker *anitaPos = new TMarker(xAnita,yAnita,23);
  TMarker *eventPos = new TMarker(xEvent,yEvent,23);
  anitaPos->SetMarkerColor(kRed);
  eventPos->SetMarkerColor(kGreen);
  eventPos->Draw("");
  anitaPos->DrawMarker(xAnita,yAnita);
  canMap->Print("eventOnContinent.eps");

}

///////////////////////////////////
void MyCorrelator::backProjectEvent(int eventNumber, double sourceLon, double sourceLat, 
				    double sourceHeight, double &thetaWaveProj,
				    double &phiWaveProj)
{
  if (!fEventTree) initialize();
  if (eventStartedFlag!=eventNumber) eventStartedFlag=startEachEvent(eventNumber);
  
  fUsefulAdu5Ptr->getThetaAndPhiWave(sourceLon, sourceLat, sourceHeight, thetaWaveProj, phiWaveProj);
  //cout<<"sourceLon: "<<sourceLon<<", sourceLat: "<<sourceLat<<", sourceHeight: "<<sourceHeight
  //  <<"thetaWave: "<<thetaWaveProj*rad2deg<<"phi wave: "<<phiWaveProj*rad2deg<<endl;
  thetaWaveProj=thetaWaveProj*rad2deg*(-1.);
  phiWaveProj=phiWaveProj*rad2deg;

}
////////////////////////////
TGraph *MyCorrelator::deconvolveWaveform(TGraph *grCoherent, int polFlag, int drawFlag)
{
  //  int upSampleFactor=1;
  Double_t deltaTInt=1./(2.6);
  //cout<<grCoherent->GetN()<<", "<<2.6*260/grCoherent->GetN()<<", "<<1/(2.6*260/grCoherent->GetN())<<endl;
  //Double_t deltaTInt=1./(2.6*double(upSampleFactor)*260/grCoherent->GetN());
  
  FILE *fp1;
  //  if (polFlag==0) fp1= fopen("/rh5stuff/64bit/src/anita/analysis/agoodhue/systemImpulseResponse/systemImpulseResponseInvertedAvgNew.txt","r");//vpol
  if (polFlag==0) fp1= fopen("${ANITA_ANALYSIS_AGOODHUE}/systemImpulseResponse/systemImpulseResponseInvertedAvgNew.txt","r");//vpol
  // else fp1= fopen("/rh5stuff/64bit/src/anita/analysis/agoodhue/systemImpulseResponse/systemImpulseResponseInvertedAvgNewHoriz.txt","r");//hpol
  else fp1= fopen("${ANITA_ANALYSIS_AGOODHUE}/systemImpulseResponse/systemImpulseResponseInvertedAvgNewHoriz.txt","r");//hpol
  //std::ifstream impulseResponse_file(filename);
  int npoints=260;  
  int ncols=0;
  double voltageTemp, timeTemp;
  double voltageImpulse[npoints*3], timeImpulse[npoints*3];
  for (int i=0;i<npoints*3;i++){
    voltageImpulse[i]=0;
    timeImpulse[i]=0;
  }
  
  //get system response from a file
  for (int i=0;i<npoints;i++)
    {
      ncols=fscanf(fp1,"%lf\t%lf",&timeTemp,&voltageTemp);
      timeImpulse[i]=timeTemp;
      voltageImpulse[i]=voltageTemp;
    }
  for (int i=0;i<npoints*3;i++){
    if (i>npoints-2) timeImpulse[i]=timeImpulse[i-1]-timeImpulse[i-2]+timeImpulse[i-1];
  }
  
  TGraph *grImpulseResponse=new TGraph(npoints*3,timeImpulse,voltageImpulse);
  double *oldYImpulse = grImpulseResponse->GetY();
  double *oldXImpulse = grImpulseResponse->GetX();
  double deltaTImpulse=oldXImpulse[1]-oldXImpulse[0];
  int lengthImpulse=grImpulseResponse->GetN();
  FFTWComplex *theFFTImpulse=FFTtools::doFFT(lengthImpulse,oldYImpulse);
  int newLengthImpulse=(lengthImpulse/2)+1;
  double deltaFImpulse=1/(deltaTImpulse*lengthImpulse); //Hz
  deltaFImpulse*=1e3; //MHz
  
  double magImpulse[newLengthImpulse], phaseImpulse[newLengthImpulse],frequencyImpulse[newLengthImpulse];
  
  for (int i=0;i<newLengthImpulse;i++){
    if (i>0) frequencyImpulse[i]=frequencyImpulse[i-1]+deltaFImpulse;
    magImpulse[i]=sqrt(theFFTImpulse[i].re*theFFTImpulse[i].re+theFFTImpulse[i].im*theFFTImpulse[i].im);
    phaseImpulse[i]=atan(theFFTImpulse[i].im/theFFTImpulse[i].re);
  }
  
  if (drawFlag==1){
    TCanvas *cimpulse=new TCanvas("cimpulse","cimpulse",800,800);
    cimpulse->Divide(1,3);
    cimpulse->cd(1);
    grImpulseResponse->Draw("al");
    TGraph *grmagImpulse=new TGraph(newLengthImpulse,frequencyImpulse,magImpulse);
    TGraph *grphaseImpulse=new TGraph(newLengthImpulse,frequencyImpulse,phaseImpulse);
    cimpulse->cd(2);
    grmagImpulse->Draw("al");
    cimpulse->cd(3);
    grphaseImpulse->Draw("al");
    
  }

  //if (printFlag==1) cout<<"deltaT response (ns): "<<deltaTImpulse<<",deltaF response (MHZ): "
  //		<<deltaFImpulse<<", npoints new response: "<<newLengthImpulse
  //		<<", npoints response: "<<lengthImpulse<<endl;
  
  //interpolate coherently summed waveform to the correct frequency binning
  TGraph *grInterp;
  grInterp=FFTtools::getInterpolatedGraph(grCoherent, deltaTInt);
  //cout<<"deltaTint: "<<deltaTInt<<endl;
  //take FFT of the coherent waveform,, zero padded to 300 ns
  int newCoherentLength=260*3;
  //  int newCoherentLength=grCoherent->GetN()*(300/100);
  double newY[newCoherentLength], newX[newCoherentLength];
  double timeThis, voltThis;
  for (int i=0;i<newCoherentLength;i++){
    newY[i]=0;
    newX[i]=0;
  }
  
  for (int i=0;i<newCoherentLength;i++){
    if (i<grInterp->GetN()){
      grInterp->GetPoint(i,timeThis,voltThis);
      // newY[i+1000]=voltThis;
      // newX[i]=timeThis-20;
      newY[i+260]=voltThis;
      newX[i]=timeThis;
    }
    else{
      //      newY[i]=0;
      newX[i]=newX[i-1]-newX[i-2]+newX[i-1];
    }
  }
  TGraph *grNewCoherent=new TGraph(newCoherentLength,newX,newY);

  double *oldY = grNewCoherent->GetY();
  double *oldX = grNewCoherent->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grNewCoherent->GetN();
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  int newLength=(length/2)+1;
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  //cout<<"deltaF (MHz): "<<deltaF<<endl;
  //cout<<"new length: "<<newLength<<endl;
  
  //define variables for FFT of coherent waveform
  double frequencyCoherent[newLength];//MHz
  double magnitudeCoherent[newLength];
  double phaseCoherent[newLength];
  frequencyCoherent[0]=0;
  
  for(int i=0;i<newLength;i++) {
    if (i>0) frequencyCoherent[i]=frequencyCoherent[i-1]+deltaF;
    magnitudeCoherent[i]=sqrt(theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im);
    phaseCoherent[i]=atan(theFFT[i].im/theFFT[i].re);
    //if (frequencyCoherent[i]<1300) 
	//cout<<"freq: "<<frequencyCoherent[i]<<", mag: "<<magnitudeCoherent[i]<<", phase: "<<phaseCoherent[i]<<endl;
  }
  //if (printFlag==1) cout<<"deltaT coherent (ns): "<<deltaT<<",deltaF coherent (MHZ): "<<deltaF<<", npoints new coherent: "<<newLength
  //		<<", npoints coherent: "<<length<<endl;
  
  if (drawFlag==1){
    TCanvas *ccoherent=new TCanvas("ccoherent","ccoherent",800,800);
  ccoherent->Divide(1,3);
  ccoherent->cd(1);
  grNewCoherent->Draw("al");
  TGraph *grmagCoherent=new TGraph(newLength,frequencyCoherent,magnitudeCoherent);
  TGraph *grphaseCoherent=new TGraph(newLength,frequencyCoherent,phaseCoherent);
  ccoherent->cd(2);
  grmagCoherent->Draw("al");
  ccoherent->cd(3);
  grphaseCoherent->Draw("al");
  
  }

  FFTWComplex *thedeconvolvedFFT=new FFTWComplex[(length/2)+1];
  //now deconvolve by applying system response
  for(int i=0;i<newLength;i++) {
    if (frequencyCoherent[i]<1200 && magImpulse[i]>0 && frequencyCoherent[i]>300){
      magnitudeCoherent[i]=magnitudeCoherent[i]/magImpulse[i];
      phaseCoherent[i]=phaseCoherent[i]-phaseImpulse[i];
      thedeconvolvedFFT[i].re=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].re+theFFTImpulse[i].im*theFFT[i].im);
      thedeconvolvedFFT[i].im=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].im-theFFTImpulse[i].im*theFFT[i].re);
    }
    else{
      magnitudeCoherent[i]=0;
      phaseCoherent[i]=0;
      thedeconvolvedFFT[i].re=0;
      thedeconvolvedFFT[i].im=0;
    }
  }

  
  double *thedeconvolvedWaveform = FFTtools::doInvFFT(length,thedeconvolvedFFT);
  TGraph *grdeconvolved=new TGraph(length,oldX,thedeconvolvedWaveform);  
  TGraph *grmagdecon=new TGraph(newLengthImpulse,frequencyCoherent,magnitudeCoherent);
  TGraph *grphasedecon=new TGraph(newLengthImpulse,frequencyCoherent,phaseCoherent);
  
  if (drawFlag==1){
    TCanvas *cdeconvolved=new TCanvas("cdeconvolved","cdeconvolved",800,800);
    cdeconvolved->Divide(1,3);
    cdeconvolved->cd(1);
    grdeconvolved->Draw("al");
    cdeconvolved->cd(2);
    grmagdecon->Draw("al");
    cdeconvolved->cd(3);
    grphasedecon->Draw("al");
  }
  
  return grdeconvolved;
}
////////////////////////////
TGraph *MyCorrelator::deconvolveWaveformUsingStephens(TGraph *grCoherent, int polFlag, int drawFlag)
{
  //  int upSampleFactor=1;
  Double_t deltaTInt=1./(2.6);
  //cout<<grCoherent->GetN()<<", "<<2.6*260/grCoherent->GetN()<<", "<<1/(2.6*260/grCoherent->GetN())<<endl;
  //Double_t deltaTInt=1./(2.6*double(upSampleFactor)*260/grCoherent->GetN());
  
  FILE *fp1;
  if (polFlag==0) fp1=fopen("/home/dailey.110/analysis/systemImpulseResponse/systemImpulseResponseUsingStephens.txt","r");
  //if (polFlag==1) fp1= fopen("/home/dailey.110/abby/agoodhue/hooverCode/Aesop/data/impulse_responses/ch1hAbby.tpl","r");
  //else fp1 = fopen("/home/dailey.110/abby/agoodhue/anitaIICode/myCorrelator/systemImpulseResponse/systemImpulseResponseNoAntenna.txt","r");
  else 
    fp1= fopen("/home/dailey.110/analysis/anitaIICode/myCorrelator/systemImpulseResponse/ch1hAbby.tpl","r");
 

 //std::ifstream impulseResponse_file(filename);
  int npoints=260;
  int npointsStephen=905;
  if (polFlag==0) npointsStephen=260; //vpol use mine
  if (polFlag==1) npointsStephen=905;
  // npointsStephen=905; //hpol use stephens

  int ncols=0;
  double voltageTemp, timeTemp;
  double voltageImpulse[npointsStephen], timeImpulse[npointsStephen];
  for (int i=0;i<npointsStephen;i++){
    voltageImpulse[i]=0;
    timeImpulse[i]=0;
  }
  
  //get system response from a file
  for (int i=0;i<npointsStephen;i++)
    {
      ncols=fscanf(fp1,"%lf\t%lf",&timeTemp,&voltageTemp);
      timeImpulse[i]=timeTemp;
      voltageImpulse[i]=voltageTemp;
    }
  if (polFlag==1){
    timeImpulse[0]=10;
    for (int i=1;i<npointsStephen;i++){
      timeImpulse[i]=0.1+timeImpulse[i-1];
    }
  }
  
  double interpImpulse[npoints*3];
  double interpTime[npoints*3];
  for (int i=0;i<npoints*3;i++) interpImpulse[i]=0;
  for (int i=0;i<npoints*3;i++) interpTime[i]=0;

  //cout<<"deltat/time: "<<deltaTInt/(timeImpulse[1]-timeImpulse[0])<<endl;
  //cout<<"time: "<<timeImpulse[1]-timeImpulse[0]<<endl;

  TGraph *grImpulseResponse=new TGraph(npointsStephen,timeImpulse,voltageImpulse);
  TGraph *grInterpImpulse;
  grInterpImpulse=FFTtools::getInterpolatedGraph(grImpulseResponse, deltaTInt);
  double *YImpulse=grInterpImpulse->GetY();
  double *XImpulse=grInterpImpulse->GetX();
  for (int i=0;i<npoints*3;i++){
    if (i<grInterpImpulse->GetN()) {
      interpTime[i]=XImpulse[i];
      if (polFlag==1) interpImpulse[i]=YImpulse[i]*deltaTInt/(timeImpulse[1]-timeImpulse[0]);
      else interpImpulse[i]=YImpulse[i]*deltaTInt/(1/2.6)*3.846;
    }
    else interpTime[i]=interpTime[i-1]+interpTime[i-1]-interpTime[i-2];
  }
  
  TGraph *grInterpImpulse2=new TGraph(npoints*3,interpTime, interpImpulse);
  double *oldYImpulse = grInterpImpulse2->GetY();
  double *oldXImpulse = grInterpImpulse2->GetX();
  double deltaTImpulse=oldXImpulse[1]-oldXImpulse[0];
  int lengthImpulse=grInterpImpulse2->GetN();
  FFTWComplex *theFFTImpulse=FFTtools::doFFT(lengthImpulse,oldYImpulse);
  int newLengthImpulse=(lengthImpulse/2)+1;
  double deltaFImpulse=1/(deltaTImpulse*lengthImpulse); //Hz
  deltaFImpulse*=1e3; //MHz
  
  double magImpulse[newLengthImpulse], phaseImpulse[newLengthImpulse],frequencyImpulse[newLengthImpulse];
  
  frequencyImpulse[0]=0;
  for (int i=0;i<newLengthImpulse;i++){
    if (i>0) frequencyImpulse[i]=frequencyImpulse[i-1]+deltaFImpulse;
    magImpulse[i]=sqrt(theFFTImpulse[i].re*theFFTImpulse[i].re+theFFTImpulse[i].im*theFFTImpulse[i].im);
    phaseImpulse[i]=atan(theFFTImpulse[i].im/theFFTImpulse[i].re);
  }
  
  if (drawFlag==1){
    TCanvas *cimpulse=new TCanvas("cimpulse","cimpulse",800,800);
    cimpulse->Divide(1,3);
    cimpulse->cd(1);
    grInterpImpulse2->GetXaxis()->SetTitle("time (ns)");
    grInterpImpulse2->GetYaxis()->SetTitle("Effective Height (m/ns)");
    grInterpImpulse2->Draw("al");
    TGraph *grmagImpulse=new TGraph(newLengthImpulse,frequencyImpulse,magImpulse);
    TGraph *grphaseImpulse=new TGraph(newLengthImpulse,frequencyImpulse,phaseImpulse);
    cimpulse->cd(2);
    grmagImpulse->GetXaxis()->SetTitle("frequency (MHz)");
    grmagImpulse->GetYaxis()->SetTitle("Fourier Transform Magnitude");
    grmagImpulse->Draw("al");
    cimpulse->cd(3);
    grphaseImpulse->GetXaxis()->SetTitle("frequency (MHz)");
    grphaseImpulse->GetYaxis()->SetTitle("Phase Angle (Rad)");
    grphaseImpulse->Draw("al");
    
  }
  //cout<<"deltaT signalchain: "<<deltaTImpulse<<", deltaF signalchain: "<<deltaFImpulse
  //  <<", new length signalchain: "<<newLengthImpulse<<", length singlachain: "<<grInterpImpulse2->GetN()<<endl;

  char filename[150];
  if (polFlag==1) sprintf(filename,"/home/dailey.110/analysis/picoSecondImpulse/antenna_hh_impulse_response.dat");
  else sprintf(filename,"/home/dailey.110/analysis/picoSecondImpulse/antenna_impulse_response.dat");
  //cout<<filename<<endl;
  std::ifstream antennaResponse_file(filename);
  int ctr=0;
  int npointsAntenna=1000;
  double timeAntenna[npointsAntenna*3], voltageAntenna[npointsAntenna*3];
  for (int i=0;i<npointsAntenna*3;i++){
    timeAntenna[i]=0;
    voltageAntenna[i]=0;
  }
  
  //get antenna response from  file
  while(antennaResponse_file >> timeTemp >> voltageTemp)
    {
      timeAntenna[ctr]=timeTemp*1e9+850-29.6; //ns
      voltageAntenna[ctr+100]=voltageTemp;//mV
      if (ctr==npointsAntenna-1) break;
      ctr++;
    }
  
  for (int i=0;i<npointsAntenna*3;i++){  
    if (i>=npointsAntenna-10) timeAntenna[i]=timeAntenna[i-1]-timeAntenna[i-2]+timeAntenna[i-1];
  }

  TGraph *grAntenna=new TGraph(npointsAntenna*3,timeAntenna,voltageAntenna);  
  //now draw this system response waveform
  TGraph *grInterpAntenna;
  grInterpAntenna=FFTtools::getInterpolatedGraph(grAntenna, deltaTInt);
  
  double *voltageInterpAntenna=grInterpAntenna->GetY();
  double *timeInterpAntenna=grInterpAntenna->GetX();
  //cout<<"deltaTint/timeant "<<deltaTInt/(timeAntenna[1]-timeAntenna[0])<<endl;
  for (int i=0;i<grInterpAntenna->GetN();i++)
    voltageInterpAntenna[i]*=deltaTInt/(timeAntenna[1]-timeAntenna[0]);

  //take the fourier transform to add it on to the system response
  FFTWComplex *theFFTAntenna=FFTtools::doFFT(grInterpAntenna->GetN(),voltageInterpAntenna);
  double deltaTAntenna=timeInterpAntenna[1]-timeInterpAntenna[0];
  double deltaFAntenna=1/(deltaTAntenna*grInterpAntenna->GetN());
  int newLengthAntenna=(grInterpAntenna->GetN()/2)+1;
  deltaFAntenna*=1e3;
  double frequencyAntenna[newLengthAntenna];
  double magnitudeAntenna[newLengthAntenna];
  double phaseAntenna[newLengthAntenna];
  frequencyAntenna[0]=0;

  //fill antenna response arrays
  for (int i=0;i<newLengthAntenna;i++){
    frequencyAntenna[i]=frequencyAntenna[i-1]+deltaFAntenna;
    magnitudeAntenna[i]=sqrt(theFFTAntenna[i].re*theFFTAntenna[i].re+theFFTAntenna[i].im*theFFTAntenna[i].im);
    phaseAntenna[i]=atan(theFFTAntenna[i].im/theFFTAntenna[i].re);
  }
  //  cout<<"deltaT antenna: "<<deltaTAntenna<<", deltaF antenna: "<<deltaFAntenna
  //  <<", new length antenna: "<<newLengthAntenna<<", length antenna: "<<grInterpAntenna->GetN()<<endl;
   
  //antenna fft graphs
  TGraph *grantennaFFTMag=new TGraph(newLengthAntenna,frequencyAntenna,magnitudeAntenna);
  TGraph *grantennaFFTPhase=new TGraph(newLengthAntenna,frequencyAntenna,phaseAntenna);
  if (drawFlag==1){
    TCanvas *cantennaFFT=new TCanvas("cantennaFFT","cantennaFFT",800,800);
    cantennaFFT->Divide(1,3);
    cantennaFFT->cd(1);
    grInterpAntenna->Draw("al");
    cantennaFFT->cd(2);
    grantennaFFTMag->Draw("al");
    cantennaFFT->cd(3);
    grantennaFFTPhase->Draw("al");
    
   }

  //////////////////NOW I'VE READ IN ANTENNA AND SIGNAL CHAIN RESPONSE INTO theFFTAntenna and theFFTImpulse   MUST TAKE OUT EACH to DECONVOLVE
 
  //interpolate waveform of interest to the correct frequency binning
  TGraph *grInterp;
  grInterp=FFTtools::getInterpolatedGraph(grCoherent, deltaTInt);
  //cout<<"deltaTint: "<<deltaTInt<<endl;
  //take FFT of the coherent waveform,, zero padded to 300 ns
  int newCoherentLength=npoints*3;
  //  int newCoherentLength=grCoherent->GetN()*(300/100);
  double newY[newCoherentLength], newX[newCoherentLength];
  double timeThis, voltThis;
  for (int i=0;i<newCoherentLength;i++){
    newY[i]=0;
    newX[i]=0;
  }
  
  for (int i=0;i<newCoherentLength;i++){
    if (i<grInterp->GetN()){
      grInterp->GetPoint(i,timeThis,voltThis);
      // newY[i+1000]=voltThis;
      // newX[i]=timeThis-20;
      newY[i+npoints]=voltThis;
      newX[i]=timeThis;
    }
    else{
      //      newY[i]=0;
      newX[i]=newX[i-1]-newX[i-2]+newX[i-1];
    }
  }
  TGraph *grNewCoherent=new TGraph(newCoherentLength,newX,newY);

  double *oldY = grNewCoherent->GetY();
  double *oldX = grNewCoherent->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grNewCoherent->GetN();
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  int newLength=(length/2)+1;
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  //cout<<"deltaF (MHz): "<<deltaF<<endl;
  //cout<<"new length: "<<newLength<<endl;
  
  //define variables for FFT of coherent waveform
  double frequencyCoherent[newLength];//MHz
  double magnitudeCoherent[newLength];
  double phaseCoherent[newLength];
  frequencyCoherent[0]=0;
  
  for(int i=0;i<newLength;i++) {
    if (i>0) frequencyCoherent[i]=frequencyCoherent[i-1]+deltaF;
    magnitudeCoherent[i]=sqrt(theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im);
    phaseCoherent[i]=atan(theFFT[i].im/theFFT[i].re);
    //if (frequencyCoherent[i]<1300) 
	//cout<<"freq: "<<frequencyCoherent[i]<<", mag: "<<magnitudeCoherent[i]<<", phase: "<<phaseCoherent[i]<<endl;
  }
  //if (printFlag==1) cout<<"deltaT coherent (ns): "<<deltaT<<",deltaF coherent (MHZ): "<<deltaF<<", npoints new coherent: "<<newLength
  //	<<", npoints coherent: "<<length<<endl;
  
  if (drawFlag==1){
    TCanvas *ccoherent=new TCanvas("ccoherent","ccoherent",800,800);
    ccoherent->Divide(1,3);
    ccoherent->cd(1);
    grNewCoherent->Draw("al");
    TGraph *grmagCoherent=new TGraph(newLength,frequencyCoherent,magnitudeCoherent);
    TGraph *grphaseCoherent=new TGraph(newLength,frequencyCoherent,phaseCoherent);
    ccoherent->cd(2);
    grmagCoherent->Draw("al");
    ccoherent->cd(3);
    grphaseCoherent->Draw("al");
    
  }
  
  ///now that I've gotten graph in question and shoved it inot theFFT, do deconvolution

  
  FFTWComplex *thedeconvolvedFFT=new FFTWComplex[(length/2)+1];
  FFTWComplex *thedeconvolvedFFTStep1=new FFTWComplex[(length/2)+1];
  //now deconvolve by applying system response
  for(int i=0;i<newLength;i++) {
    if (polFlag==0){
      if (frequencyCoherent[i]<1200 && frequencyCoherent[i]>200){
	magnitudeCoherent[i]=magnitudeCoherent[i]/magImpulse[i]/magnitudeAntenna[i];
	phaseCoherent[i]=phaseCoherent[i]-phaseImpulse[i]-phaseAntenna[i];
	thedeconvolvedFFTStep1[i].re=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].re+theFFTImpulse[i].im*theFFT[i].im)*(2*(1200-frequencyCoherent[i])/(1000.)+1);
	thedeconvolvedFFTStep1[i].im=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].im-theFFTImpulse[i].im*theFFT[i].re)*(2*(1200-frequencyCoherent[i])/(1000.)+1);
	
	thedeconvolvedFFT[i].re=1./(pow(magnitudeAntenna[i],2))*(theFFTAntenna[i].re*thedeconvolvedFFTStep1[i].re+theFFTAntenna[i].im*thedeconvolvedFFTStep1[i].im)*-1.;//*sqrt(377/50);
	thedeconvolvedFFT[i].im=1./(pow(magnitudeAntenna[i],2))*(theFFTAntenna[i].re*thedeconvolvedFFTStep1[i].im-theFFTAntenna[i].im*thedeconvolvedFFTStep1[i].re)*-1.;//*sqrt(377/50);
      }
      else{
	magnitudeCoherent[i]=0;
	phaseCoherent[i]=0;
	thedeconvolvedFFT[i].re=0;
	thedeconvolvedFFT[i].im=0;
      }
    }
    else {//hpol
      if (frequencyCoherent[i]<1200 && frequencyCoherent[i]>200){
	magnitudeCoherent[i]=magnitudeCoherent[i]/magImpulse[i]/magnitudeAntenna[i];
	phaseCoherent[i]=phaseCoherent[i]-phaseImpulse[i]-phaseAntenna[i];
	thedeconvolvedFFTStep1[i].re=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].re+theFFTImpulse[i].im*theFFT[i].im)*(2*(1200-frequencyCoherent[i])/(1000.)+1);
	thedeconvolvedFFTStep1[i].im=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].im-theFFTImpulse[i].im*theFFT[i].re)*(2*(1200-frequencyCoherent[i])/(1000.)+1);
	
	thedeconvolvedFFT[i].re=1./(pow(magnitudeAntenna[i],2))*(theFFTAntenna[i].re*thedeconvolvedFFTStep1[i].re+theFFTAntenna[i].im*thedeconvolvedFFTStep1[i].im)/2.;//*sqrt(377/50)/2.;//splitter not there
	thedeconvolvedFFT[i].im=1./(pow(magnitudeAntenna[i],2))*(theFFTAntenna[i].re*thedeconvolvedFFTStep1[i].im-theFFTAntenna[i].im*thedeconvolvedFFTStep1[i].re)/2.;//*sqrt(377/50)/2.;
      }
      else{
	magnitudeCoherent[i]=0;
	phaseCoherent[i]=0;
	thedeconvolvedFFT[i].re=0;
	thedeconvolvedFFT[i].im=0;
      }
    }
  }
 
  double *thedeconvolvedWaveform = FFTtools::doInvFFT(length,thedeconvolvedFFT);
  TGraph *grdeconvolved=new TGraph(length,oldX,thedeconvolvedWaveform);  
  TGraph *grmagdecon=new TGraph(newLengthImpulse,frequencyCoherent,magnitudeCoherent);
  TGraph *grphasedecon=new TGraph(newLengthImpulse,frequencyCoherent,phaseCoherent);
  
  if (drawFlag==1){
    TCanvas *cdeconvolved=new TCanvas("cdeconvolved","cdeconvolved",800,800);
    cdeconvolved->Divide(1,3);
    cdeconvolved->cd(1);
    grdeconvolved->Draw("al");
    cdeconvolved->cd(2);
    grmagdecon->Draw("al");
    cdeconvolved->cd(3);
    grphasedecon->Draw("al");
  }  
  return grdeconvolved;

}
////////////////////////////
TGraph *MyCorrelator::deconvolveWaveformUsingRyans(TGraph *grCoherent, int polFlag, int drawFlag)
{
  Double_t deltaTInt=1./(2.6);  
  char filenameRyan[150];
  sprintf(filenameRyan,"systemImpulseResponse/sumPicoImpulse.root");
  TFile *rootfileRyan=new TFile(filenameRyan,"READ");
  TGraph *grRyan;
  if (polFlag==0) grRyan=(TGraph*)rootfileRyan->Get("grImpRespV");
  else grRyan=(TGraph*)rootfileRyan->Get("grImpRespH");
  
  int npoints=260;
  double interpTime[npoints*3];
  double interpImpulse[npoints*3];
  for (int i=0;i<npoints*3;i++){ 
    interpTime[i]=0;
    interpImpulse[i]=0;
  }

  //interpolate ryan's to deltaTInt
  TGraph *grInterpImpulse;
  grInterpImpulse=FFTtools::getInterpolatedGraph(grRyan, deltaTInt);
  double *YImpulse=grInterpImpulse->GetY();
  double *XImpulse=grInterpImpulse->GetX();
  for (int i=0;i<npoints*3;i++){
    if (i<grInterpImpulse->GetN()) {
      interpTime[i]=XImpulse[i];
      interpImpulse[i]=YImpulse[i]*deltaTInt;
    }
    else interpTime[i]=interpTime[i-1]+interpTime[i-1]-interpTime[i-2];
  }
  
  //now we've lengthened to 3x, so do FFT
  TGraph *grInterpImpulse2=new TGraph(npoints*3,interpTime, interpImpulse);
  double *oldYImpulse = grInterpImpulse2->GetY();
  double *oldXImpulse = grInterpImpulse2->GetX();
  double deltaTImpulse=oldXImpulse[1]-oldXImpulse[0];
  int lengthImpulse=grInterpImpulse2->GetN();
  FFTWComplex *theFFTImpulse=FFTtools::doFFT(lengthImpulse,oldYImpulse);
  int newLengthImpulse=(lengthImpulse/2)+1;
  double deltaFImpulse=1/(deltaTImpulse*lengthImpulse); //Hz
  deltaFImpulse*=1e3; //MHz
  
  double magImpulse[newLengthImpulse], phaseImpulse[newLengthImpulse],frequencyImpulse[newLengthImpulse];
  
  frequencyImpulse[0]=0;
  for (int i=0;i<newLengthImpulse;i++){
    if (i>0) frequencyImpulse[i]=frequencyImpulse[i-1]+deltaFImpulse;
    magImpulse[i]=sqrt(theFFTImpulse[i].re*theFFTImpulse[i].re+theFFTImpulse[i].im*theFFTImpulse[i].im);
    phaseImpulse[i]=atan(theFFTImpulse[i].im/theFFTImpulse[i].re);
  }
  
  if (drawFlag==1){
    TCanvas *cimpulse=new TCanvas("cimpulse","cimpulse",800,800);
    cimpulse->Divide(1,3);
    cimpulse->cd(1);
    grInterpImpulse2->Draw("al");
    TGraph *grmagImpulse=new TGraph(newLengthImpulse,frequencyImpulse,magImpulse);
    TGraph *grphaseImpulse=new TGraph(newLengthImpulse,frequencyImpulse,phaseImpulse);
    cimpulse->cd(2);
    grmagImpulse->Draw("al");
    cimpulse->cd(3);
    grphaseImpulse->Draw("al");
    
  }
  //cout<<"deltaT signalchain: "<<deltaTImpulse<<", deltaF signalchain: "<<deltaFImpulse
  //  <<", new length signalchain: "<<newLengthImpulse<<", length singlachain: "<<grInterpImpulse2->GetN()<<endl;
  

  //now get antenna responses for H and V
  char filename[150];
  if (polFlag==1) sprintf(filename,"/home/dailey.110/analysis/picoSecondImpulse/antenna_hh_impulse_response.dat");
  else sprintf(filename,"/home/dailey.110/analysis/picoSecondImpulse/antenna_impulse_response.dat");
  std::ifstream antennaResponse_file(filename);
  int ctr=0;
  int npointsAntenna=1000;
  double timeAntenna[npointsAntenna*3], voltageAntenna[npointsAntenna*3];
  for (int i=0;i<npointsAntenna*3;i++){
    timeAntenna[i]=0;
    voltageAntenna[i]=0;
  }
  double timeTemp, voltageTemp;
  
  //get antenna response from  file
  while(antennaResponse_file >> timeTemp >> voltageTemp)
    {
      timeAntenna[ctr]=timeTemp*1e9+850-29.6; //ns
      voltageAntenna[ctr+100]=voltageTemp;//mV
      if (ctr==npointsAntenna-1) break;
      ctr++;
    }
  
  for (int i=0;i<npointsAntenna*3;i++){  
    if (i>=npointsAntenna-10) timeAntenna[i]=timeAntenna[i-1]-timeAntenna[i-2]+timeAntenna[i-1];
  }

  TGraph *grAntenna=new TGraph(npointsAntenna*3,timeAntenna,voltageAntenna);  
  //now draw this system response waveform
  TGraph *grInterpAntenna;
  grInterpAntenna=FFTtools::getInterpolatedGraph(grAntenna, deltaTInt);
  
  double *voltageInterpAntenna=grInterpAntenna->GetY();
  double *timeInterpAntenna=grInterpAntenna->GetX();
  
  //cout<<"deltaTint/timeant "<<deltaTInt/(timeAntenna[1]-timeAntenna[0])<<endl;
  for (int i=0;i<grInterpAntenna->GetN();i++)
    voltageInterpAntenna[i]*=deltaTInt/(timeAntenna[1]-timeAntenna[0]);
  
  //take the fourier transform 
  FFTWComplex *theFFTAntenna=FFTtools::doFFT(grInterpAntenna->GetN(),voltageInterpAntenna);
  double deltaTAntenna=timeInterpAntenna[1]-timeInterpAntenna[0];
  double deltaFAntenna=1/(deltaTAntenna*grInterpAntenna->GetN());
  int newLengthAntenna=(grInterpAntenna->GetN()/2)+1;
  deltaFAntenna*=1e3;
  double frequencyAntenna[newLengthAntenna];
  double magnitudeAntenna[newLengthAntenna];
  double phaseAntenna[newLengthAntenna];
  frequencyAntenna[0]=0;

  //fill antenna response arrays
  for (int i=0;i<newLengthAntenna;i++){
    frequencyAntenna[i]=frequencyAntenna[i-1]+deltaFAntenna;
    magnitudeAntenna[i]=sqrt(theFFTAntenna[i].re*theFFTAntenna[i].re+theFFTAntenna[i].im*theFFTAntenna[i].im);
    phaseAntenna[i]=atan(theFFTAntenna[i].im/theFFTAntenna[i].re);
  }
  //cout<<"deltaT antenna: "<<deltaTAntenna<<", deltaF antenna: "<<deltaFAntenna
  //  <<", new length antenna: "<<newLengthAntenna<<", length antenna: "<<grInterpAntenna->GetN()<<endl;
   
  //antenna fft graphs
  TGraph *grantennaFFTMag=new TGraph(newLengthAntenna,frequencyAntenna,magnitudeAntenna);
  TGraph *grantennaFFTPhase=new TGraph(newLengthAntenna,frequencyAntenna,phaseAntenna);
  if (drawFlag==1){
    TCanvas *cantennaFFT=new TCanvas("cantennaFFT","cantennaFFT",800,800);
    cantennaFFT->Divide(1,3);
    cantennaFFT->cd(1);
    grInterpAntenna->Draw("al");
    cantennaFFT->cd(2);
    grantennaFFTMag->Draw("al");
    cantennaFFT->cd(3);
    grantennaFFTPhase->Draw("al");
    
   }

  //////////////////NOW I'VE READ IN ANTENNA AND SIGNAL CHAIN RESPONSE INTO theFFTAntenna and theFFTImpulse   MUST TAKE OUT EACH to DECONVOLVE
 
  //interpolate waveform of interest to the correct frequency binning
  TGraph *grInterp;
  grInterp=FFTtools::getInterpolatedGraph(grCoherent, deltaTInt);
  //cout<<"deltaTint: "<<deltaTInt<<endl;
  //cout<<grCoherent->GetN()<<", "<<grInterp->GetN()<<endl;
  //take FFT of the coherent waveform,, zero padded to 300 ns
  int newCoherentLength=npoints*3;
  //  int newCoherentLength=grCoherent->GetN()*(300/100);
  double newY[newCoherentLength], newX[newCoherentLength];
  double timeThis, voltThis;
  for (int i=0;i<newCoherentLength;i++){
    newY[i]=0;
    newX[i]=0;
  }
  
  for (int i=0;i<newCoherentLength;i++){
    if (i<grInterp->GetN()){
      grInterp->GetPoint(i,timeThis,voltThis);
      // newY[i+1000]=voltThis;
      // newX[i]=timeThis-20;
      newY[i+npoints]=voltThis;
      newX[i]=timeThis;
    }
    else{
      //      newY[i]=0;
      newX[i]=newX[i-1]-newX[i-2]+newX[i-1];
    }
  }
  TGraph *grNewCoherent=new TGraph(newCoherentLength,newX,newY);

  double *oldY = grNewCoherent->GetY();
  double *oldX = grNewCoherent->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grNewCoherent->GetN();
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  int newLength=(length/2)+1;
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  //cout<<"deltaF (MHz): "<<deltaF<<endl;
  //cout<<"new length: "<<newLength<<endl;
  
  //define variables for FFT of coherent waveform
  double frequencyCoherent[newLength];//MHz
  double magnitudeCoherent[newLength];
  double phaseCoherent[newLength];
  frequencyCoherent[0]=0;
  
  for(int i=0;i<newLength;i++) {
    if (i>0) frequencyCoherent[i]=frequencyCoherent[i-1]+deltaF;
    magnitudeCoherent[i]=sqrt(theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im);
    phaseCoherent[i]=atan(theFFT[i].im/theFFT[i].re);
    //if (frequencyCoherent[i]<1300) 
	//cout<<"freq: "<<frequencyCoherent[i]<<", mag: "<<magnitudeCoherent[i]<<", phase: "<<phaseCoherent[i]<<endl;
  }
  //if (printFlag==1) cout<<"deltaT coherent (ns): "<<deltaT<<",deltaF coherent (MHZ): "<<deltaF<<", npoints new coherent: "<<newLength
  //	<<", npoints coherent: "<<length<<endl;
  
  if (drawFlag==1){
    TCanvas *ccoherent=new TCanvas("ccoherent","ccoherent",800,800);
    ccoherent->Divide(1,3);
    ccoherent->cd(1);
    grNewCoherent->Draw("al");
    TGraph *grmagCoherent=new TGraph(newLength,frequencyCoherent,magnitudeCoherent);
    TGraph *grphaseCoherent=new TGraph(newLength,frequencyCoherent,phaseCoherent);
    ccoherent->cd(2);
    grmagCoherent->Draw("al");
    ccoherent->cd(3);
    grphaseCoherent->Draw("al");
    
  }
  
  ///now that I've gotten graph in question and shoved it inot theFFT, do deconvolution

  
  FFTWComplex *thedeconvolvedFFT=new FFTWComplex[(length/2)+1];
  FFTWComplex *thedeconvolvedFFTStep1=new FFTWComplex[(length/2)+1];
  //now deconvolve by applying system response
  for(int i=0;i<newLength;i++) {
    if (polFlag==0){
      if (frequencyCoherent[i]<1200 && frequencyCoherent[i]>200){
	magnitudeCoherent[i]=magnitudeCoherent[i]/magImpulse[i]/magnitudeAntenna[i];
	phaseCoherent[i]=phaseCoherent[i]-phaseImpulse[i]-phaseAntenna[i];
	thedeconvolvedFFTStep1[i].re=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].re+theFFTImpulse[i].im*theFFT[i].im)*(2*(1200-frequencyCoherent[i])/(1000.)+1);
	thedeconvolvedFFTStep1[i].im=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].im-theFFTImpulse[i].im*theFFT[i].re)*(2*(1200-frequencyCoherent[i])/(1000.)+1);
	
	thedeconvolvedFFT[i].re=1./(pow(magnitudeAntenna[i],2))*(theFFTAntenna[i].re*thedeconvolvedFFTStep1[i].re+theFFTAntenna[i].im*thedeconvolvedFFTStep1[i].im)*sqrt(377/50);
	thedeconvolvedFFT[i].im=1./(pow(magnitudeAntenna[i],2))*(theFFTAntenna[i].re*thedeconvolvedFFTStep1[i].im-theFFTAntenna[i].im*thedeconvolvedFFTStep1[i].re)*sqrt(377/50);
      }
      else{
	magnitudeCoherent[i]=0;
	phaseCoherent[i]=0;
	thedeconvolvedFFT[i].re=0;
	thedeconvolvedFFT[i].im=0;
      }
    }
    else {//hpol
      if (frequencyCoherent[i]<1200 && frequencyCoherent[i]>200){
	magnitudeCoherent[i]=magnitudeCoherent[i]/magImpulse[i]/magnitudeAntenna[i];
	phaseCoherent[i]=phaseCoherent[i]-phaseImpulse[i]-phaseAntenna[i];
	thedeconvolvedFFTStep1[i].re=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].re+theFFTImpulse[i].im*theFFT[i].im)*(2*(1200-frequencyCoherent[i])/(1000.)+1);
	thedeconvolvedFFTStep1[i].im=1./(pow(magImpulse[i],2))*(theFFTImpulse[i].re*theFFT[i].im-theFFTImpulse[i].im*theFFT[i].re)*(2*(1200-frequencyCoherent[i])/(1000.)+1);
	
	thedeconvolvedFFT[i].re=1./(pow(magnitudeAntenna[i],2))*(theFFTAntenna[i].re*thedeconvolvedFFTStep1[i].re+theFFTAntenna[i].im*thedeconvolvedFFTStep1[i].im)*sqrt(377/50)/2.;//splitter not there
	thedeconvolvedFFT[i].im=1./(pow(magnitudeAntenna[i],2))*(theFFTAntenna[i].re*thedeconvolvedFFTStep1[i].im-theFFTAntenna[i].im*thedeconvolvedFFTStep1[i].re)*sqrt(377/50)/2.;
      }
      else{
	magnitudeCoherent[i]=0;
	phaseCoherent[i]=0;
	thedeconvolvedFFT[i].re=0;
	thedeconvolvedFFT[i].im=0;
      }
    }
  }
 
  double *thedeconvolvedWaveform = FFTtools::doInvFFT(length,thedeconvolvedFFT);
  TGraph *grdeconvolved=new TGraph(length,oldX,thedeconvolvedWaveform);  
  TGraph *grmagdecon=new TGraph(newLengthImpulse,frequencyCoherent,magnitudeCoherent);
  TGraph *grphasedecon=new TGraph(newLengthImpulse,frequencyCoherent,phaseCoherent);
  
  if (drawFlag==1){
    TCanvas *cdeconvolved=new TCanvas("cdeconvolved","cdeconvolved",800,800);
    cdeconvolved->Divide(1,3);
    cdeconvolved->cd(1);
    grdeconvolved->Draw("al");
    cdeconvolved->cd(2);
    grmagdecon->Draw("al");
    cdeconvolved->cd(3);
    grphasedecon->Draw("al");
  }  
  return grdeconvolved;

}



/////////////////////////////////////////////////////// END OF MAIN POINTING CODE STUFF

//////////////NOT MAIN CODE STUFF
///////////////////////////
void MyCorrelator::getCRPolContent()
{
  /*double crTheta=16.178; // from horizontal  for ev 271
  double anLat=-83.2062;
  double anLon=-162.124;
  double anAlt=35482.3;
  double crLat=-82.3907;
  double crLon=-156.344;
  */
  /*double crTheta=24.898; // from horizontal  for ev 144
  double crLat=-75.4823;
  double crLon=-150.682;
  double anAlt=34939.8;
  double anLat=-76.093;
  double anLon=-149.756;
  */
  double crTheta=31.363; // from horizontal  for ev 145
  double anLat=-76.3077;
  double anLon=-151.849;
  double anAlt=35821;
  double crLat=-75.8289;
  double crLon=-152.756;
  
  Vector thisBField;
  if (!fAntarctica) initializeAntarctica();
  thisBField=fAntarctica->getBField(crLon,crLat);
  cout<<"B field values: ";
  thisBField.Print();
  
  //now get CR direction in terms of  X,Y,Z and stuff it in a vector;
  double x1, y1, z1, x2, y2, z2;
  LatLonAlt2xyz(anLat,anLon,anAlt,x1,y1,z1);
  LatLonAlt2xyz(crLat,crLon,fAntarctica->surfaceAboveGeoid(crLon, crLat),x2,y2,z2);
  
  cout<<"Anita Pos: x1: "<<x1<<", y1: "<<y1<<", z1: "<<z1<<endl;
  cout<<"CR Pos: x2: "<<x2<<", y2: "<<y2<<", z2: "<<z2<<endl;

  Vector CRDirBeforeRotation;
  CRDirBeforeRotation.SetX(x1-x2);
  CRDirBeforeRotation.SetY(y1-y2);
  CRDirBeforeRotation.SetZ(z1-z2);
  cout<<"CR Direction before rotation: ";
  CRDirBeforeRotation.Print();

  Vector rotationVector;
  rotationVector.SetZ(sin(anLat*deg2rad));
  rotationVector.SetY(cos(anLat*deg2rad)*cos(anLon*deg2rad));
  rotationVector.SetX(cos(anLat*deg2rad)*sin(anLon*deg2rad));

  cout<<"rotation vector: ";
  rotationVector.Print();
  
  Vector CRDirAfterRotation;
  
  CRDirAfterRotation=CRDirBeforeRotation.Rotate(3.1415, rotationVector);
  cout<<"cr direction after rotation: ";
  CRDirAfterRotation.Print();
  
  CRDirAfterRotation.SetX(-1.*CRDirAfterRotation.GetX());
  CRDirAfterRotation.SetY(-1.*CRDirAfterRotation.GetY());
  CRDirAfterRotation.SetZ(-1.*CRDirAfterRotation.GetZ());

  cout<<"cr direction after rotation and flip: ";
  CRDirAfterRotation.Print();

  Vector vCrossB;
  
  vCrossB=CRDirAfterRotation.Cross(thisBField);

  cout<<"v cross b: ";
  vCrossB.Print();

  //now that I have E=vxB, I need to project on to H and V from the payload.
  double HPolContent;
  double VPolContent;

  VPolContent=vCrossB.Dot(rotationVector);
  cout<<"Before Fresnel: vpol: "<<VPolContent;

  Vector HPol;
  
  HPol=vCrossB.Cross(rotationVector);
  HPolContent=HPol.Mag()*cos(crTheta*deg2rad);
  cout<<", hpol: "<<HPolContent<<", ratio V/H: "<<VPolContent/HPolContent<<endl;
  
  //add fresnel
  crTheta=(90-crTheta)*deg2rad;
  float n_ice=1.325;//of firn from icemc
  double fresnel_pokey=2*cos(crTheta)/(cos(crTheta)+sqrt(n_ice*n_ice-sin(crTheta)*sin(crTheta)));//reflection
  double fresnel_slappy=-2*n_ice*cos(crTheta)/(n_ice*n_ice*cos(crTheta)+sqrt(n_ice*n_ice-sin(crTheta)*sin(crTheta)));//reflection

  HPolContent*=fresnel_slappy;
  VPolContent*=fresnel_pokey;
  
  cout<<"pokey: "<<fresnel_pokey<<", slappy: "<<fresnel_slappy<<endl<<endl;;

  cout<<"VPol: "<<VPolContent<<", HPol: "<<HPolContent<<", Ratio V/H Expected: "<<VPolContent/HPolContent<<endl;


}
////////////////////////////////
void MyCorrelator::calcANITAIIUsingStephens()
{
  Double_t deltaTInt=1./(2.6);
  //cout<<grCoherent->GetN()<<", "<<2.6*260/grCoherent->GetN()<<", "<<1/(2.6*260/grCoherent->GetN())<<endl;
  //Double_t deltaTInt=1./(2.6*double(upSampleFactor)*260/grCoherent->GetN());
  
  FILE *fp1;
  fp1= fopen("/home/dailey.110/analysis/anitaIICode/myCorrelator/systemImpulseResponse/ch1hAbby.tpl","r");
  //std::ifstream impulseResponse_file(filename);
  int npoints=260;
  int npointsStephen=905;
  npointsStephen=905; //hpol use stephens

  int ncols=0;
  double voltageTemp, timeTemp;
  double voltageImpulse[npointsStephen], timeImpulse[npointsStephen];
  for (int i=0;i<npointsStephen;i++){
    voltageImpulse[i]=0;
    timeImpulse[i]=0;
  }
  
  //get system response from a file
  for (int i=0;i<npointsStephen;i++)
    {
      ncols=fscanf(fp1,"%lf\t%lf",&timeTemp,&voltageTemp);
      timeImpulse[i]=timeTemp;
      //cout<<timeTemp<<", voltage: "<<voltageTemp<<endl;
      voltageImpulse[i]=voltageTemp;
    }
  timeImpulse[0]=10;
  for (int i=1;i<npointsStephen;i++){
    timeImpulse[i]=0.1+timeImpulse[i-1];
  }
  
  double interpImpulse[npoints*3];
  double interpTime[npoints*3];
  for (int i=0;i<npoints*3;i++) interpImpulse[i]=0;
  for (int i=0;i<npoints*3;i++) interpTime[i]=0;
  
  TGraph *grImpulseResponse=new TGraph(npointsStephen,timeImpulse,voltageImpulse);
  TGraph *grInterpImpulse;
  grInterpImpulse=FFTtools::getInterpolatedGraph(grImpulseResponse, deltaTInt);
  double *YImpulse=grInterpImpulse->GetY();
  double *XImpulse=grInterpImpulse->GetX();
  for (int i=0;i<npoints*3;i++){
    if (i<grInterpImpulse->GetN()) {
      interpTime[i]=XImpulse[i];
      interpImpulse[i]=YImpulse[i];
    }
    else interpTime[i]=interpTime[i-1]+interpTime[i-1]-interpTime[i-2];
  }
  
  TGraph *grInterpImpulse2=new TGraph(npoints*3,interpTime, interpImpulse);
  double *oldYImpulse = grInterpImpulse2->GetY();
  double *oldXImpulse = grInterpImpulse2->GetX();
  double deltaTImpulse=oldXImpulse[1]-oldXImpulse[0];
  int lengthImpulse=grInterpImpulse2->GetN();
  FFTWComplex *theFFTImpulse=FFTtools::doFFT(lengthImpulse,oldYImpulse);
  int newLengthImpulse=(lengthImpulse/2)+1;
  double deltaFImpulse=1/(deltaTImpulse*lengthImpulse); //Hz
  deltaFImpulse*=1e3; //MHz
  cout<<"newlengthimpulse: "<<newLengthImpulse<<", deltaf: "<<deltaFImpulse<<endl;
  
  double magImpulse[newLengthImpulse], phaseImpulse[newLengthImpulse],frequencyImpulse[newLengthImpulse];
  
  frequencyImpulse[0]=0;
  for (int i=0;i<newLengthImpulse;i++){
    if (i>0) frequencyImpulse[i]=frequencyImpulse[i-1]+deltaFImpulse;
    magImpulse[i]=sqrt(theFFTImpulse[i].re*theFFTImpulse[i].re+theFFTImpulse[i].im*theFFTImpulse[i].im);
    phaseImpulse[i]=atan(theFFTImpulse[i].im/theFFTImpulse[i].re);
  }
  
  TCanvas *cimpulse=new TCanvas("cimpulse","cimpulse",800,800);
  cimpulse->Divide(1,3);
  cimpulse->cd(1);
  grInterpImpulse2->Draw("al");
  TGraph *grmagImpulse=new TGraph(newLengthImpulse,frequencyImpulse,magImpulse);
  TGraph *grphaseImpulse=new TGraph(newLengthImpulse,frequencyImpulse,phaseImpulse);
  cimpulse->cd(2);
  grmagImpulse->Draw("al");
  cimpulse->cd(3);
  grphaseImpulse->Draw("al");

  ///////////now read in gains of rfcm anita i
  int npointsANITAI=62;
  int npointsANITAII=68;
  float freqTemp, gainTemp, nfTemp;
  float freqANITAI[npointsANITAI];
  float gainANITAI[npointsANITAI];
  float freqANITAII[npointsANITAII];
  float gainANITAII[npointsANITAII];
  
  FILE *fp2;
  fp2= fopen("/home/dailey.110/analysis/hooverCode/Aesop/data/rfcm_data/rfcm02_1.dat","r");
  for (int i=0;i<npointsANITAI;i++)
    {
      ncols=fscanf(fp2,"%f,%f,%f",&freqTemp,&gainTemp,&nfTemp);
      freqANITAI[i]=freqTemp;
      if (gainTemp<1e9) gainANITAI[i]=gainTemp;
      else gainANITAI[i]=0;
    }

  FILE *fp3;
  fp3= fopen("/home/dailey.110/analysis/hooverCode/Aesop/data/rfcm_data/rfcm01_ch2_pam01.dat","r");
  for (int i=0;i<npointsANITAII;i++)
    {
      ncols=fscanf(fp3,"%f,%f,%f",&freqTemp,&gainTemp,&nfTemp);
      freqANITAII[i]=freqTemp;
      if (gainTemp<1e9) gainANITAII[i]=gainTemp;
      else gainANITAII[i]=0;
    }
  
  
  TCanvas *cgain=new TCanvas("cgain","cgain",800,800);
  cgain->Divide(1,2);
  cgain->cd(1);
  TGraph *grgainANITAI=new TGraph(npointsANITAI, freqANITAI, gainANITAI);
  grgainANITAI->Draw("al");
  cgain->cd(2);
  TGraph *grgainANITAII=new TGraph(npointsANITAII, freqANITAII, gainANITAII);
  grgainANITAII->Draw("al"); 
  
  double magImpulseNew[newLengthImpulse], phaseImpulseNew[newLengthImpulse];
  
  int j, k;
  //now adjust theFFT for ratio of gains
  for (int i=0;i<newLengthImpulse;i++){
    for (j=0;j<npointsANITAII;j++){
      if (freqANITAII[j]/1e6>frequencyImpulse[i]) break;
    }
    for (k=0;k<npointsANITAI;k++){
      if (freqANITAI[k]/1e6>frequencyImpulse[i]) break;
    }
    //cout<<"frequencyImpulse[i]: "<<frequencyImpulse[i]<<"freqi: "<<freqANITAI[k]<<", freqii: "<<freqANITAII[j]<<endl;
    if (frequencyImpulse[i]>150 && frequencyImpulse[i]<1250){
      //cout<<"k: "<<k<<", j: "<<j<<", gain1: "<<gainANITAI[k]<<", gain 2: "<<gainANITAII[j]<<", gainII-gainI: "<<gainANITAII[j]-gainANITAI[k]<<endl;
      theFFTImpulse[i].re*=pow(10,(gainANITAII[j]-gainANITAI[k])/10.);
      theFFTImpulse[i].im*=pow(10,(gainANITAII[j]-gainANITAI[k])/10.);
      magImpulseNew[i]=sqrt(theFFTImpulse[i].re*theFFTImpulse[i].re+theFFTImpulse[i].im*theFFTImpulse[i].im);
      phaseImpulseNew[i]=atan(theFFTImpulse[i].im/theFFTImpulse[i].re);
    }
    else{
      magImpulseNew[i]=0;
      phaseImpulseNew[i]=0;
    }
  }
  
 double *thesystemresponseNew=FFTtools::doInvFFT(lengthImpulse,theFFTImpulse);
 TGraph *grSystem=new TGraph(lengthImpulse,oldXImpulse,thesystemresponseNew);

 TCanvas *newcan=new TCanvas("newcan","newcan",800,800);
 newcan->Divide(1,1);
 newcan->cd(1);
 grSystem->Draw("al");

 char filename[150];
 sprintf(filename,"systemImpulseResponseUsingStephens.txt");
 ofstream fsystemResponse(filename);
  for (int i=0;i<lengthImpulse/3;i++){
    fsystemResponse<<oldXImpulse[i]-10<<"\t"<<thesystemresponseNew[i]<<endl;
  }
  fsystemResponse.close();

 
}
/////////////////////////////
void MyCorrelator::calcSystemImpulseResponse()
{  
  //factor to upsample the data to, for interpolation
  int upSampleFactor=1;
  Double_t deltaTInt=1./(2.6*upSampleFactor);
  
  //open file with impulse input
  char filename[200];
  //sprintf(filename,"/home/dailey.110/abby/agoodhue/picoSecondImpulse/picoPulseThruTestStandSetupCh1.dat");//with 30 dB on picopulser
  sprintf(filename,"${ANITA_ANALYSIS_AGOODHUE}/picoSecondImpulse/picoPulseThruTestStandSetupCh1.dat");//with 30 dB on picopulser 
  cout<<filename<<endl;
  std::ifstream inputPulse_file(filename);
  double timeTemp, voltageTemp;
  int ctr=0;
  int npoints=1000;
  double time[npoints*3], voltage[npoints*3];
  
  //initialize to 0
  for (int i=0;i<npoints*3;i++){
    time[i]=0;
    voltage[i]=0;
  }
  
  //get input pulse from a file
  while(inputPulse_file >> timeTemp >> voltageTemp)
    {
      time[ctr]=timeTemp*1e9-80+0.1; //ns
      //if (time[ctr]<55 || time[ctr]>100) voltage[ctr]=0;
      if (ctr>500) voltage[ctr-500]=voltageTemp*1e3;//mV
      if (ctr==npoints-1) break;
      ctr++;
    }

  for (int i=0;i<npoints*3;i++){  
    if (i>=npoints-100) time[i]=time[i-1]-time[i-2]+time[i-1];
    if (time[i]<10 || time[i]>40) voltage[i]=0;
  }

  //draw input pulse
  TGraph *inputPulse=new TGraph(npoints*3,time,voltage);
  TCanvas *cImpulse2=new TCanvas("cimpulse2","cimpulse2",800,800);
  cImpulse2->cd(0);

  //interpolate the graph to the same delta T as the data will be
  TGraph *grInterpImpulse;
  grInterpImpulse=FFTtools::getInterpolatedGraph(inputPulse, deltaTInt);
  grInterpImpulse->Draw("al");
  
  //take FFT of the pulse
  double *oldY = grInterpImpulse->GetY();
  double *oldX = grInterpImpulse->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grInterpImpulse->GetN();
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  int newLength=(length/2)+1;
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz

  //define  variables for FFT of input pulse
  double frequencyInput[newLength];//MHz
  double magnitudeInput[newLength];
  double phaseInput[newLength];
  frequencyInput[0]=0;
  
  FFTWComplex *theFFTAvg=new FFTWComplex[(length/2)+1];
  
  //fill the magnitude and phase, and bring it from 30 dB to 48 dB (down by 18 dB in power)
  for(int i=0;i<newLength;i++) {
    if (i>0) frequencyInput[i]=frequencyInput[i-1]+deltaF;
    theFFT[i].re*=pow(10,-0.9);
    theFFT[i].im*=pow(10,-0.9);
    if (i>3 && i<newLength-4){
      theFFTAvg[i].re=(theFFT[i].re+theFFT[i-1].re+theFFT[i+1].re+theFFT[i+2].re+theFFT[i-2].re+theFFT[i+3].re+theFFT[i-3].re+theFFT[i+4].re+theFFT[i-4].re)/9;
      theFFTAvg[i].im=(theFFT[i].im+theFFT[i-1].im+theFFT[i+1].im+theFFT[i-2].im+theFFT[i+2].im+theFFT[i-3].im+theFFT[i+3].im+theFFT[i-4].im+theFFT[i+4].im)/9;
    }
    else{ 
      theFFTAvg[i].re=theFFT[i].re;
      theFFTAvg[i].im=theFFT[i].im;
    }
    magnitudeInput[i]=sqrt(theFFTAvg[i].re*theFFTAvg[i].re+theFFTAvg[i].im*theFFTAvg[i].im);
    phaseInput[i]=atan(theFFTAvg[i].im/theFFTAvg[i].re);
  }
  
  TGraph *grinputFFTMag=new TGraph(newLength,frequencyInput,magnitudeInput);
  TGraph *grinputFFTPhase=new TGraph(newLength,frequencyInput,phaseInput);
  TCanvas *cinputFFT=new TCanvas("cinputFFT","cinputFFT",800,800);
  cinputFFT->Divide(1,2);
  cinputFFT->cd(1);
  grinputFFTMag->Draw("al");
  cinputFFT->cd(2);
  grinputFFTPhase->Draw("al");

  
  //now start looking at the data from anita from that pulse
  int trigTime_this;
  int trigTime_prev;
  int trigTimeDelta;
  int numEventsIncluded=0;
  double timeThis, voltThis;
  double timeData[260*upSampleFactor*3], voltageData[260*upSampleFactor*3], voltageToDeconvolve[260*upSampleFactor*3];;
  
  for (int i=0;i<260*upSampleFactor*3;i++){
    timeData[i]=0;
    voltageData[i]=0;
    voltageToDeconvolve[i]=0;
  }
  
  //initialize data and loop over events
  if (!fEventTree) initialize();
  int nentries=fEventTree->GetEntries();
 
  for (int eventctr=0;eventctr<nentries;eventctr++){
    eventctr=10;
    fHeadTree->GetEntry(eventctr);
    if (eventStartedFlag!=(int)fHeadPtr->eventNumber) eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);

    //cut for events by time
    trigTime_prev=trigTime_this;
    trigTime_this=fHeadPtr->triggerTimeNs;
    trigTimeDelta=trigTime_this-trigTime_prev;
    if (trigTimeDelta<=0)trigTimeDelta+=int(1e9);
    
    //cut by time for the events
    if((fHeadPtr->triggerTimeNs%100000000)>56.980e6 && 
       (fHeadPtr->triggerTimeNs%100000000)<57.02e6 && trigTimeDelta>90000000){ 
      if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
      numEventsIncluded++;
      cout<<"eventctr: "<<eventctr<<endl;
      int ant=11;//12V
      //these are the events we are interested in, pulses.
      TGraph *gr1 = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
      TGraph *grInterp;
      grInterp=FFTtools::getInterpolatedGraph(gr1, deltaTInt);

      //fill the data array from the graph
      for (int i=0;i<260*3;i++){
	if (i<grInterp->GetN()){
	  grInterp->GetPoint(i,timeThis,voltThis);
	  if (numEventsIncluded==1) timeData[i]=timeThis-4.59;
	  if (timeData[i]<40) voltageData[i]=0;
	  else if (timeData[i]>100) voltageData[i]=0;
	  else voltageData[i]+=voltThis;
	  voltageToDeconvolve[i]+=voltThis;
	}
	else{
	  timeData[i]=(timeData[i-1]-timeData[i-2])+timeData[i-1];
	  voltageData[i]=0;
	  voltageToDeconvolve[i]=0;
	}
      }
       
      delete grInterp;
      delete gr1;
    }
    //just look at the first event
    cout<<"numEventsIncluded: "<<numEventsIncluded<<endl;
    if (numEventsIncluded==1) break;
  }

  //make the average voltage over all the events we did
  for (int i=0;i<260*upSampleFactor*3;i++){
    voltageData[i]=voltageData[i]/numEventsIncluded;
    voltageToDeconvolve[i]=voltageToDeconvolve[i]/numEventsIncluded;
  }

  //now draw this average waveform
  TGraph *grAverage=new TGraph(260*upSampleFactor*3,timeData,voltageData);
  TCanvas *cAverage=new TCanvas("cAverage","cAverage",800,800);
  cAverage->cd(0);
  grAverage->Draw("al");

  //take FFT of the data averaged pulse
  double *oldYData = grAverage->GetY();
  double *oldXData = grAverage->GetX();
  double deltaTData=oldXData[1]-oldXData[0];
  int lengthData=grAverage->GetN();
  FFTWComplex *theFFTData=FFTtools::doFFT(lengthData,oldYData);
  int newLengthData=(lengthData/2)+1;
  double deltaFData=1/(deltaTData*lengthData); //Hz
  deltaFData*=1e3; //MHz

  //define arrays for FFT of data
  double frequencyData[newLengthData];//MHz
  double magnitudeData[newLengthData];
  double phaseData[newLengthData];
  frequencyData[0]=0;
 
  FFTWComplex *theFFTDataAvg=new FFTWComplex[(lengthData/2)+1];
  
  //fill FFT of data
  for(int i=0;i<newLengthData;i++) {
    if (i>0) frequencyData[i]=frequencyData[i-1]+deltaFData;
    if (i>1 && i<newLengthData-2){
      theFFTDataAvg[i].re=(theFFTData[i].re+theFFTData[i-1].re+theFFTData[i+1].re+theFFTData[i+2].re+theFFTData[i-2].re)/5.;
      theFFTDataAvg[i].im=(theFFTData[i].im+theFFTData[i-1].im+theFFTData[i+1].im+theFFTData[i-2].im+theFFTData[i+2].im)/5.;
    }
    else{ 
      theFFTDataAvg[i].re=theFFTData[i].re;
      theFFTDataAvg[i].im=theFFTData[i].im;
    }
    magnitudeData[i]=sqrt(theFFTDataAvg[i].re*theFFTDataAvg[i].re+theFFTDataAvg[i].im*theFFTDataAvg[i].im);
    phaseData[i]=atan(theFFTDataAvg[i].im/theFFTDataAvg[i].re);
  }
  
  TGraph *grdataFFTMag=new TGraph(newLengthData,frequencyData,magnitudeData);
  TGraph *grdataFFTPhase=new TGraph(newLengthData,frequencyData,phaseData);
  TCanvas *cdataFFT=new TCanvas("cdataFFT","cdataFFT",800,800);
  cdataFFT->Divide(1,2);
  cdataFFT->cd(1);
  grdataFFTMag->Draw("al");
  cdataFFT->cd(2);
  grdataFFTPhase->Draw("al");
  
  //get the impulse response (grAverage - inputPulse)
  double magnitudeSystemResponse[newLengthData];
  double phaseSystemResponse[newLengthData];
  FFTWComplex *thesystemFFT2=new FFTWComplex[(lengthData/2)+1];
  cout<<"deltaT input (ns): "<<deltaT<<",deltaF input (MHZ): "<<deltaF<<", npoints new input: "<<newLength
      <<", npoints input: "<<length<<endl;
  cout<<"deltaT data (ns): "<<deltaTData<<",deltaF data (MHZ): "<<deltaFData<<", npoints new data: "<<newLengthData
      <<", npoints data: "<<lengthData<<endl;
  
  //calculate impulse response
  for (int i=0;i<newLengthData;i++){
    //    cout<<"frequency data: "<<frequencyData[i]<<endl;
    if (frequencyData[i]<200 || frequencyData[i]>1200) magnitudeSystemResponse[i]=0;
    else magnitudeSystemResponse[i]=magnitudeData[i]/magnitudeInput[i];
    phaseSystemResponse[i]=phaseData[i]-phaseInput[i];
    if (phaseSystemResponse[i]<-3.14159265/2) phaseSystemResponse[i]+=3.14159265;
    if (phaseSystemResponse[i]>3.14159265/2) phaseSystemResponse[i]-=3.14159265;
    if (frequencyData[i]>200 && frequencyData[i]<1200){
      thesystemFFT2[i].re=(1/pow(magnitudeInput[i],2))*(theFFTAvg[i].re*theFFTDataAvg[i].re+theFFTAvg[i].im*theFFTDataAvg[i].im);
      thesystemFFT2[i].im=(1/pow(magnitudeInput[i],2))*(theFFTAvg[i].re*theFFTDataAvg[i].im-theFFTAvg[i].im*theFFTDataAvg[i].re);
    }
    else{
      thesystemFFT2[i].re=0;
      thesystemFFT2[i].im=0;
    }
  }
  
  //add in the antenna response stephen gave me
  //open file with antenna response for vpol
  sprintf(filename,"${ANITA_ANALYSIS_AGOODHUE}/picoSecondImpulse/antenna_hh_impulse_response.dat");
  //cout<<filename<<endl;
  std::ifstream antennaResponse_file(filename);
  ctr=0;
  int npointsAntenna=1000;
  double timeAntenna[npointsAntenna*3], voltageAntenna[npointsAntenna*3];
  for (int i=0;i<npointsAntenna*3;i++){
    timeAntenna[i]=0;
    voltageAntenna[i]=0;
  }
  
  //get antenna response from  file
  while(antennaResponse_file >> timeTemp >> voltageTemp)
    {
      timeAntenna[ctr]=timeTemp*1e9+850-29.6; //ns
      voltageAntenna[ctr+100]=voltageTemp;//mV
      if (ctr==npointsAntenna-1) break;
      ctr++;
    }
  
  for (int i=0;i<npointsAntenna*3;i++){  
    if (i>=npointsAntenna-10) timeAntenna[i]=timeAntenna[i-1]-timeAntenna[i-2]+timeAntenna[i-1];
  }

  TGraph *grAntenna=new TGraph(npointsAntenna*3,timeAntenna,voltageAntenna);  
  //now draw this system response waveform
  TCanvas *cAntenna=new TCanvas("cAntenna","cAntenna",800,800);
  cAntenna->cd(0);
  TGraph *grInterpAntenna;
  grInterpAntenna=FFTtools::getInterpolatedGraph(grAntenna, deltaTInt);
  grInterpAntenna->Draw("al");
  
  double *voltageInterpAntenna=grInterpAntenna->GetY();
  double *timeInterpAntenna=grInterpAntenna->GetX();
  
  //take the fourier transform to add it on to the system response
  FFTWComplex *theFFTAntenna=FFTtools::doFFT(grInterpAntenna->GetN(),voltageInterpAntenna);
  double deltaTAntenna=timeInterpAntenna[1]-timeInterpAntenna[0];
  double deltaFAntenna=1/(deltaTAntenna*grInterpAntenna->GetN());
  int newLengthAntenna=(grInterpAntenna->GetN()/2)+1;
  deltaFAntenna*=1e3;
  double frequencyAntenna[newLengthAntenna];
  double magnitudeAntenna[newLengthAntenna];
  double phaseAntenna[newLengthAntenna];
  frequencyAntenna[0]=0;

  //fill antenna response arrays
  for (int i=0;i<newLengthAntenna;i++){
    frequencyAntenna[i]=frequencyAntenna[i-1]+deltaFAntenna;
    magnitudeAntenna[i]=sqrt(theFFTAntenna[i].re*theFFTAntenna[i].re+theFFTAntenna[i].im*theFFTAntenna[i].im);
    phaseAntenna[i]=atan(theFFTAntenna[i].im/theFFTAntenna[i].re);
  }
  cout<<"deltaT antenna: "<<deltaTAntenna<<", deltaF antenna: "<<deltaFAntenna
      <<", new length antenna: "<<newLengthAntenna<<", length antenna: "<<grInterpAntenna->GetN()<<endl;
   
  //antenna fft graphs
  TGraph *grantennaFFTMag=new TGraph(newLengthData,frequencyData,magnitudeAntenna);
  TGraph *grantennaFFTPhase=new TGraph(newLengthData,frequencyData,phaseAntenna);
  TCanvas *cantennaFFT=new TCanvas("cantennaFFT","cantennaFFT",800,800);
  cantennaFFT->Divide(1,2);
  cantennaFFT->cd(1);
  grantennaFFTMag->Draw("al");
  cantennaFFT->cd(2);
  grantennaFFTPhase->Draw("al");
  
  FFTWComplex *thesystemAntennaFFT2=new FFTWComplex[(lengthData/2)+1];
  //add the antenna to the system response
  for (int i=0;i<newLengthData;i++){
    magnitudeSystemResponse[i]=magnitudeSystemResponse[i];//ABBY *magnitudeAntenna[i];
    phaseSystemResponse[i]=phaseSystemResponse[i]; // ABBY +phaseAntenna[i];
    if (phaseSystemResponse[i]<-3.14159265/2) phaseSystemResponse[i]+=3.14159265;
    if (phaseSystemResponse[i]>3.14159265/2) phaseSystemResponse[i]-=3.14159265;
    thesystemAntennaFFT2[i].re=(thesystemFFT2[i].re*theFFTAntenna[i].re-thesystemFFT2[i].im*theFFTAntenna[i].im);
    thesystemAntennaFFT2[i].im=(thesystemFFT2[i].im*theFFTAntenna[i].re+thesystemFFT2[i].re*theFFTAntenna[i].im);
  }
  FFTWComplex *thesystemAntennaFFT2Avg=new FFTWComplex[(lengthData/2)+1];
  for (int i=0;i<newLengthData;i++){
    if (i>1 && i<newLengthData-2){
      thesystemAntennaFFT2Avg[i].re=(thesystemAntennaFFT2[i].re+thesystemAntennaFFT2[i-1].re+thesystemAntennaFFT2[i+1].re+
				     thesystemAntennaFFT2[i-2].re+thesystemAntennaFFT2[i+2].re)/5.;
      thesystemAntennaFFT2Avg[i].im=(thesystemAntennaFFT2[i].im+thesystemAntennaFFT2[i-1].im+thesystemAntennaFFT2[i+1].im+
				     thesystemAntennaFFT2[i-2].im+thesystemAntennaFFT2[i+2].im)/5.;
    }
    else {
      thesystemAntennaFFT2Avg[i].re=thesystemAntennaFFT2[i].re;
      thesystemAntennaFFT2Avg[i].im=thesystemAntennaFFT2[i].im;
    }
  }
  
  //graph the system response in fourier space
  TGraph *grsystemFFTMag=new TGraph(newLengthData,frequencyData,magnitudeSystemResponse);
  TGraph *grsystemFFTPhase=new TGraph(newLengthData,frequencyData,phaseSystemResponse);
  TCanvas *csystemFFT=new TCanvas("csystemFFT","csystemFFT",800,800);
  csystemFFT->Divide(1,2);
  csystemFFT->cd(1);
  grsystemFFTMag->Draw("al");
  csystemFFT->cd(2);
  grsystemFFTPhase->Draw("al");

  double *thesystemWaveform2=FFTtools::doInvFFT(lengthData,thesystemFFT2);//thesystemAntennaFFT2Avg ABBY
  for (int i=1;i<lengthData-1;i++) thesystemWaveform2[i]=thesystemWaveform2[i-1]+thesystemWaveform2[i+1]+thesystemWaveform2[i]/3.;
  for (int i=0;i<lengthData;i++) if (oldXData[i]<28 || oldXData[i]>40) thesystemWaveform2[i]=0; //40-60, 30-45

  double theNewSystemWaveform[lengthData/3];
  for (int i=0;i<lengthData/3;i++) theNewSystemWaveform[i]=thesystemWaveform2[i];//used to be 140

  TGraph *grSystem=new TGraph(lengthData/3,oldXData,theNewSystemWaveform);
  //now draw this system response waveform
  TCanvas *cSystem=new TCanvas("cSystem","cSystem",800,800);
  cSystem->cd(0);
  grSystem->Draw("al");
  
  //print out the system response for future use
  sprintf(filename,"systemImpulseResponse.txt");
  ofstream fsystemResponse(filename);
  for (int i=0;i<lengthData/3;i++){
    fsystemResponse<<oldXData[i]<<"\t"<<theNewSystemWaveform[i]<<endl;
  }
  fsystemResponse.close();

  
  //now test the deconvolution
  TGraph *grToDecon=new TGraph(int(grAverage->GetN()/3),timeData,voltageToDeconvolve);
  //  TGraph *grToDecon=new TGraph(int(grAverage->GetN()/3),grAverage->GetX(),grAverage->GetY());
  TGraph *grdeconvolved;
  grdeconvolved=deconvolveWaveformUsingStephens(grToDecon,0,1);
  

}
//////////////////////////////
void MyCorrelator::loopOverStartingEvents()
{
  if (!fEventTree) initialize();
  int nentries=fEventTree->GetEntries();

  for (int eventctr=0;eventctr<nentries;eventctr++){
    //get the event
    if (eventctr%(nentries/10)==0) 
      cout<<eventctr<<" neutrinos.  "<<((double)eventctr/(double)nentries)*100<< "% complete."<<endl;
    //get the event
    fHeadTree->GetEntry(eventctr);
    eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
  }


}
///////////////////////////
void MyCorrelator::calcGPSFlags()
{
  if (!fEventTree) initialize();
  int nentries=fEventTree->GetEntries();
  int eventNumber;
  int ngpsBadFlag=0;
  char filename[150];
  sprintf(filename,"rootOutputs/outputGPS%d.root",
	  fCurrentRun);
  if (printFlag==1) cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");
  TTree *tgpsout=new TTree("tgpsout","gpsbadflag");
  tgpsout->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tgpsout->Branch("gpsBadFlag",&gpsBadFlag,"gpsBadFlag/I");
  

  for (int eventctr=0;eventctr<nentries;eventctr++){
    //get the event
    //if (eventctr%(nentries/1000)==0) 
    // cout<<eventctr<<" neutrinos.  "<<((double)eventctr/(double)nentries)*100<< "% complete."<<endl;
    //get the event
    fHeadTree->GetEntry(eventctr);
    eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
    eventNumber=fHeadPtr->eventNumber;
    tgpsout->Fill();
    if (gpsBadFlag==1) ngpsBadFlag++;
  }
  cout<<"bad gps total: "<<ngpsBadFlag<<endl;
  tgpsout->Write();
  rootfile->Close("R");
}
void MyCorrelator::firstLastEachRun(int &firstEvent, int &lastEvent, int &numberOfEvents)
{
  int eventCtrStart=0;
  int eventCtrEnd=1000000;
  if (!fEventTree) initialize();
  if (eventCtrEnd>fEventTree->GetEntries()) eventCtrEnd=fEventTree->GetEntries();
  if (eventCtrStart>fEventTree->GetEntries()) eventCtrStart=fEventTree->GetEntries();
  fHeadTree->GetEntry(eventCtrStart);
  fEventTree->GetEvent(eventCtrStart);
  firstEvent=fHeadPtr->eventNumber;
  lastEvent=fHeadPtr->eventNumber+eventCtrEnd-eventCtrStart;
  numberOfEvents=lastEvent-firstEvent;
  if (printFlag==1) cout<<"first Event Number: "<<firstEvent<<", last Event Number: "<<lastEvent
			<<", number of Events: "<<lastEvent-firstEvent<<endl;

}
/////////////////////
void MyCorrelator::calcVarnerFlags()
{
  if (!fEventTree) initialize();
  int nentries=fEventTree->GetEntries();
  int varnerFlag;
  int varnerFlag2;
  int eventNumber;
  int varnerFlagCtr=0;
  int varnerFlag2Ctr=0;
  
  char filename[150];
  sprintf(filename,"rootOutputs/outputVarner%d.root",
	  fCurrentRun);
  if (printFlag==1) cout<<"outputting to file: "<<filename<<endl;
  //TFile *rootfile = new TFile(filename,"RECREATE");
  TTree *tvarner=new TTree("tvarner","varnerFlag");
  tvarner->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tvarner->Branch("varnerFlag",&varnerFlag,"varnerFlag/I");
  tvarner->Branch("varnerFlag2",&varnerFlag2,"varnerFlag2/I");
  
  for (int eventctr=0;eventctr<nentries;eventctr++){
    //get the event
    if (eventctr%(nentries/100)==0) 
      cout<<eventctr<<" neutrinos.  "<<((double)eventctr/(double)nentries)*100<< "% complete."<<endl;
    //get the event
    fHeadTree->GetEntry(eventctr);
    eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
    eventEntryGottenFlag=getEventEntry();
    varnerFlag=isVarnerEvent(fHeadPtr->eventNumber);
    varnerFlag2=isVarnerEvent2(fHeadPtr->eventNumber);
    eventNumber=fHeadPtr->eventNumber;
    tvarner->Fill();
    if (varnerFlag2!=varnerFlag) cout<<"event number: "<<eventNumber<<"varnerFlagOld: "
				      <<varnerFlag<<", varnerFlagNew: "<<varnerFlag2<<endl;
    if (varnerFlag==1) varnerFlagCtr++;
    if (varnerFlag2==1) varnerFlag2Ctr++;
    
  }
  cout<<"varner old total: "<<varnerFlagCtr<<", varner new total: "<<varnerFlag2Ctr<<endl;
  tvarner->Write();

}
////////////////////////////////////////////////////////////////////////////////////
void MyCorrelator::taylorPhiMaskandTriggerBits(int numEvents, 
					       double isbesttrigpm1_denom[NUM_PHI_SECTORS], 
					       double isbesttrig_denom[NUM_PHI_SECTORS],
					       double isbesttrigpm1orphimask_denom[NUM_PHI_SECTORS], 
					       double isbesttrigorphimask_denom[NUM_PHI_SECTORS],
					       double isbesttrigpm1_num[NUM_PHI_SECTORS],
					       double isbesttrig_num[NUM_PHI_SECTORS],
					       double isbesttrigpm1orphimask_num[NUM_PHI_SECTORS], 
					       double isbesttrigorphimask_num[NUM_PHI_SECTORS])
{
  
  if (!fEventTree){
    initialize();
  }
  
  if (numEvents==0)
    numEvents=fEventTree->GetEntries();

  int phiMaskArray[NUM_ANTS_WITH_NADIRS];
  int triggeredPhi[NUM_ANTS_WITH_NADIRS];
  int thisPhiMask;
  float eventNumber;
  int i=0;
  float phiMap;
  float deltaPhi;

  // int nantennas=int(NUM_ANTS_WITH_NADIRS);
  
  TProfile *pIsBestTriggered=new TProfile("pIsBestTriggered","If not phi masked, Is Best Triggered?",16,0,16,-1,2);
  TProfile *pIsBestTriggeredPlusMinusOne=new TProfile("pIsBestTriggeredPlusMinusOne","If not phi masked, Is Best Triggered +/- 1 ?",16,0,16,-1,2);
  TProfile *pIsBestTriggeredOrPhiMasked=new TProfile("pIsBestTriggeredOrPhiMasked","Is Best Triggered Or Phi Masked?",16,0,16,-1,2);
  TProfile *pIsBestTriggeredOrPhiMaskedPlusMinusOne=new TProfile("pIsBestTriggeredOrPhiMaskedPlusMinusOne","Is Best Triggered +/- 1  or Phi Masked?",16,0,16,-1,2);
 
  char filename[150];

  sprintf(filename,"rootOutputs/051309/output%d.root",fCurrentRun); 
  TFile *rootfile=new TFile(filename,"READ");
  TTree *ndata=(TTree*)rootfile->Get("ndata");
  cout<<"number of entries: "<<ndata->GetEntries()<<endl;
  ndata->SetBranchAddress("eventNumber",&eventNumber);
  ndata->SetBranchAddress("phiMap",&phiMap);
  ndata->SetBranchAddress("deltaPhi",&deltaPhi);

  for (int eventctr=0;eventctr<numEvents;eventctr++){
    //get the event
    fHeadTree->GetEntry(eventctr);
    eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
    int isTaylorThis=isTaylor(fHeadPtr->eventNumber);
    
    if (isTaylorThis==1){
      ndata->GetEntry(i);
      i++;
      getTriggeredL2Phi(fHeadPtr,triggeredPhi);
      thisPhiMask=isPhiMaskingOn(fHeadPtr->eventNumber,phiMaskArray); 
      int bestPhiSector=getBestPhiSector(phiMap*deg2rad);
      int bestPhiSectorLeft=bestPhiSector-1;
      int bestPhiSectorRight=bestPhiSector+1;
      if (bestPhiSectorLeft<0) bestPhiSectorLeft+=16;
      if (bestPhiSectorRight>15) bestPhiSectorLeft-=16;
      
      if (fabs(deltaPhi)<7){
	if (triggeredPhi[bestPhiSector]==0 && //phiMaskArray[bestPhiSectorRight]==0 && phiMaskArray[bestPhiSectorLeft]==0 &&
	    phiMaskArray[bestPhiSector]==0 && phiMaskArray[bestPhiSectorLeft]==0 && phiMaskArray[bestPhiSectorRight]==0 &&
	    triggeredPhi[bestPhiSectorRight]==0 && triggeredPhi[bestPhiSectorLeft]==0){
	  pIsBestTriggeredOrPhiMaskedPlusMinusOne->Fill(bestPhiSector,0,1);
	  isbesttrigpm1orphimask_denom[bestPhiSector]++;
	  
	  cout<<"Best Phi Sector is untriggered and unmasked for event: "<<fHeadPtr->eventNumber
	      <<", best phi: "<<bestPhiSector<<endl;
	  for (int phictr=0;phictr<NUM_PHI_SECTORS;phictr++){
	    cout<<"triggerbit phi: "<<phictr<<" = "<<triggeredPhi[phictr]<<endl;
	    //cout<<phiMaskArray[phictr]<<endl;
	  }
	}
	else{ 
	  pIsBestTriggeredOrPhiMaskedPlusMinusOne->Fill(bestPhiSector,1,1);
	  isbesttrigpm1orphimask_denom[bestPhiSector]++;
	  isbesttrigpm1orphimask_num[bestPhiSector]++;
	}

	if (triggeredPhi[bestPhiSector]==0 && phiMaskArray[bestPhiSector]==0){
	  pIsBestTriggeredOrPhiMasked->Fill(bestPhiSector,0,1);
	  isbesttrigorphimask_denom[bestPhiSector]++;
	}
	else{
	  pIsBestTriggeredOrPhiMasked->Fill(bestPhiSector,1,1);
	  isbesttrigorphimask_denom[bestPhiSector]++;
	  isbesttrigorphimask_num[bestPhiSector]++;
	}
	
	if (phiMaskArray[bestPhiSector]==0){
	  if (triggeredPhi[bestPhiSector]==0){
	    pIsBestTriggered->Fill(bestPhiSector,0,1);
	    isbesttrig_denom[bestPhiSector]++;
	  }
	  else{
	    pIsBestTriggered->Fill(bestPhiSector,1,1);
	    isbesttrig_denom[bestPhiSector]++;
	    isbesttrig_num[bestPhiSector]++;
	  }	  
	  if (triggeredPhi[bestPhiSector]==0 && triggeredPhi[bestPhiSectorRight]==0 && triggeredPhi[bestPhiSectorLeft]==0){
	    pIsBestTriggeredPlusMinusOne->Fill(bestPhiSector,0,1);
	    isbesttrigpm1_denom[bestPhiSector]++;
	  }
	  else{
	    pIsBestTriggeredPlusMinusOne->Fill(bestPhiSector,1,1);
	    isbesttrigpm1_denom[bestPhiSector]++;
	    isbesttrigpm1_num[bestPhiSector]++;
	  }
	}      
	
      }
      // else cout<<"Misreconstructed event: "<<fHeadPtr->eventNumber<<endl;      
    }
  }
  
  /* TCanvas *cBestTriggered=new TCanvas("cBestTriggered","cBestTriggered",800,800);
  cBestTriggered->Divide(1,2);
  cBestTriggered->cd(1);
  pIsBestTriggered->Draw();
  cBestTriggered->cd(2);
  pIsBestTriggeredPlusMinusOne->Draw();

  TCanvas *cBestTriggeredOrPhiMasked=new TCanvas("cBestTriggeredOrPhiMasked","cBestTriggeredOrPhiMasked",800,800);
  cBestTriggeredOrPhiMasked->Divide(1,2);
  cBestTriggeredOrPhiMasked->cd(1);
  pIsBestTriggeredOrPhiMasked->Draw();
  cBestTriggeredOrPhiMasked->cd(2);
  pIsBestTriggeredOrPhiMaskedPlusMinusOne->Draw();
  */

  //  tBaseline->Print();
  // rootfile->cd();
  //tBaseline->Write();
  
}
/////////////////////////////
void MyCorrelator::taylorEfficiency(int numEvents)
{
  if (!fEventTree){
    initialize();
    loadTurfTree();
  }
  if (numEvents==0)
    numEvents=fEventTree->GetEntries();
  int npoints=100000;
  double distance;
  int i=0;
  int phiMaskArray[NUM_PHI_SECTORS];
  int thisPhiMask=0;
  int prevPhiMask=0;
  int timeOfMaskedInterval[npoints];
  int beginIntervalMasked=0;
  int intervalctrMasked=0;
  int realTimeThis=0;
  int realTimePrev=0;
  int numSeconds=0;
  int triggeredPhi[NUM_PHI_SECTORS];
  int triggeredAnt[NUM_ANTS_WITH_NADIRS];
  int peakAntenna;

  double snr;
  double rmsNoise;
  double maxPeakVal;
  double deadtimeFraction;

  TTimeStamp t1;
  TTimeStamp t2;
 
  //float DecOne2008 = TTimeStamp(2008,12,01,0,0,0).GetSec();
 
  char filename[150];
  //sprintf(filename,"outputTaylorEfficiencyRun%d.root",fCurrentRun);
  sprintf(filename,"outputTaylorEfficiency/run%d.root",fCurrentRun);
  if (printFlag==1) cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");
  TNtuple *ndata = new TNtuple("ndata","stuff to plot",
			       "eventNumber:isTaylor:Day:Second:isPhiMask:distance:snr:peakAntenna:deadTimeFraction");
  //"eventNumber:isTaylor:realTime:isPhiMask:distance:snr:peakAntenna:deadTimeFraction");
  
  for (int eventctr=0;eventctr<numEvents;eventctr++){
    if (eventctr%(numEvents/10)==0) 
   cout<<eventctr<<" neutrinos.  "<<((double)eventctr/(double)numEvents)*100<< "% complete."<<endl;
    //get the event

    fHeadTree->GetEntry(eventctr);
    eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);        
    int turfEntry = fTurfRateTree->GetEntryNumberWithBestIndex(fHeadPtr->realTime);
    fTurfRateTree->GetEntry(turfEntry);
    


    prevPhiMask=thisPhiMask;  
    thisPhiMask=isPhiMaskingOn(fHeadPtr->eventNumber,phiMaskArray);
    realTimePrev=realTimeThis;
    realTimeThis=fHeadPtr->realTime;
    if ((prevPhiMask==0 && thisPhiMask==1) || (eventctr==0 && thisPhiMask==1)) beginIntervalMasked=fHeadPtr->realTime;
    if ((prevPhiMask==1 && thisPhiMask==0 && eventctr!=0) || (eventctr==numEvents-1 && thisPhiMask==1)){
      timeOfMaskedInterval[intervalctrMasked]=fHeadPtr->realTime-beginIntervalMasked;
      intervalctrMasked++;
    }

    int isTaylorThis=isTaylor(fHeadPtr->eventNumber,distance);
    int realTimeDay = (int)((float)fHeadPtr->realTime/86400.);
    float realTimeSec = ((float)fHeadPtr->realTime/86400.-realTimeDay)*86400;
    //if(realTimeDay*86400+realTimeSec-fHeadPtr->realTime!=0)cout<<setprecision(12)<<realTimeDay<<" "<<realTimeSec<<" "<<realTimeMin*60+realTimeSec-fHeadPtr->realTime<<setprecision(6)<<endl;
    //cout<<setprecision(12)<<realTimeDay<<" "<<realTimeSec<<" "<<realTimeSec+realTimeDay*86400-fHeadPtr->realTime<<setprecision(6)<<endl;

    if (realTimePrev==realTimeThis-1){//the next second happened
      deadtimeFraction=fTurfPtr->deadTime/65535.;
      //ndata->Fill(fHeadPtr->eventNumber,0,fHeadPtr->realTime,
      ndata->Fill(fHeadPtr->eventNumber,0,realTimeDay,realTimeSec,
		  isPhiMaskingOn(fHeadPtr->eventNumber,phiMaskArray),-1,-1,-1,deadtimeFraction);
      numSeconds++;

      cout<<fHeadPtr->eventNumber<<" "<<setprecision(12)<<fHeadPtr->realTime<<"\t"<<setprecision(6)<<t1.AsString()<<endl;      
    }
    
    if (isTaylorThis==1){
      maxPeakVal=0;
      getEventEntry();
      getTriggeredL2Phi(fHeadPtr,triggeredPhi);
      getTriggeredAnt(triggeredPhi,triggeredAnt);
      for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
	if (triggeredAnt[ant]==1){
	  TGraph *gr = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);	  	
	  Double_t peakTime, peakVal;
	  Double_t peakBin=FFTtools::getPeakBin(gr);
	  gr->GetPoint(int(peakBin),peakTime,peakVal);
	  if (peakVal>maxPeakVal){
	    maxPeakVal=peakVal;
	    peakAntenna=ant;
	    snr=getSNR(gr,rmsNoise);
	  }
	  delete gr;
	}
      }
      
      ndata->Fill(fHeadPtr->eventNumber,isTaylorThis,realTimeDay,realTimeSec,//fHeadPtr->realTime,
		  //ndata->Fill(fHeadPtr->eventNumber,isTaylorThis,,//(float)fHeadPtr->realTime-DecOne2008,
		  isPhiMaskingOn(fHeadPtr->eventNumber,phiMaskArray),distance,snr,peakAntenna,-1);

       cout<<fHeadPtr->eventNumber<<" "<<setprecision(12)<<fHeadPtr->realTime<<"\t"<<setprecision(6)<<t1.AsString()<<endl;
      i++;
    }
    t1 = fHeadPtr->realTime;
    // cout<<fHeadPtr->eventNumber<<" "<<setprecision(12)<<fHeadPtr->realTime<<"\t"<<setprecision(6)<<t1.AsString()<<endl;
  }
  
  
  int totalTimeMasked=0;
  for (int ctr=0;ctr<intervalctrMasked;ctr++){
    totalTimeMasked+=timeOfMaskedInterval[ctr];
  }

  double totalEfficiency=double(i)/double(numSeconds-totalTimeMasked);
  cout<<"number of pulses seen: "<<i<<", totalTimeElapsed: "<<numSeconds-totalTimeMasked
      <<", totalEfficiency: "<<totalEfficiency<<endl;
  // cout<<"total time masked: "<<totalTimeMasked<<endl;
  // cout<<"numseconds: "<<numSeconds<<endl;

  rootfile->Write("ndata");
  fEventTree->Delete();  
  fTurfRateTree->Delete();
  ndata->Delete();

  }
////////////////////////////////////////////
void MyCorrelator::trackTaylor(int numEvents){
  if (!fEventTree){ initialize();
  }
  if (numEvents==0)
    numEvents=fEventTree->GetEntries();
  double distance;
  double distance_TD;
  double maxPeakVal, maxAdjustedPeakVal;
  double fresnelCoefficient;
  //int triggeredPhi[NUM_PHI_SECTORS];
  //int triggeredAnt[NUM_ANTS_WITH_NADIRS];
  int peakAntenna;
  int i=0;
  int eventNumberArray;
  double snr;
  double rmsNoise;

  char filename[150];
  sprintf(filename,"outputtrackTaylorRun%d.root",fCurrentRun);
  if (printFlag==1) cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");
  TNtuple *ndata = new TNtuple("ndata","stuff to plot",
			       "eventNumber:distance:peakVal:peakAdjustedVal:snr:rms:fresnel");
  //get taylor dome location
  double x2,y2,z2; 
  getTDxyz(x2,y2,z2); 
  double deltaTInt=1/2.6;

  for (int eventctr=0;eventctr<numEvents;eventctr++){
    if (eventctr%(numEvents/10)==0) 
      cout<<eventctr<<" neutrinos.  "<<((double)eventctr/(double)numEvents)*100<< "% complete."<<endl;
    //get the event
    snr=0;
    maxPeakVal=0;
    maxAdjustedPeakVal=0;
    rmsNoise=0;
    fHeadTree->GetEntry(eventctr);
    if (eventStartedFlag!=(int)fHeadPtr->eventNumber) eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
    if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber){ 
      eventEntryGottenFlag=getEventEntry();
    }
    if (isTaylor(fHeadPtr->eventNumber,distance)==1 && isChannelSaturated(fHeadPtr->eventNumber)==-1){
      // cout<<"Taylor Dome event: "<<fHeadPtr->eventNumber<<", distance to payload: "<<distance<<endl;
      //now calculate the amplitude*r vs. r
      distance_TD=distance;
      eventNumberArray=fHeadPtr->eventNumber;
      // getTriggeredPhi(fHeadPtr,triggeredPhi);
      // getTriggeredAnt(triggeredPhi,triggeredAnt);
      for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
	//	if (triggeredAnt[ant]==1){
	  TGraph *gr = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);
	  TGraph *grInterp;
	  grInterp=FFTtools::getInterpolatedGraph(gr, deltaTInt);
	  TGraph *grNotch=simpleNotchFilter(grInterp,380,500);
	  TGraph *grNotch2=simpleNotchFilter(grNotch,235,285);
	  
	  TGraph *grFiltered;
	  grFiltered=FFTtools::simplePassBandFilter(grNotch2,200,1100);	
	  
	  Double_t peak2peak, peakToPeakOverTwo;

	  peak2peak=getPeak2Peak(grFiltered);
	  peakToPeakOverTwo=peak2peak/2.;
	  
	  if (peakToPeakOverTwo>maxPeakVal){
	    snr=getSNR(grFiltered,rmsNoise);
	    maxPeakVal=peakToPeakOverTwo;
	    maxAdjustedPeakVal=peakToPeakOverTwo-rmsNoise;
	    peakAntenna=ant;
	  }
	  delete gr;
	  delete grNotch;
	  delete grNotch2;
	  delete grInterp;
	  delete grFiltered;
	  //	}
      }
      
      //cout <<"i "<<i<< ", eventNumber: " << fHeadPtr->eventNumber << " distance: " << distance
      //   << "  Amplitude: " << maxPeakVal[i] << endl;

      //now calculate the angle that the ray left the surface of the ice (assume flat surface)     
      double thetaWave, phiWave;
  
      double x1,y1,z1;
      LatLonAlt2xyz(fAdu5APatPtr->latitude*(-1.),fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,x1,y1,z1);
  
      double angleFromEarthCenter=0.;
      //angleFromEarthCenter=fabs(fAdu5APatPtr->latitude*(-1.)-latTaylor)+fabs(ninetyMinusLatTD*(lonTaylor-fAdu5APatPtr->longitude)/90.); //fix this line
      
      fUsefulAdu5Ptr->getThetaAndPhiWave(lonTaylor, latTaylor, heightTaylor, thetaWave, phiWave);
      //angleFromEarthCenter=fUsefulAdu5Ptr->getAngleBetweenPayloadAndSource(lonTaylor, latTaylor, heightTaylor);
      
      double angleRayLeftIce=angleFromEarthCenter*rad2deg-thetaWave*rad2deg+90;
      //cout<<"thetaWave: "<<thetaWave*rad2deg<<", angleFromEarthCenter: "<<angleFromEarthCenter*rad2deg
      //  <<", angleRayLeftIce: "<<angleRayLeftIce<<endl;
      
      //do the fresnel coefficients
      float n_ice=1.325;//of firn from icemc
      float transmitted_angle=angleRayLeftIce*deg2rad;//angle between ray and ice normal in ice
      float incident_angle=asin(sin(transmitted_angle)/n_ice);//snelling: angle between ray and air normal in air
      float magnification=sqrt(tan(incident_angle)/tan(transmitted_angle));//magnification factor from icemc
      float fresnel_pokey=(2*cos(incident_angle))/((1./n_ice)*cos(incident_angle)+cos(transmitted_angle));
      float fresnel_factor=magnification*fresnel_pokey;//assume pokey only (totally vertically polarized signal)
      fresnelCoefficient=fresnel_factor;
      // cout<<"transmitted: "<<transmitted_angle<<", incident: "<<incident_angle<<", magnification: "
      // <<magnification<<", fresnel_pokey: "<<fresnel_pokey<<", fresnel coefficient: "<<fresnel_factor<<endl;      

      // cout<<"event number: "<<eventNumberArray<<", distance: "<<distance_TD<<", maxPeakVal: "<<maxPeakVal<<", maxAdjustedPeak: "
      //  <<maxAdjustedPeakVal<<", snr: "<<snr<<", rmsnoise: "<<rmsNoise<<", fresnel: "<<fresnelCoefficient<<endl;
      ndata->Fill(eventNumberArray,distance_TD,maxPeakVal,maxAdjustedPeakVal,snr,rmsNoise,fresnelCoefficient);
      i++;
    }
  } 
  
  rootfile->Write("ndata");
  
}
//////////////////END STUFF FOR POINTING TO TD
///////////////////// A Function for finding TD events from telemetered header files.
void MyCorrelator::trackTaylorForEventInsertion(int numEvents)
{
  char headerName[FILENAME_MAX];  
  sprintf(headerName,"%s/run%d/headFile%d.root",fCurrentBaseDir,fCurrentRun,fCurrentRun);

  fHeadFile = TFile::Open(headerName);
  if(!fHeadFile) {
    cout << "Couldn't open: " << headerName << "\n";
  }
  fHeadTree = (TTree*) fHeadFile->Get("headTree");
 
  if(!fHeadTree) {
    cout << "Couldn't get headTree from " << headerName << endl;

  }

  
  fHeadFile = TFile::Open(headerName);
  if(!fHeadFile) {
    cout << "Couldn't open: " << headerName << "\n";
  }
  fHeadTree = (TTree*) fHeadFile->Get("headTree");
  if(!fHeadTree) {
    cout << "Couldn't get headTree from " << headerName << endl;
    }
 
  fHeadIndex = (TTreeIndex*) fHeadTree->GetTreeIndex();
  std::cerr << fHeadTree->GetEntries() <<"\n";
  
  //got header.

 if (fAdu5APatPtr) delete fAdu5APatPtr;
  if (fAdu5aPatTree) delete fAdu5aPatTree;
  
  char gpsName[FILENAME_MAX];
  sprintf(gpsName,"%s/run%d/gpsFile%d.root",fCurrentBaseDir,fCurrentRun,
	  fCurrentRun);
  fGpsFile = TFile::Open(gpsName);
  if(!fGpsFile) {
    cout << "Couldn't open: " << gpsName << "\n";
  }
 
  fAdu5aPatTree = (TTree*) fGpsFile->Get("adu5PatTree");
 
  if(!fAdu5aPatTree) {
    cout << "Couldn't get adu5aPatTree\n";
  }
  else {
    fAdu5aPatTree->SetBranchAddress("pat",&fAdu5APatPtr);
  }
  fAdu5aPatEntry=0;
  fAdu5aPatTree->BuildIndex("realTime");
  
  //got gps.

  if (numEvents==0)
    numEvents=fHeadTree->GetEntries();
  cout<<"Number of events to process: "<<numEvents<<endl;
  double distance_TD;
  int eventNumberArray;
  int runNumberArray;
  int i=0;
  int runNumber;
  int runNumber_prev;
  
  char filename[150];
  sprintf(filename,"outputtrackTaylorForEventInsertionRun.root");
  if (printFlag==1) cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");
  TNtuple *ndata = new TNtuple("ndata","stuff to plot",
			       "eventNumber:distance:run");

  //get taylor dome location
  double x2,y2,z2; 
  getTDxyz(x2,y2,z2); 
  
  for (int eventctr=0;eventctr<numEvents;eventctr++){
    if (eventctr%(numEvents/10)==0) 
      cout<<eventctr<<" neutrinos.  "<<((double)eventctr/(double)numEvents)*100<< "% complete."<<endl;
    //get the event
    fHeadTree->GetEntry(eventctr);
    runNumber_prev=runNumber;
    runNumber=fHeadPtr->run;
    //if (runNumber_prev!=runNumber) cout<<"Run: "<<runNumber<<", event: "<<fHeadPtr->eventNumber<<endl;;

    int patEntry = fAdu5aPatTree->GetEntryNumberWithBestIndex(fHeadPtr->realTime);
    fAdu5aPatTree->GetEntry(patEntry);
    
    double x1,y1,z1;
    LatLonAlt2xyz(fAdu5APatPtr->latitude*(-1.),fAdu5APatPtr->longitude,fAdu5APatPtr->altitude,x1,y1,z1);
    double distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    double timeOfFlight=distance/(C_light);//in nanoseconds
    double delay=timeOfFlight-40e3;     
    int trig_time=fHeadPtr->triggerTimeNs;

    if (fabs(trig_time-delay)<500){  
      //if (printFlag==1) cout<<"Event "<<fHeadPtr->eventNumber<<" is a Taylor Dome Event."<<endl;
      distance_TD=distance;
      eventNumberArray=fHeadPtr->eventNumber;
      runNumberArray=fHeadPtr->run;
      i++;
      ndata->Fill(eventNumberArray,distance_TD,runNumberArray);
      
    }
  }
   
  rootfile->Write("ndata");
  
}

///////////////////
void MyCorrelator::createBaseline(int numEvents)
{
  if (!fEventTree) initialize();
  if (numEvents==0)
    numEvents=fEventTree->GetEntries();
  
  int ant;
  int nPointsForTree;

  int upSampleFactor=1;
  Double_t deltaTInt=1./(2.6*upSampleFactor);
  float frequencyArray[260*upSampleFactor];
 
  double vertFFTre[NUM_ANTS_WITH_NADIRS][260*upSampleFactor];
  double horizFFTre[NUM_ANTS_WITH_NADIRS][260*upSampleFactor];
  double vertFFTim[NUM_ANTS_WITH_NADIRS][260*upSampleFactor]; 
  double horizFFTim[NUM_ANTS_WITH_NADIRS][260*upSampleFactor];
  float vertFFT[260*upSampleFactor];
  float horizFFT[260*upSampleFactor];
  int eventsIncludedCtr=0;
  int newLength;
  
  char filename[150];
  sprintf(filename,"outputcreateBaselineRun%d_test.root",fCurrentRun);
  if (printFlag==1) cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");
  TTree *tBaseline = new TTree("tBaseline","stuff for baseline");
  tBaseline->Branch("ant",&ant,"ant/I");
  tBaseline->Branch("nPointsForTree",&nPointsForTree,"nPointsForTree/I");
  tBaseline->Branch("frequencyArray",
		    frequencyArray,"frequencyArray[nPointsForTree]/F");
  tBaseline->Branch("vertFFT",vertFFT,"vertFFT[nPointsForTree]/F");
  tBaseline->Branch("horizFFT",horizFFT,"horizFFT[nPointsForTree]/F");

  for (int i=0;i<NUM_ANTS_WITH_NADIRS;i++){
    for (int j=0;j<260*upSampleFactor;j++){
      vertFFTre[i][j]=0;
      horizFFTre[i][j]=0;
      vertFFTim[i][j]=0;
      horizFFTim[i][j]=0;
    }
  }

  for (int eventctr=0;eventctr<numEvents;eventctr++){
    if (eventctr%(numEvents/10)==0) 
      cout<<eventctr<<" neutrinos.  "<<((double)eventctr/(double)numEvents)*100<< "% complete."<<endl;
    //get the event
    fHeadTree->GetEntry(eventctr);
    if (eventStartedFlag!=(int)fHeadPtr->eventNumber) eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
    if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
    if (isChannelSaturated(fHeadPtr->eventNumber)==-1 
	&& (fHeadPtr->trigType&(1<<0))==1 && fHeadPtr->triggerTimeNs>1000000){ //RF triggers
      for (ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){

	TGraph *grv = fUsefulEventPtr->getGraph(ant,AnitaPol::kVertical);	
	TGraph *grh = fUsefulEventPtr->getGraph(ant,AnitaPol::kHorizontal);	
	TGraph *grhInterp;
	TGraph *grvInterp;
	grvInterp=FFTtools::getInterpolatedGraph(grv, deltaTInt);
	grhInterp=FFTtools::getInterpolatedGraph(grh, deltaTInt);
	TGraph *gr1_resized = Resizeplots(grvInterp);//force plot to 256 points
	TGraph *gr1Horiz_resized = Resizeplots(grhInterp);//force plot to 256 points
    
	double *vY = gr1_resized->GetY();
	double *vX = gr1_resized->GetX();
	double deltaT=vX[1]-vX[0];
	int length=gr1_resized->GetN();
	FFTWComplex *theVertFFT=FFTtools::doFFT(length,vY);
	
	double *hY = gr1Horiz_resized->GetY();
	FFTWComplex *theHorizFFT=FFTtools::doFFT(length,hY);

	newLength=(length/2)+1;
	double deltaF=1/(deltaT*length); //Hz
	deltaF*=1e3; //MHz	
	//	cout<<"newLength is "<<newLength<<"\n";
	//cout<<"deltaF is "<<deltaF<<"\n";
	for(int i=0;i<newLength;i++) {
	  
	  vertFFTre[ant][i]+=pow(theVertFFT[i].re,2);
	  vertFFTim[ant][i]+=pow(theVertFFT[i].im,2);
	  horizFFTre[ant][i]+=pow(theHorizFFT[i].re,2);
	  horizFFTim[ant][i]+=pow(theHorizFFT[i].im,2);
	  if (i==0) frequencyArray[i]=0;
	  if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;
	  // cout<<"freq is "<<frequencyArray[i]<<"\n";
	}
	delete gr1_resized;
	delete gr1Horiz_resized;
	delete grhInterp;
	delete grvInterp;
	delete grv;
	delete grh;
	
      }
      eventsIncludedCtr++;
    }
  }
  
  cout<<"number of events included: "<<eventsIncludedCtr<<endl;

  for (ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){
    for (int i=0;i<newLength;i++){
      vertFFTre[ant][i]=vertFFTre[ant][i]/eventsIncludedCtr;
      horizFFTre[ant][i]=horizFFTre[ant][i]/eventsIncludedCtr;
      vertFFTim[ant][i]=vertFFTim[ant][i]/eventsIncludedCtr;
      horizFFTim[ant][i]=horizFFTim[ant][i]/eventsIncludedCtr;
      
      vertFFT[i]=sqrt(vertFFTre[ant][i]+vertFFTim[ant][i]);
      horizFFT[i]=sqrt(horizFFTre[ant][i]+horizFFTim[ant][i]);
      
      if (vertFFT[i]!=0) vertFFT[i]=10*log10(vertFFT[i]/10.);
      if (horizFFT[i]!=0) horizFFT[i]=10*log10(horizFFT[i]/10.);
      
      if (frequencyArray[i]<200 || frequencyArray[i]>1200){
	vertFFT[i]=-1000;
	horizFFT[i]=-1000;
      }
    }
    
    nPointsForTree=newLength;
    tBaseline->Fill();
    
  }
  
  tBaseline->Print();
  rootfile->cd();
  tBaseline->Write();

}
Double_t fitf(Double_t *x, Double_t *par){
   
    Double_t fitval = par[1]*(x[0]/(pow(par[0],2)))*TMath::Exp(-pow(x[0],2)/(2*pow(par[0],2)));
    return fitval;
}
Double_t fitRician(Double_t *x, Double_t *par){
   
    Double_t fitval = par[2]*(x[0]/(pow(par[0],2)))*TMath::Exp(-(pow(x[0],2)+ pow(par[1],2))/(2*pow(par[0],2)))*TMath::BesselI0(x[0]*par[1]/pow(par[0],2));
    return fitval;
}

void MyCorrelator::createSigma(int numEvents, int antctr){
  
  if (!fEventTree) initialize();
  if (numEvents==0)
    numEvents=fEventTree->GetEntries();

  int pol=0;
  char filename[150];
  if(pol==0){
    sprintf(filename,"rayleighfits_run%d_%d.root",fCurrentRun,antctr);
  }
  if(pol==1){
    sprintf(filename,"rayleighfitsHoriz_run%d_%d.root",fCurrentRun,antctr);
  }
  if (printFlag==1) cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");

  TTree *data = new TTree("data","data");

  // int antctr;
  float freqctr;
  double index;
  double parameters[20][2];//Pol, Index, Sigma, Normalization
  int newLength;
  int upSampleFactor=1;
  Double_t deltaTInt=1./(2.6*upSampleFactor);
  char printer[256];
  data->Branch("ant",&antctr,"ant/I");
  data->Branch("parameters",&parameters,"parameters[20][2]/D");
  
  TF1 *rayleighfit_0;
  /* rayleighfit_0= new TF1("rayleighfit_0",fitf,0,1000,2);
  rayleighfit_0->SetLineColor(kRed);
  rayleighfit_0->SetParameter(0,5);//sigma
  rayleighfit_0->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_1;
  /*  rayleighfit_1= new TF1("rayleighfit_1",fitf,0,1000,2);
  rayleighfit_1->SetLineColor(kRed);
  rayleighfit_1->SetParameter(0,5);//sigma
  rayleighfit_1->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_2;
  /* rayleighfit_2= new TF1("rayleighfit_2",fitf,0,1000,2);
  rayleighfit_2->SetLineColor(kRed);
  rayleighfit_2->SetParameter(0,5);//sigma
  rayleighfit_2->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_3;
  /* rayleighfit_3= new TF1("rayleighfit_3",fitf,0,1000,2);
  rayleighfit_3->SetLineColor(kRed);
  rayleighfit_3->SetParameter(0,5);//sigma
  rayleighfit_3->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_4;
  /*  rayleighfit_4= new TF1("rayleighfit_4",fitf,0,1000,2);
  rayleighfit_4->SetLineColor(kRed);
  rayleighfit_4->SetParameter(0,5);//sigma
  rayleighfit_4->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_5;
  /*rayleighfit_5= new TF1("rayleighfit_5",fitf,0,1000,2);
  rayleighfit_5->SetLineColor(kRed);
  rayleighfit_5->SetParameter(0,5);//sigma
  rayleighfit_5->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_6;
  /* rayleighfit_6= new TF1("rayleighfit_6",fitf,0,1000,2);
  rayleighfit_6->SetLineColor(kRed);
  rayleighfit_6->SetParameter(0,5);//sigma
  rayleighfit_6->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_7;
  /*rayleighfit_7= new TF1("rayleighfit_7",fitf,0,1000,2);
  rayleighfit_7->SetLineColor(kRed);
  rayleighfit_7->SetParameter(0,5);//sigma
  rayleighfit_7->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_8;
  /* rayleighfit_8= new TF1("rayleighfit_8",fitf,0,1000,2);
  rayleighfit_8->SetLineColor(kRed);
  rayleighfit_8->SetParameter(0,5);//sigma
  rayleighfit_8->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_9;
  /*rayleighfit_9= new TF1("rayleighfit_9",fitf,0,1000,2);
  rayleighfit_9->SetLineColor(kRed);
  rayleighfit_9->SetParameter(0,5);//sigma
  rayleighfit_9->SetParameter(1,20);//Normalization*/
 TF1 *rayleighfit_10;
 /* rayleighfit_10= new TF1("rayleighfit_10",fitf,0,1000,2);
  rayleighfit_10->SetLineColor(kRed);
  rayleighfit_10->SetParameter(0,5);//sigma
  rayleighfit_10->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_11;
  /*  rayleighfit_11= new TF1("rayleighfit_11",fitf,0,1000,2);
  rayleighfit_11->SetLineColor(kRed);
  rayleighfit_11->SetParameter(0,5);//sigma
  rayleighfit_11->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_12;
  /*rayleighfit_12= new TF1("rayleighfit_12",fitf,0,1000,2);
  rayleighfit_12->SetLineColor(kRed);
  rayleighfit_12->SetParameter(0,5);//sigma
  rayleighfit_12->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_13;
  /*rayleighfit_13= new TF1("rayleighfit_13",fitf,0,1000,2);
  rayleighfit_13->SetLineColor(kRed);
  rayleighfit_13->SetParameter(0,5);//sigma
  rayleighfit_13->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_14;
  /*rayleighfit_14= new TF1("rayleighfit_14",fitf,0,1000,2);
  rayleighfit_14->SetLineColor(kRed);
  rayleighfit_14->SetParameter(0,5);//sigma
  rayleighfit_14->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_15;
  /* rayleighfit_15= new TF1("rayleighfit_15",fitf,0,1000,2);
  rayleighfit_15->SetLineColor(kRed);
  rayleighfit_15->SetParameter(0,5);//sigma
  rayleighfit_15->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_16;
  /* rayleighfit_16= new TF1("rayleighfit_16",fitf,0,1000,2);
  rayleighfit_16->SetLineColor(kRed);
  rayleighfit_16->SetParameter(0,5);//sigma
  rayleighfit_16->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_17;
  /*rayleighfit_17= new TF1("rayleighfit_17",fitf,0,1000,2);
  rayleighfit_17->SetLineColor(kRed);
  rayleighfit_17->SetParameter(0,5);//sigma
  rayleighfit_17->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_18;
  /* rayleighfit_18= new TF1("rayleighfit_18",fitf,0,1000,2);
  rayleighfit_18->SetLineColor(kRed);
  rayleighfit_18->SetParameter(0,5);//sigma
  rayleighfit_18->SetParameter(1,20);//Normalization*/
  TF1 *rayleighfit_19;
  /*rayleighfit_19= new TF1("rayleighfit_19",fitf,0,1000,2);
  rayleighfit_19->SetLineColor(kRed);
  rayleighfit_19->SetParameter(0,5);//sigma
  rayleighfit_19->SetParameter(1,20);//Normalization*/

  TH1F *myHist[20];
  TH1F *myHistHoriz[20];
  char *histname = new char[10];
  char *histname1 = new char[10];
  char *fitname = new char[30];
  for(int i=0;i<20;i++){
    sprintf(histname, "hist_%d",i);
    myHist[i]= new TH1F(histname,"V_{event}",1200,0.,1200.);
    
    sprintf(histname1, "histHoriz_%d",i);
    myHist[i]= new TH1F(histname1,"V_{event}",1200,0.,1200.);
    
  }
    cout<<"fHeadTree events is "<<fHeadTree->GetEntries()<<"\n";
    cout<<"numevents is "<<numEvents<<"\n";
    //for (antctr=0;antctr<3;antctr++){//NUM_ANTS_WITH_NADIRS
    
    
    rayleighfit_0= new TF1("rayleighfit_0",fitf,0,1000,2);
    rayleighfit_0->SetLineColor(kRed);
    rayleighfit_0->SetParameter(0,5);//sigma
    rayleighfit_0->SetParameter(1,1);//A0
    
    rayleighfit_1= new TF1("rayleighfit_1",fitf,0,1000,2);
    rayleighfit_1->SetLineColor(kRed);
    rayleighfit_1->SetParameter(0,5);//sigma
    rayleighfit_1->SetParameter(1,1);//A0
    
    rayleighfit_2= new TF1("rayleighfit_2",fitf,0,1000,2);
    rayleighfit_2->SetLineColor(kRed);
    rayleighfit_2->SetParameter(0,5);//sigma
    rayleighfit_2->SetParameter(1,1);//A0
    
    rayleighfit_3= new TF1("rayleighfit_3",fitf,0,1000,2);
    rayleighfit_3->SetLineColor(kRed);
    rayleighfit_3->SetParameter(0,5);//sigma
    rayleighfit_3->SetParameter(1,1);//A0
   
    rayleighfit_4= new TF1("rayleighfit_4",fitf,0,1000,2);
    rayleighfit_4->SetLineColor(kRed);
    rayleighfit_4->SetParameter(0,5);//sigma
    rayleighfit_4->SetParameter(1,1);//A0
   
    rayleighfit_5= new TF1("rayleighfit_5",fitf,0,1000,2);
    rayleighfit_5->SetLineColor(kRed);
    rayleighfit_5->SetParameter(0,5);//sigma
    rayleighfit_5->SetParameter(1,1);//A0
   
    rayleighfit_6= new TF1("rayleighfit_6",fitf,0,1000,2);
    rayleighfit_6->SetLineColor(kRed);
    rayleighfit_6->SetParameter(0,5);//sigma
    rayleighfit_6->SetParameter(1,1);//A0
   
    rayleighfit_7= new TF1("rayleighfit_7",fitf,0,1000,2);
    rayleighfit_7->SetLineColor(kRed);
    rayleighfit_7->SetParameter(0,5);//sigma
    rayleighfit_7->SetParameter(1,1);//A0
    
    rayleighfit_8= new TF1("rayleighfit_8",fitf,0,1000,2);
    rayleighfit_8->SetLineColor(kRed);
    rayleighfit_8->SetParameter(0,5);//sigma
    rayleighfit_8->SetParameter(1,1);//A0
    
    rayleighfit_9= new TF1("rayleighfit_9",fitf,0,1000,2);
    rayleighfit_9->SetLineColor(kRed);
    rayleighfit_9->SetParameter(0,5);//sigma
    rayleighfit_9->SetParameter(1,1);//A0
    
    rayleighfit_10= new TF1("rayleighfit_10",fitf,0,1000,2);
    rayleighfit_10->SetLineColor(kRed);
    rayleighfit_10->SetParameter(0,5);//sigma
    rayleighfit_10->SetParameter(1,1);//A0
   
    rayleighfit_11= new TF1("rayleighfit_11",fitf,0,1000,2);
    rayleighfit_11->SetLineColor(kRed);
    rayleighfit_11->SetParameter(0,5);//sigma
    rayleighfit_11->SetParameter(1,1);//A0
    
    rayleighfit_12= new TF1("rayleighfit_12",fitf,0,1000,2);
    rayleighfit_12->SetLineColor(kRed);
    rayleighfit_12->SetParameter(0,5);//sigma
    rayleighfit_12->SetParameter(1,1);//A0
    
    rayleighfit_13= new TF1("rayleighfit_13",fitf,0,1000,2);
    rayleighfit_13->SetLineColor(kRed);
    rayleighfit_13->SetParameter(0,5);//sigma
    rayleighfit_13->SetParameter(1,1);//A0
   
    rayleighfit_14= new TF1("rayleighfit_14",fitf,0,1000,2);
    rayleighfit_14->SetLineColor(kRed);
    rayleighfit_14->SetParameter(0,5);//sigma
    rayleighfit_14->SetParameter(1,1);//A0
    
    rayleighfit_15= new TF1("rayleighfit_15",fitf,0,1000,2);
    rayleighfit_15->SetLineColor(kRed);
    rayleighfit_15->SetParameter(0,5);//sigma
    rayleighfit_15->SetParameter(1,1);//A0
    
    rayleighfit_16= new TF1("rayleighfit_16",fitf,0,1000,2);
    rayleighfit_16->SetLineColor(kRed);
    rayleighfit_16->SetParameter(0,5);//sigma
    rayleighfit_16->SetParameter(1,1);//A0
   
    rayleighfit_17= new TF1("rayleighfit_17",fitf,0,1000,2);
    rayleighfit_17->SetLineColor(kRed);
    rayleighfit_17->SetParameter(0,5);//sigma
    rayleighfit_17->SetParameter(1,1);//A0
    
    rayleighfit_18= new TF1("rayleighfit_18",fitf,0,1000,2);
    rayleighfit_18->SetLineColor(kRed);
    rayleighfit_18->SetParameter(0,5);//sigma
    rayleighfit_18->SetParameter(1,1);//A0
    
    rayleighfit_19= new TF1("rayleighfit_19",fitf,0,1000,2);
    rayleighfit_19->SetLineColor(kRed);
    rayleighfit_19->SetParameter(0,5);//sigma
    rayleighfit_19->SetParameter(1,1);//A0
    

    cout<<"ant is "<<antctr<<"\n";
    //cout<<"numEvents is "<<numEvents<<"\n";

    for (int eventctr=0;eventctr<numEvents;eventctr++){
      if (eventctr%(numEvents/10)==0) 
	cout<<eventctr<<" neutrinos.  "<<((double)eventctr/(double)numEvents)*100<< "% complete."<<endl;
      //get the event
      fHeadTree->GetEntry(eventctr);
      if (eventStartedFlag!=(int)fHeadPtr->eventNumber) eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
      if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
      if (isChannelSaturated(fHeadPtr->eventNumber)==-1 
	&& (fHeadPtr->trigType&(1<<0))==1 && fHeadPtr->triggerTimeNs>1000000){ //RF triggers
	
	TGraph *grv = fUsefulEventPtr->getGraph(antctr,AnitaPol::kVertical);	
	TGraph *grh = fUsefulEventPtr->getGraph(antctr,AnitaPol::kHorizontal);	
	TGraph *grhInterp;
	TGraph *grvInterp;
	grvInterp=FFTtools::getInterpolatedGraph(grv, deltaTInt);
	grhInterp=FFTtools::getInterpolatedGraph(grh, deltaTInt);
	
	double *vY = grvInterp->GetY();
	double *vX = grvInterp->GetX();
	double deltaT=vX[1]-vX[0];
	int length=grvInterp->GetN();
	float vertdiff;
	float horizdiff;
	FFTWComplex *theVertFFT=FFTtools::doFFT(length,vY);
	
	double *hY = grvInterp->GetY();
	FFTWComplex *theHorizFFT=FFTtools::doFFT(length,hY);
	
	newLength=(length/2)+1;
	double deltaF=1/(deltaT*length); //Hz
	deltaF*=1e3; //MHz	
	
	
	for(int i=0;i<newLength;i++) {
	  if(i==0) freqctr=0;
	  if(i>0) freqctr= freqctr+deltaF;
	  
	  index = freqctr/50.;
	  index = index-4.;
	  index = floor(index);
	 
	  if(index>=0. && index<20.){
	    
	    vertdiff = pow(theVertFFT[i].re,2)+pow(theVertFFT[i].im,2);
	    vertdiff = sqrt(vertdiff);
	    vertdiff = vertdiff;//- vertAvg[ant][(int)index];
	    
	    horizdiff = pow(theHorizFFT[i].re,2)+pow(theHorizFFT[i].im,2);
	    horizdiff = sqrt(horizdiff);
	    horizdiff = horizdiff; //- horizAvg[ant][(int)index];
	    
	    myHist[(int)index]->Fill(horizdiff);
	    // myHistHoriz[(int)index]->Fill(horizdiff);
	  }//index
	  
	}//newlength
	
	delete grhInterp;
	delete grvInterp;
	delete grv;
	delete grh;
	delete [] theVertFFT;
	delete [] theHorizFFT;
      	
      }//if statements
     
    }//eventctr
    double_t sigma;
    double_t A0;
    
    // for(int pol=0;pol<2;pol++){
      for(int i=0;i<20;i++){
	sprintf(fitname,"rayleighfit_%d",i);
	if(pol==0){
	  if(i<=1){
	    myHist[i]->Fit(fitname,"R","",0,500);
	  }
	  else{
	    myHist[i]->Fit(fitname,"R");
	  }
	}
	if(pol==1){
	  if(i<=1){
	    myHistHoriz[i]->Fit(fitname,"R","",0,500);
	  }
	  else{
	    myHistHoriz[i]->Fit(fitname,"R");
	  }
	}
      }
      
      
      
      TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
     
      gStyle->SetCanvasColor(0);
      c1->UseCurrentStyle();
      c1->Divide(4,5);
      for(int i=0;i<20;i++){
	c1->cd(i+1);
	sprintf(fitname,"rayleighfit_%d",i);
	myHist[i]->Draw();
	if(i==0){
	 
	  rayleighfit_0->Draw("same");
	  TLegend leg0(.5,.7,.9,.9,"");
	  leg0.SetFillColor(0);
	  
	  leg0.AddEntry((TObject*)0,"200-250MHz");
	  leg0.Draw("same");
	}
	if(i==1){
	  rayleighfit_1->Draw("same");
	}
	if(i==2){
	  rayleighfit_2->Draw("same");
	}
	if(i==3){
	  rayleighfit_3->Draw("same");
	}
	if(i==4){
	  rayleighfit_4->Draw("same");
	}
	if(i==5){
	  rayleighfit_5->Draw("same");
	}
	if(i==6){
	  rayleighfit_6->Draw("same");
	}
	if(i==7){
	  rayleighfit_7->Draw("same");
	}
	if(i==8){
	  rayleighfit_8->Draw("same");
	}
	if(i==9){
	  rayleighfit_9->Draw("same");
	}
	if(i==10){
	  rayleighfit_10->Draw("same");
	}
	if(i==11){
	  rayleighfit_11->Draw("same");
	}
	if(i==12){
	  rayleighfit_12->Draw("same");
	}
	if(i==13){
	  rayleighfit_13->Draw("same");
	}
	if(i==14){
	  rayleighfit_14->Draw("same");
	}
	if(i==15){
	  rayleighfit_15->Draw("same");
	}
	if(i==16){
	  rayleighfit_16->Draw("same");
	}
	if(i==17){
	  rayleighfit_17->Draw("same");
	}
	if(i==18){
	  rayleighfit_18->Draw("same");
	}
	if(i==19){
	  rayleighfit_19->Draw("same");
	}
	
      }
      
      sprintf(printer,"OneAnt50MHzbin_%d.png",antctr);
      c1->Print(printer);
      
      for(int i=0;i<20;i++){
	if(i==0){
	  if(rayleighfit_0->GetParameter(0)<0){
	     sigma = -1*rayleighfit_0->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_0->GetParameter(0);
	  }
	  A0 = rayleighfit_0->GetParameter(1);
	}
	if(i==1){
	   if(rayleighfit_1->GetParameter(0)<0){
	     sigma = -1*rayleighfit_1->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_1->GetParameter(0);
	  }
	   A0 = rayleighfit_1->GetParameter(1);
	}
	if(i==2){
	  if(rayleighfit_2->GetParameter(0)<0){
	    sigma = -1*rayleighfit_2->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_2->GetParameter(0);
	  }
	  A0 = rayleighfit_2->GetParameter(1);
	}
	if(i==3){
	  if(rayleighfit_3->GetParameter(0)<0){
	     sigma = -1*rayleighfit_3->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_3->GetParameter(0);
	  }
	  A0 = rayleighfit_3->GetParameter(1);
	}
	if(i==4){
	  if(rayleighfit_4->GetParameter(0)<0){
	    sigma = -1*rayleighfit_4->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_4->GetParameter(0);
	  }
	  A0 = rayleighfit_4->GetParameter(1);
	}
	if(i==5){
	  if(rayleighfit_5->GetParameter(0)<0){
	    sigma = -1*rayleighfit_5->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_5->GetParameter(0);
	  }
	  A0 = rayleighfit_5->GetParameter(1);
	}
	if(i==6){
	  if(rayleighfit_6->GetParameter(0)<0){
	    sigma = -1*rayleighfit_6->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_6->GetParameter(0);
	  }
	  A0 = rayleighfit_6->GetParameter(1);
	}
	if(i==7){
	  if(rayleighfit_7->GetParameter(0)<0){
	    sigma = -1*rayleighfit_7->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_7->GetParameter(0);
	  }
	  A0 = rayleighfit_7->GetParameter(1);
	}
	if(i==8){
	  if(rayleighfit_8->GetParameter(0)<0){
	    sigma = -1*rayleighfit_8->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_8->GetParameter(0);
	  }
	  A0 = rayleighfit_8->GetParameter(1);
	}
	if(i==9){
	  if(rayleighfit_9->GetParameter(0)<0){
	     sigma = -1*rayleighfit_9->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_9->GetParameter(0);
	  }
	  A0 = rayleighfit_9->GetParameter(1);
	}
	if(i==10){
	  if(rayleighfit_10->GetParameter(0)<0){
	    sigma = -1*rayleighfit_10->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_10->GetParameter(0);
	  }
	  A0 = rayleighfit_10->GetParameter(1);
	}
	if(i==11){
	  if(rayleighfit_11->GetParameter(0)<0){
	     sigma = -1*rayleighfit_11->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_11->GetParameter(0);
	  }
	    A0 = rayleighfit_11->GetParameter(1);
	}
	if(i==12){
	  if(rayleighfit_12->GetParameter(0)<0){
	     sigma = -1*rayleighfit_12->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_12->GetParameter(0);
	  }
	  A0 = rayleighfit_12->GetParameter(1);
	}
	if(i==13){
	  if(rayleighfit_13->GetParameter(0)<0){
	     sigma = -1*rayleighfit_13->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_13->GetParameter(0);
	  }
	  A0 = rayleighfit_13->GetParameter(1);
	}
	if(i==14){
	   if(rayleighfit_14->GetParameter(0)<0){
	     sigma = -1*rayleighfit_14->GetParameter(0);
	   }
	   else{
	     sigma = rayleighfit_14->GetParameter(0);
	   }
	   A0 = rayleighfit_14->GetParameter(1);
	}
	if(i==15){
	   if(rayleighfit_15->GetParameter(0)<0){
	     sigma = -1*rayleighfit_15->GetParameter(0);
	   }
	   else{
	     sigma = rayleighfit_15->GetParameter(0);
	   }
	   A0 = rayleighfit_15->GetParameter(1);
	}
	if(i==16){
	  if(rayleighfit_16->GetParameter(0)<0){
	    sigma = -1*rayleighfit_1->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_16->GetParameter(0);
	  }
	  A0 = rayleighfit_16->GetParameter(1);
	}
	if(i==17){
	  if(rayleighfit_17->GetParameter(0)<0){
	    sigma = -1*rayleighfit_1->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_17->GetParameter(0);
	  }
	  A0 = rayleighfit_17->GetParameter(1);
	}
	if(i==18){
	  if(rayleighfit_18->GetParameter(0)<0){
	    sigma = -1*rayleighfit_18->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_18->GetParameter(0);
	  }
	  A0 = rayleighfit_18->GetParameter(1);
	}
	if(i==19){
	  if(rayleighfit_19->GetParameter(0)<0){
	    sigma = -1*rayleighfit_19->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_19->GetParameter(0);
	  }
	  A0 = rayleighfit_19->GetParameter(1);
	}
	
	cout<<"for i = "<<i<<" sigma is "<<sigma<<" and A0 is "<<A0<<"\n";
	
	
	parameters[i][0]=sigma;
	parameters[i][1]=A0;
	
	
      }
      //}
    data->Fill();

    for(int k=0;k<20;k++){
      myHist[k]->Clear();
      
    }
    //delete myHist[20];
    delete rayleighfit_0;
    delete rayleighfit_1;
    delete rayleighfit_2;
    delete rayleighfit_3;
    delete rayleighfit_4;
    delete rayleighfit_5;
    delete rayleighfit_6;
    delete rayleighfit_7;
    delete rayleighfit_8;
    delete rayleighfit_9;
    delete rayleighfit_10;
    delete rayleighfit_11;
    delete rayleighfit_12;
    delete rayleighfit_13;
    delete rayleighfit_14;
    delete rayleighfit_15;
    delete rayleighfit_16;
    delete rayleighfit_17;
    delete rayleighfit_18;
    delete rayleighfit_19;
    
    //}//antctr
  
  
  

  rootfile->Write();
  rootfile->Close();

}//createsigma
///////////////////////////////
void MyCorrelator::createdBCut(int numEvents, int antctr){
  
  if (!fEventTree) initialize();
  if (numEvents==0)
    numEvents=fEventTree->GetEntries();


   readBaselineFFTs();
  int pol=0;
  char filename[150];
  if(pol==0){
    sprintf(filename,"dBCuts_run%d_%d.root",fCurrentRun,antctr);
  }
  if(pol==1){
    sprintf(filename,"dBCutsHoriz_run%d_%d.root",fCurrentRun,antctr);
  }
  if (printFlag==1) cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");

  TTree *data = new TTree("data","data");

  // int antctr;
  float freqctr;
  double index;
  double parameters[20][2];//Pol, Index, Sigma, Normalization
  int newLength;
  int upSampleFactor=1;
  Double_t deltaTInt=1./(2.6*upSampleFactor);
  char printer[256];
  data->Branch("ant",&antctr,"ant/I");
  data->Branch("parameters",&parameters,"parameters[20][2]/D");
  
  TF1 *rayleighfit_0;
  TF1 *rayleighfit_1;
  TF1 *rayleighfit_2;
  TF1 *rayleighfit_3;
  TF1 *rayleighfit_4;
  TF1 *rayleighfit_5;
  TF1 *rayleighfit_6;
  TF1 *rayleighfit_7;
  TF1 *rayleighfit_8;
  TF1 *rayleighfit_9;
  TF1 *rayleighfit_10;
  TF1 *rayleighfit_11;
  TF1 *rayleighfit_12;
  TF1 *rayleighfit_13;
  TF1 *rayleighfit_14;
  TF1 *rayleighfit_15;
  TF1 *rayleighfit_16;
  TF1 *rayleighfit_17;
  TF1 *rayleighfit_18;
  TF1 *rayleighfit_19;
 
  TH1F *myHist[20];
  TH1F *myHistHoriz[20];
  char *histname = new char[10];
  char *histname1 = new char[10];
  char *fitname = new char[30];
  for(int i=0;i<20;i++){
    sprintf(histname, "hist_%d",i);
    myHist[i]= new TH1F(histname,"V_{event}",1200,0.,1200.);
    
    sprintf(histname1, "histHoriz_%d",i);
    myHist[i]= new TH1F(histname1,"V_{event}",1200,0.,1200.);
    
  }
    cout<<"fHeadTree events is "<<fHeadTree->GetEntries()<<"\n";
    cout<<"numevents is "<<numEvents<<"\n";
    //for (antctr=0;antctr<NUM_ANTS_WITH_NADIRS;antctr++){//NUM_ANTS_WITH_NADIRS
    
    
      rayleighfit_0= new TF1("rayleighfit_0",fitRician,0,1000,3);
      rayleighfit_0->SetLineColor(kRed);
      rayleighfit_0->SetParameter(0,5);//sigma
      rayleighfit_0->SetParameter(1,1);//A0
      rayleighfit_0->SetParameter(2,5);//signal
      
      rayleighfit_1= new TF1("rayleighfit_1",fitRician,0,1000,3);
      rayleighfit_1->SetLineColor(kRed);
      rayleighfit_1->SetParameter(0,5);//sigma
      rayleighfit_1->SetParameter(1,1);//A0
      rayleighfit_0->SetParameter(2,5);//signal
      
      rayleighfit_2= new TF1("rayleighfit_2",fitf,0,1000,2);
      rayleighfit_2->SetLineColor(kRed);
      rayleighfit_2->SetParameter(0,5);//sigma
      rayleighfit_2->SetParameter(1,1);//A0
      
      rayleighfit_3= new TF1("rayleighfit_3",fitf,0,1000,2);
      rayleighfit_3->SetLineColor(kRed);
      rayleighfit_3->SetParameter(0,5);//sigma
      rayleighfit_3->SetParameter(1,1);//A0
      
      rayleighfit_4= new TF1("rayleighfit_4",fitf,0,1000,2);
      rayleighfit_4->SetLineColor(kRed);
      rayleighfit_4->SetParameter(0,5);//sigma
      rayleighfit_4->SetParameter(1,1);//A0
      
      rayleighfit_5= new TF1("rayleighfit_5",fitf,0,1000,2);
      rayleighfit_5->SetLineColor(kRed);
      rayleighfit_5->SetParameter(0,5);//sigma
      rayleighfit_5->SetParameter(1,1);//A0
      
      rayleighfit_6= new TF1("rayleighfit_6",fitf,0,1000,2);
      rayleighfit_6->SetLineColor(kRed);
      rayleighfit_6->SetParameter(0,5);//sigma
      rayleighfit_6->SetParameter(1,1);//A0
      
      rayleighfit_7= new TF1("rayleighfit_7",fitf,0,1000,2);
      rayleighfit_7->SetLineColor(kRed);
      rayleighfit_7->SetParameter(0,5);//sigma
      rayleighfit_7->SetParameter(1,1);//A0
      
      rayleighfit_8= new TF1("rayleighfit_8",fitf,0,1000,2);
      rayleighfit_8->SetLineColor(kRed);
      rayleighfit_8->SetParameter(0,5);//sigma
      rayleighfit_8->SetParameter(1,1);//A0
      
      rayleighfit_9= new TF1("rayleighfit_9",fitf,0,1000,2);
      rayleighfit_9->SetLineColor(kRed);
      rayleighfit_9->SetParameter(0,5);//sigma
      rayleighfit_9->SetParameter(1,1);//A0
      
      rayleighfit_10= new TF1("rayleighfit_10",fitf,0,1000,2);
      rayleighfit_10->SetLineColor(kRed);
      rayleighfit_10->SetParameter(0,5);//sigma
      rayleighfit_10->SetParameter(1,1);//A0
      
      rayleighfit_11= new TF1("rayleighfit_11",fitf,0,1000,2);
      rayleighfit_11->SetLineColor(kRed);
      rayleighfit_11->SetParameter(0,5);//sigma
      rayleighfit_11->SetParameter(1,1);//A0
      
      rayleighfit_12= new TF1("rayleighfit_12",fitf,0,1000,2);
      rayleighfit_12->SetLineColor(kRed);
      rayleighfit_12->SetParameter(0,5);//sigma
      rayleighfit_12->SetParameter(1,1);//A0
      
      rayleighfit_13= new TF1("rayleighfit_13",fitf,0,1000,2);
      rayleighfit_13->SetLineColor(kRed);
      rayleighfit_13->SetParameter(0,5);//sigma
      rayleighfit_13->SetParameter(1,1);//A0
      
      rayleighfit_14= new TF1("rayleighfit_14",fitf,0,1000,2);
      rayleighfit_14->SetLineColor(kRed);
      rayleighfit_14->SetParameter(0,5);//sigma
      rayleighfit_14->SetParameter(1,1);//A0
      
      rayleighfit_15= new TF1("rayleighfit_15",fitf,0,1000,2);
      rayleighfit_15->SetLineColor(kRed);
      rayleighfit_15->SetParameter(0,5);//sigma
      rayleighfit_15->SetParameter(1,1);//A0
      
      rayleighfit_16= new TF1("rayleighfit_16",fitf,0,1000,2);
      rayleighfit_16->SetLineColor(kRed);
      rayleighfit_16->SetParameter(0,5);//sigma
      rayleighfit_16->SetParameter(1,1);//A0
      
      rayleighfit_17= new TF1("rayleighfit_17",fitf,0,1000,2);
      rayleighfit_17->SetLineColor(kRed);
      rayleighfit_17->SetParameter(0,5);//sigma
      rayleighfit_17->SetParameter(1,1);//A0
      
      rayleighfit_18= new TF1("rayleighfit_18",fitf,0,1000,2);
      rayleighfit_18->SetLineColor(kRed);
      rayleighfit_18->SetParameter(0,5);//sigma
      rayleighfit_18->SetParameter(1,1);//A0
    
      rayleighfit_19= new TF1("rayleighfit_19",fitf,0,1000,2);
      rayleighfit_19->SetLineColor(kRed);
      rayleighfit_19->SetParameter(0,5);//sigma
      rayleighfit_19->SetParameter(1,1);//A0
    

    cout<<"ant is "<<antctr<<"\n";
    //cout<<"numEvents is "<<numEvents<<"\n";

    for (int eventctr=0;eventctr<numEvents;eventctr++){
      if (eventctr%(numEvents/10)==0) 
	cout<<eventctr<<" neutrinos.  "<<((double)eventctr/(double)numEvents)*100<< "% complete."<<endl;
      //get the event
      fHeadTree->GetEntry(eventctr);
      if (eventStartedFlag!=(int)fHeadPtr->eventNumber) eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
      if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
      if (isChannelSaturated(fHeadPtr->eventNumber)==-1 
	&& (fHeadPtr->trigType&(1<<0))==1 && fHeadPtr->triggerTimeNs>1000000){ //RF triggers
	
	TGraph *grv = fUsefulEventPtr->getGraph(antctr,AnitaPol::kVertical);	
	TGraph *grh = fUsefulEventPtr->getGraph(antctr,AnitaPol::kHorizontal);	
	TGraph *grhInterp;
	TGraph *grvInterp;
	grvInterp=FFTtools::getInterpolatedGraph(grv, deltaTInt);
	grhInterp=FFTtools::getInterpolatedGraph(grh, deltaTInt);
	 TGraph *grvInterp_resized = Resizeplots(grvInterp);//force plot to 256 points
	 TGraph *grhInterp_resized = Resizeplots(grhInterp);//force plot to 256 points
    
	 delete grvInterp;
	 delete grhInterp;
    
	double *vY = grvInterp_resized->GetY();
	double *vX = grvInterp_resized->GetX();
	double deltaT=vX[1]-vX[0];
	int length=grvInterp_resized->GetN();
	float vertdiff;
	float horizdiff;
	FFTWComplex *theVertFFT=FFTtools::doFFT(length,vY);
	
	double *hY = grhInterp_resized->GetY();
	FFTWComplex *theHorizFFT=FFTtools::doFFT(length,hY);
	
	newLength=(length/2)+1;
	double deltaF=1/(deltaT*length); //Hz
	deltaF*=1e3; //MHz	

	for(int i=0;i<newLength;i++) {
	  if(i==0) freqctr=0;
	  if(i>0) freqctr= freqctr+deltaF;
	  
	  index = freqctr/50.;
	  index = index-4.;
	  index = floor(index);
	 
	  if(index>=0. && index<20.){
	    
	    vertdiff = pow(theVertFFT[i].re,2)+pow(theVertFFT[i].im,2);
	    vertdiff = sqrt(vertdiff);
	    vertdiff = vertdiff-baselinearray[antctr][i];
	    
	    horizdiff = pow(theHorizFFT[i].re,2)+pow(theHorizFFT[i].im,2);
	    horizdiff = sqrt(horizdiff);
	    horizdiff = horizdiff-baselinearrayHoriz[antctr][i];
	    
	    if(pol==0){
	      myHist[(int)index]->Fill(vertdiff);
	    }
	    if(pol==1){
	      myHist[(int)index]->Fill(horizdiff);
	    }
	    // myHistHoriz[(int)index]->Fill(horizdiff);
	  }//index
	  
	}//newlength
	
	delete grhInterp_resized;
	delete grvInterp_resized;
	delete grv;
	delete grh;
	delete [] theVertFFT;
	delete [] theHorizFFT;
      	
      }//if statements
     
    }//eventctr
    double_t sigma;
    double_t A0;
    
    // for(int pol=0;pol<2;pol++){
      for(int i=0;i<20;i++){
	sprintf(fitname,"rayleighfit_%d",i);
	if(pol==0){
	  if(i<=1){
	    myHist[i]->Fit(fitname,"R","",0,1000);
	  }
	  else{
	    myHist[i]->Fit(fitname,"R");
	  }
	}
	if(pol==1){
	  if(i<=1){
	    myHistHoriz[i]->Fit(fitname,"R","",0,1000);
	  }
	  else{
	    myHistHoriz[i]->Fit(fitname,"R");
	  }
	}
      }
      
      
      
      TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
     
      gStyle->SetCanvasColor(0);
      c1->UseCurrentStyle();
      c1->Divide(4,5);
      for(int i=0;i<20;i++){
	c1->cd(i+1);
	sprintf(fitname,"rayleighfit_%d",i);
	myHist[i]->Draw();
	if(i==0){
	 
	  rayleighfit_0->Draw("same");
	  TLegend leg0(.5,.7,.9,.9,"");
	  leg0.SetFillColor(0);
	  
	  leg0.AddEntry((TObject*)0,"200-250MHz");
	  leg0.Draw("same");
	}
	if(i==1){
	  rayleighfit_1->Draw("same");
	}
	if(i==2){
	  rayleighfit_2->Draw("same");
	}
	if(i==3){
	  rayleighfit_3->Draw("same");
	}
	if(i==4){
	  rayleighfit_4->Draw("same");
	}
	if(i==5){
	  rayleighfit_5->Draw("same");
	}
	if(i==6){
	  rayleighfit_6->Draw("same");
	}
	if(i==7){
	  rayleighfit_7->Draw("same");
	}
	if(i==8){
	  rayleighfit_8->Draw("same");
	}
	if(i==9){
	  rayleighfit_9->Draw("same");
	}
	if(i==10){
	  rayleighfit_10->Draw("same");
	}
	if(i==11){
	  rayleighfit_11->Draw("same");
	}
	if(i==12){
	  rayleighfit_12->Draw("same");
	}
	if(i==13){
	  rayleighfit_13->Draw("same");
	}
	if(i==14){
	  rayleighfit_14->Draw("same");
	}
	if(i==15){
	  rayleighfit_15->Draw("same");
	}
	if(i==16){
	  rayleighfit_16->Draw("same");
	}
	if(i==17){
	  rayleighfit_17->Draw("same");
	}
	if(i==18){
	  rayleighfit_18->Draw("same");
	}
	if(i==19){
	  rayleighfit_19->Draw("same");
	}
	
      }
      
      sprintf(printer,"OneAnt50MHzbin_%d.png",antctr);
      c1->Print(printer);
      
      for(int i=0;i<20;i++){
	if(i==0){
	  if(rayleighfit_0->GetParameter(0)<0){
	     sigma = -1*rayleighfit_0->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_0->GetParameter(0);
	  }
	  A0 = rayleighfit_0->GetParameter(1);
	}
	if(i==1){
	   if(rayleighfit_1->GetParameter(0)<0){
	     sigma = -1*rayleighfit_1->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_1->GetParameter(0);
	  }
	   A0 = rayleighfit_1->GetParameter(1);
	}
	if(i==2){
	  if(rayleighfit_2->GetParameter(0)<0){
	    sigma = -1*rayleighfit_2->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_2->GetParameter(0);
	  }
	  A0 = rayleighfit_2->GetParameter(1);
	}
	if(i==3){
	  if(rayleighfit_3->GetParameter(0)<0){
	     sigma = -1*rayleighfit_3->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_3->GetParameter(0);
	  }
	  A0 = rayleighfit_3->GetParameter(1);
	}
	if(i==4){
	  if(rayleighfit_4->GetParameter(0)<0){
	    sigma = -1*rayleighfit_4->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_4->GetParameter(0);
	  }
	  A0 = rayleighfit_4->GetParameter(1);
	}
	if(i==5){
	  if(rayleighfit_5->GetParameter(0)<0){
	    sigma = -1*rayleighfit_5->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_5->GetParameter(0);
	  }
	  A0 = rayleighfit_5->GetParameter(1);
	}
	if(i==6){
	  if(rayleighfit_6->GetParameter(0)<0){
	    sigma = -1*rayleighfit_6->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_6->GetParameter(0);
	  }
	  A0 = rayleighfit_6->GetParameter(1);
	}
	if(i==7){
	  if(rayleighfit_7->GetParameter(0)<0){
	    sigma = -1*rayleighfit_7->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_7->GetParameter(0);
	  }
	  A0 = rayleighfit_7->GetParameter(1);
	}
	if(i==8){
	  if(rayleighfit_8->GetParameter(0)<0){
	    sigma = -1*rayleighfit_8->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_8->GetParameter(0);
	  }
	  A0 = rayleighfit_8->GetParameter(1);
	}
	if(i==9){
	  if(rayleighfit_9->GetParameter(0)<0){
	     sigma = -1*rayleighfit_9->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_9->GetParameter(0);
	  }
	  A0 = rayleighfit_9->GetParameter(1);
	}
	if(i==10){
	  if(rayleighfit_10->GetParameter(0)<0){
	    sigma = -1*rayleighfit_10->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_10->GetParameter(0);
	  }
	  A0 = rayleighfit_10->GetParameter(1);
	}
	if(i==11){
	  if(rayleighfit_11->GetParameter(0)<0){
	     sigma = -1*rayleighfit_11->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_11->GetParameter(0);
	  }
	    A0 = rayleighfit_11->GetParameter(1);
	}
	if(i==12){
	  if(rayleighfit_12->GetParameter(0)<0){
	     sigma = -1*rayleighfit_12->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_12->GetParameter(0);
	  }
	  A0 = rayleighfit_12->GetParameter(1);
	}
	if(i==13){
	  if(rayleighfit_13->GetParameter(0)<0){
	     sigma = -1*rayleighfit_13->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_13->GetParameter(0);
	  }
	  A0 = rayleighfit_13->GetParameter(1);
	}
	if(i==14){
	   if(rayleighfit_14->GetParameter(0)<0){
	     sigma = -1*rayleighfit_14->GetParameter(0);
	   }
	   else{
	     sigma = rayleighfit_14->GetParameter(0);
	   }
	   A0 = rayleighfit_14->GetParameter(1);
	}
	if(i==15){
	   if(rayleighfit_15->GetParameter(0)<0){
	     sigma = -1*rayleighfit_15->GetParameter(0);
	   }
	   else{
	     sigma = rayleighfit_15->GetParameter(0);
	   }
	   A0 = rayleighfit_15->GetParameter(1);
	}
	if(i==16){
	  if(rayleighfit_16->GetParameter(0)<0){
	    sigma = -1*rayleighfit_1->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_16->GetParameter(0);
	  }
	  A0 = rayleighfit_16->GetParameter(1);
	}
	if(i==17){
	  if(rayleighfit_17->GetParameter(0)<0){
	    sigma = -1*rayleighfit_1->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_17->GetParameter(0);
	  }
	  A0 = rayleighfit_17->GetParameter(1);
	}
	if(i==18){
	  if(rayleighfit_18->GetParameter(0)<0){
	    sigma = -1*rayleighfit_18->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_18->GetParameter(0);
	  }
	  A0 = rayleighfit_18->GetParameter(1);
	}
	if(i==19){
	  if(rayleighfit_19->GetParameter(0)<0){
	    sigma = -1*rayleighfit_19->GetParameter(0);
	  }
	  else{
	    sigma = rayleighfit_19->GetParameter(0);
	  }
	  A0 = rayleighfit_19->GetParameter(1);
	}
	
	cout<<"for i = "<<i<<" sigma is "<<sigma<<" and A0 is "<<A0<<"\n";
	
	
	parameters[i][0]=sigma;
	parameters[i][1]=A0;
	
	
      }
      //}
    data->Fill();

    for(int k=0;k<20;k++){
      myHist[k]->Clear();
      
    }
    //delete myHist[20];
    delete rayleighfit_0;
    delete rayleighfit_1;
    delete rayleighfit_2;
    delete rayleighfit_3;
    delete rayleighfit_4;
    delete rayleighfit_5;
    delete rayleighfit_6;
    delete rayleighfit_7;
    delete rayleighfit_8;
    delete rayleighfit_9;
    delete rayleighfit_10;
    delete rayleighfit_11;
    delete rayleighfit_12;
    delete rayleighfit_13;
    delete rayleighfit_14;
    delete rayleighfit_15;
    delete rayleighfit_16;
    delete rayleighfit_17;
    delete rayleighfit_18;
    delete rayleighfit_19;
    
    //}//antctr
  
  
  

  rootfile->Write();
  rootfile->Close();

}//createdbCut


////////////////////////////////
//______________________________________________________________________________
///////////BEGIN GENERAL PURPOSE STUFF

MyCorrelator*  MyCorrelator::Instance()
{
   //static function
   return (fgInstance) ? (MyCorrelator*) fgInstance : new MyCorrelator();
}
int MyCorrelator::getHeaderEntry()
{
  
  if(!fEventTree){
    int eventTree=loadEventTree(); 
    eventTree=0;
  }
  
  if(fEventEntry<fHeadTree->GetEntries())
    fHeadTree->GetEntry(fEventEntry);
  else {
    std::cout << "No more entries in header tree" << endl;
    return -1;
  }
      
  return 0;
}

int MyCorrelator::getEventEntry()
{
  // cout<<"here in getEventEntry \n";
  if(!fEventTree){
    int eventTree=loadEventTree();
    eventTree=0;
  }
 
  if(fEventEntry<fEventTree->GetEntries()){
    //cout<<"fEventEntry is "<<fEventEntry<<"\n";
    //cout<<"fEventTree entries is "<<fEventTree->GetEntries()<<"\n";
    fEventTree->GetEntry(fEventEntry);
    //cout<<"made it past! \n";
    
  }
   else {
      std::cout << "No more entries in event tree" << endl;
      return -1;
   }
  
  //   fUsefulEventPtr = new UsefulAnitaEvent(fRawEventPtr,WaveCalType::kVoltageTime);  
  //fUsefulEventPtr = new UsefulAnitaEvent(fRawEventPtr,WaveCalType::kVTLabAG); 

  if(fUseCalibratedEventFile) {
    // cout<<"here in Calibrated EventFile \n";
    if(fUsefulEventPtr)
      delete fUsefulEventPtr;
    if(fCalType==WaveCalType::kDefault){
      // cout<<"here  fCalType is Default\n";
      //cout<<"fCalEventPtr is "<<fCalEventPtr<<"\n";
      fUsefulEventPtr = new UsefulAnitaEvent(fCalEventPtr);
     
    }
    else{
      
      fUsefulEventPtr = new UsefulAnitaEvent((RawAnitaEvent*)fCalEventPtr,fCalType,fHkPtr);
      
    }
  }
  else if (fUseEventFile) {
    cout<<"here 1\n";
    if(fUsefulEventPtr)
      delete fUsefulEventPtr;
    fUsefulEventPtr = new UsefulAnitaEvent(fRawEventPtr,WaveCalType::kVTFullAGCrossCorClock); 
    //fUsefulEventPtr = new UsefulAnitaEvent(fRawEventPtr,fCalType,fHeadPtr);  
  }
  // 
  //Need to make configurable at some point
  //This will also need to be modifed to make realEvent accessible outside here
  if (fUseCalibratedEventFile && fCalType==WaveCalType::kDefault) return fCalEventPtr->eventNumber;
  else if (fUseEventFile) return fRawEventPtr->eventNumber;
  else return fUsefulEventPtr->eventNumber;
}


void MyCorrelator::closeCurrentRun()
{
  
  if(fHeadFile)
    fHeadFile->Close();
  if(fEventFile)
    fEventFile->Close();
  if(fTurfRateFile)
    fTurfRateFile->Close();
  if(fSumTurfRateFile)
    fSumTurfRateFile->Close();
  if(fSurfHkFile)
    fSurfHkFile->Close();
  if(fAvgSurfHkFile)
    fAvgSurfHkFile->Close();

  fHeadFile=0;
  fEventFile=0;
  fTurfRateFile=0;
  fSumTurfRateFile=0;
  fSurfHkFile=0;
  fAvgSurfHkFile=0;

  fHeadTree=0;
  fEventTree=0;
  fTurfRateTree=0;
  fSumTurfRateTree=0;
  fSurfHkTree=0;
  fAvgSurfHkTree=0;
}


int MyCorrelator::loadEventTree()
{       
  Int_t fGotCalEventFile=0;
  Int_t fGotEventFile=0;
  char eventName[FILENAME_MAX];
  char headerName[FILENAME_MAX];
  cout<<"filename_max is "<<FILENAME_MAX<<" \n";
  sprintf(headerName,"%s/run%d/headFile%d.root",fCurrentBaseDir,fCurrentRun,fCurrentRun);
  cout<<"headername is "<<headerName<<"\n";
  cout<<"here! \n";
  if(1) {
    cout<<"fCurrentBaseDir is "<<fCurrentBaseDir<<"\n";
      //Will try and use calibrated event files  
      sprintf(eventName,"%s/run%d/calEventFile%d.root",fCurrentBaseDir,fCurrentRun,fCurrentRun);    
      fEventTree = new TChain("eventTree");
      fEventTree->Add(eventName);
      cout<<"Eventname is "<<eventName<<"\n";
      for(int extra=1;extra<100;extra++) {
	sprintf(eventName,"%s/run%d/calEventFile%d_%d.root",fCurrentBaseDir,fCurrentRun,fCurrentRun,extra);
	//TFile fpTest(eventName);
	
	TFile *fpTest = TFile::Open(eventName);
	if(!fpTest) 
	  break;
	else {
	  delete fpTest;
	  fEventTree->Add(eventName);
	}
      }
      if(fEventTree->GetEntries()>0) {
	fGotCalEventFile=1;
	fUseCalibratedEventFile=1;
	fEventTree->SetBranchAddress("event",&fCalEventPtr);
      }
      else {
	fUseCalibratedEventFile=0;
	fGotCalEventFile=0;
      }
    }
    
    if(!fGotCalEventFile) {
      sprintf(eventName,"%s/run%d/eventFile%d.root",fCurrentBaseDir,fCurrentRun,fCurrentRun);    
      fEventTree = new TChain("eventTree");
      fEventTree->Add(eventName);
      
      for(int extra=1;extra<100;extra++) {
	sprintf(eventName,"%s/run%d/eventFile%d_%d.root",fCurrentBaseDir,fCurrentRun,fCurrentRun,extra);
	TFile *fpTest = TFile::Open(eventName);
	if(!fpTest) 
	  break;
	else {
	  delete fpTest;
	  fEventTree->Add(eventName);
	}
      }
      if(fEventTree->GetEntries()>0) {
	fGotEventFile=1;
	fUseEventFile=1;
	fEventTree->SetBranchAddress("event",&fRawEventPtr);
      }
      else{
	fGotEventFile=0;
	fUseEventFile=0;
      }
    }
    if (!fGotCalEventFile && !fGotEventFile){
      cout<<"Here at importing sim data! \n";
      sprintf(eventName,"%s/icefinal.root",fCurrentBaseDir);
      sprintf(headerName,"%s/headFile.root",fCurrentBaseDir);
      cout<<"eventName is "<<eventName<<"\n";
      fEventTree = new TChain("eventTree");
      fEventTree->Add(eventName);
      
     
      sprintf(eventName,"%s/icefinal.root",fCurrentBaseDir);    
      // fEventTree = new TChain("eventTree");
      // fEventTree->Add(eventName);
      for(int extra=1;extra<3;extra++) {
	sprintf(eventName,"%s/icefinal_%d.root",fCurrentBaseDir,extra);
	TFile *fpTest = TFile::Open(eventName);
	if(!fpTest){
	   cout<<"broke at extra = "<<extra<<"\n";
	  break;
	}
	else {
	  delete fpTest;
	  fEventTree->Add(eventName);
	 
	}
      }
      
      //How to work with the pointers? 
      //fEventTree->SetBranchAddress("fUsefulEventPtr",&fUsefulEventPtr);???
      //do we need weights?
      // mc events?


      if (fUsefulEventPtr) delete fUsefulEventPtr;//work on IceMC simulated data???
      if(fEventTree->GetEntries()>0) {
	fEventTree->SetBranchAddress("UsefulAnitaEvent",&fUsefulEventPtr);  //simulated data is in useful event format  //event
	/*fEventTree->SetBranchAddress("weightEvent",&mcWeightEvent);
	fEventTree->SetBranchAddress("exponentEvent",&mcExponentEvent);
	fEventTree->SetBranchAddress("latEvent",&mcLatEvent);
	fEventTree->SetBranchAddress("lonEvent",&mcLonEvent);
	fEventTree->SetBranchAddress("altEvent",&mcAltEvent);*/
	}
    }
    
    if(fEventTree->GetEntries()<1) {
      cout << "Couldn't open: " << eventName << "\n";
      return -1;
    }
    
     if(fCurrentRun ==1){
      cout<<"loading sim head file! \n";
      cout<<"headername is "<<headerName<<"\n";
      fHeadFile = TFile::Open(headerName);
      fHeadTree = (TTree*) fHeadFile->Get("headTree");
       
      if(!fHeadTree) {
	cout << "Couldn't get headTree from " << headerName << endl;
	return -1;
      }

      }
     else{
       fHeadFile = TFile::Open(headerName);
       if(!fHeadFile) {
	 cout << "Couldn't open: " << headerName << "\n";
	 return -1;
    }
       fHeadTree = (TTree*) fHeadFile->Get("headTree");
       
       if(!fHeadTree) {
	 cout << "Couldn't get headTree from " << headerName << endl;
	 return -1;
      }
    }
  fHeadTree->SetBranchAddress("header",&fHeadPtr);
  fEventEntry=0;
  fHeadTree->BuildIndex("eventNumber");
  fHeadIndex = (TTreeIndex*) fHeadTree->GetTreeIndex();
  std::cerr << fEventTree << "\t" << fHeadTree << "\n";
  std::cerr << fHeadTree->GetEntries() << "\t"
	    << fEventTree->GetEntries() << "\n";
  return 0;
}

void MyCorrelator::loadTurfTree() 
{
  char turfName[FILENAME_MAX];
  sprintf(turfName,"%s/run%d/turfRateFile%d.root",fCurrentBaseDir,
	  fCurrentRun,fCurrentRun);
  fTurfRateFile = new TFile(turfName);
  if(!fTurfRateFile) {
    cout << "Couldn't open: " << turfName << "\n";
    return;
  }
  fTurfRateTree = (TTree*) fTurfRateFile->Get("turfRateTree");
  if(!fTurfRateTree) {
    cout << "Couldn't get turfRateTree from " << turfName << endl;
    return;
  }
  fTurfRateTree->SetBranchAddress("turf",&fTurfPtr);
  fTurfRateEntry=0;
} 


int MyCorrelator::getTurfEntry() 
{
  if(!fTurfRateTree) {
    loadTurfTree();     
  }
  if(fTurfRateEntry<fTurfRateTree->GetEntries())
    fTurfRateTree->GetEntry(fTurfRateEntry);
  else {
    std::cout << "No more entries in turf rate tree" << endl;
      return -1;
   }
   //   std::cout << fTurfRateEntry << "\t" << fTurfPtr->realTime 
   //	     << "\t" << fTurfPtr->ppsNum << std::endl;

   return 0;
}
int MyCorrelator::getSurfIDEntry() 
{
  if(!fEventTree) {
    int eventTree=loadEventTree(); 
    eventTree=0;
  }
  if(fEventEntry<fEventTree->GetEntries())
    fEventTree->GetEntry(fEventEntry);
  else {
    std::cout << "No more entries in event tree" << endl;
    return -1;
  }
//    //   std::cout << fTurfRateEntry << "\t" << fTurfPtr->realTime 
//    //	     << "\t" << fTurfPtr->ppsNum << std::endl;

  return 0;
}


void MyCorrelator::loadSurfTree()
{
 char surfName[FILENAME_MAX];
 sprintf(surfName,"%s/run%d/surfHkFile%d.root",fCurrentBaseDir,
	 fCurrentRun,fCurrentRun);
 fSurfHkFile = new TFile(surfName);
 if(!fSurfHkFile) {
   cout << "Couldn't open: " << surfName << "\n";
   return;
      }
 fSurfHkTree = (TTree*) fSurfHkFile->Get("surfHkTree");
 if(!fSurfHkTree) {
   cout << "Couldn't get surfHkTree from " << surfName << endl;
   return;
 }
 fSurfHkTree->SetBranchAddress("surf",&fSurfPtr);
 fSurfHkEntry=0;


}

void MyCorrelator::loadAcqdTree()
{
  char acqname[FILENAME_MAX];
  sprintf(acqname,"%s/run%d/auxFile%d*.root",fCurrentBaseDir,fCurrentRun,fCurrentRun);
  fAcqdFile = new TFile(acqname);
  if(!fAcqdFile) {
    cout << "Couldn't open: " << acqname << "\n";
    return;
  }
  
  fAcqTree = new TChain("acqdStartTree");
  
  if(!fAcqTree) {
    cout << "Couldn't get AcqdTree from " << acqname << endl;
   return;
  }
  
  fAcqTree->Add(acqname);
  fAcqTree->SetBranchAddress("acqd",&fAcqPtr);
  
}

int MyCorrelator::getSurfEntry() 
{
  if(!fSurfHkTree) {
    loadSurfTree();     
   }
   //   std::cerr << 
   if(fSurfHkEntry<fSurfHkTree->GetEntries())
      fSurfHkTree->GetEntry(fSurfHkEntry);
   else {
      std::cout << "No more entries in surfHkTree" << endl;
      return -1;
   }
   //   std::cout << fSurfHkEntry << "\t" << fSurfPtr->realTime 
   //	     << "\t" << fSurfPtr->ppsNum << std::endl;

   return 0;
}

void MyCorrelator::loadAvgSurfTree()
{
  char surfName[FILENAME_MAX];
  sprintf(surfName,"%s/run%d/avgSurfHkFile%d.root",fCurrentBaseDir,
	  fCurrentRun,fCurrentRun);
  fSurfHkFile = new TFile(surfName);
  if(!fSurfHkFile) {
    cout << "Couldn't open: " << surfName << "\n";
    return;
  }
  fAvgSurfHkTree = (TTree*) fSurfHkFile->Get("avgSurfHkTree");
  if(!fAvgSurfHkTree) {
    cout << "Couldn't get avgSurfHkTree from " << surfName << endl;
    return;
  }
  fAvgSurfHkTree->SetBranchAddress("avgsurf",&fAvgSurfPtr);
  fAvgSurfHkEntry=0;

}

int MyCorrelator::getAvgSurfEntry() 
{
  if(!fAvgSurfHkTree) {
    loadAvgSurfTree();
  }
   //   std::cerr << 
   if(fAvgSurfHkEntry<fAvgSurfHkTree->GetEntries())
      fAvgSurfHkTree->GetEntry(fAvgSurfHkEntry);
   else {
      std::cout << "No more entries in avgSurfHkTree" << endl;
      return -1;
   }
   //   std::cout << fAvgSurfHkEntry << "\t" << fAvgSurfPtr->realTime 
   //	     << "\t" << fAvgSurfPtr->ppsNum << std::endl;

   return 0;
}

void MyCorrelator::loadSumTurfTree()
{

  char sumTurfName[FILENAME_MAX];
  sprintf(sumTurfName,"%s/run%d/sumTurfRateFile%d.root",fCurrentBaseDir,
	  fCurrentRun,fCurrentRun);
  fSumTurfRateFile = new TFile(sumTurfName);
  if(!fSumTurfRateFile) {
    cout << "Couldn't open: " << sumTurfName << "\n";
    return;
  }
  fSumTurfRateTree = (TTree*) fSumTurfRateFile->Get("sumTurfRateTree");
  if(!fSumTurfRateTree) {
    cout << "Couldn't get sumTurfRateTree from " << sumTurfName << endl;
    return;
  }
  fSumTurfRateTree->SetBranchAddress("sumturf",&fSumTurfPtr);
  fSumTurfRateEntry=0;

}

int MyCorrelator::loadGpsTrees()
{
  if (fAdu5APatPtr) delete fAdu5APatPtr;
  if (fAdu5aPatTree) delete fAdu5aPatTree;
  
  char gpsName[FILENAME_MAX];
  //if (fCurrentRun==203) sprintf(gpsName,"%s/run%d/gpsEvent%d.root",fCurrentBaseDir,
  //			fCurrentRun,fCurrentRun);
  //else 
  //sprintf(gpsName,"%s/run%d/gpsSsFile%d.root",fCurrentBaseDir,
  //	  fCurrentRun,fCurrentRun);
  sprintf(gpsName,"/home/dailey.110/analysis/testrun/ssDataAndToolsgpsSsFile%d.root",fCurrentRun);
  fGpsFile = TFile::Open(gpsName);
  if(!fGpsFile) {
    cout << "Couldn't open: " << gpsName << "\n";

    if(fCurrentRun==1){
      sprintf(gpsName,"%s/gpsFile.root",fCurrentBaseDir);
      cout<<"gpsName is "<<gpsName<<"\n";
      fGpsFile = TFile::Open(gpsName);
    }
    else{
      sprintf(gpsName,"%s/run%d/gpsFile%d.root",fCurrentBaseDir,
	      fCurrentRun,fCurrentRun);
      fGpsFile = TFile::Open(gpsName);
    }
    if (!fGpsFile){
      cout << "Couldn't open: " << gpsName << "\n";
      return -1;
    }
    else cout<<"Opened gps: "<<gpsName<<"\n";
  }
  else cout<<"Opened gps: "<<gpsName<<"\n";

  fAdu5aPatTree = (TTree*) fGpsFile->Get("adu5PatTree");
  if(!fAdu5aPatTree) {
    cout << "Couldn't get adu5aPatTree\n";
  }
  else {
    fAdu5aPatTree->SetBranchAddress("pat",&fAdu5APatPtr);
  }
  fAdu5aPatEntry=0;
  fAdu5aPatTree->BuildIndex("realTime");
  return 0;
}



int MyCorrelator::getSumTurfEntry() 
{
  if(!fSumTurfRateTree) {
    loadSumTurfTree();
  }
   if(fSumTurfRateEntry<fSumTurfRateTree->GetEntries())
      fSumTurfRateTree->GetEntry(fSumTurfRateEntry);
   else {
      std::cout << "No more entries in sumTurf rate tree" << endl;
      return -1;
   }
   //   std::cout << fSumTurfRateEntry << "\t" << fSumTurfPtr->realTime 
   //	     << "\t" << fSumTurfPtr->ppsNum << std::endl;

   return 0;
}

UInt_t MyCorrelator::getCurrentEvent()
{
  if(fHeadPtr) return fHeadPtr->eventNumber; 
  return 0;
}




void MyCorrelator::GetPatrickEvents(int eventCtrStart, int eventCtrEnd, int thermalFlag, int whichPolarization, std::string current_dir)
{
 
  //  initializeBaseList();
  
  if (!fEventTree) initialize();
  int mcmflag,tdflag,tdreflectionflag,rfonlyflag,calpulserflag,syncslipflag;
  //int ctr_brian;
   int eventNumber;
  int eventPointedFlag, eventTracedFlag;
  if (eventCtrEnd>fEventTree->GetEntries()) eventCtrEnd=fEventTree->GetEntries();
  if (eventCtrStart>fEventTree->GetEntries()) eventCtrStart=fEventTree->GetEntries();
  fHeadTree->GetEntry(eventCtrStart);
  
  // fEventTree->GetEvent(eventCtrStart);
  int firstEvent=(int)fHeadPtr->eventNumber;
  int lastEvent=fHeadPtr->eventNumber+eventCtrEnd-eventCtrStart;
  if (printFlag==1) cout<<"first Event Number: "<<firstEvent<<", last Event Number: "<<lastEvent
			<<", number of Events: "<<lastEvent-firstEvent<<endl;

  
  char filename[150];
  char headername[150];
  char gpsname[150];

  char filename_90[150];
  char headername_90[150];
  char gpsname_90[150];
  cout<<"current_dir is "<<current_dir<<"\n";

  sprintf(filename,"%s/calEventFile%d.root",current_dir.c_str(),fCurrentRun);//change this
  sprintf(headername,"%s/headFile%d.root",current_dir.c_str(),fCurrentRun);
  sprintf(gpsname,"%s/gpsFile%d.root",current_dir.c_str(),fCurrentRun);

  cout<<"outputting to file: "<<filename<<endl;
  TFile *rootfile = new TFile(filename,"RECREATE");
  TTree *eventTree = fEventTree->CloneTree(0);

  TFile *headerfile = new TFile(headername,"RECREATE");
  TTree *headerTree = fHeadTree->CloneTree(0);
  

  TFile *gpsfile = new TFile(gpsname,"RECREATE");
  TTree *gpsTree = fAdu5aPatTree->CloneTree(0);

  sprintf(filename_90,"%s/calEventFile%d_90.root",current_dir.c_str(),fCurrentRun);//change this
  sprintf(headername_90,"%s/headFile%d_90.root",current_dir.c_str(),fCurrentRun);
  sprintf(gpsname_90,"%s/gpsFile%d_90.root",current_dir.c_str(),fCurrentRun);

  TFile *rootfile_90 = new TFile(filename_90,"RECREATE");
  TTree *eventTree_90 = fEventTree->CloneTree(0);

  TFile *headerfile_90 = new TFile(headername_90,"RECREATE");
  TTree *headerTree_90 = fHeadTree->CloneTree(0);
  

  TFile *gpsfile_90 = new TFile(gpsname_90,"RECREATE");
  TTree *gpsTree_90 = fAdu5aPatTree->CloneTree(0);

  int n=0;
  int n_counter=0;

  int ctr10=0;
  int ctr90=0;
  for (int eventctr=eventCtrStart;eventctr<eventCtrEnd;eventctr+=1){//LOOP START
    passed=0;
    percentage=0.;
    eventTracedFlag=0;
    eventPointedFlag=0;
    //get the event
    fHeadTree->GetEntry(eventctr);
    fEventTree->GetEntry(eventctr);
    //fAdu5aPatTree->GetEntry(eventctr);
    eventStartedFlag=startEachEvent(fHeadPtr->eventNumber);
    eventNumber=fHeadPtr->eventNumber;
    
    if(n_counter==10){
      n_counter=0;
      n=(int)10*(gRandom->Rndm());
    }
    /* if (thermalFlag==2){//all non-rf 
      if ((fHeadPtr->trigType&(1<<0))==0 && isCalPulser(fHeadPtr->eventNumber)==0 
	  && isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber)==0 && isTaylor(fHeadPtr->eventNumber)==0
	  && isTaylorReflection(fHeadPtr->eventNumber)==0 && isMainRFCMOn(fHeadPtr->eventNumber)){ 
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isShortTrace(fHeadPtr->eventNumber)==0 ){
    
    // if(taylorFlag ==1 || taylorFlag==2){
       if (isTaylor(fHeadPtr->eventNumber)==1 || isTaylorReflection(fHeadPtr->eventNumber)==1){
	 cout<<"eventNumber is "<<eventNumber<<" is taylor is "<<isTaylor(fHeadPtr->eventNumber)<<" reflect is "<<isTaylorReflection(fHeadPtr->eventNumber)<<"\n";
  
	if (eventEntryGottenFlag!=(int)fHeadPtr->eventNumber) eventEntryGottenFlag=getEventEntry();
	if (isChannelSaturated(fHeadPtr->eventNumber)==-1 && !isSyncSlip(fHeadPtr->eventNumber) 
	    && isMainRFCMOn(fHeadPtr->eventNumber) && isShortTrace(fHeadPtr->eventNumber)==0){
    */	 

    if (fHeadPtr->trigType&(1<<1) || fHeadPtr->trigType&(1<<2) || fHeadPtr->trigType&(1<<3) ||
	!fHeadPtr->trigType&(1<<0)){
      rfonlyflag=0;

    }
    else rfonlyflag=1;
    
    calpulserflag=isCalPulser(fHeadPtr->eventNumber);
    tdflag=isTaylor(fHeadPtr->eventNumber);
    tdreflectionflag=isTaylorReflection(fHeadPtr->eventNumber);
    mcmflag=isMcMBoreholeOrSeaveyFromList(fHeadPtr->eventNumber);
    syncslipflag=isSyncSlip(fHeadPtr->eventNumber);
    
    if (rfonlyflag!=0 && calpulserflag==0 && tdflag==0
	&& mcmflag==0 && tdreflectionflag==0 && !syncslipflag && gpsBadFlag==0){
     
      if(n==n_counter){
	eventTree->Fill();
	headerTree->Fill();
	gpsTree->Fill();
	ctr10++;
      }
      else{
	eventTree_90->Fill();
	headerTree_90->Fill();
	gpsTree_90->Fill();
	ctr90++;
      }
    }//non rf calpulser TD or mcmurdo
 
    // }
    // }
      //}//flag==2
    n_counter++;
  }//loopend
  cout<<"ctr10 is "<<ctr10<<" ctr90 is "<<ctr90<<"\n";
  rootfile = eventTree->GetCurrentFile();
  rootfile->Write();
  
  headerfile = headerTree->GetCurrentFile();
  headerfile->Write();
  
  gpsfile = gpsTree->GetCurrentFile();
  gpsfile->Write();

  rootfile_90 = eventTree_90->GetCurrentFile();
  rootfile_90->Write();
  
  headerfile_90 = headerTree_90->GetCurrentFile();
  headerfile_90->Write();
  
  gpsfile_90 = gpsTree_90->GetCurrentFile();
  gpsfile_90->Write();

  rootfile->Close();
  headerfile->Close();
  gpsfile->Close();

  rootfile_90->Close();
  headerfile_90->Close();
  gpsfile_90->Close();

}//end
//////////////////////////


TGraph* MyCorrelator::fitSineWave(TGraph *g, int bDraw,int ant){
	 
	 // get the FFT and PSD
	 // use the magnitude and phase of the fft as the initial guess
	 FFTtools f;
  	 TGraph *gfft = f.makePowerSpectrumMilliVoltsNanoSeconds(g);
	 FFTWComplex* cfft = f.doFFT(g->GetN(), g->GetY());
	 int maxbin = f.getPeakBin(gfft);
	 double peakfreq, mag;
	 gfft->GetPoint(maxbin, peakfreq, mag);

	 std::cout << "First Guess\t peakfreq: " << peakfreq << "\t sqrtmag:" << sqrt(mag) << "\t N:" << g->GetN() << "\t Mag:" << sqrt(cfft[maxbin].getAbs()) <<"\t phase:" << cfft[maxbin].getPhase() <<  std::endl;
	
	//build a sine wave fitter
	//TF1 *fit = new TF1("sineWave",sineWave,0,100,3);
	 TF1 *fit = new TF1("sineWave",sineWave,20,80,3);
	fit->SetNumberFitPoints(g->GetN());
	//double parms[3] = {peakfreq/1000.,  sqrt(cfft[maxbin].getAbs()),   cfft[maxbin].getPhase() };
	double parms[3] = {peakfreq/1000.,  sqrt(mag),   cfft[maxbin].getPhase() };
	
	//fit->SetParameters(parms);
	fit->SetParameter(0,parms[0]);
	fit->SetParameter(1,parms[1]);
	fit->SetParameter(2,parms[2]);
	//fit->SetParLimits(0, 0, 1.5);
	//fit->SetParLimits(1, 200, 1000);
	//fit->SetParLimits(2, 0, 2*TMath::Pi());
	//fit->FixParameter(0,0.538);
	std::cout << parms[0] << "\t" << parms[1] << "\t" << parms[2] << std::endl;
	g->Fit("sineWave", "QR");
	
	TGraph *gsinefit = new TGraph(g->GetN());
	double* parms2 = fit->GetParameters();
	std::cout <<"sine wave parameters are freq: "<< parms2[0] << "\t mag: " << parms2[1] << "\t phase:" << parms2[2] <<  std::endl;
	//parms2[0] =  peakfreq; parms2[1] = sqrt(cfft[maxbin].getAbs()); parms2[2] =  cfft[maxbin].getPhase();
	for(int ip=0;ip<gsinefit->GetN(); ip++){
		double t, y;
		g->GetPoint(ip, t, y);
		
		gsinefit->SetPoint(ip, t, sineWave(&t, parms2 ) );
	}

	TGraph *g2 = new TGraph(g->GetN());
	// subtract the fit
	
	for( int ip=0; ip<g->GetN(); ip++){
		double t,y, ys;
		g->GetPoint(ip, t, y);
		ys = sineWave(&t, parms2);
		//std::cout << t << "\t" << y << "\t" << ys << "\t" << y-ys << std::endl;
		g2->SetPoint(ip, t, y-ys);
	}
	/*	if(bDraw==1){
	  char printer[256];
	  TCanvas *check = new TCanvas("check","check",800,800);
	  check->Divide(1,2);
	  check->cd(1);
	  g->Draw("AL");
	  
	  check->cd(2);
	  g2->Draw("AL");
	  sprintf(printer,"/home/dailey.110/analysis/filter_study/subtraction/steph_sine_%d.png",ant);
	  check->Print(printer);
	  }*/
	

	delete gfft; 
	delete [] cfft;
	delete fit;
	delete gsinefit;

	return g2;
}

double MyCorrelator::integrateTDPower(TGraph *g){
	Double_t integral = 0.;

	//TGraph g2(g->GetN());
	Double_t oldt, oldv, t, v;
	g->GetPoint(0, oldt, oldv);	
	for(int ip=1; ip<g->GetN(); ip++){
		//double t, v;
		g->GetPoint(ip, t, v);
		//g2.SetPoint(ip, t, v*v);
		//trapezoid rule
		integral += (t-oldt) * 0.5*( v*v + oldv*oldv);
		oldt=t;
		oldv=v;
	}
	//return g2.Integral();
	return integral;
}

////////////////////////////////////////
double MyCorrelator::solveGamma_plus(double theta, double psi, double delta){
  double gamma;
  double sqrt_val;
  double sin_delta = sin(delta);
  double cos_delta = cos(delta);
  double arg = psi-theta;
  sqrt_val = 1-pow(cos_delta/cos(arg),2);
   sqrt_val = sqrt(sqrt_val);
  gamma = sin(2*arg)*(1+sqrt_val)/(2*sin_delta);
  gamma = acos(gamma);
  gamma = theta + gamma;
 
  return gamma; 
}
double MyCorrelator::solveGamma_minus(double theta, double psi, double delta){
  double gamma;
  double sqrt_val;
  double sin_delta = sin(delta);
  double cos_delta = cos(delta);
  double arg = psi-theta;
  sqrt_val = 1-pow(cos_delta/cos(arg),2);
  sqrt_val = sqrt(sqrt_val);
  gamma = sin(2*arg)*(1-sqrt_val)/(2*sin_delta);
  gamma = acos(gamma);

  gamma = theta + gamma;

  return gamma; 
}

void MyCorrelator::GeomMethod(int ant,int pol,vector<double> Freq,vector<double> bandWidth,vector<double> cutFreqs){
  double minFreq;// = Freq - bandWidth;
  double maxFreq;// = Freq + bandWidth;
  double magFFT[2000]={0};
  double magFFT_dB[2000]={0};
  double frequencyArray[2000]={0};
  double deltaT, deltaF;
  int length;
  int newLength;
  double phase_single[2000];
  double phase_single_unshifted[2000];
  double phase_new[2000];
  
  double *times;
  double *volts;
 
  
  double delta;
  
  if(pol==0){
    length = grEv[ant]->GetN();
  }
  else if(pol==1){
    length = grEvHoriz[ant]->GetN();
  }

   if(pol==0){
    times = grEv[ant]->GetX();
    volts = grEv[ant]->GetY();
   
  }
  else if(pol==1){
    times = grEvHoriz[ant]->GetX();
    volts = grEvHoriz[ant]->GetY();
  }
 
 
 
  double average_phase;
  double average_x;
  double average_y;
  int num_avg=0;

  //Get Fourier Transform

  deltaT=times[1]-times[0];
   
  newLength=(length/2)+1;
  deltaF=1/(deltaT*length); //GHz
  deltaF*=1e3; //MHz
  
  FFTWComplex *theFFT=FFTtools::doFFT(length,volts);
  
  //cout<<"Freq is "<<Freq[0]<<"\n";
  for(int i=0;i<newLength;i++) {
    if (i==0) frequencyArray[i]=0;
    if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;//make a freq Array. 
    
    
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      
      magFFT[i] = sqrt(pow(theFFT[i].re,2)+pow(theFFT[i].im,2));
      phase_single[i]=atan2(theFFT[i].im,theFFT[i].re);
      phase_single[i]=phase_single[i];//+2*TMath::Pi()*k_unwrap;
      
      phase_single_unshifted[i]=atan2(theFFT[i].im,theFFT[i].re);
      phase_new[i]=phase_single[i];
      
    }//>200 <1200
      
    
    else{
      phase_single[i]=0.;
      phase_single_unshifted[i]=0.;
      magFFT[i]=-1000;
    }
    
  }//i  
  
  for(int i=0;i<newLength;i++) {
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      magFFT_dB[i]=10*log10(sqrt(magFFT[i]/double(1))/10.);//correct RMS
    }    
    else {
      magFFT_dB[i]=-1000;
      
    }
    
  }
  
    double mean_mag;

    double val1;
    double val2;
    double val_max=TMath::Pi()/2.;
    double val_min=TMath::Pi()/2.;
    
    vector<double> lower_bounds;
    vector<double> upper_bounds;

   
    int one_flag=0;
    int two_flag=0;
    
    for(int j=0;j<(int)Freq.size();j++){
      lower_bounds.clear();
      upper_bounds.clear();
      minFreq = Freq[j]-(bandWidth[j]);//start of notch region
      maxFreq = Freq[j]+(bandWidth[j]);//end of notch region
     
      //set vectors that contain where notched regions are
      for(int k=0;k<(int)cutFreqs.size();k++){
	if(cutFreqs[k] > minFreq && cutFreqs[k]<maxFreq){
	  lower_bounds.push_back(cutFreqs[k]-deltaF);//one bin to left of step (needed to make sure we dont include step in average)
	  upper_bounds.push_back(cutFreqs[k]+deltaF);//one bin to right of step (neede to make sure we dont include step in average)
	}
      }//k

      
      lower_bounds.push_back(1200);
      upper_bounds.push_back(1200);
      int k_lower=0;
      int pass_flag=0;
     
      //start process
      for(int i=0;i<newLength;i++){//lowest freq in range to start of notch

	val_max=TMath::Pi()/2.;
	val_min=TMath::Pi()/2.;

	if(frequencyArray[i]>=minFreq && frequencyArray[i]<=maxFreq){//if inside notched region
	  one_flag=0;
	  two_flag=0;
	  pass_flag=0;
	 
	  if(frequencyArray[i] > upper_bounds[k_lower]+1){//increment to next region
	      k_lower++;
	  }
	 
	  if(frequencyArray[i] < lower_bounds[k_lower]) {//can do all freqs up to bin before step
	    pass_flag=1;
	  }
	
	  if(pass_flag==1){//calculate average for method. Want to use up to 7 bins for average, but cannot include where jump occurs so samples is not constant
	    
	    average_x = theFFT[i].re+theFFT[i-1].re + theFFT[i+1].re;
	    average_y = theFFT[i].im+theFFT[i-1].im + theFFT[i+1].im;
	    mean_mag = magFFT[i-1]+magFFT[i]+magFFT[i+1];
	    one_flag=1;
	    num_avg=3;
	  }

	  pass_flag=0;//being careful about flags?
	 
	  if(one_flag==1){
	    pass_flag=1;
	    for(int k=0;k<(int)cutFreqs.size();k++){//can we expand average out one more bin on either side?
	      if(frequencyArray[i-2] > cutFreqs[k]-1 && frequencyArray[i-2] < cutFreqs[k]+1 ){
		pass_flag=0;
	      }
	      if(frequencyArray[i+2] > cutFreqs[k]-1 && frequencyArray[i+2] < cutFreqs[k]+1 ){
		pass_flag=0;
	      }
	    }//k
	    if(pass_flag==1){
	      //yes we can
	      average_x = average_x + theFFT[i-2].re + theFFT[i+2].re;
	      average_y = average_y + theFFT[i-2].im + theFFT[i+2].im;
	      mean_mag = mean_mag + magFFT[i-2] + magFFT[i+2];
	      num_avg+=2;
	      two_flag=1;
	    }
	  }
	  pass_flag=0;
	  
	  if(two_flag==1){
	    pass_flag=1;
	    for(int k=0;k<(int)cutFreqs.size();k++){//if we expanded to 5 bins for average, can we go to 7 bins?
	      if(frequencyArray[i-3] > cutFreqs[k]-1 && frequencyArray[i-3] < cutFreqs[k]+1){
		pass_flag=0;
	      }
	      if(frequencyArray[i+3] > cutFreqs[k]-1 && frequencyArray[i+3] < cutFreqs[k]+1){
		pass_flag=0;
	      }
	    }//k
	    if(pass_flag==1){
	      //yes
	      average_x = average_x + theFFT[i-3].re + theFFT[i+3].re;
	      average_y = average_y + theFFT[i-3].im + theFFT[i+3].im;
	      mean_mag = mean_mag + magFFT[i-3] + magFFT[i+3];
	      num_avg+=2;
	    }
	  }

	  //get average mag and phase
	  average_x = average_x/num_avg;
	  average_y = average_y/num_avg;
	  mean_mag = mean_mag/num_avg;
	  
	  average_phase = atan2(average_y,average_x);
	 
	  delta = abs(average_phase - phase_single[i]);
	  //do solution of phasor equation

	  //DO I NEED BOTH EQUATIONS WITH MY SIMPLIFICATION?
	  val1 = solveGamma_plus(average_phase, phase_single[i], delta);
	  val2 = solveGamma_minus(average_phase, phase_single[i], delta);
	  if(val1 != val1 && val2 != val2){//somethign went wrong!
	   cout<<"Problem! val1,val2 = "<<val1<<" " <<val2<<" mean, phase_single, delta are "<<average_phase<<" "<<phase_single[i]<<" "<<average_x<<" "<<average_y<<"\n";
	  }
	  /////////SINCE CHANGE OF FORMULA, DONT KNOW IF WE NEED THIS ANYMORE.
	   if(cos(val1-mean_mag) > cos(val2-mean_mag)){
	     val_max += val1;
	     val_min += val2;
	   }
	   else{
	     val_max+=val2;
	     val_min+=val1;
	   }
	  
	   //DOES VAL_MIN==VAL_MAX?
	   if(val1 == val1 && val2 == val2 && one_flag==1){//being careful about NaNs
	     if(magFFT[i]>=mean_mag){
	       phase_new[i]=val_max;
	       
	     }
	     if(magFFT[i]<mean_mag){
	       phase_new[i]=val_min;
	      
	     }
	   }
	   
	}//minFreq
	

      }//i
     
    }//j=Freq.size

    //change the Fourier components
    double x;
    double y;
    for(int i=0;i<newLength;i++){
      if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
	x = magFFT[i]*cos(phase_new[i]);
	y = magFFT[i]*sin(phase_new[i]);

	theFFT[i].re = x;
	theFFT[i].im = y;
	
      }
    }
    //FFT back
    double *filteredVals = FFTtools::doInvFFT(length,theFFT);
    TGraph *grTime = new TGraph(length,times,volts);
    double *times1 = grTime->GetX();
    for(int i=0;i<length;i++){
      volts[i] = filteredVals[i];
    }
    if(pol==0){
    delete grEv[ant];
    grEv[ant] = new TGraph(length,times1,filteredVals);
    }
    if(pol==1){
      delete grEvHoriz[ant];
      grEvHoriz[ant] = new TGraph(length,times1,filteredVals);
    }
    delete [] theFFT;
    delete [] filteredVals;
    delete grTime;
   
  
}//geommethod
void MyCorrelator::DrawFreqDomain(TGraph *gr, int eventnumber,char *namer){
  
  double magFFT[2000];
  double frequencyArray[2000];
  double phaseArray[2000];
  double deltaT, deltaF;
  int length;
  int newLength;
 
  for (int i=0;i<2000;i++){
    magFFT[i]=0.;
    frequencyArray[i]=-1.;
    phaseArray[2000]=0.;
  }
  
  // for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){      
  // if (pol!=0 || ant!=1){//get rid of 2V
  double *Y = gr->GetY();
  double *X = gr->GetX();
  deltaT=X[1]-X[0];
  length=gr->GetN();
  
  newLength=(length/2)+1;
  deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz   

  FFTWComplex *theFFT=FFTtools::doFFT(length,Y);
  
  
  
  for(int i=0;i<newLength;i++) {
   
    if (i==0) frequencyArray[i]=0;
    if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;//make a freq Array. DOES THIS CHANGE FOR EVERY ANT GRAPH??
    
    
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      magFFT[i] = sqrt(pow(theFFT[i].re,2)+pow(theFFT[i].im,2));
      phaseArray[i]=atan2(theFFT[i].im,theFFT[i].re);
     
      
    }
    else{
      phaseArray[i]=0.;
    }
    
  }   
  delete [] theFFT;
  
  for(int i=0;i<newLength;i++) {
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      // magFFT[i]=10*log10(pow(magFFT[i]*deltaF,2)/.01);//correct RMS
      magFFT[i]=10*log10(magFFT[i]/10.);
    }
    
  }

  TH2F *htime_axes = new TH2F("time_axes",";Time(ns);Volts(mV)",10,0,100,10,-30,40);
  TH2F *haxes_mag = new TH2F("axes_phase1",";Frequency (MHz);Magnitude (dB)",10,0,1400,10,-40,60);
  TH2F *haxes_phase = new TH2F("axes_phase",";Frequency (MHz);Phase (Radians)",10,0,1400,10,-4,4);
  
  gr->GetXaxis()->SetTitle("Time (ns)");
  gr->GetYaxis()->SetTitle("Voltage (mV)");
  TGraph *grMag = new TGraph(newLength,frequencyArray,magFFT);
 
  TGraph  *grPhase = new TGraph(newLength,frequencyArray,phaseArray);
  
  TCanvas *c0 = new TCanvas("c0","c0",800,800);
    c0->Divide(1,3);
    c0->cd(1);
    gr->Draw("AL");
   
    c0->cd(2);
    haxes_mag->Draw();
    grMag->Draw("same");
  
    c0->cd(3);
    haxes_phase->Draw();
    grPhase->Draw("same");
    
    char printer2[256];
    sprintf(printer2,"Time_mag_phase_%i_%s.png",eventnumber,namer);
    c0->Print(printer2);
  
}//DrawFreq


void MyCorrelator::GetHealPixMap(){
  //Healpix_Map<double> hpix_map(256,RING);

  long nside = 64;
  long *pixel_num = new long;
  double theta;
  double phi;

  theta = TMath::Pi();
  phi = 0.;
  //vector<double> *vec;
  double *vec;
  // vec[0]= cos(phi)*sin(theta);
  //vec[1]= sin(phi)*sin(theta);
  //vec[2]= cos(theta);

  //ang2vec(theta,phi,vec);
 
  //vec2pix_nest(nside,vec,pixel_num);
  
  //cout<<"pixel_num is "<<pixel_num<<"\n";
  //ang2pix_nest(nside,theta,phi,pixel_num);
  //double trial = (double)(*pixel_num);
  cout<<"pixel_num is "<< *pixel_num <<"\n";
  // ang2pix_nest(nside,theta,phi,pixel_num);
  // pixel_num = hpix_map.ang2pix(vec);
 

  
}

void MyCorrelator::InterpPhase(int ant,int pol,vector<double> Freq,vector<double> bandWidth){
  double minFreq;// = Freq - bandWidth;
  double maxFreq;// = Freq + bandWidth;
  double magFFT[2000]={0};
  double frequencyArray[2000]={0};
  double deltaT, deltaF;
  int length;
  int newLength;
  double phase_single[2000];
  double phase_single_unshifted[2000];
  
  int time_offset;
  double max_voltage=-1000;
  int max_index=-1;
 
  if(pol==0){
    length = grEv[ant]->GetN();
  }
  else if(pol==1){
    length = grEvHoriz[ant]->GetN();
  }
  double SNR;
  double rmsNoise;
  
  double min_time = length/4.;
  double max_time = 3.*length/4.;
  double volt_holder_ant[length];
  double *times;
  double *volts;
  int index_start=-1;
  int index_end=-1;
  int i_pre=-1;
  int i_post=-1;
  double notchwidth=0.;
  double volts_nonshifted[length];
 
  double x;
  double y;

  double max_voltage_pre;
  double max_voltage_post;

   TH2F *haxes_phase = new TH2F("axes_phase",";Frequency (MHz);Phase (Radians)",10,0,1400,10,-4,4);
    char printer[256];


  if(pol==0){
    times = grEv[ant]->GetX();
    volts = grEv[ant]->GetY();
    
  }
  else if(pol==1){
    times = grEvHoriz[ant]->GetX();
    volts = grEvHoriz[ant]->GetY();
  }
 
  TGraph *grTime = new TGraph(length,times,volts);
  grTime->SetLineWidth(3);
  for(int k=0;k<length;k++){//grEv[ant]->GetN()
    volts_nonshifted[k]=volts[k];
    if(k>min_time && k <max_time){ 
      volt_holder_ant[k] = volts[k];
      if(abs(volts[k])>=max_voltage) {
	max_voltage = abs(volts[k]);
	max_voltage_pre = volts[k];
	max_index=k;
      }
    }
    else{
      //volt_holder_ant[k]=0.;
      volt_holder_ant[k]=volts[k];
    }
   
  }
    time_offset=max_index;
    
    //cout<<"ant is "<<ant<<" time_offset "<<time_offset<<"\n";
    for(int i=0;i<length;i++){
      if(i - time_offset <0){
	volts[length + i -time_offset] = volt_holder_ant[i];

      }
      else{
	volts[i-time_offset]=volt_holder_ant[i];
	
      }
    }  

    deltaT=times[1]-times[0];
   
    newLength=(length/2)+1;
    deltaF=1/(deltaT*length); //GHz
    deltaF*=1e3; //MHz
    
    TGraph *grTimer_shifted = new TGraph(length,times,volts);
    grTimer_shifted->SetLineWidth(3);
    grTimer_shifted->GetXaxis()->SetTitle("Time (ns)");
    grTimer_shifted->GetYaxis()->SetTitle("Volts (mV)");
    /*TCanvas *c0a = new TCanvas("c0a","c0a",800,800);
    
     c0a->cd(1);
     grTimer_shifted->Draw("AL");
    
     
     sprintf(printer,"time_shifted_%d.png",ant);
     c0a->Print(printer);
    */
    FFTWComplex *theFFT=FFTtools::doFFT(length,volts);

     for(int i=0;i<newLength;i++) {
      if (i==0) frequencyArray[i]=0;
      if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;//make a freq Array. 
      
      if(frequencyArray[i]>=200 && index_start <0){
	index_start =i;
      }
      if(frequencyArray[i]<=1200){
	index_end =i;
      }
      if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
	
	magFFT[i] = sqrt(pow(theFFT[i].re,2)+pow(theFFT[i].im,2));
	phase_single[i]=atan2(theFFT[i].im,theFFT[i].re);
	phase_single_unshifted[i]=atan2(theFFT[i].im,theFFT[i].re);
	
      }//>200 <1200
      
      
      else{
	phase_single[i]=0.;
	phase_single_unshifted[i]=0.;
	magFFT[i]=-1000;
      }
      
    }//i 

     for(int j=0;j<Freq.size();j++){
     
       minFreq = Freq[j]-bandWidth[j];
       maxFreq = Freq[j]+bandWidth[j];
       
       /////////////Move each part to same mean
       i_pre =-10;
       i_post =-10;
       notchwidth = (int)ceil(2*bandWidth[j]/deltaF);
       
       for(int i=0;i<newLength;i++){
	 if(frequencyArray[i]>=minFreq && frequencyArray[i]<=maxFreq){
	   if(i_pre <0){
	     i_pre = i-1;
	   }
   
	   phase_single[i]=phase_single[i_pre];
	   
	 }

	 
       }//i
       
     }//j


     
     for(int i=0;i<newLength;i++){
       if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
	 x = magFFT[i]*cos(phase_single[i]);
	 y = magFFT[i]*sin(phase_single[i]);
	 
	 theFFT[i].re = x;
	 theFFT[i].im = y;
	 
       }
     }
     double *filteredVals = FFTtools::doInvFFT(length,theFFT);
     
     for(int i=0;i<length;i++){
       volts[i] = filteredVals[i];
     }

     for(int i=0;i<length;i++){
      
      if(i + time_offset>length){
	volts[i-length+time_offset] = filteredVals[i];	
      }
      else{
	volts[i+time_offset]=filteredVals[i];
      }
      
    }
      TH2F *htime_axes = new TH2F("time_axes",";Time(ns);Volts(mV)",10,0,100,10,-30,40);
     TGraph *grtime_shifted = new TGraph(grEv[ant]->GetN(),grEv[ant]->GetX(),grEv[ant]->GetY());
     grtime_shifted->SetLineColor(kRed);
     grtime_shifted->SetLineWidth(2);
     delete [] theFFT;
     delete [] filteredVals;
      grTime->GetXaxis()->SetTitle("Time (ns)");
     grTime->GetYaxis()->SetTitle("Volts (mV)");
     TGraph *grPhase = new TGraph(newLength,frequencyArray,phase_single);
     grPhase->SetLineColor(kRed);
     grPhase->SetLineWidth(2);
     TGraph *grPhase_unshifted = new TGraph(newLength,frequencyArray,phase_single_unshifted);
     grPhase_unshifted->SetLineWidth(3);
     TCanvas *c1 = new TCanvas("c1","c1",800,800);
     c1->Divide(1,2);
     c1->cd(1);
     htime_axes->Draw();
     grTime->Draw("same");
     grtime_shifted->Draw("same");
     c1->cd(2);
     haxes_phase->Draw();
     grPhase_unshifted->Draw("same");
     grPhase->Draw("same");
     
     sprintf(printer,"phase_interp_%d.png",ant);
     c1->Print(printer);
     

     //////////// PLotting!
     /*
     length = grEv[ant]->GetN();
     
	
    double *times1 = grEv[ant]->GetX();
    double * volts1 = grEv[ant]->GetY();
	
     double magFFT_dB[2000];
     TGraph *grTime1 = new TGraph(length,times1,volts1);
     grTime1->SetLineWidth(3);
     grTime1->GetXaxis()->SetTitle("Time (ns)");
     grTime1->GetYaxis()->SetTitle("Volts (mV)");
     TH2F *haxes_phase1 = new TH2F("axes_phase",";Frequency (MHz);Phase (Radians)",10,0,1400,10,-4,4);
     TH2F *haxes_mag = new TH2F("axes_phase1",";Frequency (MHz);Magnitude (dB)",10,0,1400,10,-20,40);
     ////////finished time shift
     
     TGraph *grPhase;
     TGraph *grPhase_unshifted;
    
     deltaT=times1[1]-times1[0];
     cout<<"deltaT is "<<deltaT<<"\n";
     newLength=(length/2)+1;
     deltaF=1/(deltaT*length); //GHz
     deltaF*=1e3; //MHz
     
     FFTWComplex *theFFT1=FFTtools::doFFT(length,volts1);
     
     
     for(int i1=0;i1<newLength;i1++) {
       if (i1==0) frequencyArray[i1]=0;
       if (i1>0) frequencyArray[i1]=frequencyArray[i1-1]+deltaF;//make a freq Array. 
       
      
       if (frequencyArray[i1]>=200 && frequencyArray[i1]<=1200){
	 
	 magFFT[i1] = sqrt(pow(theFFT1[i1].re,2)+pow(theFFT1[i1].im,2));
	 phase_single[i1]=atan2(theFFT1[i1].im,theFFT1[i1].re);
	 cout<<"freq is "<<frequencyArray[i1]<<" phase is "<<phase_single[i1]<<"\n";
	 
       }//>200 <1200
       
       
       else{
	 phase_single[i1]=0.;
	 
	 magFFT[i1]=-1000;
       }
       
     }//i1  
     
     for(int i1=0;i1<newLength;i1++) {
       if (frequencyArray[i1]>=200 && frequencyArray[i1]<=1200){
	 magFFT_dB[i1]=10*log10(sqrt(magFFT[i1]/double(1))/10.);//correct RMS
	 
	 
       }    
       else {
	 magFFT_dB[i1]=-1000;
	 
       }
       
     }
     grTime1->GetXaxis()->SetTitle("Time (ns)");
     grTime1->GetYaxis()->SetTitle("Voltage (mV)");
     TGraph *grMag = new TGraph(newLength,frequencyArray,magFFT_dB);
     grPhase = new TGraph(newLength,frequencyArray,phase_single);
     grPhase->SetLineWidth(3);
     TCanvas *c0a = new TCanvas("c0a","c0a",800,800);
     c0a->Divide(1,3);
     c0a->cd(1);
     grTime1->Draw("AL");
     c0a->cd(2);
     haxes_mag->Draw();
     grMag->Draw("same");
     c0a->cd(3);
     haxes_phase->Draw();
     grPhase->Draw("same");

     sprintf(printer,"Time_mag_phase_interp_%i.png",ant);
     c0a->Print(printer);
     
     */

}//InterpPhase
