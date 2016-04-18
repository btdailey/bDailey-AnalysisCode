/////////////////////////////////////////////////////////////////////////////
/////  MyCorrelator.h        MyCorrelator                                /////
/////                                                                    /////
/////  Description:                                                      /////
/////      The Marvellous ANITA Graphical Interface and Classy Display   /////
/////     (Magic Display) is a simple event display for ANITA. It's not  /////
/////     nearly as all singing and dancing as Ped's display, but it     /////
/////     does provide a convenient interface to the data and it has     /////
/////     some FFT etc. capabilities.                                    ///// 
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#ifndef MAGICDISPLAY_H
#define MAGICDISPLAY_H

//Includes
#include <iomanip>

#include "TChain.h"
#include "TStyle.h"
#include "TGraph.h"
#include "AnitaConventions.h"
#include "AnitaGeomTool.h"
#include "AnitaEventCalibrator.h"
#include "Antarctica.h"
#include "FFTtools.h"
#include "TTimeStamp.h"
#include "RampdemReader.h"
#include <vector>

class TCanvas;
class TPad;
class RawAnitaHeader;
class PrettyAnitaHk;
class RawAnitaEvent;
class UsefulAnitaEvent;
class UsefulAdu5Pat;
class TurfRate;
class SummedTurfRate;
class AveragedSurfHk;
class SurfHk;
class TButton;
class TTreeIndex;
class TFile;
class AcqdStart;
class Adu5Pat;
class CalibratedAnitaEvent;
class TNtuple;

//define some global variables
const double pi=3.14159265;
const double deg2rad=pi/180.;
const double rad2deg=180./pi;
const double C_light=0.299792458;
const int NUM_PHI_SECTORS=16;
const int NUM_ANTS_WITH_NADIRS=40;
const int NUM_ANTS_NO_NADIRS=32;
const int ROUGH_BIN_SIZE=2;//2
const int NUM_BINS_ROUGH_THETA=int(180/ROUGH_BIN_SIZE);
const int NUM_BINS_ROUGH_PHI=int(360/ROUGH_BIN_SIZE);
const double FINE_BIN_SIZE=0.3;
const int NUM_BINS_FOR_TRIG_ARRAYS=360*5;
const int NUM_BINS_FINE_THETA=75;
const int NUM_BINS_FINE_PHI=150;
const int LOWEST_THETA=60;
const int HIGHEST_THETA=25;
const int LOWEST_THETA_COS=30;//30
const int HIGHEST_THETA_COS =115;//115;
const int NUM_BINS_ROUGH_COS =(HIGHEST_THETA_COS-LOWEST_THETA_COS)/ROUGH_BIN_SIZE;

//const int NUM_BINS_ROUGH_THETA = NUM_BINS_ROUGH_COS;
const double cosHighestTheta = cos(HIGHEST_THETA_COS*deg2rad);
const double cosLowestTheta = cos(LOWEST_THETA_COS*deg2rad);
const double step_size_cos_rough = (cosLowestTheta - cosHighestTheta)/ NUM_BINS_ROUGH_COS;
const int NUM_BINS_PLOT_COS = 65;
const int NUM_CHANNELS=90;
const int NUM_DEGREES_OFF_CENTER=75;
const double ADU5_FORE_PHI=22.5;
const double PHI_SECTOR_ANGLE=22.5;
const double phiStartOffset=NUM_BINS_FINE_PHI*FINE_BIN_SIZE/2;
const double thetaStartOffset=NUM_BINS_FINE_THETA*FINE_BIN_SIZE/2;
const double latTaylor=77.8803;
const double lonTaylor=158.45925;
const double heightTaylor=2260-97.;//meters
const double latWilly=77.861;
const double lonWilly=167.056;
const double heightWilly=-49.;//meters
const int writeData=1;
const int printFlag=0;
const int groupDelayFlag=1;
const int newnotchflag=1;//1 for Brian, 0 for abby
//const int phase_flag=1;//0=old phase, 1 = new (random) phase, 2 = interp phase;

//const int notchFilterFlag =2;//0 for no-fill,1 for rayleigh, 2 for wiener, 3 for interpolated, 4 for interpolated+noise
//int notchFilterFlag=0;
const int debug_flag=0;
const float xSize=750;
const float ySize=625;
const float xOffset=xSize/2.;
const float yOffset=ySize/2.;
const float scale=271.5/2.19496e+06;

class MyCorrelator 
{
 public:
  
  MyCorrelator(char *baseDir, int run, WaveCalType::WaveCalType_t calType=WaveCalType::kVoltageTime);
  MyCorrelator();
  ~MyCorrelator();

  TStyle* RootStyle();
  //Control Panel Functions
  Int_t getCurrentRun()
  {return fCurrentRun;}
  UInt_t getCurrentEvent();


  void closeCurrentRun();

  //initialization stuff
  int startEachEvent(int myEventNumber);
  void initialize();
  void initializeAntarctica();
  void initializeBaseList();
  void setupCosSinTanArray();
  void readInFirstLastEventsOfEachRun();
  void getPositionsOfEachAntenna();
  void getGraphsThisEvent( int windowWaveformFlag, double &snrPeak, double &maxSignalPeak, int &peakAnt);
  TGraph *Resizeplots(TGraph *grWave);
 
  void processEventsFromAList(int drawMaps, int rfOnlyFlag, int whichMcMFlag, int whichPolarization);
  void loopOverEvents(int eventCtrStart, int eventCtrEnd, int drawMaps, int taylorFlag, int thermalFlag, int whichPolarization, std::string current_dir, int filter_number, int phase_number);
  void loopOverThirdEvents(int whichThird, int drawMaps, int taylorFlag, int thermalFlag, int whichPolarization, std::string current_dir);
  void readBaselineFFTs();
  void GetBaselineperPhi(int pol, double *baseline, int nantennasToUse,std::vector<int>& whichAntennasToUse);
  //general tools
  double getMaximum(int n, double *array, int &index);
  double getRMS(TGraph *gr);
  int getPeakAntenna(int myEventNumber, int nantennas);
  double getPeak2Peak(TGraph *gr);
  double getSNR(TGraph *gr, double &rmsNoise);
  double getRMSOfRange(TGraph *gr, double xLow, double xHigh);
  double getPeakHilbert(TGraph *gr);
  TGraph *simpleNotchFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq);
  TGraph *complicatedNotchFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq, int ant, int pol,double *baseY);
  TGraph *nofillNotchFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq);
  TGraph *wienerFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq, int ant, int pol,double *baseY);
  TGraph *interpolatedFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq);
  void InterpPhase(int ant,int pol,std::vector<double> Freq,std::vector<double> bandWidth);
  void GeomMethod(int ant,int pol,std::vector<double> Freq,std::vector<double> bandWidth, std::vector<double> cutFreqs);
  double solveGamma_plus(double theta, double psi, double delta);
  double solveGamma_minus(double theta, double psi, double delta);
  double getCoherentPolarization(TGraph *grV, TGraph *grH, double &polarizationFraction);
  
  //flag stuff
  int isPayloadBlast(int myEventNumber);
  int isChannelSaturated(int myEventNumber);
  int isVarnerEvent(int myEventNumber);
  int isVarnerEvent2(int myEventNumber);
  int isSyncSlip(int myEventNumber);
  int isShortTrace(int eventNumber);
  int isNadirRFCMOn(int eventNumber);
  int isMainRFCMOn(int eventNumber);
  int isBigEnoughPeakToPeak(int eventNumber);
  int isDCOffsetLarge(int eventNumber);
  int isTaylor(int eventNumber);
  int isTaylorReflection(int eventNumber);
  int isTaylor(int eventNumber, double &distance);
  int isTaylor(int eventNumber, double &distance, double &deltaT);
  int isMcMBoreholeByTime(int eventNumber);
  int isMcMSeaveyByTime(int eventNumber);
  int isMcMBoreholeOrSeaveyFromList(int eventNumber);
  int isCalPulser(int eventNumber);

  //distance and coordinates
  void LatLonAlt2xyz(double lat, double lon, double alt, double &x, double &y, double &z);
  void xyz2LatLonAlt(double &lat, double &lon, double &alt, double x, double y, double z);
  double getMcMDistance(int eventNumber);
  double findDistanceBetweenTwoThings(double lat1, double lon1, double alt1, double lat2, double lon2, double alt2);
  void getTDxyz(double &x, double &y, double &z);
  void getWillyxyz(double &x, double &y, double &z);
  void getRelXYFromLatLong(double latitude, double longitude,double &x, double &y);

  //specific tools
  void getClosestNAntennas(int nantennasToUse, double peakPhi, std::vector<int>& whichAntennasToUse, int nadirFlag);
  void getGroupsofAntennas(int nantennasToUse, int nadirFlag);
  void getTriggeredPhi(RawAnitaHeader *hdPtr,int triggeredPhi[NUM_PHI_SECTORS]);
  void getTriggeredL2Phi(RawAnitaHeader *hdPtr,int triggeredPhi[NUM_PHI_SECTORS]);
  void getTriggeredAnt(int triggeredPhi[NUM_PHI_SECTORS],int triggeredAnt[NUM_ANTS_WITH_NADIRS]);
  void getTriggeredAntpm1PhiSector(int triggeredPhi[NUM_PHI_SECTORS],int triggeredAnt[NUM_ANTS_WITH_NADIRS]);
  void getTriggeredAntOf3PhiSectors(int triggeredPhi[NUM_PHI_SECTORS],int triggeredAnt[NUM_ANTS_WITH_NADIRS]);
  int allowedPhisPairOfAntennas(double &lowerAngle, double &higherAngle, double &centerTheta1, 
				double &centerTheta2, double &centerPhi1, double &centerPhi2, int ant1, int ant2);
  int isPhiMaskingOn(int eventNumber, int phiMaskedArray[NUM_PHI_SECTORS]);
  void getPhiMaskedAnts(int phiMaskArray[NUM_PHI_SECTORS],int phiMaskedAnts[NUM_ANTS_WITH_NADIRS]);
  void getPhiMaskedAntspm1(int phiMaskArray[NUM_PHI_SECTORS],int phiMaskedAnts[NUM_ANTS_WITH_NADIRS]);
  int getBestPhiSector(double phiWaveRadians);
  float getHWTriggerAngle(int eventNumber);
  int isPeakPhiTriggeredOrMasked(RawAnitaHeader *hdPtr, double peakPhiInterp);
  int isPeakPhiMasked(RawAnitaHeader *hdPtr, double peakPhiInterp);
  int isPeakTriggeredpm1(RawAnitaHeader *hdPtr, double peakPhiInterp);
  double getHeadingOfEvent(int eventNumber, double peakPhiMapDegrees);
  void GetPatrickEvents(int eventCtrStart, int eventCtrEnd, int thermalFlag, int whichPolarization, std::string current_dir);
   //filtering
  void adaptiveFilterCoherent(TGraph *grCoherent, int pol, double dBCut, int nfreq, 
			      double *frequencies, int drawFlag, double bandWidth);//, float &mean);
  void adaptiveFilter(int pol, double dBCut, int nfreq, double *frequencies, int drawFlag, double bandWidth, float &mean_freq, double *freqArray, double *FFTarray, double *baseX, double *baseY, int baseN,int myEventNumber);
  void adaptiveFilterPartialPayload(int pol, double dBCut, int nfreq, 
				    double *frequencies,double *bandwidth,double *magPeaks, int drawFlag, double bandWidth, double peakPhi,int nAntennasToUse,
				    std::vector<int>& whichAntennasToUse,float &mean_freq, int myEventNumber, int nadirFlag, int antenna_groups);
  void GetFrequenciestoCut(int antenna,std::vector< std::vector<double> > &antennaFreq,std::vector< std::vector<double> > &bandwidth,
			   std::vector< std::vector<double> > &PeakMag, std::vector<double>& uniquefreqs,std::vector<double>& uniquebandwidth, int nfreq,std::vector<double>& uniquePhase,std::vector<double>& uniquePhase_bandwidth);
  void GetSatelliteFrequenciestoCut(int antenna,std::vector<double> &uniquefreqs,std::vector<double> &uniquebandwidth, double &satellite_freq,double &satelite_bandwidth,
				    double &satelite_freq2, double &satelite_bandwidth2, int &satellite_flag);

  void GetFFTandBaseline(int nantennasToUse,std::vector<int>& whichAntennasToUse,int pol, double *frequencyArray,double *baseX, double *baseY2, double *magFFT2, int nfreq, int dBCut);
void GetFFTandBaseline_old(int nantennasToUse,std::vector<int>& whichAntennasToUse,int pol, double *frequencyArray,double *baseX, double *baseY2, double *magFFT2); 
 void applyAdaptiveFilter(double centerFrequency, double bandWidth, int polFlag);
 
  void applyAdaptiveFilter_singleAnt(double centerFrequency, double bandWidth, int polFlag,int ant, double *baseY);
  void applySatelliteFilter(double centerFrequency, double bandWidth,int ant,int pol,double *baseY);
 
  //pointing stuff
  Double_t getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave);//in degrees
  Double_t getGroupDelay(Double_t phiToAntBoresight, Double_t thetaWave);//in radians
  void doCorrelationMap(int myEventNumber, double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], 
			int onlyTriggered, int nadirFlag, int whichPolarization);
  int pointThisEvent(int eventNumber, int drawMaps, TNtuple *ndata, TNtuple *ndata2,  TNtuple *ndata3, TNtuple *ndata4, 
		     double &peakThetaFinal, double &peakPhiFinal, int whichPolarization,int &xCorPassFlag,
		     int &ratioOfPeaksPassFlag, int &elevationAnglePassFlag, int &peakCrossCorrFlag,
		     int &polFractionFlag, int &peakHilbertFlag, int &triggerflag, double &finaltheta);
  int traceBackTo0Altitude(int eventNumber, double thetaWave, double phiWave, 
			    double &sourceLon, double &sourceLat);
  int traceBackToContinent_Brian(int eventNumber, double thetaWave, double phiWave, double &sourceLon, double &sourceLat, double &sourceAlt);
  int checkIfNearAnyBase(int eventNumber,double peakThetaFinal, double peakPhiFinal, double sourceLon, 
		      double sourceLat, double sourceAlt, std::string &baseName);
  int checkIfNearSpecificBase(int eventNumber,double peakThetaFinal, double peakPhiFinal, double sourceLon, 
			      double sourceLat, double sourceAlt, std::string baseName, int baseIndex);
  void drawEventOnContinent(int eventNumber,double sourceLon,double sourceLat);
  void drawCorrelationMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], int myEventNumber,int trial);
  TGraph *makeCoherentlySummedWaveform(int myEventNumber, double peakTheta, double peakPhi,
				       int nadirFlag, int polFlag, int nantennas, int drawFlag);
  TGraph *makeCoherentlySummedDeconvolvedWaveform(int myEventNumber, double peakTheta, 
						  double peakPhi, int nadirFlag, int polFlag, 
						  int nantennas, int drawFlag);
  void findPeakOfMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], double &peakVal, 
		     double &peakTheta, double &peakPhi);
  void drawRefinedMap(double mapCorValRefined[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI], double peakTheta, double peakPhi, int myEventNumber);
  void doRefinedMap(int myEventNumber, double mapCorValRefined[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI], 
		    double peakTheta, double peakPhi, int onlyTriggered, int nadirFlag, int whichPolarization);
  void findPeakOfRefinedMap(double mapCorVal[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI], double &peakVal, 
			    int &peakThetaBin, int &peakPhiBin);
  void doInterpolationPeakFinding(double mapCorVal[NUM_BINS_FINE_THETA][NUM_BINS_FINE_PHI], double &peakVal, 
				  double &peakTheta, double &peakPhi, double &FWHMTheta, double &FWHMPhi,
				  double peakThetaRough, double peakPhiRough, int drawFlag);
  double getSNROfMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI]);
  void findPeakOfMapBin(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], double &peakVal, 
			       int &peakThetaBin, int &peakPhiBin);
  void findSecondaryPeakOfMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], double &secondaryPeakVal, 
			      double &secondaryPeakTheta, double &secondaryPeakPhi);
  void findTertiaryPeakOfMap(double mapCorVal[NUM_BINS_ROUGH_THETA][NUM_BINS_ROUGH_PHI], 
			     double &tertiaryPeakVal, 
			     double &tertiaryPeakTheta, double &tertiaryPeakPhi);
  void backProjectEvent(int eventNumber, double sourceLon, double sourceLat, 
			double sourceHeight, double &thetaWaveProj,
			double &phiWaveProj);
  TGraph *deconvolveWaveform(TGraph *grCoherent, int polFlag, int drawFlag);
  TGraph *deconvolveWaveformUsingStephens(TGraph *grCoherent, int polFlag, int drawFlag);
  TGraph *deconvolveWaveformUsingRyans(TGraph *grCoherent, int polFlag, int drawFlag);
  //int getLabChip(RawAnitaEvent *fRawEventPtr,int chanIndex);
  //nonmaincode stuff
  void getCRPolContent();
  void calcANITAIIUsingStephens();
  void calcSystemImpulseResponse();
  void firstLastEachRun(int &firstEvent, int &lastEvent, int &numberOfEvents);
  void taylorPhiMaskandTriggerBits(int numEvents, double isbesttrigpm1_denom[NUM_PHI_SECTORS], 
				   double isbesttrig_denom[NUM_PHI_SECTORS],
				   double isbesttrigpm1orphimask_denom[NUM_PHI_SECTORS], 
				   double isbesttrigorphimask_denom[NUM_PHI_SECTORS],
				   double isbesttrigpm1_num[NUM_PHI_SECTORS],
				   double isbesttrig_num[NUM_PHI_SECTORS],
				   double isbesttrigpm1orphimask_num[NUM_PHI_SECTORS], 
				   double isbesttrigorphimask_num[NUM_PHI_SECTORS]);
  void taylorEfficiency(int numEvents);
  void trackTaylor(int numEvents);
  void trackTaylorForEventInsertion(int numEvents);
  void createBaseline(int numEvents);
  void createSigma(int numEvents,int antctr);
  void createdBCut(int numEvents,int antctr);
  void calcGPSFlags();
  void calcVarnerFlags();
  void loopOverStartingEvents();
  void DrawFreqDomain(TGraph *gr, int eventnumber, char *namer);
  void GetHealPixMap();
  
 TGraph* fitSineWave(TGraph *g, int bDraw,int ant);
	 
 double integrateTDPower(TGraph *g);




 // Double_t fitf(Double_t *x, Double_t *par);
  //ryan's stuff
  int loadEventTree();
  int getEventEntry();
  int getHeaderEntry();
  void loadTurfTree();
  int getTurfEntry();
  void loadSurfTree();
  void loadAcqdTree();
  int getSurfEntry();
  int getSurfIDEntry();  
  void loadAvgSurfTree();
  int getAvgSurfEntry();
  void loadSumTurfTree();
  int getSumTurfEntry();
  int loadGpsTrees();

  //Instance generator
  static MyCorrelator*  Instance();
  
  TFile *fHeadFile;
  TFile *fEventFile;
  TFile *fTurfRateFile;
  TFile *fSumTurfRateFile;
  TFile *fSurfHkFile;
  TFile *fAvgSurfHkFile;
  TFile *fAcqdFile;
  TFile *fGpsFile; 
  
  //Here are the data managers
  TTree *fHeadTree;
  TChain *fEventTree;
  TTree *tgroundPulser;
  TTree *fPrettyHkTree;
  TTree *fTurfRateTree;
  TTree *fSurfHkTree;
  TTree *fSumTurfRateTree;
  TTree *fAvgSurfHkTree;
  TChain *fAcqTree;
  TTree *fAdu5aPatTree; 
  //And some useful info to keep track of what is where
  Long64_t fEventEntry;
  Long64_t fPrettyHkEntry;
  Long64_t fTurfRateEntry;
  Long64_t fSurfHkEntry;
  Long64_t fSurfIDEntry;
  Long64_t fSumTurfRateEntry;
  Long64_t fAvgSurfHkEntry;
  Long64_t fAdu5aPatEntry; 
  //And something to help with the indexing
  TTreeIndex *fHeadIndex;

  Int_t  fUseCalibratedEventFile; ///< Flag to determine whether or not to use TTrees of lovely CalibratedAnitaEvent objects
  Int_t  fUseEventFile;
  UInt_t fCurrentRun;
  Char_t fCurrentBaseDir[180];

  //int getLabChip(RawAnitaEvent *fRawEventPtr,int chanIndex);

 protected:
   static MyCorrelator *fgInstance;  
   // protect against multiple instances

 private:
   //   MyCorrelatorOption::MyCorrelatorOption_t fMainOption;
   
   RawAnitaHeader *fHeadPtr;
   PrettyAnitaHk *fHkPtr;
   AnitaEventCalibrator *myCally;
   AnitaGeomTool *fUPGeomTool;
   FFTtools *fFFTTools;
   Antarctica *fAntarctica;
   RawAnitaEvent *fRawEventPtr;
   UsefulAnitaEvent *fUsefulEventPtr;
   CalibratedAnitaEvent *fCalEventPtr; ///< Pointer to the raw event.
   UsefulAdu5Pat *fUsefulAdu5Ptr;
   TurfRate *fTurfPtr;
   SurfHk *fSurfPtr;
   SummedTurfRate *fSumTurfPtr;
   AveragedSurfHk *fAvgSurfPtr;
   AcqdStart *fAcqPtr;
   Adu5Pat *fAdu5APatPtr;
   RampdemReader *fRampdemReader;

   WaveCalType::WaveCalType_t fCalType; ///< The waveform calibration type.
};


#endif //MAGICDISPLAY_H
