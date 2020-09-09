#include "TString.h"
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TChain.h"
#include "TMath.h"

#include "classes/mTowerHit.h"
#include "classes/mTowerCluster.h"
#include "classes/mTowerEvent.h"
#include "classes/mTowerChip.h"

#include <iostream>
#include <fstream>
#include <string>


//conversion tables
const std::map< Int_t, Int_t > chipid2lane_lut = {
  { 0,40},{ 1,39},{ 2,42},{ 3,41},{ 4,44},{ 5,43},{ 6,46},{ 7,45},
  { 8,48},{ 9,47},{10,50},{11,49},{12,52},{13,51},{14,54},{15,53},
  {16,38},{17,55},{18,36},{19,37},{20,32},{21,35},{22,34},{23,33},
  {24,64},{25,63},{26,66},{27,65},{28,68},{29,67},{30,70},{31,69},
  {32,72},{33,71},{34,74},{35,73},{36,76},{37,75},{38,78},{39,77},
  {40,62},{41,79},{42,60},{43,61},{44,56},{45,59},{46,58},{47,57}
};

const std::map< Int_t, Int_t > lane2chipid_lut = {
  {40, 0},{39, 1},{42, 2},{41, 3},{44, 4},{43, 5},{46, 6},{45, 7},
  {48, 8},{47, 9},{50,10},{49,11},{52,12},{51,13},{54,14},{53,15},
  {38,16},{55,17},{36,18},{37,19},{32,20},{35,21},{34,22},{33,23},
  {64,24},{63,25},{66,26},{65,27},{68,28},{67,29},{70,30},{69,31},
  {72,32},{71,33},{74,34},{73,35},{76,36},{75,37},{78,38},{77,39},
  {62,40},{79,41},{60,42},{61,43},{56,44},{59,45},{58,46},{57,47}
};

const std::map< Int_t, Int_t > chipid2layer_lut = {
  { 0,22},{ 1,22},{ 2,20},{ 3,20},{ 4,18},{ 5,18},{ 6,16},{ 7,16},
  { 8,14},{ 9,14},{10,12},{11,12},{12,10},{13,10},{14, 8},{15, 8},
  {16, 6},{17, 6},{18, 4},{19, 4},{20, 0},{21, 0},{22, 2},{23, 2},
  {24,23},{25,23},{26,21},{27,21},{28,19},{29,19},{30,17},{31,17},
  {32,15},{33,15},{34,13},{35,13},{36,11},{37,11},{38, 9},{39, 9},
  {40, 7},{41, 7},{42, 5},{43, 5},{44, 1},{45, 1},{46, 3},{47, 3}
};

const std::map< Int_t, Int_t > lane2layer_lut = {
  {40,22},{39,22},{42,20},{41,20},{44,18},{43,18},{46,16},{45,16},
  {48,14},{47,14},{50,12},{49,12},{52,10},{51,10},{54, 8},{53, 8},
  {38, 6},{55, 6},{36, 4},{37, 4},{32, 0},{35, 0},{34, 2},{33, 2},
  {64,23},{63,23},{66,21},{65,21},{68,19},{67,19},{70,17},{69,17},
  {72,15},{71,15},{74,13},{73,13},{76,11},{75,11},{78, 9},{77, 9},
  {62, 7},{79, 7},{60, 5},{61, 5},{56, 1},{59, 1},{58, 3},{57, 3}
};

const std::map<int,bool> layer2isInv_lut = {
  { 0, kFALSE}, { 1, kTRUE}, { 2, kFALSE}, { 3, kTRUE}, 
  { 4, kFALSE}, { 5, kTRUE}, { 6, kFALSE}, { 7, kTRUE}, 
  { 8, kFALSE}, { 9, kTRUE}, {10, kFALSE}, {11, kTRUE}, 
  {12, kFALSE}, {13, kTRUE}, {14, kFALSE}, {15, kTRUE}, 
  {16, kFALSE}, {17, kTRUE}, {18, kFALSE}, {19, kTRUE}, 
  {20, kFALSE}, {21, kTRUE}, {22, kFALSE}, {23, kTRUE} 
};

bool IsLeftChip(int lane){
  int layerNr = lane2layer_lut.at(lane);
  bool isInv  = layer2isInv_lut.at(layerNr);
  int chipid = lane2chipid_lut.at(lane); 
  bool isOdd;
  if (chipid%2 == 1){ isOdd = kTRUE;}
  else {isOdd = kFALSE;}
  bool isLeft = (bool)(isOdd != isInv);
  return isLeft;
}

//for hit maps
int lane2padhitmap(int lane){
  int layerNr = lane2layer_lut.at(lane);
  bool isLeft = IsLeftChip(lane);
  
  int padid;
  if (layerNr < 6) {padid = layerNr+1;}
  else if (layerNr < 12) {padid = layerNr+7;}
  else if (layerNr < 18) {padid = layerNr+13;}
  else {padid = layerNr+19;}
  if (isLeft) {padid += 6;}
 
  return padid; 
}

//for hit maps
int lane2padoccupancy(int lane){
  int layerNr = lane2layer_lut.at(lane);
  int padid = layerNr+1;
  return padid; 
}

void Analyse_mTower(int run)
{
  //-------------------------------------------------------
  //macro to read the mTower TB data
  //using event classes: mTowerEvent, mTowerHit, mTowerChip, mTowerCluster
  //detailed code description can be found in the bachelor thesis of Aart van Bochove
  //authors: N. van der Kolk, A. van Bochove
  //-------------------------------------------------------
  
  //-------------------------------------------------------
  //Change this few variables to what applies for your analysis
  //-------------------------------------------------------
  
  TString fileLocationOutputTree = "./";
  TString fileLocation = "/eos/project/m/mtower/Data/Data_TB_February_2020/mTower_Data_DESY_Feb_2020_raw1/"; //The location of the selected data
  //TString fileLocation = "../";
  //TString fileLocation = "/data2/data/TB_mTower/mTower_Data_DESY_Feb_2020_raw1/"; //The location of the raw data
  TString maskingFileLocation = "./masking/"; //The location of the masking .txt files
 
  bool UND = false; //is true if an earlier made TTree is used, so if not using the original data.
  bool testing = false; //Testing outputfile, run over small amount of events.

  bool DB = false; //debug print statements
  bool HP = false; //find hot pixels
  bool CT = true; //Create TTree
  bool C1 = false; //Criterion 1. Only events where the accepted cluster (C2) is within a certain area of the first layer are accepted for analysis. NEEDS C2.
  bool C1ES = false; //C1, but even events with just single hits outside area are rejected.
  bool C2 = true; //Criterion 2. Every cluster in the first layer with hits behind it in layer 2 is accepted. Only events with 1 accepted cluster are accepted for analysis.
  bool C3 = false; //Criterion 3. Every event with clusters with clustersize > 1, which is ignored by C2, is rejected. NEEDS C2
  bool C4 = false; //Criterion 4. Only events where there are no hits outside of a certain area in the second layer centrated around the accepted cluster (C2) in the first layer are accepted for analysis. NEEDS C2
  bool C5 = false; //Criterion 5. The fraction (nHits in border of chip)/(nHits in a circle centrated around the accepted cluster (C2)) for the 3th to the nth layer has to be higher then a threshold, otherwise the event is rejected. NEEDS C2

  int nPixelSideC1 = 500; //length of the sides of the area for criterion C1. Actual row side is 1 pixel bigger if nPixelsGap is exactly odd to maintain symmetry.
  int deltaRowAreaC1 = 0; //nPixels deviation from row for criterion C1, with respect to center between 2 chips (may vary by 0.5 pixel). The area contains nPixelSideEC1/2 rows up and down from this point.
  int nPixelRadiusC2 = 10; //The search area behind a cluster is a circle centered around the average coordinate of the cluster. nPixelRadiusC2 is the radius of that circle.
  int nPixelRadiusC4 = 120; //Radius of circle for criterion C4
  int nPixelRadiusC5 = 80; //Radius of circle for criterion C5
  int nPixelBorderC5 = 170; //Width of border for criterion C5
  int nLayersC5 = 8; //Number of layers C5 takes into account
  double fractionC5 = 0.15; //Maximum allowed fraction hitsInBorder/hitsInCircle for criterion C5

  int nParts = 1; //In how many parts is this run done?
  int part = 1; //Which part is this?

  const int maxNChips = 48; //48 chips is the maximum size of mTower
  const int laneOffset = 32; 
  std::cout << "*** lane offset used is "<<laneOffset<<endl;
  //for the mTower cosmics test and 2019 and 2020 test beam the offset is 32, for the testbeam data of 2018 there is no offset
  const int rowsPerChip = 512; //512 for ALPIDE chip
  const int columnsPerChip = 1024; //1024 for ALPIDE chip
  double nPixelsGap = 5; //Width of gap between chips in units of pixels
  const int laneNumber[6] = {0,3,27,24,2,1}; //corresponds 'lanecode' with laneNumber-laneOffset. Element i is a chip in the (i/2)th layer. Lane 32, 35 in layer 0, 59, 56 in l1, 34, 33 in l2

  //-------------------------------------------------------
  //Setting up some variables automatically
  //-------------------------------------------------------

  if (!(C2)) {C1 = false; C3 = false; C4 = false; C5 = false;} //If C2 is not applied, all of these criteria cannot be applied
  TString beamEnergy; //for masking inputfile
  if (run == 1309 || run == 1310 || run == 1346 || run == 1375 || run == 1376) beamEnergy = "5R8";
  else if (run == 1250 || run == 1261 || run == 1308 || run == 1333 || run == 1339 || run == 1413) beamEnergy = "5R0";
  else if (run == 1257 || run == 1272 || run == 1274 || run == 1275 || run == 1338 || run == 1345) beamEnergy = "4R0";
  else if (run == 1335 || run == 1341 || run == 1262) beamEnergy = "3R0";
  else if (run == 1276 || run == 1337 || run == 1344) beamEnergy = "2R0";
  else if (run == 1336 || run == 1343 || run == 1263) beamEnergy = "1R0";
  else beamEnergy = "0";

  //later the rows and columns will be translated to an absolute coordinate system with x and y. This are the maximum values of those coordinates
  int maxX = 2*rowsPerChip + nPixelsGap - 1;
  int maxY = columnsPerChip - 1;

  //minimal and maximal values of acceptance area of criterion C1
  double minXC1 = (nPixelsGap-1)/2.+deltaRowAreaC1-nPixelSideC1/2+rowsPerChip;
  double maxXC1 = (nPixelsGap-1)/2.+deltaRowAreaC1+nPixelSideC1/2+rowsPerChip;
  double minYC1 = (columnsPerChip-nPixelSideC1)/2-.5;
  double maxYC1 = (columnsPerChip+nPixelSideC1)/2-.5;

  //set plain style for histograms
  gROOT->SetStyle("Plain");

  //-------------------------------------------------------
  //extra masking
  //-------------------------------------------------------

  TH3F* hMaskPtn = new TH3F("hMaskPtn","Hit mask",maxNChips, 0, maxNChips, columnsPerChip, 0, columnsPerChip, rowsPerChip, 0, rowsPerChip); //Histogram with pixels to mask: lane, column, row

  if (beamEnergy != "0")
    {
      ifstream in;
      TString maskingFileName = "%smTower_Data_";
      maskingFileName += beamEnergy;
      maskingFileName += "_GeV_mask.txt";
      in.open(Form(maskingFileName,maskingFileLocation.Data()));
      
      Float_t chip_id, nr_lane, hot_pixel_column, hot_pixel_row, pixel_entry;
      Double_t difference, average, std_dev;
      Int_t nlines = 0;
      
      while (1)
	{
	  in >> chip_id >> nr_lane >> hot_pixel_column >> hot_pixel_row >> pixel_entry >> difference >> average >> std_dev;
	  if (!in.good()) break;
	  hMaskPtn->Fill(nr_lane-laneOffset,hot_pixel_column,hot_pixel_row);
	  nlines++;
	}
      
      in.close();
    }

  //-------------------------------------------------------
  //create outputfile
  //-------------------------------------------------------

  TString baseName = "Run_";
  TString outputFileName = "";

  if (testing)
    {
      outputFileName = "results_";
      outputFileName += baseName;
      outputFileName += run;
      outputFileName += "_test.root";
    }
  else
    {
      outputFileName = "results_";
      outputFileName += baseName;
      outputFileName += run;
      if (C1 || C1ES)
	{
	  outputFileName += "_C1";
	  if (C1ES) outputFileName += "ES";
	  outputFileName += "_";
	  outputFileName += nPixelSideC1;
	  outputFileName += "cX";
	  if (nPixelsGap == int(nPixelsGap) && int(nPixelsGap)%2) outputFileName += nPixelSideC1+1; //if nPixelsGap is exactly odd
	  else outputFileName += nPixelSideC1;
	  outputFileName += "r_dr";
	  outputFileName += deltaRowAreaC1;
	}
      if (C2)
	{
	  outputFileName += "_C2_";
	  outputFileName += nPixelRadiusC2;
	  outputFileName += "p";
	}
      if (C3) outputFileName += "_C3";
      if (C4)
	{
	  outputFileName += "_C4_";
	  outputFileName += nPixelRadiusC4;
	  outputFileName += "p";
	}
      if (C5)
	{
	  outputFileName += "_C5_";
	  outputFileName += nPixelRadiusC5;
	  outputFileName += "pr_";
	  outputFileName += nPixelBorderC5;
	  outputFileName += "pb_";
	  outputFileName += nLayersC5;
	  outputFileName += "l_";
	  outputFileName += fractionC5;
	  outputFileName += "f";
	}
      if (beamEnergy == "0") outputFileName += "_noExtraMask";
      if (nParts != 1)
	{
	  outputFileName += "_part";
	  outputFileName += part;
	}
      outputFileName += ".root";
    }
  
  std::cout << "*** Name of outputfile = " << outputFileName << endl;
  
  TFile* outputFile = new TFile(outputFileName,"recreate");
  
  //-------------------------------------------------------
  //define histograms
  //-------------------------------------------------------
  //output histogram list
  TList* histogramList = new TList(); 

  //number of hits distributions
  TH1D* hHitsDistributionSelection = new TH1D("hitDistrSel","", 5000,0,5000);
  hHitsDistributionSelection -> Sumw2();
  hHitsDistributionSelection -> SetTitle("Distribution of the number of hits per event;number of hits;#");
  histogramList -> Add(hHitsDistributionSelection);

  //TH1D* hHitsPerLaneSel = new TH1D("hitsPerLaneSel", "Hits per lane", maxNChips+1,-1,maxNChips);
  //hHitsPerLaneSel -> Sumw2();
  //histogramList -> Add(hHitsPerLaneSel);

  TH1D* hHitsDistribution = new TH1D("hitsDistribution","Hits Distribution",5000,0,5000); //before selection
  hHitsDistribution -> Sumw2();
  hHitsDistribution -> SetTitle("Distribution of the number of hits per event;number of hits;#");
  histogramList -> Add(hHitsDistribution);
  
  //TH1D* hEventsPerLane = new TH1D("hEventsPerLane","hEventsPerLane",maxNChips+1,-1,maxNChips); //one extra bin for other lane numbers
  //hEventsPerLane -> Sumw2();
  //hEventsPerLane -> SetTitle("Events per lane;lane;events");
  //histogramList -> Add(hEventsPerLane);

  //TH1D* hHitsPerLaneNoSel = new TH1D("hHitsPerLaneNoSel", "Hits per lane", maxNChips+1,-1,maxNChips);
  //hHitsPerLaneNoSel -> Sumw2();
  //histogramList -> Add(hHitsPerLaneNoSel);

  //hitmap over all selected events and for a single event (before selection)
  //also mean hit position and spread
  //occupancy per chip (before selection) and clustersize
  TH2I* hHitMap[maxNChips];
  //TH2I* hHitMap_event[maxNChips];
  /*TH1D* hMeanCol[maxNChips];
  TH1D* hMeanRow[maxNChips];
  TH1D* hSpreadCol[maxNChips];
  TH1D* hSpreadRow[maxNChips];
  TH1D* hOccupancy[maxNChips];
  TH1D* hClusterSize[maxNChips];   */ 
      
  for (int lane = 0; lane < maxNChips; lane++)
    {
      TString name = "hitMap_lane";
      name += lane;
      hHitMap[lane] = new TH2I(name,"",columnsPerChip,0,columnsPerChip,rowsPerChip,0,rowsPerChip);
      hHitMap[lane] -> SetTitle("Hit Map;column;row");
      histogramList -> Add(hHitMap[lane]);
      
      //name = "hitMap_event_lane";
      //name += lane;
      //hHitMap_event[lane] = new TH2I(name,"",0.5*columnsPerChip,0,columnsPerChip,0.5*rowsPerChip,0,rowsPerChip);
      //hHitMap_event[lane] -> SetTitle("Hit Map single event;column;row");
      //histogramList -> Add(hHitMap_event[lane]);
      
      /*      TString nameC = "meanCol_lane";
      nameC += lane;
      hMeanCol[lane] = new TH1D(nameC,"",0.1*columnsPerChip,0,columnsPerChip);
      hMeanCol[lane] -> Sumw2();
      hMeanCol[lane] -> SetTitle("Mean column;column;#");
      histogramList -> Add(hMeanCol[lane]);

      TString nameR = "meanRow_lane";
      nameR += lane;
      hMeanRow[lane] = new TH1D(nameR,"",0.1*rowsPerChip,0,rowsPerChip);
      hMeanRow[lane] -> Sumw2();
      hMeanRow[lane] -> SetTitle("Mean row;row;#");
      histogramList -> Add(hMeanRow[lane]);
      
      TString nameSC = "spreadCol_lane";
      nameSC += lane;
      hSpreadCol[lane] = new TH1D(nameSC,"",0.1*columnsPerChip,0,columnsPerChip);
      hSpreadCol[lane] -> Sumw2();
      hSpreadCol[lane] -> SetTitle("Spread column;column;#");
      histogramList -> Add(hSpreadCol[lane]);
		     
      TString nameSR = "spreadRow_lane";
      nameSR += lane;
      hSpreadRow[lane] = new TH1D(nameSR,"",0.1*rowsPerChip,0,rowsPerChip);
      hSpreadRow[lane] -> Sumw2();
      hSpreadRow[lane] -> SetTitle("Spread row;row;#");
      histogramList -> Add(hSpreadRow[lane]);

      TString nameO = "occupancy_lane";
      nameO += lane;
      hOccupancy[lane] = new TH1D(nameO,"",1000,0,0.001);
      hOccupancy[lane] -> Sumw2();
      hOccupancy[lane] -> SetTitle("Mean occupancy;Occupancy;#");
      histogramList -> Add(hOccupancy[lane]);
      
      TString nameCS = "clusterSize_lane";
      nameCS += lane;
      hClusterSize[lane] = new TH1D(nameCS,"",100,0,100);
      hClusterSize[lane] -> Sumw2();
      hClusterSize[lane] -> SetTitle(";cluster size;#");
      histogramList -> Add(hClusterSize[lane]);     */
    }

  //-------------------------------------------------------
  //read the data
  //-------------------------------------------------------
  
  //open the input root file
  fileLocation += baseName;
  fileLocation += run;
  if (UND) fileLocation += "_eventsLeft_justcluster";
  else
    {
      fileLocation += "/rootout_raw/conv_Run_";
      fileLocation += run;
    }
  fileLocation +=".root";
  std::cout<<endl<<"*** Reading file: "<<fileLocation<<endl<<endl<<"Ignore the following warning if there is one."<<endl;
  TFile* inputFile = TFile::Open(fileLocation);
  std::cout<<endl;

  inputFile->cd(); //set to current directory
  
  //get the tree
  TTree* frames = (TTree*)inputFile->Get("Frames");
  //get the branches of the TTree
  int runNumber, fileNumber, eventIndex, nHits, dataSize, eventNumberOriginal;
  frames->SetBranchAddress("runNumber",&runNumber);
  frames->SetBranchAddress("fileNumber",&fileNumber);
  frames->SetBranchAddress("eventNumber",&eventIndex); //Note that this does not run from 0 to nEvents but resets to 0 sometimes
  frames->SetBranchAddress("nHits",&nHits);
  frames->SetBranchAddress("dataSize",&dataSize);
  vector<Int_t>* vlane = new vector<Int_t>();
  vector<Int_t>* vcolumn = new vector<Int_t>();
  vector<Int_t>* vrow = new vector<Int_t>();
  frames->SetBranchAddress("lane",&vlane);
  frames->SetBranchAddress("column",&vcolumn);
  frames->SetBranchAddress("row",&vrow);
  int nEvents = frames->GetEntries();

  int nEventsOriginal = nEvents;
  if (UND)
    {
      frames->SetBranchAddress("eventNumberOriginal",&eventNumberOriginal);
      frames->GetEntry(nEvents-1);
      nEventsOriginal = eventNumberOriginal;
    }

  //number of hits vs event number histogram
  outputFile->cd(); //set outputFile as current directory

  TH1I* hHitsvsEvent = new TH1I("hitsvsEvent","",nEventsOriginal,0,nEventsOriginal);
  hHitsvsEvent -> SetTitle("Number of hits per event;event number;number of hits");
  histogramList -> Add(hHitsvsEvent);

  inputFile->cd(); //set inputFile as current directory

  //-------------------------------------------------------
  // Create TTree with selected events
  //-------------------------------------------------------

  //creat the file and the tree
  TFile* outputTree;
  int eventNumberNew;
  int eventNumberOld;
  fileLocationOutputTree += baseName;
  fileLocationOutputTree += run;
  fileLocationOutputTree += "_eventsLeft_justcluster.root";
  if (!(CT)) fileLocationOutputTree = "tobedeleted"; //Otherwise a previously made file might be deleted later
  outputTree = new TFile(fileLocationOutputTree,"recreate");
  outputTree->cd();
  TTree Frames("Frames", "Focal frames");

  if (CT) //if CT, create the branches of the tree
    {
      eventNumberNew = -1;
      //Frames.Branch("runNumber",&runNumber,"r/I");
      //Frames.Branch("fileNumber",&fileNumber,"f/I");
      //Frames.Branch("eventNumber",&eventNumberNew,"e/I");
      //Frames.Branch("eventNumberOriginal",&eventNumberOld,"eo/I"); //This number does not reset to 0 sometimes, it goes from 0 to nEvents
      Frames.Branch("nHits",&nHits,"n/I");
      //Frames.Branch("dataSize",&dataSize,"d/I");
      //Frames.Branch("lane",&vlane);
      //Frames.Branch("column",&vcolumn);
      //Frames.Branch("row",&vrow);
    }
  else //else delete the file again. This making and deleting is needed because Frames can not be created in an if statement.
    {
      remove(fileLocationOutputTree);
    }

  inputFile->cd();

  //------------------------------------------------------- 
  // Loop over events 
  //-------------------------------------------------------
  
  std::cout<<"*** Loop over events in input file"<<endl;
  int nReadEvents = 0; //count how many events the loop has read
  int nEmptyEvents = 0; //count the number of "events" that have less than 1 hits
  //int maxCountedChips = 0; //count the number of active chips
  //int nFailsPerLane[maxNChips] = {0}; //check how often a lane has the wrong number of hits
  
  //The events which are looped over
  int minEvent = nEvents/nParts*(part-1);
  int maxEvent = nEvents/nParts*part;  
  if (testing)
    {
      minEvent = 0;
      maxEvent = 100;
    }

  for (int event = minEvent; event < maxEvent; event++)
    {
      frames->GetEntry(event);

      if (CT)
	{
	  if (UND) eventNumberOld = eventNumberOriginal;
	  else eventNumberOld = event;
	}

      //create mTowerEvent
      mTowerEvent* currentEvent = new mTowerEvent(runNumber,eventIndex);
      currentEvent->setNHits(nHits);
      currentEvent->setNChips(maxNChips);
      //currentEvent->setNLayers(maxNChips/2);
      TObjArray* hitList = currentEvent->getHits();
      if (DB)
	{
	  std::cout<<endl<<"(DB) Run: "<<runNumber<<", event: "<<event<<"/"<<nEvents-1<<", number of hits: "<<nHits;
	  if (UND) std::cout<<", original event number: "<<eventNumberOriginal;
	  std::cout<<endl<<"(DB) Loop over hits in event "<<event<<" to add hits to hitlist and apply extra mask"<<endl;
	}

      int nHitsMasked = 0; //number of extra hits masked
      for (int hit=0; hit<nHits; hit++)
	{
	  mTowerHit* currentHit = new mTowerHit();
	  Int_t lane = vlane->at(hit);
	  Int_t col = vcolumn->at(hit);
	  Int_t row = vrow->at(hit);
	  if (hMaskPtn->GetBinContent(lane-laneOffset+1,col+1,row+1) == 0) //extra masking
	    {
 	      currentHit->setCoordinates(lane, col, row);
	      hitList->Add(currentHit);
	    }
	  else nHitsMasked++;
	}
      if (DB) std::cout<<"(DB) "<<nHitsMasked<<" extra hits masked in this event"<<endl;
     
      //------------------------------------------------------- 
      //Read the event and do analysis
      //------------------------------------------------------- 

      nReadEvents++;
      
      int entries = hitList->GetEntries();
      
      if (entries < 1) nEmptyEvents++;
      else
      	{
	  //Initialising variables
	  double meanCol[maxNChips] = {0.0};
	  double mean2Col[maxNChips] = {0.0};
	  double meanRow[maxNChips] = {0.0};
	  double mean2Row[maxNChips] = {0.0};
	  int nHitsPerLane[maxNChips]={0};
	  mTowerChip* hitsInChip[maxNChips];
	  vector<vector<int>> hitsInLayer(maxNChips/2, vector<int>{});
	  for (int l = 0; l<maxNChips ; l++)
	    {
	      hitsInChip[l] = new mTowerChip(l);
	      hitsInChip[l]->setLane(l+laneOffset);
	    }

	  //loop over all entries
	  if (DB) std::cout<<"(DB) Loop over hits in event "<<event<<" to get properties"<<endl;
	  for (int hit = 0; hit < entries; hit++)
	    {
	      mTowerHit* currentHit = (mTowerHit*)hitList->At(hit);
	      int lane = currentHit->getLane();
	      int row = currentHit->getRow();
	      int column = currentHit->getColumn();

	      if ((lane-laneOffset) <maxNChips &&(lane-laneOffset)>-1 )
		{	  
		  hHitMap[lane-laneOffset]->Fill(column,row,1);
		  nHitsPerLane[lane-laneOffset]++;
		  //hHitsPerLaneNoSel->Fill(lane-laneOffset);
		  meanCol[lane-laneOffset] += column;
		  mean2Col[lane-laneOffset] += column*column;
		  meanRow[lane-laneOffset] += row;
		  mean2Row[lane-laneOffset] += row*row;
		  hitsInLayer[lane2layer_lut.at(lane)].push_back(hit);
		  hitsInChip[lane-laneOffset]->AddHit(currentHit);
		}
	      else
		{
		  std::cout<<"lane number of hit "<<hit<<" out of range: "<<lane<<endl;
		}
	    } //loop over entries
	  
	  //calculate occupancy per chip
	  /*for (int l = 0; l<maxNChips ; l++)
	    {
	      if (hitsInChip[l]->getNHits() > 0)
		{
		  //hEventsPerLane->Fill(l);
		  hOccupancy[l] -> Fill(double(hitsInChip[l]->getNHits())/double(rowsPerChip*columnsPerChip));
		}
		}*/
	  
	  //Fill raw data histogram
	  hHitsDistribution->Fill(entries);

	  //------------------------------------------------------- 
	  //Clustering
	  //-------------------------------------------------------

	  if (DB) std::cout<<"(DB) CLUSTERING"<<endl;
	  vector<vector<int>> vClusters(4, vector<int> {}); //Vector with for every cluster: lanecode, meanX, meanY, clustersize
	  int nClusters = 0;
	  	      
	  for (int lanecode = 0; lanecode < 2; lanecode++) //loop over chips in first layer to find all clusters
	    {
	      int l = laneNumber[lanecode];
	      if (hitsInChip[l]->getNHits()>0)
		{
		  //get the array of clusters
		  hitsInChip[l]->Clusterize();
		  TObjArray* clusterlist = hitsInChip[l]->getClusters();
		  for (int c = 0;c<clusterlist->GetEntries();c++) //loop over clusters
		    {
		      mTowerCluster* cluster = (mTowerCluster*) clusterlist->At(c);
		      if (cluster)
			{
			  nClusters++;

			  //get properties of cluster
			  if (DB) std::cout<<"(DB) cluster id = "<<cluster->getId()<<", and lane = "<<cluster->getLane()<<", has "<<cluster->getNHits()<<" hits. This is cluster nr "<<nClusters-1<<endl;
			  vClusters[0].push_back(lanecode); //lanecode
			  
			  //get mean coordinate of cluster.
			  int meanRow = round(cluster -> getMeanRow());
			  int meanColumn = round(cluster -> getMeanColumn());
			  
			  //change chip dependent coordinates to absolute coordinates   
			  int meanX = meanRow;
			  int meanY = meanColumn;
			  if (IsLeftChip(l+laneOffset)) //cluster is in the left chip
			    {
			      meanX = -meanX - 1 - nPixelsGap;
			      meanY = columnsPerChip - 1 - meanY;
			    }
			  meanX += rowsPerChip + nPixelsGap;

			  vClusters[1].push_back(meanX); //meanX
			  vClusters[2].push_back(meanY); //meanY

			  vClusters[3].push_back(cluster->getNHits()); //clustersize
			  //hClusterSize[l]->Fill(cluster->getNHits());
			}
		    }
		}
	    }

	  //------------------------------------------------------- 
	  //Event selection
	  //-------------------------------------------------------	      	      

	  if (DB) std::cout<<"SELECTION"<<endl;

	  bool eventRejected = false; //variable for event criteria
	  int acceptedCluster = -1; //the index of the accepted cluster - for now no cluster accepted
	  int meanXAC; //the mean x of the accepted cluster
	  int meanYAC; //the mean y of the accepted cluster

	  //------------------------------------------------------- 
	  //Criterion C3 part 1/2: (check if there is more then 1 cluster with a clustersize bigger then one)
	  //-------------------------------------------------------	      	      

	  int nBigClusters = 0; //number of clusters with a clustersize bigger then 1
	  if (C3)
	    {
	      for (int c = 0; c < nClusters; c++) //loop over clusters
		{
		  if (vClusters[0][c] < 2 && vClusters[3][c] > 1) //if cluster is in first layer and clustersize > 1
		    {
		      nBigClusters++;
		      if (nBigClusters > 1)
			{
			  eventRejected = true;
			  if (DB) std::cout<<"(DB) event rejected by C3 criterion"<<endl;
			  break; 
			}
		    }
		}
	    } //if C3

	  //------------------------------------------------------- 
	  //Criteria C1 and C2
	  //-------------------------------------------------------	      	      

	  if (C2)
	    {
	      for (int c = 0; c < nClusters; c++) //loop over clusters for criteria C1 and C2
		{
		  if (eventRejected) break;
		  if (vClusters[0][c] < 2) //if cluster is in first layer
		    {
		      int meanXC = vClusters[1][c]; //mean x of this cluster
		      int meanYC = vClusters[2][c]; //mean y of this cluster

		      //Criterion C2
		      for (int hitl2 = 0; hitl2 < hitsInLayer[1].size(); hitl2++) //loop over hits in second layer for criterion C2
			{
			  int hit = hitsInLayer[1][hitl2];

			  //get properties of hit
			  mTowerHit* currentHit = (mTowerHit*)hitList->At(hit);
			  int lane = currentHit->getLane() - laneOffset;
			  int column = currentHit->getColumn();
			  int row = currentHit->getRow();
			  
			  //change chip dependent coordinates to absolute coordinates   		      
			  int x = row;
			  int y = column;
			  if (IsLeftChip(lane+laneOffset)) //hit is in the left chip
			    {
			      x = -x - 1 - nPixelsGap;
			      y = columnsPerChip - 1 - y;
			    }
			  x += rowsPerChip + nPixelsGap;

			  if (pow(pow(x-meanXC,2)+pow(y-meanYC,2),0.5)<nPixelRadiusC2) //if hit is in search area
			    {
			      if (DB) std::cout<<"(DB) cluster nr "<<c<<" accepted by criterion C2"<<endl;
			      if (acceptedCluster != -1) //If there is already an accepted cluster
				{
				  eventRejected = true;
				  if (DB) std::cout<<"(DB) more then 1 accepted cluster in the first layer, so event rejected by criterion C2"<<endl;
				}
			      acceptedCluster = c;
			      meanXAC = meanXC;
			      meanYAC = meanYC;
			      break; //stop loop over hits
			    } //if hit in search area
			} //loop over hits in second layer for criterion C2
		      if (DB && c != acceptedCluster) std::cout<<"(DB) cluster nr "<<c<<" ignored by criterion C2"<<endl;
		      
		      //Criterion C1
		      if (C1)
			{
			  if (c == acceptedCluster && !(meanXC <= maxXC1 && meanXC >= minXC1 && meanYC <= maxYC1 && meanYC >= minYC1)) //If cluster is accepted by criterion C2, but x and/or y are not inside of search area, reject event.
			    {
			      eventRejected = true;
			      if (DB) std::cout<<"(DB) event rejected by C1 criterion"<<endl;			  
			    }
			  else if (DB) std::cout<<"(DB) cluster accepted by C1 criterion"<<endl;
			} //if C1

		    }//if cluster in first layer
		} //loop over clusters for criterion C2 (and C1)

	      if (acceptedCluster == -1) //If there is no cluster accepted
		{
		  eventRejected = true;
		  if (DB) std::cout<<"(DB) 0 accepted clusters in the first layer, so event rejected"<<endl;
		}

	    } //if C2

	  //------------------------------------------------------- 
	  //Criterion C3 (part 2/2: second check needed if the accepted cluster has only one hit)
	  //-------------------------------------------------------	      	      

	  if (C3 && !(eventRejected))
	    {
	      if (vClusters[3][acceptedCluster] == 1 && nBigClusters == 1) //if the clustersize of the accepted cluster is one but there is a cluster with more then one hit
		{
		  eventRejected = true;
		  if (DB) std::cout<<"(DB) event rejected by C3 criterion"<<endl;
		}
	      else if (DB) std::cout<<"(DB) event accepted by C3 criterion"<<endl;
	    } //if C3
	  
	  //------------------------------------------------------- 
	  //Criteria C1ES, C4 and C5
	  //-------------------------------------------------------	      	      

	  if ((C1ES || C4 || C5) && !(eventRejected))
	    {
	      int nHitsInCircleC5 = 0; //number of hits in circle for criterion C5
	      int nHitsInBorderC5 = 0; //number of hits in border for criterion C5

	      for (int layer = 0; layer < 2 + nLayersC5; layer++) //loop over layers
		{
		  if (layer == 0 && !(C1ES)) continue; //Only C1ES needs hits in the first layer
		  if (layer == 1 && !(C4)) continue; //Only C4 needs hits in the second layer
		  if (layer == 2 && !(C5)) break; //Only C5 needs hits in the later layers

		  for (int hitl2 = 0; hitl2 < hitsInLayer[layer].size(); hitl2++) //loop over hits in layers
		    {
		      int hit = hitsInLayer[layer][hitl2];
		      
		      //get properties of hit
		      mTowerHit* currentHit = (mTowerHit*)hitList->At(hit);
		      int lane = currentHit->getLane() - laneOffset;
		  
		      int column = currentHit->getColumn();
		      int row = currentHit->getRow();
		  
		      //change chip dependent coordinates to absolute coordinates   		      
		      int x = row;
		      int y = column;
		      if (IsLeftChip(lane+laneOffset)) //cluster is in the left chip
			{
			  x = -x - 1 - nPixelsGap;
			  y = columnsPerChip - 1 - y;
			}
		      x += rowsPerChip + nPixelsGap;
		      		      
		      //Criterion C1ES
		      if (C1ES && layer == 0 && !(y <= maxYC1 && y >= minYC1 && x <= maxXC1 && x >= minXC1)) //if C1ES is applied, hit is in first layer, but y and/or x are not inside of acceptance area, reject event.
			{
			  eventRejected = true;
			  if (DB) std::cout<<"(DB) event rejected by C1ES criterion"<<endl;
			  break;
			} //if C1ES
		      
		      //Criterion C4
		      if (C4 && layer == 1 && pow(pow(x-meanXAC,2)+pow(y-meanYAC,2),0.5) > nPixelRadiusC4) //if C4 is applied, hit is in second layer and hit is outside of search area
			{
			  eventRejected = true;
			  if (DB) std::cout<<"(DB) event rejected by C4 criterion"<<endl;
			  break;
			} //if C4
		      
		      //Criterion C5 (part 1/2: checking where the hit is)
		      if (C5 && layer > 1) //if C5 is applied and the hit is not in the first two layers
			{
			  if (pow(pow(x-meanXAC,2)+pow(y-meanYAC,2),0.5) < nPixelRadiusC5) nHitsInCircleC5++; //If the hit is in the circle centered behind the accepted cluster
			  if (x < nPixelBorderC5 || x > maxX - nPixelBorderC5) nHitsInBorderC5++; //If the hit is in the border close to the minimal or maximal x value
			  else if (y < nPixelBorderC5 || y > maxY - nPixelBorderC5) nHitsInBorderC5++; //If the hit is in the border close to the minimal or maximal y value
			}

		    } //loop over hits in layers
		} //loop over layers
		 	      
	      //Criterion C5 (part 2/2: rejecting if required)
	      if (C5 && (nHitsInCircleC5 == 0 || (nHitsInCircleC5 != 0 && double(nHitsInBorderC5)/nHitsInCircleC5 > fractionC5))) //if C5 is applied and there are either no hits in the circle or more hits in border then allowed
		{
		  eventRejected = true;
		  if (DB) std::cout<<"(DB) event rejected by C5 criterion"<<endl;
		}

	      if (DB && !(eventRejected))
		{
		  if (C1ES) std::cout<<"(DB) event accepted by C1ES criterion"<<endl; 
		  if (C4) std::cout<<"(DB) event accepted by C4 criterion"<<endl;
		  if (C5) std::cout<<"(DB) event accepted by C5 criterion"<<endl;
		}
	    }//if C1ES, C4 or C5

	  //------------------------------------------------------- 
	  //Fill histograms and fill TTree
	  //-------------------------------------------------------
	  if (!(eventRejected))
	    {
	      if (DB) std::cout << "Filling histograms" << endl;

	      //for (int l = 0; l<maxNChips ; l++)
	      //{
	      //  hHitsPerLaneSel->Fill(l,nHitsPerLane[l]);
	      //}
	     
	      if (UND) hHitsvsEvent->Fill(eventNumberOriginal,entries);
	      else hHitsvsEvent->Fill(event,entries);
	      hHitsDistributionSelection -> Fill(entries);
     
	      for (int l = 0; l<maxNChips ; l++)
		{
		  if (nHitsPerLane[l] > 0)
		    {
		      //calculate mean
		      meanCol[l] = meanCol[l]/nHitsPerLane[l];
		      //hMeanCol[l]->Fill(meanCol[l]);
		      meanRow[l] = meanRow[l]/nHitsPerLane[l];
		      //hMeanRow[l]->Fill(meanRow[l]);
		      //calculate spread
		      mean2Col[l] = mean2Col[l]/nHitsPerLane[l];
		      mean2Row[l] = mean2Row[l]/nHitsPerLane[l];
		      double spreadCol =  TMath::Sqrt( mean2Col[l] - ( meanCol[l]* meanCol[l]) );
		      double spreadRow =  TMath::Sqrt( mean2Row[l] - ( meanRow[l]* meanRow[l]) );
		      //hSpreadCol[l]->Fill(spreadCol);
		      //hSpreadRow[l]->Fill(spreadRow);
		    }
		  delete hitsInChip[l];
		}		      

	      if (CT) //Fill TTree
		{
		  outputTree->cd();
		  eventNumberNew++;
		  Frames.Fill();
		  inputFile->cd();
		}

	    } //if the event is not rejected
	} //end of event selection
      delete currentEvent;
    } //end of loop over events

  //(Optional) Find hot pixels and print them
  if (HP)
    {
      std::cout<<"*** Find hot pixels"<<endl;
      for (int l = 0;l<maxNChips;l++)
	{
	  int nHitsLane = hHitMap[l]->GetEntries();
	  if (nHitsLane > 0)
	    {
	      for (int c = 1; c < columnsPerChip+1; c++) //shift by 1 because bins start counting at 1
		{
		  for (int r = 1;r < rowsPerChip+1;r++)
		    {
		      int nHitsBin = hHitMap[l]->GetBinContent(c,r);
		      double fractionBin = (double)nHitsBin/(double)nHitsLane;
		      if (fractionBin > 0.001) std::cout<<"Lane "<<l<<": Bin ("<<c<<", "<<r<<") has "<<nHitsBin<<" entries, "<<fractionBin<<endl;
		    }
		}
	    }
	}
    } //if HP

  //std::cout<<"*** The number of active chips is "<<maxCountedChips<<endl;
  std::cout<<"*** There were "<<nReadEvents<<" events read by the event loop"<<endl;
  std::cout<<"*** There were "<<nEmptyEvents<<" events with less than 1 hits"<<endl;
  //std::cout<<"Fails per lane:"<<endl;
  //for (int l = 0;l<maxNChips;l++)
  //  {
  //    std::cout<<"   lane "<<l<<": "<<nFailsPerLane[l]<<endl;
  //  }
  
  //----------------------------------------------------------------------------------
 
  //Write the output to disk
  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  
  std::cout<<" I managed to get all the way to saving the tree!";

  if (CT)
    {
      outputTree->cd();
      outputTree->Write();
      outputTree->Close();
    }
}
