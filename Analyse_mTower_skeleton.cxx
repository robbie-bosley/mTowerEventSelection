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
#include "subprocessors/EventSelection.h"

#include <iostream>
#include <fstream>
#include <string>


//conversion tables (Now included in EventSelection Processor)

const std::map< Int_t, Int_t > chipid2lane_LUT = {
  { 0,40},{ 1,39},{ 2,42},{ 3,41},{ 4,44},{ 5,43},{ 6,46},{ 7,45},
  { 8,48},{ 9,47},{10,50},{11,49},{12,52},{13,51},{14,54},{15,53},
  {16,38},{17,55},{18,36},{19,37},{20,32},{21,35},{22,34},{23,33},
  {24,64},{25,63},{26,66},{27,65},{28,68},{29,67},{30,70},{31,69},
  {32,72},{33,71},{34,74},{35,73},{36,76},{37,75},{38,78},{39,77},
  {40,62},{41,79},{42,60},{43,61},{44,56},{45,59},{46,58},{47,57}
};

const std::map< Int_t, Int_t > lane2chipid_LUT = {
  {40, 0},{39, 1},{42, 2},{41, 3},{44, 4},{43, 5},{46, 6},{45, 7},
  {48, 8},{47, 9},{50,10},{49,11},{52,12},{51,13},{54,14},{53,15},
  {38,16},{55,17},{36,18},{37,19},{32,20},{35,21},{34,22},{33,23},
  {64,24},{63,25},{66,26},{65,27},{68,28},{67,29},{70,30},{69,31},
  {72,32},{71,33},{74,34},{73,35},{76,36},{75,37},{78,38},{77,39},
  {62,40},{79,41},{60,42},{61,43},{56,44},{59,45},{58,46},{57,47}
};

const std::map< Int_t, Int_t > chipid2layer_LUT = {
  { 0,22},{ 1,22},{ 2,20},{ 3,20},{ 4,18},{ 5,18},{ 6,16},{ 7,16},
  { 8,14},{ 9,14},{10,12},{11,12},{12,10},{13,10},{14, 8},{15, 8},
  {16, 6},{17, 6},{18, 4},{19, 4},{20, 0},{21, 0},{22, 2},{23, 2},
  {24,23},{25,23},{26,21},{27,21},{28,19},{29,19},{30,17},{31,17},
  {32,15},{33,15},{34,13},{35,13},{36,11},{37,11},{38, 9},{39, 9},
  {40, 7},{41, 7},{42, 5},{43, 5},{44, 1},{45, 1},{46, 3},{47, 3}
};

const std::map< Int_t, Int_t > lane2layer_LUT = {
  {40,22},{39,22},{42,20},{41,20},{44,18},{43,18},{46,16},{45,16},
  {48,14},{47,14},{50,12},{49,12},{52,10},{51,10},{54, 8},{53, 8},
  {38, 6},{55, 6},{36, 4},{37, 4},{32, 0},{35, 0},{34, 2},{33, 2},
  {64,23},{63,23},{66,21},{65,21},{68,19},{67,19},{70,17},{69,17},
  {72,15},{71,15},{74,13},{73,13},{76,11},{75,11},{78, 9},{77, 9},
  {62, 7},{79, 7},{60, 5},{61, 5},{56, 1},{59, 1},{58, 3},{57, 3}
};

const std::map<int,bool> layer2isInv_LUT = {
  { 0, kFALSE}, { 1, kTRUE}, { 2, kFALSE}, { 3, kTRUE}, 
  { 4, kFALSE}, { 5, kTRUE}, { 6, kFALSE}, { 7, kTRUE}, 
  { 8, kFALSE}, { 9, kTRUE}, {10, kFALSE}, {11, kTRUE}, 
  {12, kFALSE}, {13, kTRUE}, {14, kFALSE}, {15, kTRUE}, 
  {16, kFALSE}, {17, kTRUE}, {18, kFALSE}, {19, kTRUE}, 
  {20, kFALSE}, {21, kTRUE}, {22, kFALSE}, {23, kTRUE} 
};

//for hit maps
int lane2padhitmap(int lane){
  int layerNr = lane2layer_LUT.at(lane);
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
  int layerNr = lane2layer_LUT.at(lane);
  int padid = layerNr+1;
  return padid; 
}


//MAIN PROCESSOR

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

  //INITIAL PARAMETERS
  bool UND = false; //is true if an earlier made TTree is used, so if not using the original data.
  bool testing = false; //Testing outputfile, run over small amount of events.

  bool DB = false; //debug print statements
  bool HP = false; //find hot pixels
  bool CT = true; //Create TTree (In addition to normal output file)
  bool C1 = false; //Criterion 1. Only events where the accepted cluster (C2) is within a certain area of the first layer are accepted for analysis. NEEDS C2.
  bool C1ES = false; //C1, but even events with just single hits outside area are rejected.
  bool C2 = true; //Criterion 2. Every cluster in the first layer with hits behind it in layer 2 is accepted. Only events with 1 accepted cluster are accepted for analysis.
  bool C3 = false; //Criterion 3. Every event with clusters with clustersize > 1, which is ignored by C2, is rejected. NEEDS C2
  bool C4 = false; //Criterion 4. Only events where there are no hits outside of a certain area in the second layer centrated around the accepted cluster (C2) in the first layer are accepted for analysis. NEEDS C2
  bool C5 = false; //Criterion 5. The fraction (nHits in border of chip)/(nHits in a circle centrated around the accepted cluster (C2)) for the 3th to the nth layer has to be higher then a threshold, otherwise the event is rejected. NEEDS C2


  //Some additional parameters. You probably won't want to change these significantly.
  
  TString fileLocationOutputTree = "./";
  TString fileLocation = "/eos/project/m/mtower/Data/Data_TB_February_2020/mTower_Data_DESY_Feb_2020_raw1/"; //The location of the selected data
  //TString fileLocation = "../";
  //TString fileLocation = "/data2/data/TB_mTower/mTower_Data_DESY_Feb_2020_raw1/"; //The location of the raw data
  TString maskingFileLocation = "./masking/"; //The location of the masking .txt files
  int nPixelSideC1 = 125; //length of the sides of the area for criterion C1. Actual row side is 1 pixel bigger if nPixelsGap is exactly odd to maintain symmetry.
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
  //DEFINE HISTOGRAMS
  //-------------------------------------------------------

  //output histogram list
  TList* histogramList = new TList(); 

  //number of hits distributions
  TH1D* hHitsDistribution = new TH1D("hitsDistribution","Hits Distribution",5000,0,5000); //before selection
  hHitsDistribution -> Sumw2();
  hHitsDistribution -> SetTitle("Distribution of the number of hits per event;number of hits;#");
  histogramList -> Add(hHitsDistribution);

  TH1D* hHitsDistributionSelection = new TH1D("hitsDistributionSelection","Hits Distribution After Event Selection",5000,0,5000); //after selection
  hHitsDistributionSelection -> Sumw2();
  hHitsDistributionSelection -> SetTitle("Distribution of the number of hits per event;number of hits;#");
  histogramList -> Add(hHitsDistributionSelection);
  
  //hitmap over all selected events and for a single event (before selection)
  TH2I* hHitMap[maxNChips];
      
  for (int lane = 0; lane < maxNChips; lane++)
    {
      TString name = "hitMap_lane";
      name += lane;
      hHitMap[lane] = new TH2I(name,"",columnsPerChip,0,columnsPerChip,rowsPerChip,0,rowsPerChip);
      hHitMap[lane] -> SetTitle("Hit Map;column;row");
      histogramList -> Add(hHitMap[lane]);
      
    }

  //-------------------------------------------------------
  //read the data
  //-------------------------------------------------------
  
  //open the input root file
  fileLocation += baseName;
  fileLocation += run;
  if (UND) fileLocation += "_eventsLeft";
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
  fileLocationOutputTree += "_eventsLeft.root";
  if (!(CT)) fileLocationOutputTree = "tobedeleted"; //Otherwise a previously made file might be deleted later
  outputTree = new TFile(fileLocationOutputTree,"recreate");
  outputTree->cd();
  TTree Frames("Frames", "Focal frames");

  if (CT) //if CT, create the branches of the tree
    {
      eventNumberNew = -1;
      Frames.Branch("runNumber",&runNumber,"r/I");
      Frames.Branch("fileNumber",&fileNumber,"f/I");
      Frames.Branch("eventNumber",&eventNumberNew,"e/I");
      Frames.Branch("eventNumberOriginal",&eventNumberOld,"eo/I"); //This number does not reset to 0 sometimes, it goes from 0 to nEvents
      Frames.Branch("nHits",&nHits,"n/I");
      Frames.Branch("dataSize",&dataSize,"d/I");
      Frames.Branch("lane",&vlane);
      Frames.Branch("column",&vcolumn);
      Frames.Branch("row",&vrow);
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
		  meanCol[lane-laneOffset] += column;
		  mean2Col[lane-laneOffset] += column*column;
		  meanRow[lane-laneOffset] += row;
		  mean2Row[lane-laneOffset] += row*row;
		  hitsInLayer[lane2layer_LUT.at(lane)].push_back(hit);
		  hitsInChip[lane-laneOffset]->AddHit(currentHit);
		}
	      else
		{
		  std::cout<<"lane number of hit "<<hit<<" out of range: "<<lane<<endl;
		}
	    } //loop over entries
	  
	  //Fill raw data histogram
	  hHitsDistribution->Fill(entries);

	  // EVENT SELECTION PROCESSOR

	  if (EventSelection(DB,C1ES,C1,C2,C3,C4,C5,hitList,hitsInChip,hitsInLayer,laneNumber,laneOffset,columnsPerChip,rowsPerChip,nPixelsGap,nPixelRadiusC2,nPixelRadiusC4,nPixelRadiusC5,nLayersC5,nPixelBorderC5,fractionC5,minXC1,maxXC1,minYC1,maxYC1,maxX,maxY)==true) {

	    //------------------------------------------------------- 
	    //Fill histograms and fill TTree
	    //-------------------------------------------------------

	    if (DB) std::cout << "Filling histograms" << endl;
	    hHitsDistributionSelection -> Fill(entries);
	    
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

  std::cout<<"*** There were "<<nReadEvents<<" events read by the event loop"<<endl;
  std::cout<<"*** There were "<<nEmptyEvents<<" events with less than 1 hits"<<endl;
 
  //Write the output to disk
  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  
  std::cout<<" Saving the tree...";

  if (CT)
    {
      outputTree->cd();
      outputTree->Write();
      outputTree->Close();
    }
  std::cout<<" DONE! ";
}
