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
//const std::map< Int_t, Int_t > lane2layer_lut;                                                                                                                                                                                                
                                                                                                                                                                                                                                              
//const std::map<int,bool> layer2isInv_lut;*/

 
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


//------------------------------------------------------- 
//Clustering
//-------------------------------------------------------

bool EventSelection(bool DB, bool C1ES, bool C1, bool C2, bool C3, bool C4, bool C5, TObjArray* hitList, mTowerChip* hitsInChip[], vector<vector<int>> hitsInLayer, const int laneNumber[], const int laneOffset, const int columnsPerChip, const int rowsPerChip, double nPixelsGap, int nPixelRadiusC2, int nPixelRadiusC4, int nPixelRadiusC5, int nLayersC5, int nPixelBorderC5, double fractionC5, double minXC1, double maxXC1, double minYC1, double maxYC1, int maxX, int maxY)
{
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

  if (eventRejected){
    return false;
  } else {
    return true;
  }
  
}
  
