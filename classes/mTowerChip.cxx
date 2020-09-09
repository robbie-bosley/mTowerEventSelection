#include "mTowerChip.h"
#include "mTowerCluster.h"
#include "mTowerHit.h"
#include "TObjArray.h"
#include <iostream>

using namespace std;


ClassImp(mTowerChip)

//default constructor
mTowerChip::mTowerChip()
{
  hits = new TObjArray();
  clusters = new TObjArray();
  lane = -1;
}

mTowerChip::mTowerChip(int l)
{
  hits = new TObjArray();
  clusters = new TObjArray();
  lane = l;
}


//destructor
mTowerChip::~mTowerChip()
{
  delete hits;
  delete clusters;
}

//adding a hit to the chip
void mTowerChip::AddHit(mTowerHit* hit)
{
  hits->Add(hit);
}


//reset the status of all hits in the chip
void mTowerChip::ResetHitStatus()
{
  for (int h = 0;h < hits->GetEntries(); h++)
    {
      mTowerHit* hit = (mTowerHit*)hits->At(h);
      hit->setStatus(-1);
    }
}

//find the neighbours of a hit and add them to a TObjArray
int mTowerChip::findNeighbours(mTowerHit* currentHit, TObjArray* neighbours)
{
  int nNeighbours = neighbours->GetEntries();
  int dist = 2;
  
  for (int h = 0;h < hits->GetEntries(); h++)
    {
      mTowerHit* hit = (mTowerHit*)hits->At(h);
      if (hit->IsEqual(currentHit)) continue; //not for itself, Default equal comparison (objects are equal if they have the same address in memory).
      if ( TMath::Abs( hit->getRow() - currentHit->getRow() ) < dist && TMath::Abs( hit->getColumn() - currentHit->getColumn()) < dist )
	{
	  //add hit to collection
	  if (hit->getStatus() < 0) //do not add status of 0 or 1
	    {
	      neighbours->Add(hit);
	      hit->setStatus(0);
	      nNeighbours++;
	    }
	}
      
    }
  //cout<<"Found neighbours: "<<nNeighbours<<endl;
  return nNeighbours;
}


//cluster all hits in this chip
 int mTowerChip::Clusterize()
 {
   // use DBSCAN method of clustering
   // https://en.wikipedia.org/wiki/DBSCAN
   
   //assign a status and clusternumber to each hit
   //status:
   //   -1 : unvisited hit
   //   0  : neighbouring hit to current hit
   //   1  : hit assigned to a cluster
   //   -2 : noise hit

   ResetHitStatus(); //set all hits to -1
   
   int minClusterSize = 0;
   TObjArray* neighbours = new TObjArray();
   int currentClusterNumber = 0;

   for (int h = 0;h < hits->GetEntries(); h++)
     {
       mTowerHit* currentHit = (mTowerHit*)hits->At(h);
       //cout<<"status of hit "<<h<<" is "<<currentHit->getStatus()<<endl;
       if (currentHit->getStatus() == 1) continue; //hit has already been assigned
       //find all neighboring hits
       neighbours->Clear(); //start with a fresh list
       int nNeighbours = findNeighbours(currentHit,neighbours);
       if (nNeighbours < minClusterSize)
	 {
	   //label hit as noise
	   //cout<<"hit "<<h<<" labeled as noise"<<endl;
	   currentHit->setStatus(-2);
	   continue;
	 }
       //cout<<"Raising the cluster number by 1"<<endl;
       currentClusterNumber += 1; //add one to the cluster number
       //cout<<"cluster number = "<<currentClusterNumber<<endl;
       currentHit->setStatus(1);
       currentHit->setCluster(currentClusterNumber);
       //loop over all neighbouring hits
       //cout<<"Loop over all neighbours of hit "<<h<<endl;
       for (int n = 0;n < nNeighbours; n++)
	 {
	   mTowerHit* neighbourHit = (mTowerHit*)neighbours->At(n);
	   //cout<<"status of neighbour hit "<<n<<" is "<<neighbourHit->getStatus()<<endl;
	   if (neighbourHit->getStatus() == -2)
	     {
	       //cout<<"Neighbour was noise, now added to cluster"<<endl;
	       neighbourHit->setStatus(1);
	       neighbourHit->setCluster(currentClusterNumber);
	     }
	   if (neighbourHit->getStatus() == 1) continue; //hit has already been assigned

	   neighbourHit->setStatus(1);
	   neighbourHit->setCluster(currentClusterNumber);
	   nNeighbours = findNeighbours(neighbourHit,neighbours); 
	   
	 }
     }

   int nClusters = currentClusterNumber;
   //after all hits are processed and assigned to a cluster, make the clusters
   if (nClusters > 0)
     {
       for (int c = 1;c<nClusters+1;c++)
	 {
	   mTowerCluster* currentCluster = new mTowerCluster(c);
	   currentCluster->setLane(lane);
	   for (int h = 0;h < hits->GetEntries(); h++)
	     {
	       mTowerHit* currentHit = (mTowerHit*)hits->At(h);
	       int cluster = currentHit->getCluster();
	       if (cluster == c)
		 {
		   currentCluster->AddHit(currentHit);
		 }
	     }
	   clusters->Add(currentCluster->Clone()); //need to clone array (= duplicate hits), otherwise hits are deleted in the next line
	   delete currentCluster; 
	 }

       /*for (int c = 1;c<nClusters+1;c++)
	 {
	   mTowerCluster* currentCluster = (mTowerCluster*) clusters->At(c-1);
	   int nHitsInCluster = currentCluster->getNHits();
	   int clusterId = currentCluster->getId();
	   cout<<"Cluster "<<clusterId<<" has "<<nHitsInCluster<<" hits"<<endl;
	   }*/
     }
     
   delete neighbours; 
   return nClusters;
 }
  
