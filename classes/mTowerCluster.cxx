#include "mTowerCluster.h"
#include "mTowerHit.h"
#include "TObjArray.h"
#include <iostream>

using namespace std;

ClassImp(mTowerCluster)

//default constructor
mTowerCluster::mTowerCluster()
{
  id = 0;
  hits = new TObjArray();
}

mTowerCluster::mTowerCluster(int i)
{
  id = i;
  hits = new TObjArray();
}


//destructor
mTowerCluster::~mTowerCluster()
{
  delete hits;
}

//adding a hit to the cluster
void mTowerCluster::AddHit(mTowerHit* hit)
{
  hits->Add(hit);
}

//merging another cluster with this cluster
void mTowerCluster::AddCluster(mTowerCluster* cluster)
{
  //NOT TESTED YET
  cout<<"Merging clusters: "<<getNHits()<<" hits and "<<cluster->getHits()<<endl;
  //TObjArray* hitsToMerge = cluster->getHits();
  for (int i = 0;i < cluster->getNHits();i++)
    {
      mTowerHit* hit = (mTowerHit*)(cluster->getHits())->At(i);
      hits->Add(hit);
    }

}

//Calculating the mean position of the cluster
double mTowerCluster::getMeanRow()
{
  double meanRow = 0.0;
  int nHits = getNHits();
  if (hits && nHits != 0)
    {
      for (int i = 0;i < nHits;i++)
	{
	  meanRow += ((mTowerHit*)hits->At(i))->getRow();
	}
      meanRow = double(meanRow)/nHits;
    }
  return meanRow;
}

double mTowerCluster::getMeanColumn()
{
  double meanColumn = 0.0;
  int nHits = getNHits();
  if (hits && nHits != 0)
    {
      for (int i = 0;i < nHits;i++)
	{
	  meanColumn += ((mTowerHit*)hits->At(i))->getColumn();
	}
      meanColumn = double(meanColumn)/nHits;
    }
  return meanColumn;
}


