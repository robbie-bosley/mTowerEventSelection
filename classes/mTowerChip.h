#ifndef MTOWERCHIP_H
#define MTOWERCHIP_H

#include "TObjArray.h"
#include "mTowerHit.h"
#include "mTowerCluster.h"

//Class chip for the mTower (collection of hits in one chip)

class mTowerChip : public TObject
{

 protected:
  int lane; //lane of the chip
  TObjArray* hits; //of type mTowerHit
  TObjArray* clusters; //of type mTowerCluster
  

 public:
  mTowerChip();
  mTowerChip(int l);
  virtual ~mTowerChip();

  int getNHits() {return hits->GetEntries();}
  TObjArray* getHits() {return hits;}
  TObjArray* getClusters() {return clusters;}
  void AddHit(mTowerHit* hit);
  void ResetHitStatus();
  int findNeighbours(mTowerHit* currentHit, TObjArray* neighbours); //return the number of neighbouring hits
  int Clusterize(); //return the number of clusters and fill the clusters data member
  void setLane(int l) {lane = l;}
  int getLane() {return lane;}
  
  ClassDef(mTowerChip,1) //mTower chip
    
};

#endif

