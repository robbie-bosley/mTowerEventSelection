#ifndef MTOWERCLUSTER_H
#define MTOWERCLUSTER_H

#include "TObjArray.h"
#include "mTowerHit.h"

//Class cluster for the mTower

class mTowerCluster : public TObject
{

 protected:
  int id; //cluster id number
  int lane; //lane of the chip
  TObjArray* hits; //of type mTowerHit
  

 public:
  mTowerCluster();
  mTowerCluster(int i);
  virtual ~mTowerCluster();

  int getNHits() {return hits->GetEntries();}
  TObjArray* getHits() {return hits;}
  void AddHit(mTowerHit* hit);
  void AddCluster(mTowerCluster* cluster);
  void setId(int i) {id = i;}
  int getId() {return id;}
  void setLane(int l) {lane = l;}
  int getLane() {return lane;}
  double getMeanRow();
  double getMeanColumn();

  ClassDef(mTowerCluster,1) //mTower cluster
    
};

#endif

