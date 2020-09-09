#ifndef _MTOWERHIT_
#define _MTOWERHIT_

#include "TObject.h"


class mTowerHit : public TObject {

 private:
  int lane;
  int row;
  int column;
  int status; //for clustering
  int cluster; //for clustering

 public:
  mTowerHit();
  mTowerHit(int l, int r, int c);
  virtual ~mTowerHit();
  int getLane() {return lane;}
  int getRow() {return row;}
  int getColumn() {return column;}
  int getStatus() {return status;}
  int getCluster() {return cluster;}
  void setLane(int l) {lane = l;}
  void setRow(int r) {row = r;}
  void setColumn(int c) {column = c;}
  void setCoordinates(int l, int c, int r) {lane = l; column = c; row = r;}
  void setStatus(int s) {status = s;}
  void setCluster(int c) {cluster = c;}
  

  ClassDef(mTowerHit,2) //a simple hit for mTower data
};

#endif
