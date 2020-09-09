#ifndef MTOWEREVENT_H
#define MTOWEREVENT_H

#include "TObjArray.h"
#include "mTowerHit.h"

// Class event for the mTower

class mTowerEvent : public TObject
{

 protected:
  int runNumber;
  int eventNumber;
  int nChips;
  int nHits;
  TObjArray* hits; //of type mTowerHit
  
 public:
  mTowerEvent();
  mTowerEvent(int r, int ev);
  virtual ~mTowerEvent();
  
  int getRunNumber() {return runNumber;}
  void setRunNumber(int r) {runNumber = r;}
  int getEventNumber() {return eventNumber;}
  void setEventNumber(int ev) {eventNumber = ev;}
  int getNHits() {return nHits;}
  void setNHits(int nh) {nHits = nh;}
  TObjArray* getHits() {return hits;}
  int getNChips();
  void setNChips(int nc) {nChips = nc;}
  
  ClassDef(mTowerEvent,1) //mTower event
    
};

#endif
