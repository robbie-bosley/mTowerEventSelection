#include "mTowerEvent.h"
#include "mTowerHit.h"
#include "TObjArray.h"

ClassImp(mTowerEvent)

//default constructor
mTowerEvent::mTowerEvent()
{
  eventNumber = 0;
  runNumber = 0;
  nChips = 0;
  nHits = 0;
  hits = new TObjArray();
}

//constructor
mTowerEvent::mTowerEvent(int r, int ev)
{
  eventNumber = ev;
  runNumber = r;
  nChips = 0;
  nHits = 0;
  hits = new TObjArray();
}

//destructor
mTowerEvent::~mTowerEvent()
{
  hits->Delete(); //deletes contents of hits, only done in this class!
  delete hits; //deletes hits object and frees memory
}

//calculate and return the number of active chips
int mTowerEvent::getNChips()
{
  const int maxNChips = 48;
  int nHitsChip[maxNChips];
  for (int c = 0;c<maxNChips;c++)
    {
      nHitsChip[c]=0;
    }
  int laneOffset = 0;
  
  if (nChips == 0)
    {
      
      //calculate the number of active chips
      for (int i = 0; i<nHits; i++)
	{
	  mTowerHit* currentHit = (mTowerHit*)hits->At(i);
	  int lane = currentHit->getLane();
	  nHitsChip[lane-laneOffset]++;
	}
      
      for (int c = 0;c<maxNChips;c++)
	{
	  if (nHitsChip[c] > 0)
	    {
	      nChips++;
	    }
	}
      
      return nChips;
    }
  else
    {
      return nChips;
    }
}
