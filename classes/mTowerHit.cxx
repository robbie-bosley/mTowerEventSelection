#include "mTowerHit.h"

ClassImp(mTowerHit)

mTowerHit::mTowerHit() : TObject()
{
  //hit default constructor
  lane = 0;
  column = 0;
  row = 0;
}

mTowerHit::mTowerHit(int l, int c, int r) : TObject()
{
  //hit normal constructor
  lane = l;
  column = c;
  row = r;
}

mTowerHit::~mTowerHit()
{
  //hit default destructor
}
