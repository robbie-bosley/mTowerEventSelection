#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

// Event Selection Processor for the mTower

/*const std::map< Int_t, Int_t > chipid2lane_lut;
const std::map< Int_t, Int_t > lane2chipid_lut;
const std::map< Int_t, Int_t > chipid2layer_lut;
const std::map< Int_t, Int_t > lane2layer_lut;
const std::map<int,bool> layer2isInv_lut;*/

/*const std::map< Int_t, Int_t > chipid2lane_lut;

const std::map< Int_t, Int_t > lane2chipid_lut;

const std::map< Int_t, Int_t > chipid2layer_lut;

const std::map< Int_t, Int_t > lane2layer_lut;

const std::map<int,bool> layer2isInv_lut;*/

bool IsLeftChip(int lane);

bool EventSelection(bool DB, bool C1ES, bool C1, bool C2, bool C3, bool C4, bool C5, TObjArray* hitList, mTowerChip* hitsInChip[], vector<vector<int>> hitsInLayer, const int lanenumber[], const int laneoffset, const int columnsPerChip, const int rowsPerChip, double nPixelsGap, int nPixelRadiusC2, int nPixelRadiusC4, int nPixelRadiusC5, int nLayersC5, int nPixelBorderC5, double fractionC5, double minXC1, double maxXC1, double minYC1, double maxYC1, int maxX, int maxY);

#endif
