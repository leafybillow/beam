/* #include <vector> */
/* #include <utility> */
/* #include <iostream> */
/* #include <stdio.h> */

#ifndef BeamTypes_h
#define BeamTypes_h

using namespace std;

struct ACluster{
  pair<int, int> pRange;
  int fSplit;  // if 0, no split is detected
  vector<int> valley; // valley position;
};

struct AHit{
  double fPosition; // Hit position in this projection, unit:mm
  double fCharge; // amount of charge integrated over this hit, unit: adc
  double fRes;  // spatial resolution of this hit, unit um
  int fWidth; // A single hit width, unit: # of strips, an integer
  pair<int, int> pRange;
};

struct ATrack{
  vector<double> x;
  vector<double> y;
  vector<double> z;
  double fSlope_zx;
  double fSlope_zy;
  double fIntercept_x;
  double fIntercept_y;
  double fChi2;
  vector<int> myPattern;
};
#endif
