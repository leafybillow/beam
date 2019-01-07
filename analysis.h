#define N_HITS_MAX 10; // Buffer foro hits allowed 

/* It is just easier to use C structures to group data in a compact format.  */
/* C structures are chosen for simplicity.  */
/* Alternatively, cpp classes can be used. */
/* A user-defined class branch requires a function to retrieve entries from a tree. */
/* Padding/memory alignment could be an issue in some platforms. Call experts. */

struct GEMPlane{ 

  Double_t pos_x[N_HITS_MAX];    // position of hits (mm) in X direction: 10 cm side
  Double_t pos_y[N_HITS_MAX];    // pos in Y direction : 20 cm side
  Double_t charge_x[N_HITS_MAX]; // deposited charge of hits in X readout strips 
  Double_t charge_y[N_HITS_MAX]; // (unit): integrated ADC channel)
  Double_t res_x[N_HITS_MAX]; // 2D position resolution (um);
  Double_t res_y[N_HITS_MAX];

  Double_t correlation[N_HITS_MAX];
  Double_t cmn[6]; // Common Mode noises on each apv

  Int_t mpl_x[N_HITS_MAX];       // Multiplicty: Number of stripes with amplitude above THRESHOLD*sigma
  Int_t mpl_y[N_HITS_MAX]; 
  Int_t nhits;           // Number of hits is found

  Bool_t isSplit;
};

struct Track{
  Double_t slope_x[N_HITS_MAX];
  Double_t slope_y[N_HITS_MAX];

  Double_t theta[N_HITS_MAX]; //polar angle: degree
  Double_t phi[N_HITS_MAX];  // azimuthal angle: degree
 
  Double_t pos_det_x[N_HITS_MAX]; // extrapolated hit position on detector
  Double_t pos_det_y[N_HITS_MAX];

  Int_t nTracks; // Number of tracks is found
  Bool_t isGood; //A Good Track
};

struct QDC{
  Double_t us_lo; // low/high range channel on v965 QDC
  Double_t us_hi;
  Double_t ds_lo; // for a downstream detector , if it exists, e.g. PREX-II tandem
  Double_t ds_hi;
};
TString leaflist_gem ="pos_x[N_HITS_MAX]/D:pos_y[N_HITS_MAX]:\
charge_x[N_HITS_MAX]:charge_y[N_HITS_MAX]:\
res_x[N_HITS_MAX]:res_y[N_HITS_MAX]:\
corelation[N_HITS_MAX]:\
cmn[6]/D:\
mpl_x[N_HITS_MAX]/I:mpl_y[N_HITS_MAX]:\
nhits/I:isSplit/O";

TString leaflist_track ="slope_x[N_HITS_MAX]/D:slope_y[N_HITS_MAX]:\
theta[N_HITS_MAX]:phi[N_HITS_MAX]:\
pos_det_x[N_HITS_MAX]:pos_det_y[N_HITS_MAX]:\
nTracks/I:isGood/O";

TString leaflist_qdc ="us_lo/D:us_hi:ds_lo:ds_hi";
