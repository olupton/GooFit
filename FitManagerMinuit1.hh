#ifndef FITMANAGER_MINUIT1_HH
#define FITMANAGER_MINUIT1_HH

#include "TMinuit.h" 
extern PdfBase* pdfPointer; 
extern int numPars; 
#ifdef OMP_ON
#pragma omp threadprivate (numPars)
#pragma omp threadprivate (pdfPointer)
#endif

void FitFun(int &npar, double *gin, double &fun, double *fp, int iflag); 

class FitManager { 
public:
  enum FitAlgorithm { MIGRAD, SEEK };
  FitManager (PdfBase* dat);
  ~FitManager ();
  void setMaxCalls (double mxc) {overrideCallLimit = mxc;}
  void setupMinuit ();
  void runSeek();
  void runMigrad (); 
  void fit (FitAlgorithm algo = MIGRAD, Int_t strategy = 2);
  TMinuit* getMinuitObject () {return minuit;} 
  void getMinuitValues () const;
  void getMinosErrors() const;
  void getMinuitStatus(double& fmin, double& fedm, double& errdef, int& npari, int& nparx, int& istat) const;
  TMinuit* minuit;
  PdfBase *getPdfPointer();
private:
  double overrideCallLimit; 
};

#endif 
