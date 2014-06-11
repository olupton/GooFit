#ifndef DALITZCOHERENCE_PDF_HH
#define DALITZCOHERENCE_PDF_HH

#include "GooPdf.hh" 
#include "DalitzPlotPdf.hh" 
#include "devcomplex.hh"
#include <vector>

class SpecialResonanceCoherenceIntegrator;
  
class DalitzPlotCoherencePdf : public GooPdf {
public:
  DalitzPlotCoherencePdf(std::string n, Variable* m12, Variable* m13, DalitzPlotPdf *pdfa, DalitzPlotPdf *pdfb, GooPdf *eff);
  // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the 
  // coherent sum. The caching method requires that it be done this way or the ProdPdf
  // normalisation will get *really* confused and give wrong answers. 

  __host__ virtual fptype normalise () const;
protected:

private:
  typedef std::vector<DalitzPlotPdf*> DPPvec;
  DPPvec pdfab;
  DalitzPlotPdf *pdfa, *pdfb;
  Variable *_m12, *_m13;
  fptype* dalitzNormRange; 

  // Following variables are useful if masses and widths, involved in difficult BW calculation, 
  // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
  devcomplex<fptype>*** integrals; // Caches the integrals of the BW waves for each combination of resonances. 

  bool** redoIntegral;
  SpecialResonanceCoherenceIntegrator*** integrators;
};

class SpecialResonanceCoherenceIntegrator : public SpecialResonanceIntegrator
{
public:
  // Class used to calculate integrals of terms BW_i * BW_j^*. 
  SpecialResonanceCoherenceIntegrator (int pIdx, unsigned int ri, unsigned int rj);
  EXEC_TARGET devcomplex<fptype> devicefunction(fptype m12, fptype m13, int res_i, int res_j, fptype* p, unsigned int* indices) const;
}; 

#endif
