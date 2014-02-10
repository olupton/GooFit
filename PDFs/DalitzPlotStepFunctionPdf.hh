#ifndef DALITZPLOT_STEPFUNCTION_HH
#define DALITZPLOT_STEPFUNCTION_HH

#include "GooPdf.hh" 
#include "DalitzPlotHelpers.hh" 
#include "devcomplex.hh"
#include <vector>

class DalitzPlotStepFunctionPdf : public GooPdf {
public:
  DalitzPlotStepFunctionPdf (std::string n, Variable* m12, Variable* m13, DecayInfo* decay, fptype maxAbsCosine = 1.5);
private:
  Variable *_m12, *_m13;
  fptype _maxAbsCosine;
};

#endif
