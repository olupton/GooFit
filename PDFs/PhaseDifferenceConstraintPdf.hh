#ifndef PHASEDIFF_PDF_HH
#define PHASEDIFF_PDF_HH

#include "GooPdf.hh" 

// takes two functions a(x) and b(x)
// returns a gaussian constraint of (a(x)-b(x)) to 'mean' with error 'sigma',
// accounting for 2pi phase redundancies
class PhaseDifferenceConstraintPdf : public GooPdf
{
  public:
    PhaseDifferenceConstraintPdf(const std::string &n, PdfBase *pdfa, PdfBase *pdfb, Variable *mean, Variable *sigma);
    __host__ virtual fptype normalise() const;
};

#endif
