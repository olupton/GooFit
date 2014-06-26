#include "PhaseDifferenceConstraintPdf.hh"

EXEC_TARGET fptype device_PhaseDiffConstraint (fptype* evt, fptype* p, unsigned int* indices) {
  unsigned int
    PDFaFcnIndex(indices[1]),
    PDFaParIndex(indices[2]),
    PDFbFcnIndex(indices[3]),
    PDFbParIndex(indices[4]);
  fptype
    mean(p[indices[5]]),
    sigma(p[indices[6]]),
    PDFaValue(callFunction(evt, PDFaFcnIndex, PDFaParIndex)),
    PDFbValue(callFunction(evt, PDFbFcnIndex, PDFbParIndex)),
    diff(PDFaValue - PDFbValue - mean);

  // We want to know how far (PDFaValue - PDFbValue) is from 'mean'
  // in units of 'sigma'
  if(diff > Pi)
    diff -= 2.0*Pi;
  if(diff < -Pi)
    diff += 2.0*Pi;
  diff /= sigma;

  fptype ret(EXP(-diff*diff*0.5) * invRootPi / (sigma * root2));
  //printf("PhaseDifferenceConstraintPdf concluded %f sigma from expected (a,b,m,s) = (%f,%f,%f,%f)\n",
  //    FABS(diff), PDFaValue, PDFbValue, mean, sigma);
  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_PhaseDiffConstraint = device_PhaseDiffConstraint; 

__host__ PhaseDifferenceConstraintPdf::PhaseDifferenceConstraintPdf(const std::string &n,
    PdfBase *pdfa,
    PdfBase *pdfb,
    Variable *mean,
    Variable *sigma
    )
  : GooPdf(0, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(pdfa->getFunctionIndex());
  pindices.push_back(pdfa->getParameterIndex());
  pindices.push_back(pdfb->getFunctionIndex());
  pindices.push_back(pdfb->getParameterIndex());
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma));

  // Add as components so that observables and parameters will be registered. 
  components.push_back(pdfa);
  components.push_back(pdfb);

  GET_FUNCTION_ADDR(ptr_to_PhaseDiffConstraint);
  initialise(pindices); 
}

__host__ fptype PhaseDifferenceConstraintPdf::normalise () const
{
  host_normalisation[parameters] = 1.0;
  for(std::vector<PdfBase*>::const_iterator c_iter = components.begin(); c_iter != components.end(); c_iter++)
    (*c_iter)->normalise();
  return 1.0;
}
