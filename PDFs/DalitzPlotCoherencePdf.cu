#include "DalitzPlotCoherencePdf.hh"
#include <complex>

EXEC_TARGET devcomplex<fptype> device_DalitzCoherence_calcIntegrals (fptype m12, fptype m13, int res_i, int res_j, fptype* p, unsigned int* indices) {
  // Calculates BW_i(m12, m13) * BW_j^*(m12, m13). 
  // This calculation is in a separate function so
  // it can be cached. Note that this function expects
  // to be called on a normalisation grid, not on 
  // observed points, that's why it doesn't use 
  // cResonances. No need to cache the values at individual
  // grid points - we only care about totals. 
  unsigned int
    pdfa_index(indices[1]),
    pdfb_index(indices[2]);

  // paramIndices + pdfa_index = what PDF 'a' would get as 'indices' in this function
  unsigned int
    *pdfa_indices(paramIndices + pdfa_index),
    *pdfb_indices(paramIndices + pdfb_index);

  fptype
    motherMass_a(functorConstants[pdfa_indices[1] + 0]),
    motherMass_b(functorConstants[pdfb_indices[1] + 0]),
    daug1Mass_a(functorConstants[pdfa_indices[1] + 1]),
    daug1Mass_b(functorConstants[pdfb_indices[1] + 1]),
    daug2Mass_a(functorConstants[pdfa_indices[1] + 2]),
    daug2Mass_b(functorConstants[pdfb_indices[1] + 2]),
    daug3Mass_a(functorConstants[pdfa_indices[1] + 3]),
    daug3Mass_b(functorConstants[pdfb_indices[1] + 3]),
    eps(1e-4);
  if( FABS(motherMass_a-motherMass_b)/motherMass_a > eps
      || FABS(daug1Mass_a-daug1Mass_b)/daug1Mass_a > eps
      || FABS(daug2Mass_a-daug2Mass_b)/daug2Mass_a > eps
      || FABS(daug3Mass_a-daug3Mass_b)/daug3Mass_a > eps)
  {
    // These ought to be the same...
    printf("motherMass = (%.4f, %.4f), daug1Mass = (%.4f, %.4f), daug2Mass = (%.4f, %.4f), daug3Mass = (%.4f, %.4f)\n",
        motherMass_a, motherMass_b, daug1Mass_a, daug1Mass_b, daug2Mass_a, daug2Mass_b, daug3Mass_a, daug3Mass_b);
  }

  devcomplex<fptype> ret; 
  if(!inDalitz(m12, m13, motherMass_a, daug1Mass_a, daug2Mass_a, daug3Mass_a))
    return ret;
  fptype m23 = motherMass_a*motherMass_a + daug1Mass_a*daug1Mass_a + daug2Mass_a*daug2Mass_a + daug3Mass_a*daug3Mass_a - m12 - m13; 

  int parameter_i = parIndexFromResIndex_DP(res_i);
  unsigned int functn_i = pdfa_indices[parameter_i+2];
  unsigned int params_i = pdfa_indices[parameter_i+3];
  ret = getResonanceAmplitude(m12, m13, m23, functn_i, params_i);

  // And now from the other DalitzPlotPdf
  int parameter_j = parIndexFromResIndex_DP(res_j);
  unsigned int functn_j = pdfb_indices[parameter_j+2];
  unsigned int params_j = pdfb_indices[parameter_j+3];
  ret *= getResonanceAmplitude(m12, m13, m23, functn_j, params_j);

  return ret; 
}

EXEC_TARGET fptype device_DalitzPlotCoherence (fptype* evt, fptype* p, unsigned int* indices) {
  // This is the actual PDF function. It needs to return a Gaussian PDF constraining the 
  // measured coherence values to the calculated values
  
  // These should be what pdfa and pdfb get as 'indices'
  //
  unsigned int
    pdfa_index(indices[1]),
    pdfb_index(indices[2]),
    *pdfa_indices(paramIndices + pdfa_index),
    *pdfb_indices(paramIndices + pdfb_index);

  devcomplex<fptype> totalAmp(0, 0);
  unsigned int numResonances = indices[2]; 
  unsigned int cacheToUse    = indices[3]; 

  for (int i = 0; i < numResonances; ++i) {
    int paramIndex  = parIndexFromResIndex_DP(i);
    //fptype amp_real = p[indices[paramIndex+0]];
    //fptype amp_imag = p[indices[paramIndex+1]];
    
    // behaviour depends on compile flag, between cartesian and polar interpretation of parameters
    devcomplex<fptype> amp(makedevcomplex(p[indices[paramIndex+0]], p[indices[paramIndex+1]]));
    
    devcomplex<fptype> matrixelement;
    //((cResonances[cacheToUse][evtNum*numResonances + i]).real,
	//			     (cResonances[cacheToUse][evtNum*numResonances + i]).imag); 
    matrixelement *= amp;//.multiply(amp_real, amp_imag);
    totalAmp += matrixelement; 
  } 
   
  fptype ret = norm2(totalAmp); 
  int effFunctionIdx = parIndexFromResIndex_DP(numResonances); 
  fptype eff = callFunction(evt, indices[effFunctionIdx], indices[effFunctionIdx + 1]); 
  ret *= eff;

  //printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real, totalAmp.imag); 

  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_DalitzPlotCoherencePdf = device_DalitzPlotCoherence; 

__host__ DalitzPlotCoherencePdf::DalitzPlotCoherencePdf (
    std::string n, 
    Variable* m12, 
    Variable* m13, 
    DalitzPlotPdf *pdfa_,
    DalitzPlotPdf *pdfb_,
    GooPdf* efficiency)
  : GooPdf(0, n) 
  , _m12(m12)
  , _m13(m13)
  , pdfa(pdfa_)
  , pdfb(pdfb_)
  , pdfab(2, NULL)
  , dalitzNormRange(0)
  , integrals(0)
  , integrators(0)
{
  registerObservable(_m12);
  registerObservable(_m13);
  pdfab[0] = pdfa;
  pdfab[1] = pdfb;

  std::vector<unsigned int> pindices;
  // pindices.push_back(registerConstants(6));

  //pindices.push_back(decayInfo->resonances.size()); 
  //pindices.push_back(cacheToUse); 

  //for (std::vector<ResonancePdf*>::iterator res = decayInfo->resonances.begin(); res != decayInfo->resonances.end(); ++res) {
  //  pindices.push_back(registerParameter((*res)->amp_real));
  //  pindices.push_back(registerParameter((*res)->amp_imag));
  //  pindices.push_back((*res)->getFunctionIndex());
  //  pindices.push_back((*res)->getParameterIndex());
  //  (*res)->setConstantIndex(cIndex); 
  //  components.push_back(*res);
  //}

  pindices.push_back(pdfa->getParameterIndex()); // need to know where to look for the DalitzPlotPdf parameters
  pindices.push_back(pdfb->getParameterIndex());
  pindices.push_back(efficiency->getFunctionIndex()); // we need our own "efficiency" because it's how we'll define e.g. the K*(892) mass region
  pindices.push_back(efficiency->getParameterIndex());
  components.push_back(efficiency); // I guess this is needed...

  GET_FUNCTION_ADDR(ptr_to_DalitzPlotCoherencePdf);
  initialise(pindices);

  // We need to keep track with which resonances are changing in the 2 PDFs
  redoIntegral = new bool*[pdfab.size()]; // Flag which of the terms need re-calculating

  for(std::size_t i = 0; i < pdfab.size(); ++i)
  {
    redoIntegral[i] = new bool[pdfab[i]->getDecayInfo()->resonances.size()];
    for(std::size_t j = 0; j < pdfab[i]->getDecayInfo()->resonances.size(); ++j)
      redoIntegral[i][j] = true;
  }

  integrals    = new devcomplex<fptype>**[pdfa->getDecayInfo()->resonances.size()];
  integrators  = new SpecialResonanceCoherenceIntegrator**[pdfa->getDecayInfo()->resonances.size()];

  for (int i = 0; i < pdfa->getDecayInfo()->resonances.size(); ++i)
  {
    integrators[i]  = new SpecialResonanceCoherenceIntegrator*[pdfb->getDecayInfo()->resonances.size()];
    integrals[i]    = new devcomplex<fptype>*[pdfb->getDecayInfo()->resonances.size()];
    
    for (int j = 0; j < pdfb->getDecayInfo()->resonances.size(); ++j) {
      integrals[i][j]   = new devcomplex<fptype>(0, 0); 
      integrators[i][j] = new SpecialResonanceCoherenceIntegrator(parameters, i, j); 
    }
  }

  // I think this stops GooFit trying to be too clever.
  addSpecialMask(PdfBase::ForceSeparateNorm); 
}

__host__ fptype DalitzPlotCoherencePdf::normalise () const {
  recursiveSetNormalisation(1); // Not going to normalise efficiency, 
  // so set normalisation factor to 1 so it doesn't get multiplied by zero. 
  // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency, 
  // don't get zeroes through multiplying by the normFactor. 
  MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice); 

  int totalBins = _m12->numbins * _m13->numbins;
  if (!dalitzNormRange) {
    gooMalloc((void**) &dalitzNormRange, 6*sizeof(fptype));
  
    fptype* host_norms = new fptype[6];
    host_norms[0] = _m12->lowerlimit;
    host_norms[1] = _m12->upperlimit;
    host_norms[2] = _m12->numbins;
    host_norms[3] = _m13->lowerlimit;
    host_norms[4] = _m13->upperlimit;
    host_norms[5] = _m13->numbins;
    MEMCPY(dalitzNormRange, host_norms, 6*sizeof(fptype), cudaMemcpyHostToDevice);
    delete[] host_norms; 
  }

  for(std::size_t ab = 0; ab < pdfab.size(); ++ab)
  {
    DecayInfo *decayInfo(pdfab[ab]->getDecayInfo());
    for (unsigned int i = 0; i < decayInfo->resonances.size(); ++i)
    {
      //redoIntegral[ab][i] = forceRedoIntegrals; 
      if(!(decayInfo->resonances[i]->parametersChanged()))
        continue;
      redoIntegral[ab][i] = true; 
      decayInfo->resonances[i]->storeParameters();
    }
  }

  //forceRedoIntegrals = false; 

  // Only do this bit if masses or widths have changed.  
  thrust::constant_iterator<fptype*> arrayAddress(dalitzNormRange); 
  thrust::counting_iterator<int> binIndex(0); 

  // NB, SpecialResonanceCalculator assumes that fit is unbinned! 
  // And it needs to know the total event size, not just observables
  // for this particular PDF component. 
  //thrust::constant_iterator<fptype*> dataArray(dev_event_array); 
  //thrust::constant_iterator<int> eventSize(totalEventSize);
  //thrust::counting_iterator<int> eventIndex(0); 

  for (int i = 0; i < pdfa->getDecayInfo()->resonances.size(); ++i)
  {
    //if (redoIntegral[i]) {
    //  thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
//			thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
//			strided_range<DEVICE_VECTOR<devcomplex<fptype> >::iterator>(cachedWaves->begin() + i, 
//										    cachedWaves->end(), 
//										    decayInfo->resonances.size()).begin(), 
//			*(calculators[i]));
  //  }
    
    // Possibly this can be done more efficiently by exploiting symmetry? 
    for (int j = 0; j < pdfb->getDecayInfo()->resonances.size(); ++j) {
      if ((!redoIntegral[0][i]) && (!redoIntegral[1][j])) continue; 
      devcomplex<fptype> dummy(0, 0);
      thrust::plus<devcomplex<fptype> > complexSum; 
      (*(integrals[i][j])) = thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress)),
						      thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress)),
						      *(integrators[i][j]), 
						      dummy, 
						      complexSum); 
    }
  }      

  // End of time-consuming integrals. 
  complex<fptype> sumIntegral(0, 0);
  for (unsigned int i = 0; i < pdfa->getDecayInfo()->resonances.size(); ++i) {
    int param_i = parameters + resonanceOffset_DP + resonanceSize*i; 
    complex<fptype> amplitude_i(makecomplex(host_params[host_indices[param_i]], host_params[host_indices[param_i + 1]]));
    for (unsigned int j = 0; j < pdfb->getDecayInfo()->resonances.size(); ++j) {
      int param_j = parameters + resonanceOffset_DP + resonanceSize*j; 
      complex<fptype> amplitude_j(conj(makecomplex(host_params[host_indices[param_j]], host_params[host_indices[param_j + 1]])));
      // Notice complex conjugation

      sumIntegral += (amplitude_i * amplitude_j * complex<fptype>((*(integrals[i][j])).real, (*(integrals[i][j])).imag)); 
    }
  }

  fptype ret = real(sumIntegral); // That complex number is a square, so it's fully real
  double binSizeFactor = 1;
  binSizeFactor *= ((_m12->upperlimit - _m12->lowerlimit) / _m12->numbins);
  binSizeFactor *= ((_m13->upperlimit - _m13->lowerlimit) / _m13->numbins);
  ret *= binSizeFactor;

  host_normalisation[parameters] = 1.0/ret;
  return (fptype) ret; 
}

SpecialResonanceCoherenceIntegrator::SpecialResonanceCoherenceIntegrator (int pIdx, unsigned int ri, unsigned int rj) 
  : SpecialResonanceIntegrator(pIdx, ri, rj)
{
}

EXEC_TARGET devcomplex<fptype> SpecialResonanceCoherenceIntegrator::devicefunction(fptype m12, fptype m13, int res_i, int res_j, fptype* p, unsigned int* indices) const
{
  return device_DalitzCoherence_calcIntegrals(m12, m13, res_i, res_j, p, indices);
}
