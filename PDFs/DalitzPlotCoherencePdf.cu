#include "DalitzPlotCoherencePdf.hh"
#include <complex>

MEM_DEVICE devcomplex<fptype>* device_integrals[10];

EXEC_TARGET unsigned int get1Dindex(unsigned int i, unsigned int j, unsigned int nResA, unsigned int nResB)
{
  return i*nResB + j;
}

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
    // 2 is the 'number of resonances' which is always zero so we can reuse DalitzPlotPdf code
    pdfb_index(indices[3]),
    // 4 is the efficiency function index
    // 5 is the efficiency function parameter index
    //cacheToUse(indices[6]),
    *pdfa_indices(paramIndices + pdfa_index),
    *pdfb_indices(paramIndices + pdfb_index);;
  // paramIndices + pdfa_index = what PDF 'a' would get as 'indices' in this function

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
    &pdfa_index(indices[1]),
    &pdfb_index(indices[3]),
    *pdfa_indices(paramIndices + pdfa_index),
    *pdfb_indices(paramIndices + pdfb_index),
    &cacheToUse(indices[6]),
    &nResA(pdfa_indices[2]),
    &nResB(pdfb_indices[2]),
    &numflagpairs(indices[9]),
    evtnum(evt[indices[indices[0] + 2 + 2]] + 0.5);
  fptype
    &coherence_constraint(p[indices[7]]),
    &coherence_error(p[indices[8]]);

  devcomplex<fptype> coherence;
  for(unsigned int i = 0; i < nResA; ++i)
  {
    for(unsigned int j = 0; j < nResB; ++j)
    {
      unsigned int
        paramIndex_i(parIndexFromResIndex_DP(i)),
        paramIndex_j(parIndexFromResIndex_DP(j)),
        integral_index(get1Dindex(i, j, nResA, nResB));
      fptype
        amp_real_i(p[pdfa_indices[paramIndex_i+0]]), // these aren't necessarily literally
        amp_imag_i(p[pdfa_indices[paramIndex_i+1]]), // real and imaginary parts
        amp_real_j(p[pdfb_indices[paramIndex_j+0]]),
        amp_imag_j(p[pdfb_indices[paramIndex_j+1]]);
      devcomplex<fptype>
        amp_i(makedevcomplex(amp_real_i, amp_imag_i)),
        amp_j(makedevcomplex(amp_real_j, amp_imag_j)),
        cached_integral(device_integrals[cacheToUse][integral_index]);
      coherence += amp_i * amp_j * cached_integral;
    }
  }

  // now need to divide by the square root of the product of the integrals over the 
  // two separate Dalitz plot amp2 (over the region where the efficiency is nonzero)
  // should be able to get this as the square root of the product of device-side
  // normalisation factors
  // Doing this relies on ALL the normalisation factors being copied to the device before
  // ANY of the PDFs are evaluated on the data points.
  coherence *= SQRT(normalisationFactors[pdfa_index] * normalisationFactors[pdfb_index]);
 
  fptype ret(1.0);
  bool stilllooking(true);
  for(unsigned int nflag = 0; nflag < numflagpairs && stilllooking; ++nflag)
  {
    unsigned int
      &pair_evtnum(indices[10 + nflag*2 + 0]),
      &pair_flagval(indices[10+ nflag*2 + 1]);
    if(evtnum == pair_evtnum)
    {
      if(pair_flagval == DalitzPlotCoherencePdf::GAUSSIAN_AMPLITUDE_CONSTRAINT)
      {
        fptype coherence_norm(norm(coherence));
        ret = EXP(-(coherence_norm - coherence_constraint)*(coherence_norm - coherence_constraint)/(2.0*coherence_error*coherence_error))*invRootPi/(root2*coherence_error);
      }
      else if(pair_flagval == DalitzPlotCoherencePdf::RAW_AMPLITUDE_VALUE)
      {
        ret = norm(coherence);
      }
      else if(pair_flagval == DalitzPlotCoherencePdf::RAW_PHASE_VALUE)
      {
        ret = coherence.arg();
      }
      else
      {
        printf("DalitzPlotCoherencePdf: don't know how to interpret flag %d for event %d (from %f)\n", pair_flagval, pair_evtnum, evt[indices[indices[0] + 2 + 2]]);
      }
      stilllooking = false;
    }
  }

  if(stilllooking)
  {
    printf("DalitzPlotCoherencePdf: couldn't figure out what we're supposed to be doing for event %d (from %f)\n", evtnum, evt[indices[indices[0] + 2 + 2]]);
  }

  printf("DalitzPlotCoherencePdf: calculated coherence = (%f, %f). nResA = %d, nResB = %d\n", norm(coherence), coherence.arg(), nResA, nResB);

  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_DalitzPlotCoherencePdf = device_DalitzPlotCoherence; 

__host__ DalitzPlotCoherencePdf::DalitzPlotCoherencePdf (
    std::string n, 
    Variable* m12, 
    Variable* m13, 
    Variable* evtnum,
    Variable* coherence_constraint,
    Variable* coherence_error,
    DalitzPlotPdf *pdfa_,
    DalitzPlotPdf *pdfb_,
    GooPdf* efficiency,
    const evtNumFlagPairVec &flags)
  : GooPdf(0, n) 
  , _m12(m12)
  , _m13(m13)
  , pdfa(pdfa_)
  , pdfb(pdfb_)
  , pdfab(2, NULL)
  , dalitzNormRange(0)
  , integrals(0)
  , integrators(0)
  , cacheToUse(0)
  , nResA(pdfa->getDecayInfo()->resonances.size())
  , nResB(pdfb->getDecayInfo()->resonances.size())
  , host_integrals(0)
{
  registerObservable(_m12);
  registerObservable(_m13);
  registerObservable(evtnum);
  pdfab[0] = pdfa;
  pdfab[1] = pdfb;

  std::vector<unsigned int> pindices;
  static int cacheCount = 0;
  cacheToUse = cacheCount++;

  pindices.push_back(pdfa->getParameterIndex()); // need to know where to look for the DalitzPlotPdf parameters (this would be the constants index)
  pindices.push_back(0); // number of resonances, dummy to align us with DalitzPlotPdf
  pindices.push_back(pdfb->getParameterIndex()); // this would be 'cacheToUse' 
  // the next one has to stay as the 4th, and the previous-but-one one has to be zero, so SpecialResonanceIntegrator doesn't get confused
  pindices.push_back(efficiency->getFunctionIndex()); // we need our own "efficiency" because it's how we'll define e.g. the K*(892) mass region
  pindices.push_back(efficiency->getParameterIndex());
  pindices.push_back(cacheToUse);
  pindices.push_back(registerParameter(coherence_constraint));
  pindices.push_back(registerParameter(coherence_error));
  pindices.push_back(flags.size());
  for(evtNumFlagPairVec::const_iterator flag_iter = flags.begin(); flag_iter != flags.end(); flag_iter++)
  {
    pindices.push_back(flag_iter->first); // event number
    pindices.push_back(flag_iter->second); // flag; implicit enum -> integer conversion
  }
  components.push_back(efficiency); // I guess this is needed...

  GET_FUNCTION_ADDR(ptr_to_DalitzPlotCoherencePdf);
  initialise(pindices);

  // We need to keep track with which resonances are changing in the 2 PDFs
  redoIntegral = new bool*[pdfab.size()]; // Flag which of the terms need re-calculating

  for(std::size_t i = 0; i < pdfab.size(); ++i)
  {
    redoIntegral[i] = new bool[nResA];
    for(std::size_t j = 0; j < nResB; ++j)
      redoIntegral[i][j] = true;
  }

  integrals    = new DEVICE_VECTOR<devcomplex<fptype> >(nResA * nResB);
  host_integrals=new devcomplex<fptype>*[nResA];
  integrators  = new SpecialResonanceCoherenceIntegrator**[nResA];

  for (int i = 0; i < nResB; ++i)
  {
    integrators[i]    = new SpecialResonanceCoherenceIntegrator*[nResB];
    host_integrals[i] = new devcomplex<fptype>[nResB];
    for (int j = 0; j < nResB; ++j)
    {
      integrators[i][j] = new SpecialResonanceCoherenceIntegrator(parameters, i, j); 
      host_integrals[i][j].real = host_integrals[i][j].imag = 0.0;
    }
  }

  void *dummy(thrust::raw_pointer_cast(integrals->data()));
  MEMCPY_TO_SYMBOL(device_integrals, &dummy, sizeof(devcomplex<fptype>*), cacheToUse*sizeof(devcomplex<fptype>*), cudaMemcpyHostToDevice);

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

  // store this factor as part of the device-side integrals
  fptype binSizeFactor(
      ((_m12->upperlimit - _m12->lowerlimit) / _m12->numbins)
      * ((_m13->upperlimit - _m13->lowerlimit) / _m13->numbins));

  for(int i = 0; i < nResA; ++i)
  {
    for (int j = 0; j < nResB; ++j)
    {
      if((!redoIntegral[0][i]) && (!redoIntegral[1][j]))
        continue; 
      devcomplex<fptype> dummy(0, 0);
      thrust::plus<devcomplex<fptype> > complexSum;
      std::size_t index(i*pdfb->getDecayInfo()->resonances.size() + j);
      (*integrals)[index] = thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress)),
						      thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress)),
						      *(integrators[i][j]), 
						      dummy, 
						      complexSum) * binSizeFactor;
    }
  }

  host_normalisation[parameters] = 1.0;
  return 1.0;
}

__host__ void DalitzPlotCoherencePdf::copyIntegralsToHost ()
{
  for(int i = 0; i < nResA; ++i)
    for (int j = 0; j < nResB; ++j)
      host_integrals[i][j] = (*integrals)[i*nResB + j];
}

SpecialResonanceCoherenceIntegrator::SpecialResonanceCoherenceIntegrator (int pIdx, unsigned int ri, unsigned int rj) 
  : SpecialResonanceIntegrator(pIdx, ri, rj)
{
}

EXEC_TARGET devcomplex<fptype> SpecialResonanceCoherenceIntegrator::devicefunction(fptype m12, fptype m13, int res_i, int res_j, fptype* p, unsigned int* indices) const
{
  return device_DalitzCoherence_calcIntegrals(m12, m13, res_i, res_j, p, indices);
}
