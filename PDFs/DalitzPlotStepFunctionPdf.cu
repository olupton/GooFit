#include "DalitzPlotStepFunctionPdf.hh"


EXEC_TARGET fptype device_getHelicityCosine(fptype s12, fptype s13, fptype motherMass, fptype daug1Mass, fptype daug2Mass, fptype daug3Mass)
{
  // m12 and m13 are actually s12, s13 (i.e. they have dimension mass^2
  fptype
    m12(SQRT(s12)),
    E3(((motherMass*motherMass)-s12-(daug3Mass*daug3Mass))/(2.0*m12)),
    //E2((s12-(daug1Mass*daug1Mass)+(daug2Mass*daug2Mass))/(2.0*m12)),
    E1((s12-(daug2Mass*daug2Mass)+(daug1Mass*daug1Mass))/(2.0*m12)),
    p(SQRT(fmax(0.0,(E1*E1)-(daug1Mass*daug1Mass)))),
    q(SQRT(fmax(0.0,(E3*E3)-(daug3Mass*daug3Mass))));
  
  // I'm sure this can be done more neatly
  // cos = (2.0*E1*E3 - (m13*m13) + (m1*m1) + (m3*m3)) / (2.0*p*q);
  // this is cosine of the angle between 1,3 in the 1,2 rest frame
  fptype cos((2.0*E1*E3 - s13 + (daug1Mass*daug1Mass) + (daug3Mass*daug3Mass))/(2.0*p*q));
  if(cos > 1.0 || cos < -1.0)
  {
    printf("device_getHelicityCosine(%f, %f, %f, %f, %f, %f) produced cos = %f\n", s12, s13, motherMass, daug1Mass,
           daug2Mass, daug3Mass, cos);
    cos = 2.0;
  }
  return cos;
}

EXEC_TARGET fptype device_DalitzPlotStepFunction(fptype* evt, fptype* p, unsigned int* indices)
{
  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3];
  fptype maxAbsCosine=functorConstants[indices[1] + 4];

  fptype m12 = evt[indices[2 + indices[0] + 0]];
  fptype m13 = evt[indices[2 + indices[0] + 1]];
  fptype ret = inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass);
  
  if(ret)
  {
    fptype
      m23(motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13),
      maxCos(fmax(fmax(FABS(device_getHelicityCosine(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)),
                       FABS(device_getHelicityCosine(m13, m23, motherMass, daug3Mass, daug1Mass, daug2Mass))),
                  FABS(device_getHelicityCosine(m23, m12, motherMass, daug2Mass, daug3Mass, daug1Mass))));
    if(maxCos > 1.0)
      printf("m12 = %f, m13 = %f, maxCos = %f\n", m12, m13, maxCos);
    
    if(maxCos > maxAbsCosine)
      ret = 0.0;
  }

  //printf("m12 = %f, m13 = %f, ret = %f\n", m12, m13, ret);
  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_DalitzPlotStepFunction = device_DalitzPlotStepFunction;

__host__ DalitzPlotStepFunctionPdf::DalitzPlotStepFunctionPdf(std::string n,
                                                              Variable *m12,
                                                              Variable *m13,
                                                              DecayInfo *decayInfo,
                                                              fptype maxAbsCosine)
  : GooPdf(0, n)
  , _m12(m12)
  , _m13(m13)
  , _maxAbsCosine(maxAbsCosine)
{
  registerObservable(_m12);
  registerObservable(_m13);

  fptype decayConstants[5];
  std::vector<unsigned int> pindices;
  pindices.push_back(registerConstants(5));
  decayConstants[0] = decayInfo->motherMass;
  decayConstants[1] = decayInfo->daug1Mass;
  decayConstants[2] = decayInfo->daug2Mass;
  decayConstants[3] = decayInfo->daug3Mass;
  decayConstants[4] = maxAbsCosine;

  MEMCPY_TO_SYMBOL(functorConstants, decayConstants, 5*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);
  GET_FUNCTION_ADDR(ptr_to_DalitzPlotStepFunction);
  initialise(pindices);

  //addSpecialMask(PdfBase::ForceSeparateNorm);
}
