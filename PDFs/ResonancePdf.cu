#include "ResonancePdf.hh" 

__device__ fptype twoBodyCMmomSq(fptype rMassSq, fptype d1m, fptype d2m)
{
  // Define something to return the momentum squared so we can defer
  // deciding what to do with unphysical momenta...
  return 0.25 * rMassSq * (1.0 - ((d1m+d2m)*(d1m+d2m) / rMassSq)) * (1.0 - ((d1m - d2m)*(d1m - d2m)/rMassSq));
}

__device__ fptype twoBodyCMmom (fptype rMassSq, fptype d1m, fptype d2m) {
  // For A -> B + C, calculate momentum of B and C in rest frame of A. 
  // PDG 38.16.

  //fptype kin1 = 1 - POW(d1m+d2m, 2) / rMassSq;
  //if (kin1 > 0)
  //  kin1 = SQRT(kin1);
  //else
  //  kin1 = 1;
  
  //fptype kin2 = 1 - POW(d1m-d2m, 2) / rMassSq;
  //if (kin2 > 0)
  //  kin2 = SQRT(kin2);
  //else
  //  kin2 = 1;

  //return 0.5*SQRT(rMassSq)*kin1*kin2; 
  return SQRT(twoBodyCMmomSq(rMassSq, d1m, d2m));
}

__device__ fptype bachelorMomSq(fptype otherMass, fptype motherMass, fptype bachelorMass)
{
  fptype
    motherMassSq(motherMass*motherMass),
    massSumSq((otherMass + bachelorMass)*(otherMass + bachelorMass)),
    massDiffSq((otherMass- bachelorMass)*(otherMass - bachelorMass));
  return 0.25 * (motherMassSq - massSumSq) * (motherMassSq - massDiffSq) / (otherMass * otherMass);
}

// For D -> (R -> AB)C calculate momentum of C (and D) in rest frame of R (A+B)
// otherMass is m_R, motherMass is m_D and bachelorMass is m_C
__device__ fptype bachelorMom(fptype otherMass, fptype motherMass, fptype bachelorMass)
{
  // This will sometimes get called with masses which are kinematically forbidden (i.e. m_R + m_C > m_D)
  // In these cases we just want to return something of the right order of magnitude.
  // The way in which we do this is just copied from the twoBodyCMmom() function above.
  //fptype
  //kin1(1.0 - POW((otherMass + bachelorMass)/motherMass, 2.0)),
  //  kin2(1.0 - POW((otherMass - bachelorMass)/motherMass, 2.0));
  //if(kin1 >= 0.0)
  //  kin1 = SQRT(kin1);
  //else
  //  kin1 = SQRT(-kin1);//kin1 = 1.0;
  //kin1 = SQRT(FMAX(0.0, kin1));//(FABS(kin1));
  
  // BIG FAT WARNING
  // try and emulate what MINT does (did) here, by using the absolute value of 'kin2' if it is negative
  // all choices here are quite arbitrary....
  //if(kin2 >= 0.0)
  //  kin2 = SQRT(kin2);
  //else
  //  kin2 = SQRT(-kin2);//1.0;
  //kin2 = SQRT(FMAX(0.0, kin2));
  
  //return FMAX(1e-6, 0.5 * POW(motherMass, 2.0) * kin1 * kin2 / otherMass);
  return SQRT(bachelorMomSq(otherMass, motherMass, bachelorMass));
}

__device__ fptype dampingFactorSquare (fptype cmmom, int spin, fptype mRadius)
{
  fptype square = mRadius*mRadius*cmmom*cmmom;
  fptype dfsq = 1 + square; // This accounts for spin 1
  if (2 == spin)
    dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.

  // Spin 3 and up not accounted for. 
  return dfsq; 
}

__device__ fptype unNormalisedDampingFactorSquare(fptype measured_momentum, int spin, fptype meson_radius)
{
  printf("WARNING: you are using an incomplete, untested damping factor implementation!\n");
  if(spin == 0)
    return 1.0; 
  fptype
    square(measured_momentum * measured_momentum * meson_radius * meson_radius),
    ret(0.0);
  if(spin == 1)
    ret = 2.0 * square / (1.0 + square);
  if(spin == 2)
    ret = 13.0 * square * square / ( ((square - 3.0)*(square - 3.0)) + 9.0*square);
  return ret;
}

__device__ fptype dampingFactorRatioSquare(fptype nummom, fptype denmom, int spin, fptype meson_radius, bool useUnNormalisedFactor)
{
  if(useUnNormalisedFactor)
    return unNormalisedDampingFactorSquare(denmom, spin, meson_radius);
  else
    return dampingFactorSquare(nummom, spin, meson_radius) / dampingFactorSquare(denmom, spin, meson_radius);
}

// For spin 1:
// (p_D + p_C)_mu( -g_munu + P_mu P_nu * massFactor)(p_B - p_A)_nu where P = p_A + p_B
// = -(p_A + p_C + p_B + p_C).(p_B + p_C - p_A - p_C) + massFactor*(p_D + p_C).(p_D - p_C)*(p_A + p_B).(p_B - p_A)
// = -(s_BC - s_AC) + massFactor*(m_D^2 - m_C^2)*(m_B^2 - m_A^2)
// = -{ (s_BC - s_AC) + massFactor*(m_D^2 - m_C^2)*(m_A^2 - m_B^2) }
//
// For spin 2:
// (p_D + p_C)_mu(p_D + p_C)_nu T_mu,nu,alpha,beta (p_B - p_A)_alpha(p_B - p_A)_beta
// where T_mu,nu,alpha,beta = (1/2)*(T_mu,alpha*T_nu,beta + T_mu,beta*T_nu,alpha) - (1/3)T_mu,nu*T_alpha,beta
// and T_mu,nu = -g_mu,nu + massFactor*P_mu,P_nu
__device__ fptype spinFactorABC(unsigned int spin, fptype motherMass, fptype _mA, fptype _mB, fptype _mC, fptype _mAB, fptype _mAC, fptype _mBC, fptype massFactor)
{
  if(spin == 0)
    return 1.0;
  fptype sFactor(-1.0);
  sFactor *= ((_mBC - _mAC) + (massFactor*(motherMass*motherMass - _mC*_mC)*(_mA*_mA-_mB*_mB)));
  if (2 == spin) {
    sFactor *= sFactor;
    fptype extraterm = ((_mAB-(2*motherMass*motherMass)-(2*_mC*_mC))+massFactor*pow((motherMass*motherMass-_mC*_mC),2));
    extraterm *= ((_mAB-(2*_mA*_mA)-(2*_mB*_mB))+massFactor*pow((_mA*_mA-_mB*_mB),2));
    extraterm /= 3;
    sFactor -= extraterm;
  }
  return sFactor;
}

__device__ fptype spinFactor (unsigned int spin, fptype motherMass, fptype daug1Mass, fptype daug2Mass, fptype daug3Mass, fptype m12, fptype m13, fptype m23, unsigned int cyclic_index, fptype resMassSq, bool usenominal = true)
{
  if (0 == spin)
    return 1.0; // Should not cause branching since every thread evaluates the same resonance at the same time. 
  /*
  // Copied from BdkDMixDalitzAmp
   
  fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass)); 
  fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug3Mass)); 
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass)); 
    
  fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m12 : m12)); 
  fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m23 : m13)); 
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 

  // The above, collapsed into single tests where possible. 
  fptype _mA = (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass);
  fptype _mB = (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass); 
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass)); 

  fptype _mAC = (PAIR_23 == cyclic_index ? m13 : m23);
  fptype _mBC = (PAIR_12 == cyclic_index ? m13 : m12);
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 
  */

  // Copied from EvtDalitzReso, with assumption that pairAng convention matches pipipi0 from EvtD0mixDalitz.
  // Again, all threads should get the same branch. 
  fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass));
  fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
  fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m23 : m12)); 
  fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m12 : m13)); 
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  return spinFactorABC(spin, motherMass, _mA, _mB, _mC, _mAB, _mAC, _mBC, usenominal ? 1.0/resMassSq : 1.0/_mAB);
}

__device__ devcomplex<fptype> plainBW (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  fptype motherMass             = functorConstants[indices[1]+0];
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];
  fptype meson_radius           = functorConstants[indices[1]+4];
  fptype mother_meson_radius    = functorConstants[indices[1]+5];
 
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  const unsigned int &spin                = indices[4];
  const unsigned int &cyclic_index        = indices[5];
  const unsigned int &use_nominal_resmass = indices[6];
  const unsigned int &use_unnorm_dampingf = indices[7];
  //printf("%.3f %.3f %.3f %.3f %.3f %.3f\n", motherMass, daug1Mass, daug2Mass, daug3Mass, meson_radius, mother_meson_radius);
  //printf("%.3f %.3f %d %d %d %d\n", resmass, reswidth, spin, cyclic_index, use_nominal_resmass, use_unnorm_dampingf);

  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype frFactor(1.0);

  resmass *= resmass; 
  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2). 
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, 
					    (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), 
					    (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, 
					    (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), 
              (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));
  
  

  if (0 != spin)
    frFactor = dampingFactorRatioSquare(nominalDaughterMoms, measureDaughterMoms, spin, meson_radius, use_unnorm_dampingf);
 
  // RBW evaluation
  fptype A = (resmass - rMassSq); 
  fptype prat(measureDaughterMoms / nominalDaughterMoms);
  fptype B = resmass*reswidth * frFactor / SQRT(rMassSq);
  for(unsigned int i = 0; i < ((2*spin) + 1); ++i)
    B *= prat;
  fptype C = 1.0 / (A*A + B*B); 
  devcomplex<fptype> ret(A*C, B*C);
  // A + iB / (A^2 + B^2)
  
  // Don't want the F_D penetration factor in the mass-dependent width
  if(0 != spin)
  {
    fptype bachelorMass(PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
    frFactor *= dampingFactorRatioSquare(
        SQRT(FABS(bachelorMomSq(SQRT(resmass), motherMass, bachelorMass))),
        SQRT(FABS(bachelorMomSq(SQRT(rMassSq), motherMass, bachelorMass))),
        spin, mother_meson_radius, use_unnorm_dampingf); // using nominal/measured
  }

  ret *= SQRT(frFactor);
  fptype spinF = spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index, resmass, use_nominal_resmass); 
  ret *= spinF; 

  //if(SQRT(resmass) < 0.900)
  //{
  //  printf("%.3f %.3f %.3f\n", frFactor, spinF, prat);
  // }

  return ret; 
}

// LASS shape for including a nonresonant effective-range component
// this is parameterised like BaBar arXiv:1004.5053
// http://arxiv.org/pdf/1004.5053v3.pdf
//
// but with some modifications to correct for the variations in normalisation.
// Their BW propagator has -- as most do -- dimensions of inverse mass squared, but their
// LASS propagator is dimensionless so we need to add the appropriate dimensioned
// constants to make these match when we sum them together
//
// The shape is supposed to be:
//
// sin(delta_R)*exp(i(delta_R + phi_R + 2delta_F + 2phi_F)) + F*sin(delta_F + phi_F)*exp(i(delta_F + phi_F))
//
// which we divide by exp(i*phi_R) and replace sin(delta_R)*exp(i*delta_R) by [BW], as evaluated by that implementation
// this latter part is what causes the dimensional problems, but I think it's the simplest way to implement this
//
// This leaves us with:
//
// [BW]*exp(2*i*(delta_F + phi_F)) + F*sin(delta_F + phi_F)*exp(i*(delta_F + phi_F - phi_R))
//
// with some dimensional constants multiplying the second term
// the use of the letter F is inheritied, B for background and N for nonresonant might both make sense...
//
// Update 20/06/2014 to try and match MINT more closely
//
// MINT seems to return (Lass2.C)
//
// (m_Kpi/q) * F*sin(phi_F + delta_F)*exp(i*(phi_F + delta_F)) + R*sin(delta_R) * exp(i*delta_R + 2*i*(delta_F + phi_F) + i*phi_R)
//
//
__device__ devcomplex<fptype> lass(fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  // use the BW implementation for the resonant part...
  devcomplex<fptype> BW_part(plainBW(m12, m13, m23, indices));
  
  // need these to calculate the K,pi momentum in the Kpi frame
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];

  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  //const unsigned int &spin                = indices[4];
  const unsigned int &cyclic_index        = indices[5];
  //const unsigned int &use_nominal_resmass = indices[6];
  //const unsigned int &use_unnorm_dampingf = indices[7];

  // extra LASS parameters
  fptype lass_a                 = cudaArray[indices[8]];
  fptype lass_r                 = cudaArray[indices[9]];
  fptype lass_phi_f             = cudaArray[indices[10]];
  fptype lass_phi_r             = cudaArray[indices[11]];
  fptype lass_F                 = cudaArray[indices[12]];

  //printf("%.3f %.3f %.3f %.3f %.3f %d %.3f %.3f %.3f %.3f %.3f\n",
  //    daug1Mass, daug2Mass, daug3Mass, resmass, reswidth,
  //    cyclic_index, lass_a, lass_r, lass_phi_f, lass_phi_r, lass_F);
 
  fptype rMassSq((PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)));
  fptype q(twoBodyCMmom(rMassSq,
                        (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass),
                        (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass)));
  fptype q0(twoBodyCMmom(resmass*resmass,
                         (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass),
                         (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass)));

  fptype delta_f(atan2(2.0*lass_a*q, 2.0 + (lass_a*lass_r*q*q)));
  fptype angle_combo(delta_f + lass_phi_f);
  
  BW_part *= devcomplex<fptype>(cos((2.0*angle_combo) + lass_phi_r), sin(2.0*(angle_combo) + lass_phi_r));
  BW_part *= reswidth * resmass * resmass / q0;
  
  devcomplex<fptype> nonres_part(lass_F * sin(angle_combo), 0);
  nonres_part *= devcomplex<fptype>(cos(angle_combo), sin(angle_combo));
  //fptype rMass(SQRT((PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23))));
  nonres_part *= SQRT(rMassSq) / q;
  // now correct for the way the other shapes are normalised
  // for whatever reason, Dalitz fitters seem to like dropping constants so that BW shapes have dimensions of inverse mass^2
  // but what we've just calculated is dimensionless
  // the magic quantity is, I think, this...
  //nonres_part *= q0 / (reswidth * resmass * resmass);
  
  return BW_part + nonres_part;
}

__device__ devcomplex<fptype> polylass(fptype m12, fptype m13, fptype m23, unsigned int *indices)
{
  // need these to calculate the K,pi momentum in the Kpi frame
  fptype motherMass             = functorConstants[indices[1]+0];
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];
  
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int cyclic_index     = indices[5];
  // extra LASS parameters
  fptype lass_a                 = cudaArray[indices[6]];
  fptype lass_r                 = cudaArray[indices[7]];
  unsigned int num_poly_coeffs  = indices[8];
  unsigned int formfactor_type  = indices[9];
  // these are stored in cudaArray[indices[9]] .. cudaArray[indices[8 + num_poly_coeffs]]
  
  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype q(twoBodyCMmom(rMassSq,
                        (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass),
                        (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass)));
  fptype q0(twoBodyCMmom(resmass*resmass,
                         (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass),
                         (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass)));
  
  // This is the LASS phase shift
  fptype delta_f(atan2(2.0*lass_a*q, 2.0 + lass_a*lass_r*q*q));
  
  // Calculate the running width
  fptype resrunwidth(reswidth * q * resmass / (q0 * SQRT(rMassSq)));
  // This is the regular BW phase shift
  fptype delta_r(atan2(resrunwidth * resmass, resmass * resmass - rMassSq));

  // To match the plainBW function we would have to return:
  //   (q0 * SQRT(rMassSq) / (reswidth * resmass * resmass * q)) * sin(delta_r) * exp(i*delta_r)
  // we want to do this but with delta_r -> delta_r + delta_f, and then multiply the whole thing by a polynomial
  
  devcomplex<fptype> ret(cos(delta_r + delta_f), sin(delta_r + delta_f));
  ret *= sin(delta_r + delta_f) * /*q0*/ SQRT(rMassSq) / q;//(reswidth * resmass * resmass * q);
  
  fptype poly(1.0);
  
  if(formfactor_type == ResonancePdf::RECURSIVEPOLY)
  {
    fptype coefffornext(1.0);
    fptype expansion_parameter(SQRT(rMassSq) / resmass);
    if(num_poly_coeffs == 0)
    {
      // If we don't have any floating parameters just return 1
    }
    else
    {
      poly = 0.0;

      // This should be a0*x*x + (1-a0)(a1*x + (1-a1))
      for(unsigned int poly_index = 0; poly_index <= num_poly_coeffs; ++poly_index)
      {
        fptype coeff(poly_index == num_poly_coeffs ? 1.0 : cudaArray[indices[10 + poly_index]]);
        poly += pow(expansion_parameter, int(num_poly_coeffs - poly_index)) * coeff * coefffornext;
        coefffornext *= (1.0 - fabs(coeff));
     }
    }
  }
  else if((formfactor_type == ResonancePdf::NORMPOLY) || (formfactor_type == ResonancePdf::NORMEXPPOLY)
      || (formfactor_type == ResonancePdf::CENTRENORMEXPPOLY) || (formfactor_type == ResonancePdf::CENTRENORMEXPPOLYRECURSIVE)
      || (formfactor_type == ResonancePdf::CENTRENORMEXPPOLYDEEPRECURSIVE) || (formfactor_type == ResonancePdf::CENTRENORMEXPPOLYDEEPRECURSIVEALT))
  {
    fptype expansion_parameter(SQRT(rMassSq) / resmass);
    fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass));
    fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
    fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
    fptype centre_parameter(0.5*(motherMass - _mC + _mA + _mB) / resmass);
    fptype norm(1.0), centrenorm(1.0);
    poly = 1.0;
    for(unsigned int poly_index = 0; poly_index < num_poly_coeffs; ++poly_index)
    {
      fptype coeff(cudaArray[indices[10 + poly_index]]);
      if(poly_index > 0 && formfactor_type == ResonancePdf::CENTRENORMEXPPOLYRECURSIVE)
        coeff *= cudaArray[indices[10 + 0]]; // a,b,c -> (a, ab, ac)
      else if((formfactor_type == ResonancePdf::CENTRENORMEXPPOLYDEEPRECURSIVEALT)
          or (formfactor_type == ResonancePdf::CENTRENORMEXPPOLYDEEPRECURSIVE))
      {
        if(num_poly_coeffs != 3)
          printf("Wrong number of parameters for CENTRENORMEXPPOLYDEEPRECURSIVE[ALT]\n");
        fptype
          &c1(cudaArray[indices[10 + 0]]),
          &c2(cudaArray[indices[10 + 1]]),
          &c3(cudaArray[indices[10 + 2]]);
        bool alt(formfactor_type == ResonancePdf::CENTRENORMEXPPOLYDEEPRECURSIVEALT);
        if(poly_index == 0)
        {
          if(alt)
            // neut
            coeff = -0.452*c1 -0.676*c2 - 0.582*c3;
          else
            coeff = -0.460*c1 + 0.702*c2 - 0.543*c3;
        }
        else if(poly_index == 1)
        {
          if(alt)
            // neut
            coeff = 0.776*c1 + 0.0243*c2 - 0.631*c3;
          else
            coeff = 0.776*c1 + 0.197*c2 - 0.631*c3;
        }
        else
        {
          if(alt)
            // neut
            coeff = -0.440*c1 + 0.737*c2 - 0.513*c3;
          else
            coeff = -0.433*c1 - 0.711*c2 - 0.554*c3;
        }
      }
      poly += pow(expansion_parameter, int(poly_index+1)) * coeff;
      norm += coeff;
      centrenorm += pow(centre_parameter, int(poly_index+1)) * coeff;
    }

    if(formfactor_type == ResonancePdf::NORMEXPPOLY)
    {
      poly = EXP(poly - norm); // exponential of the polynomial, divided by the exponential of when the parameter is 1
    }
    else if(formfactor_type == ResonancePdf::CENTRENORMEXPPOLY
        or formfactor_type == ResonancePdf::CENTRENORMEXPPOLYRECURSIVE
        or formfactor_type == ResonancePdf::CENTRENORMEXPPOLYDEEPRECURSIVE
        or formfactor_type == ResonancePdf::CENTRENORMEXPPOLYDEEPRECURSIVEALT)
    {
      poly = EXP(poly - centrenorm);
    }
    else
    {
      poly /= norm;
    }
  }
  else if(formfactor_type == ResonancePdf::POLY)
  {
    fptype expansion_parameter(SQRT(rMassSq) / resmass);
    poly = pow(expansion_parameter, int(num_poly_coeffs)); // 2 coefficents: a + bx + xx
    for(unsigned int poly_index = 0; poly_index < num_poly_coeffs; ++poly_index)
      poly += pow(expansion_parameter, int(poly_index)) * cudaArray[indices[10 + poly_index]];
  }
  else if(formfactor_type == ResonancePdf::SENSIBLEPOLY)
  {
    // this expected to have as many coefficients as terms, and that the user will remember to fix one of them
    poly = 0.0;
    fptype expansion_parameter(SQRT(rMassSq) / resmass);
    for(unsigned int poly_index = 0; poly_index < num_poly_coeffs; ++poly_index)
      poly += cudaArray[indices[10 + poly_index]] * pow(expansion_parameter, int(poly_index));
  }
  else if(formfactor_type == ResonancePdf::EXPPOLY)
  {
    // form factor f(x) = exp(g(x)) so f(x) > 0
    // g(x) expressed in Chebyshev polynomials
    // with x = SQRT(rMassSq) scaled to [-1,+1]
    fptype
      rMassMin(PAIR_12 == cyclic_index ? daug1Mass + daug2Mass : (PAIR_13 == cyclic_index ? daug1Mass + daug3Mass : daug2Mass + daug3Mass)),
      rMassMax(PAIR_12 == cyclic_index ? motherMass - daug3Mass: (PAIR_13 == cyclic_index ? motherMass - daug2Mass: motherMass - daug1Mass)),
      x(-1.0 + 2.0*(SQRT(rMassSq) - rMassMin)/(rMassMax - rMassMin)),
      norm(2.0);
    for(unsigned int poly_index = 1; poly_index <= num_poly_coeffs; ++poly_index)
    {
      fptype coeff(cudaArray[indices[9 + poly_index]]);
      if(poly_index == 1)
        poly += coeff * x;
      else if(poly_index == 2)
      {
        poly += coeff * (2.0*x*x - 1.0);
        norm -= coeff * 2.0 / 3.0;
      }
      else if(poly_index == 3)
      {
        poly += coeff * (4.0*x*x*x - 3.0*x);
      }
      else
        printf("Too high an order requested from PolynomialLASS\n");  
    }
  }
  else
  {
    printf("Unknown form factor type requested\n");
  }
      
  // poly += pow(expansion_parameter, int(poly_index)) * cudaArray[indices[8 + poly_index]];
  
  ret *= poly;
  return ret;
}

// this is the general case
__device__ devcomplex<fptype> flatte_rhohelper(fptype ma, fptype mb, fptype rMassSq)
{
  fptype tmp(1.0 - (ma - mb)*(ma - mb)/rMassSq);
  tmp *= (1.0 - (ma + mb)*(ma + mb)/rMassSq);
  if(tmp > 0.0)
    return devcomplex<fptype>(SQRT(tmp), 0);
  return devcomplex<fptype>(0, SQRT(-tmp));
}

// this is the specialised case when both masses are the same
__device__ devcomplex<fptype> flatte_rhohelper(fptype m, fptype rMassSq)
{
  fptype tmp(1.0 - (4.0*m*m/rMassSq));
  if(tmp > 0.0)
    return devcomplex<fptype>(SQRT(tmp), 0);
  return devcomplex<fptype>(0, SQRT(-tmp));
}

// coupled-channel lineshape for the a(0)(980) and f(0)(980) resonances
// quite which form it takes depends on which of the above we're dealing with, and also the charge
// a(0)(980) neutral: eta pi0, K+ K-, K0 K0
// a(0)(980) charged: eta pi+, K0 K+,
// f(0)(980) neutral: pi+ pi-, pi0 pi0, K+ K-, K0 K0
// f(0)(980) charged: pi+ pi0, K0 K+
//
// a(0)(980) has parameters g_EtaPi and g_KK
// f(0)(980) has parameters g_PiPi and g_KK (sometimes called g_Pi and g_K)
//
// In some limit of strong coupling, only the ratio of these is accessible, so deal with
// a(0)(980) g_EtaPi and g_KK/g_EtaPi
// f(0)(980) g_PiPi and g_KK/g_PiPi
// rather thean the two parameters separately
//
// Note the two g_KK are not the same for the two mesons, but these parameters are assumed to be
// the same for the neutral and charged versions
//
// These also seems to be disagreement over the convention for these couplings
// Either they can occur squared, or they can be multiplied by the nominal resonance mass
// to also produce a mass-dimension-2 quantity
//
// It's not clear which is better!
//
// We don't bother with any of the spin stuff because these two mesons are both scalars.
//
// Input parameters we need:
//   meson mass
//   cyclic_index
//   g_1
//   g_KK/g_1 -- interpreted as this squared depending on the other flags
//   whether we're squaring g_1 and g_KK or multiplying by resmass
//   whether we're talking about the a(0)(980) or the f(0)(980)
//   whether it's neutral or charged
//
// In the literature I've looked at, the first two flags aren't both neccessary, because we square the
// g_i in the a(0)(980) case but not the f(0)(980) case, but I don't understand why this must be the case
__device__ devcomplex<fptype> flatte(fptype m12, fptype m13, fptype m23, unsigned int* indices)
{
  fptype resmass                = cudaArray[indices[2]];
  unsigned int cyclic_index     = indices[3];
  unsigned int square_couplings = indices[4];
  unsigned int a_meson          = indices[5];
  unsigned int charged_meson    = indices[6];
  fptype g_1                    = cudaArray[indices[7]];
  fptype g_KK_over_g_1          = cudaArray[indices[8]];
  
  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  const fptype mK0(0.497614), mK(0.493677), mPi(0.13957018), mPi0(0.1349766), mEta(0.547862);
  
  devcomplex<fptype> rho_1, rho_KK;
  
  if(a_meson)
  {
    if(charged_meson)
    {
      // charged a(0)(980)
      rho_1 = flatte_rhohelper(mEta, mPi, rMassSq); // eta pi+
      rho_KK = flatte_rhohelper(mK, mK0, rMassSq); // K0 K+
    }
    else
    {
      // neutral a(0)(980)
      rho_1 = flatte_rhohelper(mEta, mPi0, rMassSq); // eta pi0
      rho_KK = 0.5*(flatte_rhohelper(mK, rMassSq) + flatte_rhohelper(mK0, rMassSq)); // K+ K- and K0 K0
    }
  }
  else
  {
    if(charged_meson)
    {
      // charged f(0)(980)
      rho_1 = flatte_rhohelper(mPi, mPi0, rMassSq); // pi0 pi+
      rho_KK = flatte_rhohelper(mK, mK0, rMassSq); // K0 K+
    }
    else
    {
      // neutral f(0)(980)
      rho_1 = (1.0/3.0)*(2.0 * flatte_rhohelper(mPi, rMassSq) + flatte_rhohelper(mPi0, rMassSq)); // factor of 2 from isospin conservation
      rho_KK = 0.5*(flatte_rhohelper(mK, rMassSq) + flatte_rhohelper(mK0, rMassSq)); // K+ K- and K0 K0;
    }
  }
  
  // make g_i mass dimension 2 in the appropriate way
  if(square_couplings)
  {
    g_1 *= g_1;
    // in this case we interpret g_KK_over_g_1 as already being squared
    //g_KK_over_g_1 *= g_KK_over_g_1;
  }
  else
  {
    g_1 *= resmass;
    g_KK_over_g_1 *= resmass;
  }
  
  devcomplex<fptype> wid(rho_KK.real, rho_KK.imag); // not exactly a width, this has dimensions of mass^2
  wid *= g_KK_over_g_1;
  wid += rho_1;
  wid *= g_1;
  devcomplex<fptype> ret(resmass*resmass - rMassSq + wid.imag, wid.real);
  ret /= ret.abs2();
  
  return ret;
}

__device__ devcomplex<fptype> gaussian (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  // indices[1] is unused constant index, for consistency with other function types. 
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int cyclic_index     = indices[4]; 

  // Notice sqrt - this function uses mass, not mass-squared like the other resonance types. 
  fptype massToUse = SQRT(PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  massToUse -= resmass;
  massToUse /= reswidth;
  massToUse *= massToUse;
  fptype ret = EXP(-0.5*massToUse); 

  // Ignore factor 1/sqrt(2pi). 
  ret /= reswidth;

  return devcomplex<fptype>(ret, 0); 
}

// Gou-Sak h()
__device__ fptype hFun (double s, double daug2Mass, double daug3Mass) {
  // Last helper function
  const fptype _pi = 3.14159265359;
  double sm   = daug2Mass + daug3Mass;
  double SQRTs = sqrt(s);
  double k_s = twoBodyCMmom(s, daug2Mass, daug3Mass);

  double val = ((2/_pi) * (k_s/SQRTs) * log( (SQRTs + 2*k_s)/(sm)));

  return val;
}

// Gou-Sak h'()
// This is strictly correct (the derivative of h(s) w.r.t 's') only when the daughter masses are equal
__device__ fptype dh_dsFun (double s, double daug2Mass, double daug3Mass) {
  // Yet another helper function
  const fptype _pi = 3.14159265359;
  double k_s = twoBodyCMmom(s, daug2Mass, daug3Mass);
  
  double val = (hFun(s, daug2Mass, daug3Mass) * (1.0/(8.0*pow(k_s, 2)) - 1.0/(2.0 * s)) + 1.0/(2.0* _pi*s));
  return val;
}


__device__ fptype dFun (double s, double daug2Mass, double daug3Mass) {
  // Helper function used in Gronau-Sakurai
  const fptype _pi = 3.14159265359;
  double sm   = daug2Mass + daug3Mass;
  double sm24 = sm*sm/4.0; // average daughter mass squared
  double m    = sqrt(s); // nominal resonance mass
  double k_m2 = twoBodyCMmom(s, daug2Mass, daug3Mass); // this is the nominal momentum
  double val = 3.0/_pi * sm24/pow(k_m2, 2) * log((m + 2*k_m2)/sm) + m/(2*_pi*k_m2) - sm24*m/(_pi * pow(k_m2, 3));
  return val;
}

__device__ fptype fsFun (double s, double m2, double gam, double daug2Mass, double daug3Mass) {
  // Another G-S helper function
   
  double k_s   = twoBodyCMmom(s,  daug2Mass, daug3Mass);
  double k_Am2 = twoBodyCMmom(m2, daug2Mass, daug3Mass);
   
  double f     = gam * m2 / POW(k_Am2, 3);
  f           *= (POW(k_s, 2) * (hFun(s, daug2Mass, daug3Mass) - hFun(m2, daug2Mass, daug3Mass)) + (m2 - s) * pow(k_Am2, 2) * dh_dsFun(m2, daug2Mass, daug3Mass));
  // dh_dsFun is the h' function of the original paper

  return f;
}

__device__ devcomplex<fptype> gouSak (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  fptype motherMass             = functorConstants[indices[1]+0];
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];
  fptype meson_radius           = functorConstants[indices[1]+4];
  fptype mother_meson_radius    = functorConstants[indices[1]+5];

  fptype resmass                = cudaArray[indices[2]];
  fptype resmassSq              = resmass*resmass;
  fptype reswidth               = cudaArray[indices[3]];
  const unsigned int &spin                = indices[4];
  const unsigned int &cyclic_index        = indices[5]; 
  const unsigned int &use_nominal_mass    = indices[6];
  const unsigned int &use_unnorm_dampingf = indices[7];

  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype daugAMass = (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass); // These are the two daughters
  fptype daugBMass = (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass); // of the resonance
  fptype frFactor = 1;

  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2). 
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, daugAMass, daugBMass);
  fptype nominalDaughterMoms = twoBodyCMmom(resmassSq, daugAMass, daugBMass);

  if (0 != spin)
  {
    frFactor = dampingFactorRatioSquare(
        nominalDaughterMoms,
        measureDaughterMoms,
        spin, meson_radius,
        use_unnorm_dampingf);
  }

  fptype runningwidth = reswidth * frFactor * (resmass/SQRT(rMassSq)) * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin+1.0);
  
  // Implement Gou-Sak:
  fptype D = (1.0 + dFun(resmassSq, daugAMass, daugBMass)*reswidth/resmass);
  fptype E = resmassSq - rMassSq + fsFun(rMassSq, resmassSq, reswidth, daugAMass, daugBMass);
  fptype F = resmass * runningwidth;
  //SQRT(resmass) * reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor;

  D /= (E*E + F*F);
  devcomplex<fptype> retur(D*E, D*F);

  // Didn't want the F_D penetration factor in the mass-dependent width, but we do want it in the overall spin factor
  if(0 != spin)
  {
    fptype bachelorMass(PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
    frFactor *= dampingFactorRatioSquare(
        bachelorMom(resmass, motherMass, bachelorMass),
        bachelorMom(SQRT(rMassSq), motherMass, bachelorMass),
        spin, mother_meson_radius,
        use_unnorm_dampingf);
  }

  retur *= SQRT(frFactor);
  retur *= spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index, resmassSq, use_nominal_mass);
  return retur; 
}

__device__ devcomplex<fptype> nonres (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  return devcomplex<fptype>(1, 0); 
}


__device__ void getAmplitudeCoefficients (devcomplex<fptype> a1, devcomplex<fptype> a2, fptype& a1sq, fptype& a2sq, fptype& a1a2real, fptype& a1a2imag) {
  // Returns A_1^2, A_2^2, real and imaginary parts of A_1A_2^*
  a1sq = a1.abs2();
  a2sq = a2.abs2();
  a1 *= conj(a2);
  a1a2real = a1.real;
  a1a2imag = a1.imag; 
}

__device__ resonance_function_ptr ptr_to_RBW = plainBW;
__device__ resonance_function_ptr ptr_to_FLATTE = flatte;
__device__ resonance_function_ptr ptr_to_LASS = lass;
__device__ resonance_function_ptr ptr_to_polynomialLASS = polylass;
__device__ resonance_function_ptr ptr_to_GOUSAK = gouSak; 
__device__ resonance_function_ptr ptr_to_GAUSSIAN = gaussian;
__device__ resonance_function_ptr ptr_to_NONRES = nonres;


ResonancePdf::ResonancePdf (string name,
    const AmplitudeInfo &amp_,
    Variable* mass, 
    Variable* width, 
    unsigned int sp, 
    unsigned int cyc,
    bool useNominalMass,
    bool useUnNormalisedDampingFactors) 
  : GooPdf(0, name), amp(amp_)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Making room for index of decay-related constants. Assumption:
  // These are mother mass and three daughter masses in that order.
  // They will be registered by the object that uses this resonance,
  // which will tell this object where to find them by calling setConstantIndex. 

  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width)); 
  pindices.push_back(sp);
  pindices.push_back(cyc); 
  pindices.push_back(useNominalMass);
  pindices.push_back(useUnNormalisedDampingFactors);

  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_RBW, sizeof(void*));
  initialise(pindices); 
}

ResonancePdf::ResonancePdf (string name,
    const AmplitudeInfo &amp_,
    Variable* mass,
    Variable* g_1,
    Variable* g_KK_over_g_1,
    unsigned int cyc,
    CouplingTreatment square_couplings,
    WhichMeson a_meson,
    MesonCharge charged_meson)
: GooPdf(0, name), amp(amp_)
{
  vector<unsigned int> pindices;
  pindices.push_back(0);
  // Making room for index of decay-related constants. Assumption:
  // These are mother mass and three daughter masses in that order.
  // They will be registered by the object that uses this resonance,
  // which will tell this object where to find them by calling setConstantIndex.

  pindices.push_back(registerParameter(mass));
  pindices.push_back(cyc);
  pindices.push_back(square_couplings);
  pindices.push_back(a_meson);
  pindices.push_back(charged_meson);
  pindices.push_back(registerParameter(g_1));
  pindices.push_back(registerParameter(g_KK_over_g_1));
  
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_FLATTE, sizeof(void*));
  initialise(pindices); 
}

ResonancePdf::ResonancePdf (string name,
    const AmplitudeInfo &amp_,
    Variable* mass,
    Variable* width,
    Variable* lass_a,
    Variable* lass_r,
    Variable* lass_phi_f,
    Variable* lass_phi_r,
    Variable* lass_F,
    unsigned int sp,
    unsigned int cyc)
: GooPdf(0, name), amp(amp_)
{
  vector<unsigned int> pindices;
  pindices.push_back(0);
  // Making room for index of decay-related constants. Assumption:
  // These are mother mass and three daughter masses in that order.
  // They will be registered by the object that uses this resonance,
  // which will tell this object where to find them by calling setConstantIndex.
  
  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width));
  pindices.push_back(sp);
  pindices.push_back(cyc);
  pindices.push_back(0);// useNominalMass -- this lets us defer to plainBW()
  pindices.push_back(0);// useUnNormalisedDampingFactors -- this lets us defer to plainBW()
  pindices.push_back(registerParameter(lass_a));
  pindices.push_back(registerParameter(lass_r));
  pindices.push_back(registerParameter(lass_phi_f));
  pindices.push_back(registerParameter(lass_phi_r));
  pindices.push_back(registerParameter(lass_F));

  
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_LASS, sizeof(void*));
  initialise(pindices);
  
  if(sp)
    std::cout << "WARNING from LASS shape constructor: spin " << sp << " requested, at present this is ignored" << std::endl;
}

ResonancePdf::ResonancePdf(std::string name,
    const AmplitudeInfo &amp_,
    Variable* mass,
    Variable* width,
    Variable* lass_a,
    Variable* lass_r,
    const std::vector<Variable*> &poly_coeffs,
    unsigned int sp,
    unsigned int cyc,
    FormFactorType fftype)
: GooPdf(0, name), amp(amp_)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(0); // copying the rest of the constructors...
  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width));
  pindices.push_back(sp);
  pindices.push_back(cyc);
  pindices.push_back(registerParameter(lass_a));
  pindices.push_back(registerParameter(lass_r)); // might as well match the other LASS shape up to this point
  pindices.push_back(poly_coeffs.size());
  pindices.push_back(fftype);
  for(std::vector<Variable*>::const_iterator poly_iter = poly_coeffs.begin(); poly_iter != poly_coeffs.end(); poly_iter++)
    pindices.push_back(registerParameter(*poly_iter));
  
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_polynomialLASS, sizeof(void*));
  initialise(pindices);
  
  if(sp)
    std::cout << "WARNING from Polynomial LASS shape constructor: spin " << sp << " requested, at present this is ignored" << std::endl;
}

ResonancePdf::ResonancePdf (string name,
    const AmplitudeInfo &amp_,
    unsigned int sp, 
    Variable* mass, 
    Variable* width, 
    unsigned int cyc,
    bool useNominalMass,
    bool useUnNormalisedDampingFactors) 
: GooPdf(0, name), amp(amp_)
{
  // Same as BW except for function pointed to. 
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width)); 
  pindices.push_back(sp);
  pindices.push_back(cyc);
  pindices.push_back(useNominalMass);
  pindices.push_back(useUnNormalisedDampingFactors);

  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_GOUSAK, sizeof(void*));
  initialise(pindices); 
}

ResonancePdf::ResonancePdf (string name, const AmplitudeInfo &amp_) 
  : GooPdf(0, name), amp(amp_)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Dummy index for constants - won't use it, but calling 
  // functions can't know that and will call setConstantIndex anyway. 
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_NONRES, sizeof(void*));
  initialise(pindices); 
}

ResonancePdf::ResonancePdf (string name,
    const AmplitudeInfo &amp_,
    Variable* mean, 
    Variable* sigma,
    unsigned int cyc) 
  : GooPdf(0, name), amp(amp_)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Dummy index for constants - won't use it, but calling 
  // functions can't know that and will call setConstantIndex anyway. 
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma)); 
  pindices.push_back(cyc); 

  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_GAUSSIAN, sizeof(void*));
  initialise(pindices); 	
}
