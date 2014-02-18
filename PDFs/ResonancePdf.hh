#ifndef RESONANCE_PDF_HH
#define RESONANCE_PDF_HH

#include "GooPdf.hh" 
#include "devcomplex.hh" 
typedef devcomplex<fptype> (*resonance_function_ptr) (fptype, fptype, fptype, unsigned int*); 

class ResonancePdf : public GooPdf {
  friend class TddpPdf;
  friend class DalitzPlotPdf;
  friend class IncoherentSumPdf; 
public:
  // Constructor for regular BW 
  ResonancePdf (string name, 
			  Variable* ar, 
			  Variable* ai, 
			  Variable* mass, 
			  Variable* width, 
			  unsigned int sp, 
			  unsigned int cyc);
  
  enum CouplingTreatment { MULTIPLY_BY_NOMINAL_MASS = 0, SQUARE };
  enum WhichMeson { F = 0, A };
  enum MesonCharge { CHARGED, NEUTRAL };
  
  // Constructor for Flatte lineshape
  ResonancePdf (string name,
                          Variable* ar,
                          Variable* ai,
                          Variable* mass,
                          Variable* g_1,
                          Variable* g_KK_over_g_1,
                          unsigned int cyc,
                          CouplingTreatment square_couplings,
                          WhichMeson a_meson,
                          MesonCharge charged_meson);
  // Old Flatte
  //ResonancePdf (string name,
  //      Variable* ar,
  //      Variable* ai,
  //      Variable* mass,
  //      Variable* g_pi,
  //      Variable* g_k,
  //      unsigned int sp,
  //      unsigned int cyc);

  // Constructor for LASS lineshape
  ResonancePdf (string name,
        Variable* ar,
        Variable* ai,
        Variable* mass,
        Variable* width,
        Variable* lass_a,
        Variable* lass_r,
        Variable* lass_phi_f,
        Variable* lass_phi_r,
        Variable* lass_F,
        unsigned int sp,
        unsigned int cyc);
  
  // Constructor for Brian's LASS lineshape
  ResonancePdf (string name,
                Variable* ar,
                Variable* ai,
                Variable* mass,
                Variable* width,
                Variable* lass_a,
                Variable* lass_r,
                const std::vector<Variable*> &poly_coeffs,
                unsigned int sp,
                unsigned int cyc);
  
  // Gounaris-Sakurai
  ResonancePdf (string name,
			  Variable* ar,
			  Variable* ai,
			  unsigned int sp,
			  Variable* mass,
			  Variable* width,
			  unsigned int cyc);

  // Nonresonant constructor
  ResonancePdf (string name, 
			  Variable* ar, 
			  Variable* ai);  

  // Gaussian constructor
  ResonancePdf (string name,
			  Variable* ar, 
			  Variable* ai,
			  Variable* mean, 
			  Variable* sigma,
			  unsigned int cyc);

private:
  void setConstantIndex (unsigned int idx) {host_indices[parameters + 1] = idx;}

  Variable* amp_real;
  Variable* amp_imag;
  /*
  Variable* mass;
  Variable* width;
  unsigned int spin;
  unsigned int cyclic_index;
  unsigned int eval_type;
  unsigned int resonance_type; 
  */ 
};

#endif
