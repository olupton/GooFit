#ifndef GAUSSIAN_NO_OBS_THRUST_FUNCTOR_HH
#define GAUSSIAN_NO_OBS_THRUST_FUNCTOR_HH

#include "GooPdf.hh" 

class GaussianNoObsPdf : public GooPdf {
public:
  GaussianNoObsPdf (std::string n, Variable* _x,Variable* m, Variable* s);
  __host__ fptype integrate (fptype lo, fptype hi) const;
  __host__ virtual bool hasAnalyticIntegral () const {return true;}
private:
};

#endif
