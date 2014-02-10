#include "StaircasePdf.hh"

__device__ fptype device_Staircase (fptype* evt, fptype* p, unsigned int* indices)
{
  fptype x = evt[indices[2 + indices[0]]];
  unsigned int x_int = FLOOR(0.5 + x);
  unsigned int ret = indices[1];
  // indices[1] should be the number of step points we have
  for(unsigned int i = 0; i < indices[1]; i++)
  {
    if(x_int < p[indices[2 + i]])
    {
      ret = i;
      break;
    }
  }
  
  //printf("Staircase returning %u based on %u\n", ret, x_int);
  return ret;
}

__device__ device_function_ptr ptr_to_Staircase = device_Staircase;
device_function_ptr hptr_to_Staircase = device_Staircase;

__host__ StaircasePdf::StaircasePdf(std::string n, Variable* _x, const std::vector<Variable*> &x0list)
  : GooPdf(_x, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(x0list.size());
  for(std::vector<Variable*>::const_iterator x0 = x0list.begin(); x0 != x0list.end(); x0++)
    pindices.push_back(registerParameter(*x0));
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_Staircase, sizeof(void*));
  initialise(pindices);
  std::cout << "StaircasePdf::StaircasePdf(" << n << ", ...)" << std::endl;
}

//__host__ fptype StaircasePdf::integrate (fptype lo, fptype hi) const {
//  unsigned int* indices = host_indices+parameters;
  
  // the part where the function is zero is 
  
  //fptype x0 = host_params[indices[1]];
  //return (hi - x0);
//}

