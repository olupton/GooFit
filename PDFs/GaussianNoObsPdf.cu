#include "GaussianNoObsPdf.hh"

__device__ fptype device_GaussianNoObs(fptype* evt, fptype* p, unsigned int* indices) {
  //fptype x = evt[indices[2 + indices[0]]];
  fptype x = p[indices[1]];
  fptype mean = p[indices[2]];
  fptype sigma = p[indices[3]];
  fptype ret = -0.5*(x-mean)*(x-mean)/(sigma*sigma);
  //printf("GaussianNoObs(%f, %f, %f) = EXP(%f ~ %f)\n", x, mean, sigma, ret, LOG(EXP(ret)));
  return EXP(ret);
}

__device__ device_function_ptr ptr_to_GaussianNoObs = device_GaussianNoObs;

__host__ GaussianNoObsPdf::GaussianNoObsPdf (std::string n, Variable* _x, Variable* mean, Variable* sigma)
  : GooPdf(NULL, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(_x));
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma));
  cudaMemcpyFromSymbol((void**) &host_fcn_ptr, ptr_to_GaussianNoObs, sizeof(void*));
  initialise(pindices); 
}

__host__ fptype GaussianNoObsPdf::integrate (fptype lo, fptype hi) const
{
  //std::cout << "GaussianNoObsPdf::integrate(" << lo << ", " << hi << ")" << std::endl;
  printf("GaussianNoObsPdf::integrate(%f, %f)\n", lo, hi);
  //static const fptype root2 = sqrt(2.);
  static const fptype rootPi = sqrt(atan2(0.0,-1.0));
  //static const fptype rootPiBy2 = rootPi / root2;
  
  unsigned int* indices = host_indices+parameters;
  //fptype xscale = root2*host_params[indices[2]];

  /*
  std::cout << "Gaussian integral: " 
	    << xscale << " "
	    << host_params[indices[1]] << " "
	    << host_params[indices[2]] << " "
	    << ERF((hi-host_params[indices[1]])/xscale) << " "
	    << ERF((lo-host_params[indices[1]])/xscale) << " "
	    << rootPiBy2*host_params[indices[2]]*(ERF((hi-host_params[indices[1]])/xscale) -
						  ERF((lo-host_params[indices[1]])/xscale)) 
	    << std::endl; 
  */
  //return rootPiBy2*host_params[indices[2]]*(ERF((hi-host_params[indices[1]])/xscale) - 
  //					    ERF((lo-host_params[indices[1]])/xscale));

  // Integral over all R. 
  fptype sigma = host_params[indices[3]];
  sigma *= root2*rootPi; // root2 is a #define
  return sigma;
}

