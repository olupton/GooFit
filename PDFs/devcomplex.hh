#ifndef DEVCOMPLEX_HH
#define DEVCOMPLEX_HH

#include "GlobalCudaDefines.hh"
#include <complex>
using std::complex;

//#define USE_POLAR_AMPLITUDES 1
#define USE_PHASES_IN_DEGREES 1
#ifdef USE_PHASES_IN_DEGREES
// this is pi/180
#define PHASE_CONVERSION_FACTOR 0.01745329251
#else
#define PHASE_CONVERSION_FACTOR 1.0
#endif

#if CUDA_VERSION > 5000
#define HOST_STRING __host__
#else
#define HOST_STRING
#endif


template <typename T> struct devcomplex {
  T real;
  T imag; 

  HOST_STRING EXEC_TARGET devcomplex<T> () : real(0), imag(0) {}
  HOST_STRING EXEC_TARGET devcomplex<T> (T r, T i) : real(r), imag(i) {}

  HOST_STRING EXEC_TARGET devcomplex<T>& operator= (T other) {
    real = other;
    imag = 0; 
    return *this; 
  }

   HOST_STRING EXEC_TARGET devcomplex<T>& operator+= (const devcomplex<T>& other) {
    real += other.real;
    imag += other.imag;
    return *this;
  }

  HOST_STRING  EXEC_TARGET devcomplex<T>& operator-= (const devcomplex<T>& other) {
    real -= other.real;
    imag -= other.imag;
    return *this;
  }

 HOST_STRING EXEC_TARGET devcomplex<T>& operator*= (T mult) {
    real *= mult;
    imag *= mult;
    return *this;
  }

 HOST_STRING EXEC_TARGET devcomplex<T>& operator/= (T mult) {
    real /= mult;
    imag /= mult;
    return *this;
  }

 HOST_STRING EXEC_TARGET devcomplex<T>& operator*= (const devcomplex<T>& other) {
    multiply(other.real, other.imag);
    return *this;
  }

 HOST_STRING EXEC_TARGET devcomplex<T>& multiply (const T other_real, const T other_imag) {
    T nreal = real * other_real - imag * other_imag;
    T nimag = real * other_imag + imag * other_real; 
    real = nreal;
    imag = nimag; 
    return *this; 
  }

 HOST_STRING EXEC_TARGET T abs2 () {return real*real + imag*imag;}

 HOST_STRING EXEC_TARGET T arg () { return atan2(imag, real); }
};

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> operator+ (const devcomplex<T>& one, const devcomplex<T>& other) {
  return devcomplex<T>(one.real+other.real, one.imag+other.imag); 
}

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> operator- (const devcomplex<T>& one, const devcomplex<T>& other) {
  return devcomplex<T>(one.real-other.real, one.imag-other.imag); 
}

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> operator+ (const devcomplex<T>& one, T other) {
  return devcomplex<T>(one.real+other, one.imag); 
}

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> operator/ (const devcomplex<T>& one, const devcomplex<T>& other) {
  T inverse(1);
  inverse /= (other.real*other.real + other.imag*other.imag);
  return devcomplex<T>(inverse*(one.real*other.real+one.imag*other.imag), inverse*(one.imag*other.real - one.real*other.imag)); 
}

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> operator* (const devcomplex<T>& one, const devcomplex<T>& other) {
  return devcomplex<T>(one.real*other.real - one.imag*other.imag, one.real*other.imag + one.imag*other.real); 
}

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> operator* (const devcomplex<T>& one, const int& other) {
  return devcomplex<T>(one.real*other, one.imag*other); 
}

template<typename T> HOST_STRING EXEC_TARGET devcomplex<T> operator/ (const T& a, const devcomplex<T>& b) {
  T inverse(1);
  inverse /= (b.real*b.real + b.imag*b.imag);
  return devcomplex<T>(inverse*a*b.real, -inverse*a*b.imag); 
}

template<typename T> HOST_STRING EXEC_TARGET devcomplex<T> operator* (T a, const devcomplex<T>& b) {
  return devcomplex<T>(a*b.real, a*b.imag); 
}

template <typename T> HOST_STRING EXEC_TARGET T norm (const devcomplex<T>& z) {
  return SQRT(z.real*z.real + z.imag*z.imag); 
}

template <typename T> HOST_STRING EXEC_TARGET T norm2 (const devcomplex<T>& z) {
  return (z.real*z.real + z.imag*z.imag); 
}

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> exp (const devcomplex<T>& z) {
  T mult = EXP(z.real); 
  return devcomplex<T>(mult*COS(z.imag), mult*SIN(z.imag)); 
}

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> conj (const devcomplex<T>& z) {
  return devcomplex<T>(z.real, -(z.imag)); 
}

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> makedevcomplex(T a, T b)
{
#ifdef USE_POLAR_AMPLITUDES
  return devcomplex<T>(a*COS(b*PHASE_CONVERSION_FACTOR), a*SIN(b*PHASE_CONVERSION_FACTOR));
#else
  return devcomplex<T>(a, b);
#endif
}

template <typename T> HOST_STRING EXEC_TARGET devcomplex<T> makedevcomplex(T a, T b, T a_del, T b_del, bool add)
{
  if(add)
    return makedevcomplex(a+a_del, b+b_del);
  else
    return makedevcomplex(a-a_del, b-b_del);
}

template <typename T> __host__ complex<T> makecomplex(T a, T b)
{
#ifdef USE_POLAR_AMPLITUDES
    return complex<T>(a*COS(b*PHASE_CONVERSION_FACTOR), a*SIN(b*PHASE_CONVERSION_FACTOR));
#else
      return complex<T>(a, b);
#endif
}

typedef devcomplex<fptype> fpcomplex;

#endif
