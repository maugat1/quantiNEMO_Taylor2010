/** @file random.h
*
*   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
*
*   quantiNEMO:
*   quantiNEMO is an individual-based, genetically explicit stochastic
*   simulation program. It was developed to investigate the effects of
*   selection, mutation, recombination, and drift on quantitative traits
*   with varying architectures in structured populations connected by
*   migration and located in a heterogeneous habitat.
*
*   quantiNEMO is built on the evolutionary and population genetics
*   programming framework NEMO (Guillaume and Rougemont, 2006, Bioinformatics).
*
*
*   Licensing:
*   This file is part of quantiNEMO.
*
*   quantiNEMO is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   quantiNEMO is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with quantiNEMO.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef uniformH
#define uniformH

#include "mtrand.h"
#include "output.h"

/**Random number generation class, uses various types of random generation depending on the implementation.*/
class RAND: public MTRand {
private:
	typedef unsigned long uint32;  // unsigned integer type, at least 32 bits

public:

	RAND(){}          // do nothing
	RAND(unsigned long seed) : MTRand(seed){}
	RAND(const unsigned long* seed, int size) : MTRand(seed, size){}

	void set_seed(unsigned long s){seed(s);}
	void set_seed(const unsigned long* s, int size){seed(s, size);}
	void set_seed_default(){seed_default();}
	void set_seed_time(){
		unsigned long seeds[] = {time(NULL), clock()};
		seed(seeds, 2);
	}

//	RAND( uint32 *const bigSeed, uint32 const seedLength = N ){MTRand::MTRand(bigSeed,seedLength);}  // or an array
//	RAND(){MTRand::MTRand();}  // auto-initialize with /dev/urandom or time() and clock()
//	RAND( const uint32& s ){MTRand::MTRand(s);}  // initialize with a simple uint32

	/** Returns a random number from [0.0, 1.0[ uniformly distributed. **/
	double Uniform () {return next();}

	/** Returns a uniformly distributed entire random number from [0, max[.*/
	unsigned int Uniform (unsigned int max) {return (unsigned int)(max*next());}

	/**Returns a random boolean number.*/
	bool Bool() {return (next() < 0.5);}


  // -----------------------------------------------------------------------------
  /** Returns a deviate distributed as a gamma distribution of integer order ia, i.e., a waiting time
    * to the iath event in a Poisson process of unit mean, using ran1(idum) as the source of
    * uniform deviates. Return as a double floating-point number an integer value that is a random
    * deviate drawn from a binomial distribution of n trials each of
    * probability p, using Randoms () as a source of uniform random deviates.
    * Adapted from the algorithm described in the book Numerical Recipes by Press et al.
    */
  double Gamma(int ia){
    int j;
    double am,e,s,v1,v2,x,y;
    if (ia < 1) fatal("The shape in the function Gamma() must be postiive!\n");
    if (ia < 6) {               // Use direct method, adding waiting times.
      x=1.0; 
      for (j=1;j<=ia;j++) x *= Uniform();
      x = -log(x);
    }
    else {  // Use rejection method.
    do {
      do {
        do {                                      // These four lines generate the 
          v1=Uniform();                           // tangent of a random angle, i.e.,
          v2=2.0*Uniform()-1.0;                   // they are equivalent to y = tan(π * ran1(idum)).
        } while (v1*v1+v2*v2 > 1.0);
        y=v2/v1;
        am=ia-1;
        s=sqrt(2.0*am+1.0);
        x=s*y+am;                                 // We decide whether to reject x:
      } while (x <= 0.0);                         // Reject in region of zero probability.
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);          // Ratio of prob. fn. to comparison fn.
    } while (Uniform() > e);                      // Reject on basis of a second uniform deviate.
  }
  return x;
}

  // -----------------------------------------------------------------------------
  /**From the Numerical Recieps.
   A function to return the natural log of the Gamma function of x.
   Adapted from the algorithm described in the book Numerical Recipes by
   Press et al.
   It can be used to compute factorials since ln(n!) = lnGamma(n + 1)
  */
	double gammln (double xx) {
    double x,y,tmp,ser=1.000000000190015;
    static double cof[6]={76.18009172947146,-86.50532032941677,
      24.01409824083091,-1.231739572450155,
      0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
		tmp=x+5.5;
		tmp -= (x+0.5)*log(tmp);
		for (j = 0; j < 6; ++j) ser += cof[j]/++y;

    return -tmp+log(2.5066282746310005*ser/x);
  }

	// -----------------------------------------------------------------------------
	/**From the Numerical Recieps.*/
	unsigned int Poisson (double mean) {
    static double sq,alxm,g,oldm=(-1.0);
    double em,t,y;

    if (mean < 12.0){
      if (mean != oldm){
        oldm=mean;
        g=exp(-mean);
      }
      em = -1;
      t=1.0;
      do {
        ++em;
        t *= Uniform();
      } while (t > g);
    }
    else{
      if (mean != oldm)
        {
        oldm=mean;
        sq=sqrt(2.0*mean);
        alxm=log(mean);
        g=mean*alxm-gammln(mean+1.0);
        }
      do {
        do {
          y=tan(M_PI*Uniform());
          em=sq*y+mean;
        } while (em < 0.0);
        em=floor(em);
        t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (Uniform() > t);
    }
    return (unsigned int)em;
  }

// -----------------------------------------------------------------------------
/*
   BinomialRandom // added by Samuel 27-10-2006

   Return as a double floating-point number an integer value that is a random
   deviate drawn from a binomial distribution of n trials each of
   probability p, using Randoms () as a source of uniform random deviates.
   Adapted from the algorithm described in the book Numerical Recipes by
   Press et al.
*/

	double Binomial (double p, long N)
  {
    long j;
    static long iOldN = -1;
    double dAngle, dDeviate, dMean, dPtemp, dSqrt, dTangent, dTemp1, dTemp2;
    static double dLnFactN, dPold = -1, dLnP, dQ, dLnQ;

    if (p < 0 || p > 1 || N < 0) {
      fatal("Negative or null parameters for a binomial variate - Exiting\n\n");
    }

    dPtemp = ( p <= 0.5 ? p : 1 - p);
    dMean = N * dPtemp;  /* mean of the deviate to be produced. */

    /* Use the direct method if N is not too large.
       This can require up to 25 calls to random */

    if (N < 25) {
      dDeviate = 0;
      for (j = 0; j < N; j++)
        if (Uniform () < dPtemp)
          dDeviate = dDeviate + 1;
    }
    else {
      if (dMean < 1) {
        /* if less than one event is expected out of 25 or more trials,then the
           distribution is quite accurately Poisson. Use direct method. */
        dTemp1 = exp(-dMean);
        dTemp2 = 1.0;
        for (j = 0; j <= N; j++) {
          dTemp2 = dTemp2 * Uniform ();
          if (dTemp2 < dTemp1) break;
        }

        dDeviate = (j <= N ? j : N);
      }
      else { /* Use rejection */

        if (N != iOldN) {
          /* if N has changed or it's the first call, initialize */
          dLnFactN = gammln((double) N + 1);
          iOldN = N;
        }

        if (dPtemp != dPold) {
          /* if dPtemp has changed or it's the first call, initialize. */
          dPold = dPtemp;
          dQ = 1 - dPtemp;
          dLnP = log(dPtemp);
          dLnQ = log(dQ);
        } /* if */

        dSqrt = sqrt(2 * dMean * dQ);

        /* Rejection method with a Lorentzian comparison function. */

        do {
          do {
            dAngle = M_PI * Uniform ();
            dTangent = tan(dAngle);
            dTemp1 = dSqrt * dTangent + dMean;
          } while (dTemp1 < 0 || dTemp1 >= (N + 1)); /* Reject */

          dTemp1 = floor(dTemp1); /* discrete distribution */

          dTemp2 = 1.2 * dSqrt * (1 + dTangent * dTangent) *
                   exp(dLnFactN - gammln(dTemp1 + 1) -
                       gammln(N - dTemp1 + 1) +
                       dTemp1 * dLnP + (N - dTemp1) * dLnQ);

        } while (Uniform () > dTemp2);

        /* Reject on average about 1.5 time per deviate */

        dDeviate = dTemp1;

      } /* else */  /* end of rejection */
    }
    
    if (dPtemp != p)
      dDeviate = N - dDeviate; /* undo the symmetry tranformation */

    return (dDeviate);

  } /* BinomialRandom */


  //----------------------------------------------------------------------------
  /* return a random position of the array following the distribution
   * by default a non-cumulative distribution is assumed!!!
   * if a non-cumulative array is passed a temp array is used to get a
   * cumulative distribution -> original array is not modified!
   * following fl (Fred Hospital)
   * array:    array with the values
   * size:     number of elements of the array
   * isCumul: -1: array is not yet set reverse cumulative             -> random smallest value
   *           0: array is not yet set cumulative (default)           -> random biggest value
   *           1: array is already cumulative or reverse cumulative
	 */
  template <typename T>
	int AfterDistribution(T* array, int size, int isCumul){
    T* a;

		if(size<=1) return 0;

    switch(isCumul){
			case -1: // make array reverse cumulative  (1/n) (attention if fitenss is 0!)
						a = new T[size];
						a[0] = array[0] ? 1.0/array[0] : 1e9;        // first element
            for(int i = 1; i<size; ++i){
							a[i] = array[i] ? (a[i-1] + 1.0/array[i]) : (a[i-1] + 1e9);
						}
            break;
      case  0: // make array cumulative
            a = new T[size];
            a[0] = array[0];            // first element;
            for(int i=1; i<size; ++i){  // start at the second position
              a[i] = a[i-1] + array[i];
            }
            break;
      case  1: // alredy cumulative or reverse cumulative array: make nothing
						return AfterDistribution(array, size);
		}

		// only passsed here if the cuulative array had to be created
		int pos = AfterDistribution(a, size);
		delete[] a;      // delete the temp array if it was created
		return pos;
	}

	//----------------------------------------------------------------------------
	/** array must be cumulative */
	template <typename T>
	int AfterDistribution(T* array, int size){
		if(size<=1) return 0;

		double value, test;
		int start, end, middle;

		// get randomly a position based on the array values
		if (array[size-1] == 0)  return Uniform(size);    // if all values in the array have value 0
		value = array[size - 1]*Uniform(); // get a random uniform number
		if (array[0] >= value)  return 0;
		start = 0;
		end = size - 1;
		while (end - start > 1) {
			middle = (start + end) / 2;
			test = array[middle];
			if (value < test) 		  end = middle;
			else if (value > test)	start = middle;
			else					          return middle;
		}
		return end;
	}

	//----------------------------------------------------------------------------
	/* Returns a Normal random variate based on a unit variate,
		using a random generator as a source of uniform deviates.
		Adapted from the algorithm described in the book Numerical Recipes by
		Press et al.
	*/
	double Normal (const double& dMean = 0, const double& dStdDev = 1){
		double w, x1, x2;

		do {
			x1 = 2. * next() - 1.;
			x2 = 2. * next() - 1.;
			w = x1*x1 + x2*x2;
		} while (w >= 1. || w < 1E-30);

		w = sqrt((-2.*log(w))/w);
		x1 *= w;
		return (x1 * dStdDev + dMean);
	} /* NormalRandom */

	//----------------------------------------------------------------------------
	/* returns a variate such that the log of the variate is normally distributed.*/
	double LogNormal (const double& dMean, const double& dStdDev){
		if(dMean<=0) fatal("LogNormalRandom: mean must be positive!\n");
		if(dStdDev<=0) fatal("LogNormalRandom: sd must be postive!\n");

		return exp (Normal(log (dMean), log (dStdDev)));
	}

	/* -----------------------------------------------------------------------------
    cumulative Normal distribution

		Returns the integral until z of a Normal distribution N(0, 1)
		Taken from http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
  */
	double IntegralNormal (double x){
		const double b1 =  0.319381530;
    const double b2 = -0.356563782;
    const double b3 =  1.781477937;
    const double b4 = -1.821255978;
		const double b5 =  1.330274429;
    const double p  =  0.2316419;
    const double c  =  0.39894228;

    if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
			return (1.0 - c * exp( -x * x / 2.0 ) * t *
			( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
		}
		else {
			double t = 1.0 / ( 1.0 - p * x );
			return ( c * exp( -x * x / 2.0 ) * t *
			( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
		}
	}

	/* -----------------------------------------------------------------------------
		cumulative Normal distribution (High Precision Approximation)
		A high (double) precision approximation is described in the article Better
		Approximations to Cumulative Normal Fuctiona by Graeme West. The algorithm
		is based on Hart's algorithm 5666 (1968).

		Returns the integral until x of a Normal distribution N(0, 1)
		Taken from http://www4.ncsu.edu/unity/users/p/pfackler/www/ECG790C/accuratecumnorm.pdf
	*/
	double Cumnorm(double x){
		const double a = 7.07106781186547;
		const double b = 3.52624965998911E-02;
		const double c = 0.700383064443688;
		const double d = 6.37396220353165;
		const double e = 33.912866078383;
		const double f = 112.079291497871;
		const double g = 221.213596169931;
		const double h = 220.206867912376;
		const double i = 8.83883476483184E-02;
		const double j = 1.75566716318264;
		const double k = 16.064177579207;
		const double l = 86.7807322029461;
		const double m = 296.564248779674;
		const double n = 637.333633378831;
		const double o = 793.826512519948;
		const double p = 440.413735824752;
		const double q = 2.506628274631;

		double XAbs = x<0 ? -x : x;     // make it absolute
		double Cumnorm;

		if(XAbs > 37){
			Cumnorm = 0;
		}
		else{
			if(XAbs < a){
				Cumnorm = ((b * XAbs + c) * (XAbs + d) * (XAbs + e) * (XAbs + f)
								* (XAbs + g) * (XAbs + h) * exp(-XAbs*XAbs/2.))
								/ ((i * XAbs + j) * (XAbs + k) * (XAbs + l) * (XAbs + m)
								* (XAbs + n) * (XAbs + o) * (XAbs + p));
			}
			else{
				Cumnorm = exp(-XAbs*XAbs/2.)
								/ (XAbs + 1. / (XAbs + 2. / (XAbs + 3. / (XAbs + 4. / (XAbs + 0.65)))))
								/ q;
			}
		}
		if(x > 0) return (1-Cumnorm);
		return Cumnorm;
	}


	// -----------------------------------------------------------------------------
	double ProbDensNorm(double x, double mu=0, double sd=1){
		const double PI = 3.1415926535897932384626433832795028841972;
		double m = (x-mu)/sd;
		return exp(-m*m/2.)/(sqrt(2*PI)*sd);
	}


};
#endif

