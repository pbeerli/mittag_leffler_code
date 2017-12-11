/*-----------------------------------------------------------------
  Bayesian inference of population genetic forces: drift, migraiton, divergence
  allowing for the n-coalescent, the f-coalescent, and the BSC-coalescent
 
  Peter Beerli
  Department of Scientific Computing
  Florida State University
  Tallahassee FL 32306-4120
  beerli@fsu.edu
 
  Copyright 2017 Peter Beerli, Tallahassee FL

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject
 to the following conditions:
 
 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*-----------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sighandler.h"
#include "romberg.h"

// ROMBERG code is from Wikipedia (downloaded 2017)
// https://en.wikipedia.org/wiki/Romberg%27s_method
// Code from wikipedia
void dump_row(size_t i, double *R);
double romberg(double (*f)(double, complex double *, long), complex double *args, long arglen,
	       double a, double b, size_t max_steps, double acc);

double romberg(double (*f)(double, complex double *, long), complex double *args, long arglen,
	       double a, double b, size_t max_steps, double acc)
{
  //double R1[max_steps], R2[max_steps]; //buffers
  //double *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
  double *R1;
  double *R2;
  R1 = (double *) mycalloc(max_steps,sizeof(double));
  R2 = (double *) mycalloc(max_steps,sizeof(double));
  double *Rp = R1;
  double *Rc = R2;
  double h = (b-a); //step size
  Rp[0] = (f(a,args, arglen) + f(b, args, arglen))*h*.5; //first trapezoidal step
   size_t i;
   size_t j;
   for(i = 1; i < max_steps; ++i){
      h /= 2.;
      double c = 0;
      size_t ep = 1 << (i-1); //2^(n-1)
      for(j = 1; j <= ep; ++j){
	c += f(a+(2*j-1)*h, args, arglen);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(j = 1; j <= i; ++j){
	double n_k = ( 1 << (2*j));
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

#ifdef ROMBERGERTESTER
      //Dump ith column of R, R[i,i] is the best estimate so far
      dump_row(i, Rc);
#endif
      if(isinf(Rc[i]))
	{
	  //return Rc[i-1];
	  double rcval = Rc[i-1];
	  myfree(Rc);
	  myfree(Rp);
	  return rcval;
	}
      if(i > 1 && fabs(Rp[i-1]-Rc[i]) < acc){
	{
	  //return Rc[i-1];
	  double rcval = Rc[i-1];
	  myfree(Rc);
	  myfree(Rp);
	  return rcval;
	}
      }
      //swap Rn and Rc as we only need the last row
      double *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
   //return Rp[max_steps-1]; //return our best guess
   double rpval = Rp[max_steps-1];
   myfree(Rc);
   myfree(Rp);
   return rpval;
}

#ifdef ROMBERGERTESTER
double f(double x, double *args, long arglen)
{
  double v = args[0] * exp(-(x*x) * args[0]);
  return v;
}


void dump_row(size_t i, double *R){
   printf("R[%2zu] = ", i);
   for(size_t j = 0; j <= i; ++j){
      printf("%f ", R[j]);
   }
   printf("\n");
}



int main(int argc,char **argv)
{
  double *args;
  long order = atol(argv[1]);
  args = calloc(3,sizeof(double));
  args[0] = 1.0;
  printf("romberg(l*exp(-u*l),a,b,3)=%f\n",romberg(&f,args,3, 0.0,1.0,order,0.0000001));
  return 0;
}
 
#endif
