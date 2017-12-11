#ifndef MIGRATION_RANDOM
#define MIGRATION_RANDOM
/* -------------------------------------------------------
   R A N D O M   G E N E R A T O R   R O U T I N E S 
 
   creates options structures,
   reads options from parmfile if present
 
   prints options,
   and finally helps to destroy itself.
 
   Peter Beerli 1996, Seattle
   beerli@fsu.edu
   
(c) Peter Beerli 2013 Tallahassee FL
 
Addition of Quasi random module Hongmei Chi (2005)

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

 
$Id: random.h 2157 2013-04-16 22:27:40Z beerli $
   ------------------------------------------------------- */
extern long *seed;
extern long *newseed;
extern char *generator;
//#ifdef MERSENNE_TWISTER
//#include "SFMT.h"
//#endif
/*-----------------------------------------------------*
 * calculates the start seed for randum
 * picks up the global variable iseed
 * PB 94
 *-----------------------------------------------------*/
extern void getseed (option_fmt * options);

/*-----------------------------------------------------*
 * RANDUM 
 * generates an uniform random number between o..1
 * using seed (generated by getseed)
 * JF <93 (6bit) and Mary K. Kuhner 1997 (32bit)
 * RANDINT
 * generates an uniform random integer number
 * between a..b
 * using seed (generated by getseed)
 *-----------------------------------------------------*/
#ifdef PTHREADS
extern MYREAL randum_thread (void);
extern long random_integer(long low, long high);
#define RANDUM randum_thread
#ifndef QUASIRANDOM
#define UNIF_RANDUM randum_thread     
#endif
#define RANDINT   random_integer
#define RANDDOUBLE(a,b)  (MYREAL) ((a) + (randum_thread() * ((b) - (a))))  
#else
extern MYREAL randum (void);
extern long random_integer(long low, long high);
#define RANDUM randum
#ifndef QUASIRANDOM
#define UNIF_RANDUM randum     
#endif
#define RANDINT   random_integer
#define RANDDOUBLE(a,b)  (MYREAL) ((a) + (randum() * ((b) - (a))))  
#endif
#endif
#ifdef QUASIRANDOM
extern double get_quasi();
#define UNIF_RANDUM  get_quasi

#endif
extern MYREAL random_beta(MYREAL a, MYREAL b);
extern MYREAL gamma_rand(MYREAL a, MYREAL b);
extern MYREAL trunc_gamma_rand(MYREAL alpha, MYREAL beta, MYREAL lower, MYREAL upper);
extern void assign_random_startsites(long **random_startsites, long fulllength, long shortsites, long rrepeats);