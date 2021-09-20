/*********************************************************************
*                           L I B _ R N D                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
**********************************************************************
*                                                                    *
* Random number generator.                                           *
*                                                                    *
**********************************************************************
* (c) 1997, 1999 Makoto Matsumoto and Takuji Nishimura               *
* This library comes with its own artistic license appended below.   *
*********************************************************************/

/* A C-program for MT19937: Real number version (1999/10/28)         */
/* genrand() generates one pseudorandom real number (double)         */
/* which is uniformly distributed on [0,1]-interval, for each        */
/* call. sgenrand(seed) sets initial values to the working area      */
/* of 624 words. Before genrand(), sgenrand(seed) should be          */
/* called once. (seed is any 32-bit integer.)                        */
/* Integer generator is obtained by modifying two lines.             */
/* Coded by Takuji Nishimura, considering the suggestions by         */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.                 */
/* This library is free software under the artistic license          */
/* appended to the source code below.                                */
/* Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura.   */
/* Any feedback is very welcome. For any question, comments,         */
/* see http://www.math.keio.ac.jp/matumoto/emt.html or email         */
/* matumoto@math.keio.ac.jp                                          */
/* REFERENCE                                                         */
/* M. Matsumoto and T. Nishimura,                                    */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform    */
/* Pseudo-Random Number Generator",                                  */
/* ACM Transactions on Modeling and Computer Simulation,             */
/* Vol. 8, No. 1, January 1998, pp 3--30.                            */


#include <stdio.h>
#include "lib_rnd.h"

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

/* global variables */
static unsigned long g_mt[N]; /* the array for the state vector  */
static int g_mti = N + 1; /* g_mti==N+1 means g_mt[N] is not initialized */

/*===================================================================*/
/* Initialization of array with a seed. Theoretically,               */
/* there are 2^19937-1 possible states as an intial state.           */
/* This function allows to choose any of 2^19937-1 ones.             */
/* Essential bits in "seed_array[]" is following 19937 bits:         */
/* (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1].  */
/* (seed_array[0]&LOWER_MASK) is discarded.                          */
/* Theoretically,                                                    */
/* (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1]   */
/* can take any values except all zeros.                             */
void sgenrand(unsigned long seed)
{
  int i;

  for (i = 0; i < N; i++) {
    g_mt[i] = seed & 0xffff0000;
    seed = 69069 * seed + 1;
    g_mt[i] |= (seed & 0xffff0000) >> 16;
    seed = 69069 * seed + 1;
  }
  g_mti = N;
}

/*===================================================================*/
/* real random number generator */
double genrand()
{
  unsigned long y;
  static unsigned long mag01[2] = {0x0, MATRIX_A};

  if (g_mti >= N) { /* generate N words at one time */
    int kk;
    if (g_mti == N + 1) /* if sgenrand() has not been called, */
      sgenrand(4357); /* a default initial seed is used   */

    for (kk = 0; kk < N - M; kk++) {
      y = (g_mt[kk] & UPPER_MASK) | (g_mt[kk + 1] & LOWER_MASK);
      g_mt[kk] = g_mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (; kk < N - 1; kk++) {
      y = (g_mt[kk] & UPPER_MASK) | (g_mt[kk + 1] & LOWER_MASK);
      g_mt[kk] = g_mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (g_mt[N - 1] & UPPER_MASK) | (g_mt[0] & LOWER_MASK);
    g_mt[N - 1] = g_mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];
    g_mti = 0;
  }

  y = g_mt[g_mti++];
  y ^= TEMPERING_SHIFT_U(y);
  y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
  y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
  y ^= TEMPERING_SHIFT_L(y);

  return ((double)y * 2.3283064370807974e-10);
}

/*                         The "Artistic License"

                                Preamble

The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some
semblance of artistic control over the development of the package,
while giving the users of the package the right to use and distribute
the Package in a more-or-less customary fashion, plus the right to make
reasonable modifications.

Definitions:

    - "Package" refers to the collection of files distributed by the
      Copyright Holder, and derivatives of that collection of files
      created through textual modification.

    - "Standard Version" refers to such a Package if it has not been
      modified, or has been modified in accordance with the wishes of
      the Copyright Holder as specified below.

    - "Copyright Holder" is whoever is named in the copyright or
      copyrights for the package.

    - "You" is you, if you're thinking about copying or distributing
      this Package.

    - "Reasonable copying fee" is whatever you can justify on the basis
      of media cost, duplication charges, time of people involved, and
      so on. (You will not be required to justify it to the Copyright
      Holder, but only to the computing community at large as a market
      that must bear the fee.)

    - "Freely Available" means that no fee is charged for the item
      itself, though there may be fees involved in handling the item.
      It also means that recipients of the item may redistribute it
      under the same conditions they received it.

1.  You may make and give away verbatim copies of the source form of the
    Standard Version of this Package without restriction, provided that
    you duplicate all of the original copyright notices and associated
    disclaimers.

2.  You may apply bug fixes, portability fixes and other modifications
    derived from the Public Domain or from the Copyright Holder. A
    Package modified in such a way shall still be considered the
    Standard Version.

3.  You may otherwise modify your copy of this Package in any way,
    provided that you insert a prominent notice in each changed file
    stating how and when you changed that file, and provided that you do
    at least ONE of the following:

    a) place your modifications in the Public Domain or otherwise make
       them Freely Available, such as by posting said modifications to
       Usenet or an equivalent medium, or placing the modifications on a
       major archive site such as uunet.uu.net, or by allowing the
       Copyright Holder to include your modifications in the Standard
       Version of the Package.

    b) use the modified Package only within your corporation or
       organization.

    c) rename any non-standard executables so the names do not conflict
       with standard executables, which must also be provided, and
       provide a separate manual page for each non-standard executable
       that clearly documents how it differs from the Standard Version.

    d) make other distribution arrangements with the Copyright Holder.

4.  You may distribute the programs of this Package in object code or
    executable form, provided that you do at least ONE of the following:

    a) distribute a Standard Version of the executables and library
       files, together with instructions (in the manual page or
       equivalent) on where to get the Standard Version.

    b) accompany the distribution with the machine-readable source of
       the Package with your modifications.

    c) give non-standard executables non-standard names, and clearly
       document the differences in manual pages (or equivalent),
       together with instructions on where to get the Standard Version.

    d) make other distribution arrangements with the Copyright Holder.

5.  You may charge a reasonable copying fee for any distribution of this
    Package. You may charge any fee you choose for support of this
    Package. You may not charge a fee for this Package itself. However,
    you may distribute this Package in aggregate with other (possibly
    commercial) programs as part of a larger (possibly commercial)
    s.ftware distribution provided that you do not advertise this
    Package as a product of your own.

6.  The scripts and library files supplied as input to or produced as
    output from the programs of this Package do not automatically fall
    under the copyright of this Package, but belong to whomever
    generated them, and may be sold commercially, and may be aggregated
    with this Package.

7.  C subroutines (or comparably compiled subroutines in other
    languages) supplied by you and linked into this Package shall not be
    considered part of this Package.

8.  Aggregation of this Package with a commercial distribution is always
    permitted provided that the use of this Package is embedded; that
    is, when no overt attempt is made to make this Package's interfaces
    visible to the end user of the commercial distribution. Such use
    shall not be construed as a distribution of this Package.

9.  The name of the Copyright Holder may not be used to endorse or
    promote products derived from this s.ftware without specific prior
    written permission.

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
    WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
    MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.

                                The End

*/

