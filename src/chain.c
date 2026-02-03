/* This file is part of AMC, a library for Affine Motion Calculation
 *
 * Copyright (C) 2024 Michael Carley
 *
 * AMC is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. AMC is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AMC.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>

#include "amc.h"

#include "amc-private.h"

/**
 *
 * @addtogroup chain
 *
 * @{
 * 
 */

/** 
 * Allocate a transform chain, a sequence of transforms to be applied
 * successively
 * 
 * @param nt maximum number of transforms in chain.
 * 
 * @return newly allocated ::amc_transform_chain_t.
 */

amc_transform_chain_t *amc_transform_chain_alloc(int nt)

{
  amc_transform_chain_t *C ;
  
  C = (amc_transform_chain_t *)g_malloc0(sizeof(amc_transform_chain_t)) ;

  memset(C, 0, sizeof(amc_transform_chain_t)) ;

  amc_transform_chain_transform_number(C)     = 0 ;
  amc_transform_chain_transform_number_max(C) = nt ;

  C->T = (amc_transform_t **)g_malloc0(nt*sizeof(amc_transform_t *)) ;
  memset(C->T, 0, nt*sizeof(amc_transform_t *)) ;

  return C ;
}

/** 
 * Add a transform to a chain of transforms
 * 
 * @param C an allocated ::amc_transform_chain_t;
 * @param T ::amc_transform_t to be added to chain.
 * 
 * @return 0 on success, or 1 if the maximum number of transforms in
 * \a C has been reached.
 */

int amc_transform_chain_transform_add(amc_transform_chain_t *C,
				      amc_transform_t *T)

{
  if ( amc_transform_chain_transform_number(C) >=
       amc_transform_chain_transform_number_max(C) ) {
    fprintf(stderr, "%s: not enough space allocated for %d transforms\n",
	    __FUNCTION__, amc_transform_chain_transform_number(C) + 1) ;
    return 1 ;
  }

  amc_transform_chain_transform(C,amc_transform_chain_transform_number(C)) = T ;
  amc_transform_chain_transform_number(C) ++ ;
  
  return 0 ;
}

static gint amc2d_transform_chain_evaluate(amc_transform_chain_t *C, gint order,
					   amc_transform_t *T)

{
  gint i, j ;
  gdouble *A, *B, *Bdot, tmp[9] ;
  amc_transform_t *S ;
  
  amc_transform_matrix_identity(T, 0) ;

  A = amc_transform_matrix(T, 0) ;

  /*order 0*/
  for ( i = 0 ; i < amc_transform_chain_transform_number(C) ; i ++ ) {
    S = amc_transform_chain_transform(C, i) ;
    B = amc_transform_matrix(S, 0) ;
    amc2d_matrix_matrix_mul(1.0, B, A, 0, A) ;
  }

  if ( order == 0 ) return 0 ;

  amc_transform_matrix_zero(T, 1) ;

  A = amc_transform_matrix(T, 1) ;

  for ( i = 0 ; i < amc_transform_chain_transform_number(C) ; i ++ ) {
    amc2d_matrix_identity(tmp) ;
    for ( j = 0 ; j < i ; j ++ ) {
      S = amc_transform_chain_transform(C, j) ;
      B = amc_transform_matrix(S, 0) ;
      amc2d_matrix_matrix_mul(1.0, B, tmp, 0.0, tmp) ;
    }
    S = amc_transform_chain_transform(C, i) ;
    Bdot = amc_transform_matrix(S, 1) ;
    amc2d_matrix_matrix_mul(1.0, Bdot, tmp, 0.0, tmp) ;
    for ( j = i+1 ; j < amc_transform_chain_transform_number(C) ; j ++ ) {
      S = amc_transform_chain_transform(C, j) ;
      B = amc_transform_matrix(S, 0) ;
      amc2d_matrix_matrix_mul(1.0, B, tmp, 0.0, tmp) ;
    }
    for ( j = 0 ; j < 9 ; j ++ ) A[j] += tmp[j] ;
  }

  if ( order == 1 ) return 0 ;

  fprintf(stderr,
	  "%s: shouldn't get here (higher derivatives not implemented yet\n",
	  __FUNCTION__) ;
  
  return 0 ;
}

static gint amc3d_transform_chain_evaluate(amc_transform_chain_t *C, gint order,
					  amc_transform_t *T)

{
  gint i, j ;
  gdouble *A, *B, *Bdot, tmp[16] ;
  amc_transform_t *S ;
  
  amc_transform_matrix_identity(T, 0) ;

  A = amc_transform_matrix(T, 0) ;

  /*order 0*/
  for ( i = 0 ; i < amc_transform_chain_transform_number(C) ; i ++ ) {
    S = amc_transform_chain_transform(C, i) ;
    B = amc_transform_matrix(S, 0) ;
    amc3d_matrix_matrix_mul(1.0, B, A, 0, A) ;
  }

  if ( order == 0 ) return 0 ;

  amc_transform_matrix_zero(T, 1) ;

  A = amc_transform_matrix(T, 1) ;

  for ( i = 0 ; i < amc_transform_chain_transform_number(C) ; i ++ ) {
    amc3d_matrix_identity(tmp) ;
    for ( j = 0 ; j < i ; j ++ ) {
      S = amc_transform_chain_transform(C, j) ;
      B = amc_transform_matrix(S, 0) ;
      amc3d_matrix_matrix_mul(1.0, B, tmp, 0.0, tmp) ;
    }
    S = amc_transform_chain_transform(C, i) ;
    Bdot = amc_transform_matrix(S, 1) ;
    amc3d_matrix_matrix_mul(1.0, Bdot, tmp, 0.0, tmp) ;
    for ( j = i+1 ; j < amc_transform_chain_transform_number(C) ; j ++ ) {
      S = amc_transform_chain_transform(C, j) ;
      B = amc_transform_matrix(S, 0) ;
      amc3d_matrix_matrix_mul(1.0, B, tmp, 0.0, tmp) ;
    }
    for ( j = 0 ; j < 16 ; j ++ ) A[j] += tmp[j] ;
  }

  if ( order == 1 ) return 0 ;

  fprintf(stderr,
	  "%s: shouldn't get here (higher derivatives not implemented yet\n",
	  __FUNCTION__) ;
  
  return 0 ;
}

/** 
 * Evaluate a sequence of transforms in a chain, including time derivatives
 * 
 * @param C chain of transforms, which have been evaluated at required time; 
 * @param order maximum derivative to evaluate;
 * @param T on exit, contains sequence of transforms, including time 
 * derivatives, which can be applied to a point.
 * 
 * @return 0 on success.
 */

int amc_transform_chain_evaluate(amc_transform_chain_t *C, gint order,
				 amc_transform_t *T)

{
  amc_transform_check(T) ;

  if ( amc_transform_dimension(T) == 2 ) {
    return amc2d_transform_chain_evaluate(C, order, T) ;
  }
  
  return amc3d_transform_chain_evaluate(C, order, T) ;

  return 0 ;
}

/**
 *
 * @}
 * 
 */
