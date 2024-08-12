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

#include "tinyexpr.h"

#include "amc.h"

#ifdef HAVE_LIBMATHEVAL
#include <matheval.h>
#endif /*HAVE_LIBMATHEVAL*/

/**
 * @file   amc.c
 * @author  <michael@michael.paraffinalia.co.uk>
 * @date   Fri May 17 09:08:48 2024
 * 
 * @brief  
 * 
 * 
 */

/** 
 * Allocate an AMC transform of dimension \a dim and maximum
 * derivative order \a order
 * 
 * @param dim dimension of transform (must be two or three);
 * @param order highest order of time derivative of transform to be evaluated.
 * 
 * @return newly allocated ::amc_transform_t
 */

amc_transform_t *amc_transform_alloc(int dim, gint order)

{
  amc_transform_t *T ;

  if ( dim != 2 && dim != 3 ) {
    g_error("%s: dimension (%d) must be 2 or 3\n", __FUNCTION__, dim) ;
  }
  
  T = (amc_transform_t *)g_malloc0(sizeof(amc_transform_t)) ;
  memset(T, 0, sizeof(amc_transform_t)) ;

  amc_transform_dimension(T)       = dim ;
  amc_transform_matrix_size(T)     = ( dim == 2 ? 9 : 16 ) ;
  amc_transform_order_max(T)       = order ;
  amc_transform_order(T)           = -1 ;
  amc_transform_variable_number(T) = 0 ;
  amc_transform_tag(T)             = 0 ;
  amc_transform_compiled(T)        = 0 ;
  /*matrices have msize defined entries*/
  T->matrix_values  =
    (gdouble *)g_malloc0((T->msize)*(order+1)*sizeof(gdouble)) ;
  T->matrix_strings =
    (char **) g_malloc0((T->msize)*(order+1)*sizeof(char *)) ;
  T->matrix_expr    = g_malloc0((T->msize)*(order+1)*sizeof(gpointer)) ;
    /* (te_expr **) g_malloc0((T->msize)*(order+1)*sizeof(te_expr **)) ; */
  T->vars           = g_malloc0(AMC_TRANSFORM_VARIABLE_NUMBER*
			     sizeof(te_variable)) ;
  /*set everything to 0 or NULL for checking later*/
  memset(T->matrix_values , 0, (T->msize)*(order+1)*sizeof(gdouble)) ;
  memset(T->matrix_strings, 0, (T->msize)*(order+1)*sizeof(char *)) ;
  memset(T->matrix_expr   , 0, (T->msize)*(order+1)*sizeof(te_expr *)) ;

  amc_transform_variable_add(T, "t", &(T->t)) ;
  
  return T ;
}

/** 
 * Set an element of a transform as symbolic expression or numerical
 * constant
 * 
 * @param T ::amc_transform_t
 * @param order derivative order whose element is to be set
 * @param i row index of element (less than three for two-dimensional 
 * transform, less than four for three-dimensional);
 * @param j column index of element (less than three for two-dimensional 
 * transform, less than four for three-dimensional);
 * @param val numerical value of element;
 * @param str if not NULL, symbolic expression for element, which will be
 * evaluated using parser at each time, overwriting \a val, and which can
 * be symbolically differentiated.
 * 
 * @return 0 on success.
 */

int amc_transform_entry_set(amc_transform_t *T, gint order,
			    gint i, gint j, gdouble val, char *str)

{
  amc_transform_check(T) ;

  if ( order > amc_transform_order(T) ) {
    fprintf(stderr, "%s: order (%d) greater than transform order (%d)\n",
	    __FUNCTION__, order, amc_transform_order(T)) ;
    return 1 ;
  }

  if ( i > amc_transform_dimension(T) ) {
    fprintf(stderr,
	    "%s: row index (%d) must be less than or equal to "
	    "dimension (%d)\n",
	    __FUNCTION__, i, amc_transform_dimension(T)) ;
    return 1 ;
  }

  if ( j > amc_transform_dimension(T) ) {
    fprintf(stderr,
	    "%s: column index (%d) must be less than or equal to "
	    "dimension (%d)\n",
	    __FUNCTION__, i, amc_transform_dimension(T)) ;
    return 1 ;
  }

  if ( amc_transform_expr(T, order, i, j) != NULL ) {
    te_free(amc_transform_expr(T, order, i, j)) ;
  }

  amc_transform_expr(T, order, i, j) = NULL ;
  amc_transform_value(T, order, i, j) = val ;

  if ( str == NULL ) {
    amc_transform_string(T, order, i, j) = NULL ;
    return 0 ;
  }
  
  amc_transform_string(T, order, i, j) = strdup(str) ;
  
  return 0 ;
}

/** 
 * Write elements of the derivative of a transform matrix to file
 * 
 * @param f file stream for output;
 * @param T ::amc_transform_t whose entries are to be written;
 * @param order derivative order whose entries are to be written.
 * 
 * @return 0 on success.
 */

int amc_transform_matrix_write(FILE *f, amc_transform_t *T, gint order)

{
  gint i, j ;

  amc_transform_check(T) ;

  if ( order > amc_transform_order(T) ) {
    fprintf(stderr, "%s: order (%d) greater than transform order (%d)\n",
	    __FUNCTION__, order, amc_transform_order(T)) ;
    return 1 ;
  }

  for ( i = 0 ; i < amc_transform_dimension(T) + 1 ; i ++ ) {
    for ( j = 0 ; j < amc_transform_dimension(T) + 1 ; j ++ ) {
      if ( amc_transform_string(T,order,i,j) != NULL ) {
	fprintf(f, "\"%s\" ", amc_transform_string(T,order,i,j)) ;
      } else {
	fprintf(f, "%lg ", amc_transform_value(T,order,i,j)) ;
      }
    }
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

/** 
 * Compile symbolic expressions in a transform for future
 * evaluation. This must be done after symbolic expressions are
 * defined and before evaluating or applying the transform.
 * 
 * @param T an ::amc_transform_t;
 * 
 * @return 0 on success.
 */

int amc_transform_expressions_compile(amc_transform_t *T)

{
  gint i, error ;

  amc_transform_check(T) ;

  for ( i = 0 ;
	i < amc_transform_matrix_size(T)*(amc_transform_order(T)+1) ;
	i ++ ) {
    if ( T->matrix_expr[i] != NULL ) te_free(T->matrix_expr[i]) ;
    if ( T->matrix_strings[i] != NULL ) {
      T->matrix_expr[i] = te_compile(T->matrix_strings[i],
				     T->vars, T->nvars, &error) ;
      if ( error != 0 ) {
	fprintf(stderr, "%s: cannot compile expression \"%s\"\n",
		__FUNCTION__, T->matrix_strings[i]) ;
      }
    }
  }

  amc_transform_compiled(T) = 1 ;
  
  return 0 ;
}

static gint amc2d_matrix_vector_mul(gdouble *y, gdouble *A, gdouble *x)

/*
 * y := A*x (safe to perform in place with x == y)
 */
  
{
  gdouble tmp ;
  
  tmp  = A[0]*x[0] + A[1]*x[1] + A[2] ;
  y[1] = A[3]*x[0] + A[4]*x[1] + A[5] ;
  y[0] = tmp ;

  return 0 ;
}

static void amc2d_matrix_matrix_mul(gdouble al, gdouble *A, gdouble *B,
				    gdouble bt, gdouble *C)

/*
 * C := bt*C + al*A*B (safe to perform in place with C == A or B)
 */
  
{
  gdouble Ctmp[9] ;
  gint i ;

  Ctmp[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6] ;
  Ctmp[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7] ;
  Ctmp[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8] ;
  Ctmp[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6] ;
  Ctmp[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7] ;
  Ctmp[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8] ;
  Ctmp[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6] ;
  Ctmp[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7] ;
  Ctmp[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8] ;

  for ( i = 0 ; i < 9 ; i ++ ) C[i] = bt*C[i] + al*Ctmp[i] ;

  return ;
}

static gint amc3d_matrix_vector_mul(gdouble *y, gdouble *A, gdouble *x)

/*
 * y := A*x (safe to perform in place with x == y)
 */
  
{
  gdouble tmp[2] ;
  
  tmp[0]  = A[ 0]*x[0] + A[ 1]*x[1] + A[ 2]*x[2] + A[ 3] ;
  tmp[1]  = A[ 4]*x[0] + A[ 5]*x[1] + A[ 6]*x[2] + A[ 7] ;
  y[2]    = A[ 8]*x[0] + A[ 9]*x[1] + A[10]*x[2] + A[11] ;
  y[0] = tmp[0] ; y[1] = tmp[1] ;

  return 0 ;
}

static void amc3d_matrix_matrix_mul(gdouble al, gdouble *A, gdouble *B,
				    gdouble bt, gdouble *C)

/*
 * C := bt*C + al*A*B (safe to perform in place with C == A or B)
 */
  
{
  gdouble Ctmp[16] ;
  gint i ;

  Ctmp[ 0] = A[ 0]*B[ 0] + A[ 1]*B[ 4] + A[ 2]*B[ 8] + A[ 3]*B[12] ;
  Ctmp[ 1] = A[ 0]*B[ 1] + A[ 1]*B[ 5] + A[ 2]*B[ 9] + A[ 3]*B[13] ;
  Ctmp[ 2] = A[ 0]*B[ 2] + A[ 1]*B[ 6] + A[ 2]*B[10] + A[ 3]*B[14] ;
  Ctmp[ 3] = A[ 0]*B[ 3] + A[ 1]*B[ 7] + A[ 2]*B[11] + A[ 3]*B[15] ;

  Ctmp[ 4] = A[ 4]*B[ 0] + A[ 5]*B[ 4] + A[ 6]*B[ 8] + A[ 7]*B[12] ;
  Ctmp[ 5] = A[ 4]*B[ 1] + A[ 5]*B[ 5] + A[ 6]*B[ 9] + A[ 7]*B[13] ;
  Ctmp[ 6] = A[ 4]*B[ 2] + A[ 5]*B[ 6] + A[ 6]*B[10] + A[ 7]*B[14] ;
  Ctmp[ 7] = A[ 4]*B[ 3] + A[ 5]*B[ 7] + A[ 6]*B[11] + A[ 7]*B[15] ;

  Ctmp[ 8] = A[ 8]*B[ 0] + A[ 9]*B[ 4] + A[10]*B[ 8] + A[11]*B[12] ;
  Ctmp[ 9] = A[ 8]*B[ 1] + A[ 9]*B[ 5] + A[10]*B[ 9] + A[11]*B[13] ;
  Ctmp[10] = A[ 8]*B[ 2] + A[ 9]*B[ 6] + A[10]*B[10] + A[11]*B[14] ;
  Ctmp[11] = A[ 8]*B[ 3] + A[ 9]*B[ 7] + A[10]*B[11] + A[11]*B[15] ;
  
  Ctmp[12] = A[12]*B[ 0] + A[13]*B[ 4] + A[14]*B[ 8] + A[15]*B[12] ;
  Ctmp[13] = A[12]*B[ 1] + A[13]*B[ 5] + A[14]*B[ 9] + A[15]*B[13] ;
  Ctmp[14] = A[12]*B[ 2] + A[13]*B[ 6] + A[14]*B[10] + A[15]*B[14] ;
  Ctmp[15] = A[12]*B[ 3] + A[13]*B[ 7] + A[14]*B[11] + A[15]*B[15] ;
  
  for ( i = 0 ; i < 16 ; i ++ ) C[i] = bt*C[i] + al*Ctmp[i] ;

  return ;
}

/** 
 * Apply a transform to a point
 * 
 * @param T ::amc_transform_t to apply;
 * @param order order of derivative of \a T;
 * @param xin input point;
 * @param xout output point (can be the same as \a xin).
 * 
 * @return 0 on success.
 */

int amc_transform_matrix_apply(amc_transform_t *T, gint order,
			       gdouble *xin, gdouble *xout)

{
  gdouble *A ;

  amc_transform_check(T) ;
  
  A = amc_transform_matrix(T, order) ;

  if ( amc_transform_dimension(T) == 2 ) {
    return amc2d_matrix_vector_mul(xout, A, xin) ;
  }
  
  return amc3d_matrix_vector_mul(xout, A, xin) ;

  return 0 ;
}

/** 
 * Apply a transform to a vector (such as a unit normal)
 * 
 * @param T ::amc_transform_t to apply;
 * @param order order of derivative of \a T;
 * @param xin input vector;
 * @param xout output vector (can be the same as \a xin).
 * 
 * @return 0 on success.
 */

int amc_transform_matrix_apply_vector(amc_transform_t *T, gint order,
				      gdouble *xin, gdouble *xout)

{
  gdouble *A, zero[3]={0,0,0}, dx0[3], dx1[3] ;

  amc_transform_check(T) ;
  
  A = amc_transform_matrix(T, order) ;

  if ( amc_transform_dimension(T) == 2 ) {
    amc2d_matrix_vector_mul(dx0, A, zero) ;
    amc2d_matrix_vector_mul(dx1, A, xin) ;
    xout[0] = dx1[0] - dx0[0] ; 
    xout[1] = dx1[1] - dx0[1] ; 
    return 0 ;
  }
  
  amc3d_matrix_vector_mul(dx0, A, zero) ;
  amc3d_matrix_vector_mul(dx1, A, xin) ;
  xout[0] = dx1[0] - dx0[0] ; 
  xout[1] = dx1[1] - dx0[1] ; 
  xout[2] = dx1[2] - dx0[2] ; 

  return 0 ;
}

/** 
 * Evaluate numerical values of transform matrices for all derivatives
 * up to maximum defined order of \a T
 * 
 * @param T an ::amc_transform_t;
 * @param t time for evaluation.
 * 
 * @return 0 on success.
 */

int amc_transform_matrices_evaluate(amc_transform_t *T, gdouble t)

{
  gint i ;

  T->t = t ;
  for ( i = 0 ;
	i < amc_transform_matrix_size(T)*(amc_transform_order(T)+1) ;
	i ++ ) {
    if ( T->matrix_expr[i] != NULL )
      T->matrix_values[i] = te_eval(T->matrix_expr[i]) ;
  }

  return 0 ;
}

/** 
 * Add a symbolic variable to a transform.
 * 
 * @param T an ::amc_transform_t;
 * @param str name of variable;
 * @param adr pointer to address of variable for compilation and evaluation.
 * 
 * @return 0 on success.
 */

int amc_transform_variable_add(amc_transform_t *T, char *str, gdouble *adr)
  
{
  gint i ;
  te_variable *vars = (te_variable *)(T->vars) ;
  
  if ( amc_transform_variable_number(T) >= AMC_TRANSFORM_VARIABLE_NUMBER ) {
    fprintf(stderr, "%s: not enough space for %d variablesn",	
	    __FUNCTION__, amc_transform_variable_number(T)+1) ;	
    return 1 ;
  }									

  for ( i = 0 ; i < amc_transform_variable_number(T) ; i ++ ) {
    if ( strcmp(str, vars[i].name) == 0 ) {
      fprintf(stderr, "%s: variable \"%s\" already in transform\n",
	      __FUNCTION__, str) ;
      return 1 ;
    }
  }    

  vars[amc_transform_variable_number(T)].name = strdup(str) ;
  vars[amc_transform_variable_number(T)].address = adr ;
  amc_transform_variable_number(T) ++ ;

  return 0 ;
}

/** 
 * Write variables defined in a transform to file.
 * 
 * @param f output file stream;
 * @param T an ::amc_transform_t whose variables are to be listed.
 * 
 * @return 0 on success.
 */

int amc_transform_variables_write(FILE *f, amc_transform_t *T)
  
{									
  gint i ;								
  te_variable *vars = (te_variable *)(T->vars) ;

  for ( i = 0 ; i < amc_transform_variable_number(T) ; i ++ ) {	
    fprintf(f, "%s = %lg\n", vars[i].name, *((gdouble *)(vars[i].address))) ;
  }

  return 0 ;
}

/** 
 * Allocate a transform chain, a sequence of transforms to be applied
 * successively
 * 
 * @param ntrans maximum number of transforms in chain.
 * 
 * @return newly allocated ::amc_transform_chain_t.
 */

amc_transform_chain_t *amc_transform_chain_alloc(int ntrans)

{
  amc_transform_chain_t *C ;
  
  C = (amc_transform_chain_t *)g_malloc0(sizeof(amc_transform_chain_t)) ;

  memset(C, 0, sizeof(amc_transform_chain_t)) ;

  amc_transform_chain_transform_number(C)     = 0 ;
  amc_transform_chain_transform_number_max(C) = ntrans ;

  C->T = (amc_transform_t **)g_malloc0(ntrans*sizeof(amc_transform_t *)) ;
  memset(C->T, 0, ntrans*sizeof(amc_transform_t *)) ;

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

static gint amc2d_matrix_identity(gdouble *A)

{
  A[0] = 1 ; A[1] = 0 ; A[2] = 0 ;
  A[3] = 0 ; A[4] = 1 ; A[5] = 0 ;
  A[6] = 0 ; A[7] = 0 ; A[8] = 1 ;

  return 0 ;
}

static gint amc3d_matrix_identity(gdouble *A)

{
  A[ 0] = 1 ; A[ 1] = 0 ; A[ 2] = 0 ; A[ 3] = 0 ;
  A[ 4] = 0 ; A[ 5] = 1 ; A[ 6] = 0 ; A[ 7] = 0 ;
  A[ 8] = 0 ; A[ 9] = 0 ; A[10] = 1 ; A[11] = 0 ;
  A[12] = 0 ; A[13] = 0 ; A[14] = 0 ; A[15] = 1 ;

  return 0 ;
}

/** 
 * Set a transform to the identity matrix
 * 
 * @param T an ::amc_transform_t;
 * @param order derivative whose matrix is to be set to identity.
 * 
 * @return 0 on success.
 */

int amc_transform_matrix_identity(amc_transform_t *T, gint order)

{
  gdouble *A ;
  
  amc_transform_check(T) ;

  A = amc_transform_matrix(T, order) ;

  if ( amc_transform_dimension(T) == 2 ) {
    return amc2d_matrix_identity(A) ;
  }
  
  return amc3d_matrix_identity(A) ;
}

static gint amc2d_matrix_zero(gdouble *A)

{
  A[0] = 0 ; A[1] = 0 ; A[2] = 0 ;
  A[3] = 0 ; A[4] = 0 ; A[5] = 0 ;
  A[6] = 0 ; A[7] = 0 ; A[8] = 0 ;
  
  return 0 ;
}

static gint amc3d_matrix_zero(gdouble *A)

{
  A[ 0] = 0 ; A[ 1] = 0 ; A[ 2] = 0 ; A[ 3] = 0 ;
  A[ 4] = 0 ; A[ 5] = 0 ; A[ 6] = 0 ; A[ 7] = 0 ;
  A[ 8] = 0 ; A[ 9] = 0 ; A[10] = 0 ; A[11] = 0 ;
  A[12] = 0 ; A[13] = 0 ; A[14] = 0 ; A[15] = 0 ;
  
  return 0 ;
}

/** 
 * Set a transform matrix to all zeros.
 * 
 * @param T an ::amc_transform_t
 * @param order derivative whose matrix entries are to be set to zero.
 * 
 * @return 0 on success.
 */

int amc_transform_matrix_zero(amc_transform_t *T, gint order)

{
  gdouble *A ;
  
  amc_transform_check(T) ;

  A = amc_transform_matrix(T, order) ;
  if ( amc_transform_dimension(T) == 2 ) {
    return amc2d_matrix_zero(A) ;
  }

  return amc3d_matrix_zero(A) ;
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
 * derivatives. 
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
 * Numerically estimate derivative of transformed  point
 * 
 * @param C an ::amc_transform_chain_t;
 * @param order order of derivatives to evaluate;
 * @param S workspace for transform;
 * @param t time for evaluation;
 * @param x input point;
 * @param dx on output contains estimated time derivative of transform of \a x.
 * 
 * @return 0 on success.
 */

int amc_transform_chain_derivative(amc_transform_chain_t *C, gint order,
				   amc_transform_t *S, gdouble t,
				   gdouble *x, gdouble *dx)
{
  gint i ;
  gdouble dt, x0[3], x1[3] ;

  dt = 1e-6 ;
  for ( i = 0 ; i < amc_transform_chain_transform_number(C) ; i ++ ) {
    amc_transform_matrices_evaluate(amc_transform_chain_transform(C,i),
				    t-dt/2) ;
  }
  
  amc_transform_chain_evaluate(C, order, S) ;
  amc_transform_matrix_apply(S, order, x, x0) ;

  for ( i = 0 ; i < amc_transform_chain_transform_number(C) ; i ++ ) {
    amc_transform_matrices_evaluate(amc_transform_chain_transform(C,i),
				    t+dt/2) ;
  }
  
  amc_transform_chain_evaluate(C, order, S) ;
  amc_transform_matrix_apply(S, order, x, x1) ;

  dx[0] = (x1[0] - x0[0])/dt ;
  dx[1] = (x1[1] - x0[1])/dt ;

  if ( amc_transform_dimension(S) == 2 ) return 0 ;

  dx[2] = (x1[2] - x0[2])/dt ;
  
  return 0 ;
}

/** 
 * Set transform to translation by possibly time-dependent displacement
 * 
 * @param T ::amc_transform_t to set;
 * @param dx \f$x\f$ displacement;
 * @param xstr if not NULL, symbolic expression for \a dx;
 * @param dy \f$y\f$ displacement;
 * @param ystr if not NULL, symbolic expression for \a dy;
 * @param dz \f$z\f$ displacement (ignored for two-dimensional transform);
 * @param zstr if not NULL, symbolic expression for \a dz (ignored for
 * two-dimensional transform);
 * @param order maximum order of derivatives to evaluate.
 * 
 * @return 0 on success.
 */

int amc_transform_translation(amc_transform_t *T,
			      gdouble dx, char *xstr,
			      gdouble dy, char *ystr,
			      gdouble dz, char *zstr,
			      gint order)

/*
 * translation by (dx,dy,dx) or (xstr, ystr, zstr), z component
 * ignored for two-dimensional transformations
 */

{
  gint i ;
  gdouble d[3] = {dx, dy, dz} ;
  char *dstr[] = {xstr, ystr, zstr} ;
  
  amc_transform_check(T) ;

  if ( order > amc_transform_order_max(T) ) {
    fprintf(stderr, "%s: order (%d) greater than maximum "
	    "transform order (%d)\n",
	    __FUNCTION__, order, amc_transform_order_max(T)) ;
    return 1 ;
  }

  amc_transform_matrix_identity(T, 0) ;
  for ( i = 1 ; i <= amc_transform_order_max(T) ; i ++ ) {
    amc_transform_matrix_zero(T, i) ;
  }

  for ( i = 0 ; i < amc_transform_dimension(T) ; i ++ ) {
    amc_transform_entry_set(T, 0, i, amc_transform_dimension(T),
			    d[i], dstr[i]) ;
  }
  
  amc_transform_derivatives_evaluate(T, order) ;

  return 0 ;
}

/** 
 * Set two-dimensional transform to rotation in \f$x\f$-\f$y\f$ plane
 * 
 * @param T ::amc_transform_t to set;
 * @param th rotation angle \f$\theta\f$;
 * @param str if not NULL, symbolic expression for \f$\theta\f$, which 
 * overrides numerical value \a th;
 * @param order maximum order of derivatives to evaluate.
 * 
 * @return 0 on success.
 */

int amc_transform_rotation(amc_transform_t *T, gdouble th, char *str,
			   gint order)

{
  gint i ;
  char buf[256] ;
  
  if ( amc_transform_dimension(T) != 2 ) {
    fprintf(stderr, "%s: only defined for two-dimensional transform\n",
	    __FUNCTION__) ;
    return 1 ;
  }
  
  if ( order > amc_transform_order(T) ) {
    fprintf(stderr, "%s: order (%d) greater than transform order (%d)\n",
	    __FUNCTION__, order, amc_transform_order(T)) ;
    return 1 ;
  }
  if ( order > amc_transform_order_max(T) ) {
    fprintf(stderr, "%s: order (%d) greater than maximum "
	    "transform order (%d)\n",
	    __FUNCTION__, order, amc_transform_order_max(T)) ;
    return 1 ;
  }

  amc_transform_matrix_identity(T, 0) ;
  for ( i = 1 ; i <= amc_transform_order_max(T) ; i ++ ) {
    amc_transform_matrix_zero(T, i) ;
  }

  if ( str == NULL ) {
    amc_transform_entry_set(T, 0, 0, 0,  cos(th), NULL) ;
    amc_transform_entry_set(T, 0, 0, 1, -sin(th), NULL) ;
    amc_transform_entry_set(T, 0, 1, 0,  sin(th), NULL) ;
    amc_transform_entry_set(T, 0, 1, 1,  cos(th), NULL) ;
    
    return 0 ;
  }

  sprintf(buf, "cos(%s)", str) ;
  amc_transform_entry_set(T, 0, 0, 0, 0, buf) ;
  amc_transform_entry_set(T, 0, 1, 1, 0, buf) ;
  sprintf(buf, "-sin(%s)", str) ;
  amc_transform_entry_set(T, 0, 1, 0, 0, buf) ;
  sprintf(buf, "sin(%s)", str) ;
  amc_transform_entry_set(T, 0, 0, 1, 0, buf) ;
    
  amc_transform_derivatives_evaluate(T, order) ;

  return 0 ;
}

/** 
 * Set three-dimensional transform to rotation about \f$x\f$ axis
 * 
 * @param T ::amc_transform_t to set;
 * @param th rotation angle \f$\theta\f$ about \f$x\f$ axis;
 * @param str if not NULL, symbolic expression for \f$\theta\f$, which 
 * overrides numerical value \a th;
 * @param order maximum order of derivatives to evaluate.
 * 
 * @return 0 on success.
 */

int amc_transform_rotation_x(amc_transform_t *T, gdouble th, char *str,
			     gint order)

{
  gint i ;
  char buf[256] ;
  
  if ( order > amc_transform_order_max(T) ) {
    fprintf(stderr, "%s: order (%d) greater than maximum "
	    "transform order (%d)\n",
	    __FUNCTION__, order, amc_transform_order_max(T)) ;
    return 1 ;
  }

  if ( amc_transform_dimension(T) != 3 ) {
    fprintf(stderr, "%s: only defined for three-dimensional transform\n",
	    __FUNCTION__) ;
    return 1 ;
  }

  amc_transform_matrix_identity(T, 0) ;
  for ( i = 1 ; i <= amc_transform_order_max(T) ; i ++ ) {
    amc_transform_matrix_zero(T, i) ;
  }

  if ( str == NULL ) {
    amc_transform_entry_set(T, 0, 1, 1,  cos(th), NULL) ;
    amc_transform_entry_set(T, 0, 1, 2, -sin(th), NULL) ;
    amc_transform_entry_set(T, 0, 2, 1,  sin(th), NULL) ;
    amc_transform_entry_set(T, 0, 2, 2,  cos(th), NULL) ;
  
    return 0 ;
  }

  sprintf(buf, "cos(%s)", str) ;
  amc_transform_entry_set(T, 0, 1, 1, 0, buf) ;
  amc_transform_entry_set(T, 0, 2, 2, 0, buf) ;
  sprintf(buf, "-sin(%s)", str) ;
  amc_transform_entry_set(T, 0, 1, 2, 0, buf) ;
  sprintf(buf, "sin(%s)", str) ;
  amc_transform_entry_set(T, 0, 2, 1, 0, buf) ;

  amc_transform_derivatives_evaluate(T, order) ;
  
  return 0 ;
}

/** 
 * Set three-dimensional transform to rotation about \f$y\f$ axis
 * 
 * @param T ::amc_transform_t to set;
 * @param th rotation angle \f$\theta\f$ about \f$y\f$ axis;
 * @param str if not NULL, symbolic expression for \f$\theta\f$, which 
 * overrides numerical value \a th;
 * @param order maximum order of derivatives to evaluate.
 * 
 * @return 0 on success.
 */

int amc_transform_rotation_y(amc_transform_t *T, gdouble th, char *str,
			     gint order)			     

{
  gint i ;
  char buf[256] ;
  
  if ( order > amc_transform_order_max(T) ) {
    fprintf(stderr, "%s: order (%d) greater than transform order (%d)\n",
	    __FUNCTION__, order, amc_transform_order_max(T)) ;
    return 1 ;
  }

  if ( amc_transform_dimension(T) != 3 ) {
    fprintf(stderr, "%s: only defined for three-dimensional transform\n",
	    __FUNCTION__) ;
    return 1 ;
  }

  amc_transform_matrix_identity(T, 0) ;
  for ( i = 1 ; i <= amc_transform_order_max(T) ; i ++ ) {
    amc_transform_matrix_zero(T, i) ;
  }

  if ( str == NULL ) {
    amc_transform_entry_set(T, 0, 0, 0,  cos(th), NULL) ;
    amc_transform_entry_set(T, 0, 0, 2,  sin(th), NULL) ;
    amc_transform_entry_set(T, 0, 2, 0, -sin(th), NULL) ;
    amc_transform_entry_set(T, 0, 2, 2,  cos(th), NULL) ;
  
    return 0 ;
  }

  sprintf(buf, "cos(%s)", str) ;
  amc_transform_entry_set(T, 0, 0, 0, 0, buf) ;
  amc_transform_entry_set(T, 0, 2, 2, 0, buf) ;
  sprintf(buf, "-sin(%s)", str) ;
  amc_transform_entry_set(T, 0, 2, 0, 0, buf) ;
  sprintf(buf, "sin(%s)", str) ;
  amc_transform_entry_set(T, 0, 0, 2, 0, buf) ;
    
  amc_transform_derivatives_evaluate(T, order) ;

  return 0 ;
}

/** 
 * Set three-dimensional transform to rotation about \f$z\f$ axis
 * 
 * @param T ::amc_transform_t to set;
 * @param th rotation angle \f$\theta\f$ about \f$z\f$ axis;
 * @param str if not NULL, symbolic expression for \f$\theta\f$, which 
 * overrides numerical value \a th;
 * @param order maximum order of derivatives to evaluate.
 * 
 * @return 0 on success.
 */

int amc_transform_rotation_z(amc_transform_t *T, gdouble th, char *str,
			     gint order)
			     
{
  gint i ;
  char buf[256] ;
  
  if ( order > amc_transform_order_max(T) ) {
    fprintf(stderr, "%s: order (%d) greater than transform order (%d)\n",
	    __FUNCTION__, order, amc_transform_order_max(T)) ;
    return 1 ;
  }

  if ( amc_transform_dimension(T) != 3 ) {
    fprintf(stderr, "%s: only defined for three-dimensional transform\n",
	    __FUNCTION__) ;
    return 1 ;
  }

  amc_transform_matrix_identity(T, 0) ;
  for ( i = 1 ; i <= amc_transform_order_max(T) ; i ++ ) {
    amc_transform_matrix_zero(T, i) ;
  }

  if ( str == NULL ) {
    amc_transform_entry_set(T, 0, 0, 0,  cos(th), NULL) ;
    amc_transform_entry_set(T, 0, 0, 1, -sin(th), NULL) ;
    amc_transform_entry_set(T, 0, 1, 0,  sin(th), NULL) ;
    amc_transform_entry_set(T, 0, 1, 1,  cos(th), NULL) ;
  
    return 0 ;
  }

  sprintf(buf, "cos(%s)", str) ;
  amc_transform_entry_set(T, 0, 0, 0, 0, buf) ;
  amc_transform_entry_set(T, 0, 1, 1, 0, buf) ;
  sprintf(buf, "-sin(%s)", str) ;
  amc_transform_entry_set(T, 0, 0, 1, 0, buf) ;
  sprintf(buf, "sin(%s)", str) ;
  amc_transform_entry_set(T, 0, 1, 0, 0, buf) ;
    
  amc_transform_derivatives_evaluate(T, order) ;

  return 0 ;
}

/** 
 * Evaluate derivatives of a transform up to some specified order by
 * symbolic differentiation of expressions for transform matrix. This
 * requires that the libmatheval library be installed.
 * 
 * @param T an ::amc_transform_t;
 * @param order maximum order of derivative to evaluate.
 * 
 * @return 0 on success.
 */

int amc_transform_derivatives_evaluate(amc_transform_t *T, gint order)

{
  if ( order > amc_transform_order_max(T) ) {
    fprintf(stderr, "%s: order (%d) greater than maximum available (%d)\n",
	    __FUNCTION__, order, amc_transform_order_max(T)) ;
    return 1 ;
  }

#ifndef HAVE_LIBMATHEVAL
  fprintf(stderr, "%s: libmatheval symbolic differentiation not available\n",
	  __FUNCTION__) ;
  return 0 ;
#else /*HAVE_LIBMATHEVAL*/

  void *e0, *e1 ;
  char *estr ;
  gint i, j ;

  for ( i = 0 ; i < amc_transform_matrix_size(T) ; i ++ ) {
    if ( T->matrix_strings[i] != NULL ) {
      e0 = evaluator_create(T->matrix_strings[i]) ;
      if ( e0 == NULL ) {
	fprintf(stderr, "%s: cannot parse expression \"%s\"\n",
		__FUNCTION__, T->matrix_strings[i]) ;
	return 1 ;
      }
      for ( j = 1 ; j <= order ; j ++ ) {
	e1 = evaluator_derivative(e0, "t") ;
	estr = evaluator_get_string(e1) ;
	T->matrix_strings[j*amc_transform_matrix_size(T)+i] = strdup(estr) ;
	evaluator_destroy(e0) ;
	e0 = e1 ;
      }
    } else {
      for ( j = 1 ; j <= order ; j ++ ) {
	T->matrix_values[j*amc_transform_matrix_size(T)+i] = 0 ;
      }	
    }
  }
  
#endif /*HAVE_LIBMATHEVAL*/
  
  return 0 ;
}

/** 
 * Check symbolically evaluated derivatives of transform against
 * numerical estimate
 * 
 * @param T an ::amc_transform_t;
 * @param order order of derivative to check;
 * @param t evaluation time;
 * @param dt \f$\Delta t\f$ for finite difference estimate of derivative; 
 * @param imax on exit contains index of entry with largest difference
 * between numerical estimate and symbolically evaluated derivative;
 * @param emax on exit contains largest difference
 * between numerical estimate and symbolically evaluated derivative;
 * 
 * @return 0 on success
 */

int amc_transform_derivative_check(amc_transform_t *T, gint order,
				   gdouble t, gdouble dt,
				   gint *imax, gdouble *emax)
{
  gint i, n ;
  gdouble df, d ;
    
  if ( order < 1 ) {
    fprintf(stderr, "%s: cannot evaluate derivative of order %d\n",
	    __FUNCTION__, order) ;
  }

  if ( order > amc_transform_order_max(T) ) {
    fprintf(stderr, "%s: order (%d) greater than maximum available (%d)\n",
	    __FUNCTION__, order, amc_transform_order_max(T)) ;
    return -1 ;
  }

  *emax = 0 ; *imax = 0 ;
  n = amc_transform_matrix_size(T) ;

  for ( i = 0 ; i < n ; i ++ ) {
    if ( T->matrix_expr[order*n+i] != NULL ) {
      /*numerical estimate of derivative using lower order matrix*/
      T->t = t + 0.5*dt ;
      df = te_eval(T->matrix_expr[(order-1)*n+i]) ;
      T->t = t - 0.5*dt ;
      df = (df - te_eval(T->matrix_expr[(order-1)*n+i]))/dt ;
      T->t = t ;
      d = df - te_eval(T->matrix_expr[order*n+i]) ;
      if ( fabs(d) > *emax ) {
	*emax = fabs(d) ; *imax = i ;
      }
    }
  }
  
  return 0 ;
}

char *amc_sin_derivative(int i)

{
  char *diff[] = {"sin", "cos", "-sin", "-cos"} ;
  if ( i < 0 ) {
    fprintf(stderr, "%s: negative derivative orders (%d) not implemented\n",
	    __FUNCTION__, i) ;
    return NULL ;
  }
  
  return (diff[i % 4]) ;
}

char *amc_cos_derivative(int i)

{
  char *diff[] = {"cos", "-sin", "-cos", "sin"} ;
  if ( i < 0 ) {
    fprintf(stderr, "%s: negative derivative orders (%d) not implemented\n",
	    __FUNCTION__, i) ;
    return NULL ;
  }
  
  return (diff[i % 4]) ;
}

/** 
 * Convenience function to evaluate matrices in a chain
 * 
 * @param C ::amc_transform_chain_t containing transforms;
 * @param t evaluation time;
 * @param T if not NULL, is set to transform found by applying whole chain;
 * @param order order of derivatives to generate in \a T.
 * 
 * @return 0 on success
 */

int amc_transform_chain_matrices_evaluate(amc_transform_chain_t *C, gdouble t,
					  amc_transform_t *T, gint order)

{
  gint i ;

  for ( i = 0 ; i < amc_transform_chain_transform_number(C) ; i ++ ) {
    amc_transform_matrices_evaluate(amc_transform_chain_transform(C,i), t) ;
  }

  if ( T == NULL ) return 0 ;
  
  amc_transform_chain_evaluate(C, order, T) ;
  
  return 0 ;
}

amc_transform_definition_t amc_transform_definition_parse(char *str)

{
  amc_transform_definition_t def[] = {
    AMC_TRANSFORM_DEFINITION_UNKNOWN,
    AMC_TRANSFORM_DEFINITION_MATRIX,
    AMC_TRANSFORM_DEFINITION_ROTATION_X,
    AMC_TRANSFORM_DEFINITION_ROTATION_Y,
    AMC_TRANSFORM_DEFINITION_ROTATION_Z,
    AMC_TRANSFORM_DEFINITION_UNKNOWN    
  } ;
  char *f[] = {
    "unknown",
    "matrix",
    "rotation_x",
    "rotation_y",
    "rotation_z",
    NULL} ;
  gint i ;

  for ( i = 0 ; f[i] != NULL ; i ++ ) {
    if ( strcmp(str, f[i]) == 0 ) return def[i] ;
  }
  
  return AMC_TRANSFORM_DEFINITION_UNKNOWN ;
}

gint amc_transform_parse(amc_transform_t *T, char *str, gdouble val,
			 gint order)

{
  switch ( amc_transform_definition(T) ) {
  default: g_assert_not_reached() ; break ;
  case AMC_TRANSFORM_DEFINITION_ROTATION_X:
    amc_transform_rotation_x(T, val, str, order) ;
    break ;
  case AMC_TRANSFORM_DEFINITION_ROTATION_Y:
    amc_transform_rotation_y(T, val, str, order) ;
    break ;
  case AMC_TRANSFORM_DEFINITION_ROTATION_Z:
    amc_transform_rotation_z(T, val, str, order) ;
    break ;
  }
  
  return 0 ;
}
