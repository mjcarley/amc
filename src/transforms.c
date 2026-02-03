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

/**
 *
 * @{
 * 
 * @ingroup transform
 */

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
  case AMC_TRANSFORM_DEFINITION_ROTATION:
    amc_transform_rotation(T, val, str, order) ;
    break ;
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

/**
 *
 * @}
 * 
 */
