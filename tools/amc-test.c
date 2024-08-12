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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <glib.h>

#include "amc.h"

char *progname ;
char *tests[] = {"2d", "3d", "rotor", "transition", "derivative", "sincos",
  "vector_3d", "frame_3d", "frame_read",
  NULL} ;

static gint parse_test(char *test)

{
  gint i ;  

  for ( i = 0 ; tests[i] != NULL ; i ++ ) {
    if ( strcmp(tests[i], test) == 0 ) return i+1 ;
  }
  
  return -1 ;
}

static void list_tests(FILE *f)

{
  gint i ;  

  for ( i = 0 ; tests[i] != NULL ; i ++ ) {
    fprintf(f, "%s\n", tests[i]) ;
  }
  
  return ;
}

static gint test_amc2d(void)

{
  amc_transform_t *R, *T, *N, *S ;
  amc_transform_chain_t *C ;
  gint order, dim, ntrans ;
  gdouble Om, U, V, x[2], y[2], t ;

  fprintf(stderr, "two-dimensional motion test\n") ;
  fprintf(stderr, "===========================\n") ;
  
  order = 3 ; dim = 2 ; ntrans = 8 ;
  R = amc_transform_alloc(dim, order) ; R->order = 1 ;
  S = amc_transform_alloc(dim, order) ; S->order = 1 ;
  T = amc_transform_alloc(dim, order) ; T->order = 1 ;
  N = amc_transform_alloc(dim, order) ; N->order = 1 ;

  C = amc_transform_chain_alloc(ntrans) ;
  
  Om = 2*M_PI*64 ; U = 1.4 ; V = 0.3 ; t = 1.3 ;
  x[0] = 1.3 ; x[1] = 0.7 ;

  fprintf(stderr, "variable: Omega = %lg\n", Om) ;
  fprintf(stderr, "variable:     U = %lg\n", U) ;
  fprintf(stderr, "variable:     V = %lg\n", V) ;
  fprintf(stderr, "variable:     t = %lg\n", t) ;

  fprintf(stderr, "transform R(x) rotation about origin at velocity Omega\n") ;
  fprintf(stderr, "transform T(x) translation at velocity (U,V)\n") ;
  fprintf(stderr, "transform N(x) combined rotation and translation\n") ;
  fprintf(stderr, "transform S(x) = T*R(x), by analytical evaluation\n") ;
  
  amc_transform_variable_add(R, "Omega", &Om) ;
  amc_transform_variable_add(R, "U", &U) ;
  amc_transform_variable_add(R, "V", &V) ;
  amc_transform_variable_add(T, "U", &U) ;
  amc_transform_variable_add(T, "V", &V) ;
  amc_transform_variable_add(N, "Omega", &Om) ;
  amc_transform_variable_add(N, "U", &U) ;
  amc_transform_variable_add(N, "V", &V) ;

  /*rotation matrix*/
  amc_transform_entry_set(R, 0, 0, 0, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(R, 0, 0, 1, 0, "-sin(Omega*t)") ;
  amc_transform_entry_set(R, 0, 1, 0, 0, " sin(Omega*t)") ;
  amc_transform_entry_set(R, 0, 1, 1, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(R, 0, 2, 2, 1, NULL) ;

  /*derivatives*/
  amc_transform_entry_set(R, 1, 0, 0, 0, "-Omega*sin(Omega*t)") ;
  amc_transform_entry_set(R, 1, 0, 1, 0, "-Omega*cos(Omega*t)") ;
  amc_transform_entry_set(R, 1, 1, 0, 0, " Omega*cos(Omega*t)") ;
  amc_transform_entry_set(R, 1, 1, 1, 0, "-Omega*sin(Omega*t)") ;
  amc_transform_entry_set(R, 1, 2, 2, 0, NULL) ;

  /*translation matrix*/
  amc_transform_entry_set(T, 0, 0, 0, 1, NULL) ;
  amc_transform_entry_set(T, 0, 1, 1, 1, NULL) ;
  amc_transform_entry_set(T, 0, 0, 2, 0, " U*t") ;
  amc_transform_entry_set(T, 0, 1, 2, 0, " V*t") ;
  amc_transform_entry_set(T, 0, 2, 2, 1, NULL) ;
  /*derivatives*/
  amc_transform_entry_set(T, 1, 0, 2, 0, " U") ;
  amc_transform_entry_set(T, 1, 1, 2, 0, " V") ;
  amc_transform_entry_set(T, 1, 2, 2, 0, NULL) ;

  /*net transformation, for checking*/
  amc_transform_entry_set(N, 0, 0, 0, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(N, 0, 0, 1, 0, "-sin(Omega*t)") ;
  amc_transform_entry_set(N, 0, 1, 0, 0, " sin(Omega*t)") ;
  amc_transform_entry_set(N, 0, 1, 1, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(N, 0, 0, 2, 0, " U*t") ;
  amc_transform_entry_set(N, 0, 1, 2, 0, " V*t") ;

  /*derivatives*/
  amc_transform_entry_set(N, 1, 0, 0, 0, "-Omega*sin(Omega*t)") ;
  amc_transform_entry_set(N, 1, 0, 1, 0, "-Omega*cos(Omega*t)") ;
  amc_transform_entry_set(N, 1, 1, 0, 0, " Omega*cos(Omega*t)") ;
  amc_transform_entry_set(N, 1, 1, 1, 0, "-Omega*sin(Omega*t)") ;
  amc_transform_entry_set(N, 1, 0, 2, 0, " U") ;
  amc_transform_entry_set(N, 1, 1, 2, 0, " V") ;

  amc_transform_expressions_compile(R) ;
  amc_transform_expressions_compile(T) ;
  amc_transform_expressions_compile(N) ;
  
  amc_transform_matrices_evaluate(R, t) ;
  amc_transform_matrices_evaluate(T, t) ;
  amc_transform_matrices_evaluate(N, t) ;

  amc_transform_chain_transform_add(C, R) ;
  amc_transform_chain_transform_add(C, T) ;

  amc_transform_chain_evaluate(C, 1, S) ;

  fprintf(stderr, "\nrotation matrix variables\n") ;
  amc_transform_variables_write(stderr, R) ;

  fprintf(stderr, "\n") ;
  fprintf(stderr, "point x = (%lg,%lg)\n", x[0], x[1]) ;

  order = 1 ;

  fprintf(stderr, "derivative order = %d\n", order) ;
    
  amc_transform_matrix_apply(R, order, x, y) ;
  fprintf(stderr, "R(x): %lg %lg -> %lg %lg\n", x[0], x[1], y[0], y[1]) ;

  amc_transform_matrix_apply(T, order, x, y) ;
  fprintf(stderr, "T(x): %lg %lg -> %lg %lg\n", x[0], x[1], y[0], y[1]) ;

  amc_transform_matrix_apply(N, order, x, y) ;
  fprintf(stderr, "N: %lg %lg -> %lg %lg\n", x[0], x[1], y[0], y[1]) ;

  amc_transform_matrix_apply(S, order, x, y) ;
  fprintf(stderr, "S: %lg %lg -> %lg %lg\n", x[0], x[1], y[0], y[1]) ;

  fprintf(stderr, "\nanalytical matrix N = \n") ;
  amc_transform_matrix_write(stderr, N, order) ;
  fprintf(stderr, "\nnumerical matrix S = \n") ;
  amc_transform_matrix_write(stderr, S, order) ;

  return 0 ;
}

static gint test_amc3d(void)

{
  amc_transform_t *R, *T, *N, *S ;
  amc_transform_chain_t *C ;
  gint order, dim, ntrans ;
  gdouble Om, U[3], x[3], y[3], t ;
  
  fprintf(stderr, "three-dimensional motion test\n") ;
  fprintf(stderr, "=============================\n") ;

  order = 3 ; dim = 3 ; ntrans = 8 ;
  R = amc_transform_alloc(dim, order) ; R->order = 1 ;
  S = amc_transform_alloc(dim, order) ; S->order = 1 ;
  T = amc_transform_alloc(dim, order) ; T->order = 1 ;
  N = amc_transform_alloc(dim, order) ; N->order = 1 ;

  C = amc_transform_chain_alloc(ntrans) ;
  
  Om = 2*M_PI*64 ; U[0] = 1.4 ; U[1] = 0 ; U[2] = 0 ;
  t = 1.3 ;
  
  fprintf(stderr, "variable: Omega = %lg\n", Om) ;
  fprintf(stderr, "variable:     U = %lg\n", U[0]) ;
  fprintf(stderr, "variable:     V = %lg\n", U[1]) ;
  fprintf(stderr, "variable:     W = %lg\n", U[2]) ;
  fprintf(stderr, "variable:     t = %lg\n", t) ;
  
  amc_transform_variable_add(R, "Omega", &Om) ;
  amc_transform_variable_add(R, "U", &(U[0])) ;
  amc_transform_variable_add(R, "V", &(U[1])) ;
  amc_transform_variable_add(R, "W", &(U[2])) ;

  amc_transform_variable_add(T, "U", &(U[0])) ;
  amc_transform_variable_add(T, "V", &(U[1])) ;
  amc_transform_variable_add(T, "W", &(U[2])) ;

  amc_transform_variable_add(N, "Omega", &Om) ;
  amc_transform_variable_add(N, "U", &(U[0])) ;
  amc_transform_variable_add(N, "V", &(U[1])) ;
  amc_transform_variable_add(N, "W", &(U[2])) ;

  fprintf(stderr, "transform R(x) rotation about z axis at velocity Omega\n") ;
  fprintf(stderr, "transform T(x) translation at velocity (U,V,W)\n") ;
  fprintf(stderr, "transform N(x) combined rotation and translation\n") ;
  fprintf(stderr, "transform S(x) = T*R(x), by analytical evaluation\n") ;

  /*rotation matrix*/
  amc_transform_entry_set(R, 0, 0, 0, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(R, 0, 0, 1, 0, "-sin(Omega*t)") ;
  amc_transform_entry_set(R, 0, 1, 0, 0, " sin(Omega*t)") ;
  amc_transform_entry_set(R, 0, 1, 1, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(R, 0, 2, 2, 1, NULL) ;
  amc_transform_entry_set(R, 0, 3, 3, 1, NULL) ;

  /*derivatives*/
  amc_transform_entry_set(R, 1, 0, 0, 0, "-Omega*sin(Omega*t)") ;
  amc_transform_entry_set(R, 1, 0, 1, 0, "-Omega*cos(Omega*t)") ;
  amc_transform_entry_set(R, 1, 1, 0, 0, " Omega*cos(Omega*t)") ;
  amc_transform_entry_set(R, 1, 1, 1, 0, "-Omega*sin(Omega*t)") ;
  amc_transform_entry_set(R, 1, 2, 2, 0, NULL) ;
  amc_transform_entry_set(R, 1, 3, 3, 0, NULL) ;

  /*translation matrix*/
  amc_transform_entry_set(T, 0, 0, 0, 1, NULL) ;
  amc_transform_entry_set(T, 0, 1, 1, 1, NULL) ;
  amc_transform_entry_set(T, 0, 0, 3, 0, " U*t") ;
  amc_transform_entry_set(T, 0, 1, 3, 0, " V*t") ;
  amc_transform_entry_set(T, 0, 2, 3, 0, " W*t") ;
  amc_transform_entry_set(T, 0, 2, 2, 1, NULL) ;
  amc_transform_entry_set(T, 0, 3, 3, 1, NULL) ;
  /*derivatives*/
  amc_transform_entry_set(T, 1, 0, 3, 0, " U") ;
  amc_transform_entry_set(T, 1, 1, 3, 0, " V") ;
  amc_transform_entry_set(T, 1, 2, 3, 0, " W") ;
  amc_transform_entry_set(T, 1, 3, 3, 0, NULL) ;

  /*net transformation, for checking*/
  amc_transform_entry_set(N, 0, 0, 0, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(N, 0, 0, 1, 0, "-sin(Omega*t)") ;
  amc_transform_entry_set(N, 0, 1, 0, 0, " sin(Omega*t)") ;
  amc_transform_entry_set(N, 0, 1, 1, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(N, 0, 2, 2, 1, NULL) ;
  amc_transform_entry_set(N, 0, 3, 3, 1, NULL) ;
  amc_transform_entry_set(N, 0, 0, 3, 0, " U*t") ;
  amc_transform_entry_set(N, 0, 1, 3, 0, " V*t") ;
  amc_transform_entry_set(N, 0, 2, 3, 0, " W*t") ;

  /*derivatives*/
  amc_transform_entry_set(N, 1, 0, 0, 0, "-Omega*sin(Omega*t)") ;
  amc_transform_entry_set(N, 1, 0, 1, 0, "-Omega*cos(Omega*t)") ;
  amc_transform_entry_set(N, 1, 1, 0, 0, " Omega*cos(Omega*t)") ;
  amc_transform_entry_set(N, 1, 1, 1, 0, "-Omega*sin(Omega*t)") ;
  amc_transform_entry_set(N, 1, 0, 3, 0, " U") ;
  amc_transform_entry_set(N, 1, 1, 3, 0, " V") ;
  amc_transform_entry_set(N, 1, 2, 3, 0, " W") ;
  amc_transform_entry_set(N, 1, 3, 3, 0, NULL) ;

  amc_transform_expressions_compile(R) ;
  amc_transform_expressions_compile(T) ;
  amc_transform_expressions_compile(N) ;
  
  amc_transform_matrices_evaluate(R, t) ;
  amc_transform_matrices_evaluate(T, t) ;
  amc_transform_matrices_evaluate(N, t) ;

  amc_transform_chain_transform_add(C, R) ;
  amc_transform_chain_transform_add(C, T) ;

  amc_transform_chain_evaluate(C, 1, S) ;
  
  fprintf(stderr, "\nrotation matrix variables\n") ;
  amc_transform_variables_write(stderr, R) ;
  
  x[0] = 1.3 ; x[1] = 0.7 ; x[2] = -0.8 ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "point x = (%lg,%lg,%lg)\n", x[0], x[1], x[2]) ;

  order = 1 ;
  fprintf(stderr, "derivative order = %d\n", order) ;
  amc_transform_matrix_apply(R, order, x, y) ;
  fprintf(stderr, "  R(x): %lg %lg %lg -> %lg %lg %lg\n",
	  x[0], x[1], x[2], y[0], y[1], y[2]) ;

  amc_transform_matrix_apply(T, order, x, y) ;
  fprintf(stderr, "  T(x): %lg %lg %lg -> %lg %lg %lg\n",
	  x[0], x[1], x[2], y[0], y[1], y[2]) ;
  /* fprintf(stderr, "T*R:        -> %lg %lg %lg\n", y[0], y[1], y[2]) ; */

  amc_transform_matrix_apply(N, order, x, y) ;
  fprintf(stderr, "  N(x): %lg %lg %lg -> %lg %lg %lg\n",
	  x[0], x[1], x[2], y[0], y[1], y[2]) ;

  amc_transform_matrix_apply(S, order, x, y) ;
  fprintf(stderr, "  S(x): %lg %lg %lg -> %lg %lg %lg\n",
	  x[0], x[1], x[2], y[0], y[1], y[2]) ;
  
  fprintf(stderr, "\nanalytical matrix N = \n") ;
  amc_transform_matrix_write(stderr, N, order) ;
  fprintf(stderr, "\nnumerical matrix S = \n") ;
  amc_transform_matrix_write(stderr, S, order) ;

  return 0 ;
}

static gint rotor_test_3d(void)

{
  amc_transform_t *R, *T, *N, *S, *Y ;
  amc_transform_chain_t *C ;
  gint order, dim, ntrans ;
  gdouble Om, turn, r, x[3], y[3], t ;

  /*
   * rotor blade turning at speed Omega on an aircraft making a turn
   * of radius r and turn rate turn rad/s
   */
  
  fprintf(stderr, "rotating point source on rotating aircraft\n") ;
  fprintf(stderr, "==========================================\n") ;

  dim = 3 ; order = 1 ; ntrans = 8 ;
  
  R = amc_transform_alloc(dim, order) ; R->order = 1 ;
  T = amc_transform_alloc(dim, order) ; T->order = 1 ;
  N = amc_transform_alloc(dim, order) ; N->order = 1 ;
  Y = amc_transform_alloc(dim, order) ; Y->order = 1 ;
  S = amc_transform_alloc(dim, order) ; S->order = 1 ;

  C = amc_transform_chain_alloc(ntrans) ;
  
  Om = 2*M_PI*1000/60 ;
  r = 50 ;
  turn = 2*M_PI*5/60 ;

  fprintf(stderr, "source angular velocity:  Omega = %lg\n", Om) ;
  fprintf(stderr, "radius of aircraft rotation:  r = %lg\n", r) ;
  fprintf(stderr, "aircraft angular velocity: turn = %lg\n", turn) ;
	  
  amc_transform_variable_add(R, "Omega", &Om) ;

  /*rotation about x axis*/
  amc_transform_rotation_x(R, 0, "Omega*t", order) ;
  amc_transform_derivatives_evaluate(R, 1) ;

  /*"nacelle" matrix (shift to position on "wing")*/
  amc_transform_variable_add(N, "r", &r) ;
  amc_transform_translation(N, 0, NULL, 0, "r", 0, NULL, order) ;

  /*"roll" matrix*/
  amc_transform_rotation_x(Y, 10.0*M_PI/180, NULL, order) ;
  
  /*turn matrix*/
  amc_transform_variable_add(T, "turn", &turn) ;
  amc_transform_rotation_z(T, 0, "turn*t", order) ;
  amc_transform_derivatives_evaluate(T, 1) ;

  amc_transform_chain_transform_add(C, R) ;
  amc_transform_chain_transform_add(C, N) ;
  amc_transform_chain_transform_add(C, Y) ;
  amc_transform_chain_transform_add(C, T) ;

  amc_transform_expressions_compile(R) ;
  amc_transform_expressions_compile(N) ;
  amc_transform_expressions_compile(Y) ;
  amc_transform_expressions_compile(T) ;

  x[0] = 0.1 ; x[1] = 0.5 ; x[2] = 0.3 ;

  fprintf(stderr, "writing trajectory to stdout\n") ;
  fprintf(stderr, "  t x y z u v w\n") ;
  for ( t = 0 ; t < 1 ; t += 0.001 ) {
    amc_transform_matrices_evaluate(R, t) ;
    amc_transform_matrices_evaluate(N, t) ;
    amc_transform_matrices_evaluate(Y, t) ;
    amc_transform_matrices_evaluate(T, t) ;

    fprintf(stdout, "%1.16e ", t) ;
    
    amc_transform_chain_evaluate(C, 1, S) ;
    amc_transform_matrix_apply(S, 0, x, y) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e ", y[0], y[1], y[2]) ;
    amc_transform_matrix_apply(S, 1, x, y) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e\n", y[0], y[1], y[2]) ;
    /* amc_transform_chain_derivative(C, 0, S, t, x, y) ; */
    /* fprintf(stdout, "%1.16e %1.16e %1.16e\n", y[0], y[1], y[2]) ; */
  }
  
  return 0 ;
}

static gint transition_test_3d(void)

{
  amc_transform_t *R, *T, *N, *S, *Y ;
  amc_transform_chain_t *C ;
  gint order, dim, ntrans ;
  gdouble Om, turn, r, x[3], y[3], t ;

  /*
   * rotor blade turning at speed Omega making transition to hover
   */
  
  fprintf(stderr, "rotating point source in transition to hover\n") ;
  fprintf(stderr, "============================================\n") ;

  dim = 3 ; order = 1 ; ntrans = 8 ;
  
  R = amc_transform_alloc(dim, order) ; R->order = 1 ;
  T = amc_transform_alloc(dim, order) ; T->order = 1 ;
  N = amc_transform_alloc(dim, order) ; N->order = 1 ;
  Y = amc_transform_alloc(dim, order) ; Y->order = 1 ;
  S = amc_transform_alloc(dim, order) ; S->order = 1 ;

  C = amc_transform_chain_alloc(ntrans) ;
  
  Om = 2*M_PI*500/60 ;
  r = 0.5 ;
  turn = 2*M_PI*0.25 ;

  fprintf(stderr, "source angular velocity:    Omega = %lg\n", Om) ;
  fprintf(stderr, "radius of rotation:             r = %lg\n", r) ;
  fprintf(stderr, "transition angular velocity: turn = %lg\n", turn) ;

  amc_transform_variable_add(R, "Omega", &Om) ;

  /*rotation about x axis*/
  amc_transform_rotation_x(R, 0, "Omega*t", order) ;
  /* amc_transform_derivatives_evaluate(R, 1) ; */

  /*"nacelle" matrix (shift to position on "nacelle")*/
  amc_transform_variable_add(N, "r", &r) ;
  amc_transform_translation(N, 0, "r", 0, NULL, 0, NULL, order) ;

  /*turn matrix*/
  amc_transform_variable_add(T, "turn", &turn) ;
  amc_transform_rotation_y(T, 0, "turn*t", order) ;
  /* amc_transform_derivatives_evaluate(T, 1) ; */

  amc_transform_chain_transform_add(C, R) ;
  amc_transform_chain_transform_add(C, N) ;
  amc_transform_chain_transform_add(C, T) ;

  amc_transform_expressions_compile(R) ;
  amc_transform_expressions_compile(N) ;
  amc_transform_expressions_compile(T) ;

  x[0] = 0.0 ; x[1] = 0.25 ; x[2] = 0.3 ;
  
  fprintf(stderr, "writing trajectory to stdout\n") ;
  fprintf(stderr, "  t x y z u v w\n") ;
  for ( t = 0 ; t < 1 ; t += 0.001 ) {
    amc_transform_matrices_evaluate(R, t) ;
    amc_transform_matrices_evaluate(N, t) ;
    amc_transform_matrices_evaluate(T, t) ;

    fprintf(stdout, "%1.16e ", t) ;
    
    amc_transform_chain_evaluate(C, 1, S) ;
    amc_transform_matrix_apply(S, 0, x, y) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e ", y[0], y[1], y[2]) ;
    amc_transform_matrix_apply(S, 1, x, y) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e\n", y[0], y[1], y[2]) ;
    /* amc_transform_chain_derivative(C, 0, S, t, x, y) ; */
    /* fprintf(stdout, "%1.16e %1.16e %1.16e\n", y[0], y[1], y[2]) ; */
  }
  
  return 0 ;
}

static gint derivative_test(void)

{
  amc_transform_t *T ;
  gint order, dim, i, imax ;
  gdouble Om, tol, t, dt, emax ;

  /*
   * rotor blade turning at speed Omega making transition to hover
   */
  
  dim = 3 ; order = 4 ; Om = 4 ;
  
  T = amc_transform_alloc(dim, order) ; T->order = order ;
  amc_transform_variable_add(T, "Omega", &Om) ;
  amc_transform_rotation_x(T, 0, "Omega*t", order) ;
  amc_transform_derivatives_evaluate(T, order) ;
  amc_transform_expressions_compile(T) ;

  t = 1.5 ; dt = 1e-6 ; tol = 1e-7 ;

  for ( i = 1 ; i <= order ; i ++ ) {
    amc_transform_derivative_check(T, i, t, dt, &imax, &emax) ;
    if ( emax > tol ) {
      fprintf(stderr, "FAIL: order %d, entry %d, error %lg\n", i, imax, emax) ;
    } else {
      fprintf(stderr, "PASS: order %d, maximum error %lg, entry %d\n",
	      i, emax, imax) ;
    }
  }
  
  /* amc_transform_variables_write(stderr, T) ; */
  /* amc_transform_matrix_write(stderr, T, 0) ; */
  /* amc_transform_matrix_write(stderr, T, 1) ; */
  
  return 0 ;
}

static void sin_cos_test(void)

{
  gint i ;

  for ( i = 0 ; i <= 6 ; i ++ ) {
    fprintf(stderr, "%d: %s(Om*t)*(Om)^%d %s(Om*t)*(Om)^%d\n",
	    i, amc_sin_derivative(i), i, amc_cos_derivative(i), i) ;
  }
  
  return ;
}

static gint vector_test_3d(void)

{
  amc_transform_t *R, *T, *S ;
  amc_transform_chain_t *C ;
  gint order, dim, ntrans ;
  gdouble Om, U[3], x0[3], x1[3], y0[3], y1[3], n[3], t ;
  
  fprintf(stderr, "vector transformation test\n") ;
  fprintf(stderr, "==========================\n") ;

  order = 3 ; dim = 3 ; ntrans = 8 ;
  R = amc_transform_alloc(dim, order) ; R->order = order ;
  S = amc_transform_alloc(dim, order) ; S->order = order ;
  T = amc_transform_alloc(dim, order) ; T->order = order ;

  C = amc_transform_chain_alloc(ntrans) ;
  
  Om = 2*M_PI*64 ; U[0] = 1.4 ; U[1] = 0 ; U[2] = 0 ;
  t = 1.3 ;
  
  fprintf(stderr, "variable: Omega = %lg\n", Om) ;
  fprintf(stderr, "variable:     U = %lg\n", U[0]) ;
  fprintf(stderr, "variable:     V = %lg\n", U[1]) ;
  fprintf(stderr, "variable:     W = %lg\n", U[2]) ;
  fprintf(stderr, "variable:     t = %lg\n", t) ;
  
  amc_transform_variable_add(R, "Omega", &Om) ;
  amc_transform_variable_add(R, "U", &(U[0])) ;
  amc_transform_variable_add(R, "V", &(U[1])) ;
  amc_transform_variable_add(R, "W", &(U[2])) ;

  amc_transform_variable_add(T, "U", &(U[0])) ;
  amc_transform_variable_add(T, "V", &(U[1])) ;
  amc_transform_variable_add(T, "W", &(U[2])) ;

  fprintf(stderr, "transform R(x) rotation about z axis at velocity Omega\n") ;
  fprintf(stderr, "transform T(x) translation at velocity (U,V,W)\n") ;
  fprintf(stderr, "transform S(x) = T*R(x), by analytical evaluation\n") ;

  /*rotation matrix*/
  amc_transform_entry_set(R, 0, 0, 0, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(R, 0, 0, 1, 0, "-sin(Omega*t)") ;
  amc_transform_entry_set(R, 0, 1, 0, 0, " sin(Omega*t)") ;
  amc_transform_entry_set(R, 0, 1, 1, 0, " cos(Omega*t)") ;
  amc_transform_entry_set(R, 0, 2, 2, 1, NULL) ;
  amc_transform_entry_set(R, 0, 3, 3, 1, NULL) ;
  amc_transform_derivatives_evaluate(R, order) ;

  /*translation matrix*/
  amc_transform_entry_set(T, 0, 0, 0, 1, NULL) ;
  amc_transform_entry_set(T, 0, 1, 1, 1, NULL) ;
  amc_transform_entry_set(T, 0, 0, 3, 0, " U*t") ;
  amc_transform_entry_set(T, 0, 1, 3, 0, " V*t") ;
  amc_transform_entry_set(T, 0, 2, 3, 0, " W*t") ;
  amc_transform_entry_set(T, 0, 2, 2, 1, NULL) ;
  amc_transform_entry_set(T, 0, 3, 3, 1, NULL) ;
  amc_transform_derivatives_evaluate(T, order) ;

  amc_transform_expressions_compile(R) ;
  amc_transform_expressions_compile(T) ;
  
  amc_transform_matrices_evaluate(R, t) ;
  amc_transform_matrices_evaluate(T, t) ;

  amc_transform_chain_transform_add(C, R) ;
  amc_transform_chain_transform_add(C, T) ;

  amc_transform_chain_evaluate(C, order, S) ;
  
  fprintf(stderr, "\nrotation matrix variables\n") ;
  amc_transform_variables_write(stderr, R) ;
  
  x0[0] = 1.3 ; x0[1] = 0.7 ; x0[2] = -0.8 ;
  n[0] = 1.7 ; n[1] = -0.4 ; n[2] = 1.9 ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "point  x = (%lg,%lg,%lg)\n", x0[0], x0[1], x0[2]) ;
  fprintf(stderr, "vector n = (%lg,%lg,%lg)\n", n[0], n[1], n[2]) ;

  x1[0] = x0[0] + n[0] ; x1[1] = x0[1] + n[1] ; x1[2] = x0[2] + n[2] ; 
  
  order = 1 ;
  fprintf(stderr, "derivative order = %d\n", order) ;

  amc_transform_matrix_apply(S, order, x0, y0) ;
  amc_transform_matrix_apply(S, order, x1, y1) ;
  amc_transform_matrix_apply_vector(S, order, n, n) ;

  fprintf(stderr, "%lg %lg\n", n[0], y1[0] - y0[0]) ;
  fprintf(stderr, "%lg %lg\n", n[1], y1[1] - y0[1]) ;
  fprintf(stderr, "%lg %lg\n", n[2], y1[2] - y0[2]) ;

  return 0 ;
}

static void frame_test_3d(void)

{
  amc_frame_t *F ;
  gint nc, nt, nv, i ;
  amc_transform_t *R, *T, *T1, *T2, *S ;
  amc_transform_chain_t *C1, *C2 ;
  gint order, dim, ntrans ;
  gdouble Om, U[3], x[3], y[3], t ;

  nt = 64 ; nc = 8 ; nv = 128 ;
  dim = 3 ; order = 3 ; ntrans = 8 ;
  
  fprintf(stderr, "three-dimensional frame motion test\n") ;
  fprintf(stderr, "===================================\n") ;

  F = amc_frame_alloc(nt, nc, nv) ;

  R = amc_transform_alloc(dim, order)  ; R->order = 1 ;
  T = amc_transform_alloc(dim, order)  ; T->order = 1 ;
  T1 = amc_transform_alloc(dim, order) ; T1->order = 1 ;
  T2 = amc_transform_alloc(dim, order) ; T2->order = 1 ;
  S = amc_transform_alloc(dim, order) ; S->order = order ;
  
  Om = 2*M_PI*1000/60 ;
  U[0] = 0.0 ; U[1] = 0 ; U[2] = 10 ;

  /*both sources rotate about z*/
  amc_transform_variable_add(R, "Omega", &Om) ;
  amc_transform_rotation_z(R, 0, "Omega*t", order) ;

  amc_transform_variable_add(T, "U", &(U[0])) ;
  amc_transform_variable_add(T, "V", &(U[1])) ;
  amc_transform_variable_add(T, "W", &(U[2])) ;
  amc_transform_translation(T, 0, "U*t", 0, "V*t", 0, "W*t", order) ;
  amc_transform_translation(T1, 0, NULL, -3, NULL, 0, NULL, order) ;
  amc_transform_translation(T2, 0, NULL,  3, NULL, 0, NULL, order) ;

  amc_frame_transform_add(F, R, "rotation") ;  
  amc_frame_transform_add(F, T, "translation") ;  
  amc_frame_transform_add(F, T1, "left") ;
  amc_frame_transform_add(F, T2, "right") ;

  C1 = amc_transform_chain_alloc(ntrans) ;
  C2 = amc_transform_chain_alloc(ntrans) ;

  /*left hand rotor*/
  amc_transform_chain_transform_add(C1, R ) ;
  amc_transform_chain_transform_add(C1, T1) ;
  amc_transform_chain_transform_add(C1, T ) ;

  /*right hand rotor*/
  amc_transform_chain_transform_add(C2, R ) ;
  amc_transform_chain_transform_add(C2, T2) ;
  amc_transform_chain_transform_add(C2, T ) ;

  amc_frame_transform_chain_add(F, C1, "left") ;
  amc_frame_transform_chain_add(F, C2, "right") ;  

  amc_frame_write(stderr, F) ;
  
  amc_frame_initialise(F) ;

  x[0] = 0.5 ; x[1] = 0.0 ; x[2] = 0.0 ;

  for ( t = 0 ; t < 1.0/8 ; t += 1.0/2048 ) {
    amc_frame_evaluate(F, t) ;
    for ( i = 0 ; i < amc_frame_transform_chain_number(F) ; i ++ ) {
      amc_transform_chain_evaluate(amc_frame_transform_chain(F,i), 0, S) ;
      amc_transform_matrix_apply(S, 0, x, y) ;
      fprintf(stdout, "%e %e %e ", y[0], y[1], y[2]) ;
    }
    fprintf(stdout, "\n") ;
  }
  
  return ;
}

static void frame_read_test(char *file)

{
  gint nc, nt, nv ;
  amc_frame_t *F ;

  nt = 64 ; nc = 8 ; nv = 128 ;

  F = amc_frame_alloc(nt, nc, nv) ;

  amc_frame_read(F, file) ;
  
  return ;
}

gint main(gint argc, char **argv)

{
  gint test ;
  char ch, *file ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  test = 0 ; file = NULL ;
  while ( (ch = getopt(argc, argv, "lf:t:")) != EOF ) {
    switch (ch ) {
    default: g_assert_not_reached() ; break ;
    case 'f': file = g_strdup(optarg) ; break ;
    case 'l':
      list_tests(stderr) ;
      return 0 ;
      break ;
    case 't':
      if ( (test = parse_test(optarg)) == -1 ) {
	fprintf(stderr, "%s: unrecognised test \"%s\"\n", progname, optarg) ;
	return 1 ;
      }
      break ;
    }
  }

  if ( test == 0 ) {
    fprintf(stderr, "%s: test must be specified (use option -t)\n",
	    progname) ;
    return 1 ;
  }
  
  if ( test == 1 ) {
    test_amc2d() ;

    return 0 ;
  }

  if ( test == 2 ) {
    test_amc3d() ;

    return 0 ;
  }

  if ( test == 3 ) {
    rotor_test_3d() ;

    return 0 ;
  }

  if ( test == 4 ) {
    transition_test_3d() ;

    return 0 ;
  }
    
  if ( test == 5 ) {
    derivative_test() ;

    return 0 ;
  }
    
  if ( test == 6 ) {
    sin_cos_test() ;

    return 0 ;
  }

  if ( test == 7 ) {
    vector_test_3d() ;

    return 0 ;
  }

  if ( test == 8 ) {
    frame_test_3d() ;

    return 0 ;
  }

  if ( test == 9 ) {
    if ( file == NULL ) {
      fprintf(stderr, "%s: this test requires that a filename be specified\n",
	      progname) ;
      return 1 ;
    }
    frame_read_test(file) ;

    return 0 ;
  }
  
  return 0 ;
}
