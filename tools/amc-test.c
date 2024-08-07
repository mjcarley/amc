#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "tinyexpr.h"

#include "amc.h"

int test_amc2d(void)

{
  amc_transform_t *R, *T, *N, *S ;
  amc_transform_chain_t *C ;
  int order, dim, ntrans ;
  double Om, U, V, x[2], y[2], t ;
  
  order = 3 ; dim = 2 ; ntrans = 8 ;
  R = amc_transform_alloc(dim, order) ; R->order = 1 ;
  S = amc_transform_alloc(dim, order) ; S->order = 1 ;
  T = amc_transform_alloc(dim, order) ; T->order = 1 ;
  N = amc_transform_alloc(dim, order) ; N->order = 1 ;

  C = amc_transform_chain_alloc(ntrans) ;
  
  Om = 2*M_PI*64 ; U = 1.4 ; V = 0.3 ;
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
  
  t = 1.3 ;
  amc_transform_matrices_evaluate(R, t) ;
  amc_transform_matrices_evaluate(T, t) ;
  amc_transform_matrices_evaluate(N, t) ;

  amc_transform_chain_transform_add(C, R) ;
  amc_transform_chain_transform_add(C, T) ;

  amc_transform_chain_evaluate(C, 1, S) ;
  
  amc_transform_variables_write(stderr, R) ;
  
  x[0] = 1.3 ; x[1] = 0.7 ;
  order = 1 ;
  amc_transform_matrix_apply(R, order, x, y) ;
  fprintf(stderr, "R: %lg %lg -> %lg %lg\n", x[0], x[1], y[0], y[1]) ;

  amc_transform_matrix_apply(T, order, y, y) ;
  fprintf(stderr, "T*R:        -> %lg %lg\n", y[0], y[1]) ;

  amc_transform_matrix_apply(N, order, x, y) ;
  fprintf(stderr, "N: %lg %lg -> %lg %lg\n", x[0], x[1], y[0], y[1]) ;

  amc_transform_matrix_apply(S, order, x, y) ;
  fprintf(stderr, "S: %lg %lg -> %lg %lg\n", x[0], x[1], y[0], y[1]) ;
  
  amc_transform_matrix_write(stderr, N, order) ;
  amc_transform_matrix_write(stderr, S, order) ;

  return 0 ;
}

int test_amc3d(void)

{
  amc_transform_t *R, *T, *N, *S ;
  amc_transform_chain_t *C ;
  int order, dim, ntrans ;
  double Om, U[3], x[3], y[3], t ;
  
  order = 3 ; dim = 3 ; ntrans = 8 ;
  R = amc_transform_alloc(dim, order) ; R->order = 1 ;
  S = amc_transform_alloc(dim, order) ; S->order = 1 ;
  T = amc_transform_alloc(dim, order) ; T->order = 1 ;
  N = amc_transform_alloc(dim, order) ; N->order = 1 ;

  C = amc_transform_chain_alloc(ntrans) ;
  
  Om = 2*M_PI*64 ; U[0] = 1.4 ; U[1] = 0 ; U[2] = 0 ;
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
  
  t = 1.3 ;
  amc_transform_matrices_evaluate(R, t) ;
  amc_transform_matrices_evaluate(T, t) ;
  amc_transform_matrices_evaluate(N, t) ;

  amc_transform_chain_transform_add(C, R) ;
  amc_transform_chain_transform_add(C, T) ;

  amc_transform_chain_evaluate(C, 1, S) ;
  
  amc_transform_variables_write(stderr, R) ;
  
  x[0] = 1.3 ; x[1] = 0.7 ; x[2] = -0.8 ;
  order = 1 ;
  amc_transform_matrix_apply(R, order, x, y) ;
  fprintf(stderr, "  R: %lg %lg %lg -> %lg %lg %lg\n",
	  x[0], x[1], x[2], y[0], y[1], y[2]) ;

  amc_transform_matrix_apply(T, order, y, y) ;
  fprintf(stderr, "T*R:        -> %lg %lg %lg\n", y[0], y[1], y[2]) ;

  amc_transform_matrix_apply(N, order, x, y) ;
  fprintf(stderr, "N: %lg %lg %lg -> %lg %lg %lg\n",
	  x[0], x[1], x[2], y[0], y[1], y[2]) ;

  amc_transform_matrix_apply(S, order, x, y) ;
  fprintf(stderr, "S: %lg %lg %lg -> %lg %lg %lg\n",
	  x[0], x[1], x[2], y[0], y[1], y[2]) ;
  
  amc_transform_matrix_write(stderr, N, order) ;
  amc_transform_matrix_write(stderr, S, order) ;

  return 0 ;
}

int rotor_test_3d(void)

{
  amc_transform_t *R, *T, *N, *S, *Y ;
  amc_transform_chain_t *C ;
  int order, dim, ntrans ;
  double Om, turn, r, x[3], y[3], t ;

  /*
   * rotor blade turning at speed Omega on an aircraft making a turn
   * of radius r and turn rate turn rad/s
   */
  
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
    fprintf(stdout, "%1.16e %1.16e %1.16e ", y[0], y[1], y[2]) ;
    amc_transform_chain_derivative(C, 0, S, t, x, y) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e\n", y[0], y[1], y[2]) ;
  }
  
  return 0 ;
}

int transition_test_3d(void)

{
  amc_transform_t *R, *T, *N, *S, *Y ;
  amc_transform_chain_t *C ;
  int order, dim, ntrans ;
  double Om, turn, r, x[3], y[3], t ;

  /*
   * rotor blade turning at speed Omega making transition to hover
   */
  
  dim = 3 ; order = 1 ; ntrans = 8 ;
  
  R = amc_transform_alloc(dim, order) ; R->order = 1 ;
  T = amc_transform_alloc(dim, order) ; T->order = 1 ;
  N = amc_transform_alloc(dim, order) ; N->order = 1 ;
  Y = amc_transform_alloc(dim, order) ; Y->order = 1 ;
  S = amc_transform_alloc(dim, order) ; S->order = 1 ;

  C = amc_transform_chain_alloc(ntrans) ;
  
  Om = 2*M_PI*5000/60 ;
  r = 0.5 ;
  turn = 2*M_PI*0.25 ;
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
  
  for ( t = 0 ; t < 1 ; t += 0.001 ) {
    amc_transform_matrices_evaluate(R, t) ;
    amc_transform_matrices_evaluate(N, t) ;
    amc_transform_matrices_evaluate(T, t) ;

    fprintf(stdout, "%1.16e ", t) ;
    
    amc_transform_chain_evaluate(C, 1, S) ;
    amc_transform_matrix_apply(S, 0, x, y) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e ", y[0], y[1], y[2]) ;
    amc_transform_matrix_apply(S, 1, x, y) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e ", y[0], y[1], y[2]) ;
    amc_transform_chain_derivative(C, 0, S, t, x, y) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e\n", y[0], y[1], y[2]) ;
  }
  
  return 0 ;
}

int derivative_test(void)

{
  amc_transform_t *T ;
  int order, dim, i, imax ;
  double Om, tol, t, dt, emax ;

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
  int i ;

  for ( i = 0 ; i <= 6 ; i ++ ) {
    fprintf(stderr, "%d: %s(Om*t)*(Om)^%d %s(Om*t)*(Om)^%d\n",
	    i, amc_sin_derivative(i), i, amc_cos_derivative(i), i) ;
  }
  
  return ;
}

int main(int argc, char **argv)

{
  /* test_amc2d() ; */
  /* test_amc3d() ; */

  /* rotor_test_3d() ; */
  /* transition_test_3d() ; */

  /* derivative_test() ; */

  sin_cos_test() ;
  
  return 0 ;
}
