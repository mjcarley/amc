#ifndef __AMC_H_INCLUDED__
#define __AMC_H_INCLUDED__

#define AMC_TRANSFORM_VARIABLE_NUMBER 64

/**
 * @file   amc.h
 * @author  <michael@michael.paraffinalia.co.uk>
 * @date   Fri May 17 09:13:57 2024
 * 
 * @brief  
 * 
 * 
 */

#ifdef DOXYGEN
/** 
 * @typedef amc_transform_t
 * 
 * Basic type holding affine transform matrix data (two- or three-dimensional) 
 */
typedef amc_transform_t ;

/**
 * @brief Order of derivatives to be evaluated in \a T
 */
#define amc_transform_order(T)
/**
 * @brief Maximum order of derivatives which can be evaluated in \a T
 */
#define amc_transform_order_max(T)
/**
 * @brief Number of variables defined in \a T
 */
#define amc_transform_variable_number(T)
/**
 * @brief Dimension of transform \a T (two or three)
 */
#define amc_transform_dimension(T)      
/**
 * @brief Size of affine matrix in \a T (9 or 16)
 */
#define amc_transform_matrix_size(T)    
/**
 * @brief Pointer to \f$i\f$th transform matrix in \a T
 */
#define amc_transform_matrix(T,i)	
/**
 * @brief Value of entry \f$(i,j)\f$ in matrix of derivative order 
 * \f$o\f$ in \a T 
 */
#define amc_transform_value(T,o,i,j)
/**
 * @brief Symbolic expression for entry \f$(i,j)\f$ in matrix of derivative 
 * order \f$o\f$ in \a T 
 */
#define amc_transform_string(T,o,i,j)
/**
 * @brief Compiled (tinyexpr) expression for entry \f$(i,j)\f$ in matrix of
 *  derivative order \f$o\f$ in \a T 
 */
#define amc_transform_expr(T,o,i,j)
/**
 * @brief Name of tinyexpr variable \f$i\f$ in \a T
 */
#define amc_transform_variable_name(T,i)
/**
 * @brief Address of value of tinyexpr variable \f$i\f$ in \a T
 */
#define amc_transform_variable_address(T,i)
/**
 * @brief TRUE if \a T is tagged
 */
#define amc_transform_tag(T)
/**
 * @brief TRUE if expressions in \a T have been compiled
 */
#define amc_transform_compiled(T)

#else /*DOXYGEN*/
typedef struct _amc_transform_t amc_transform_t ;

struct _amc_transform_t
{
  int dim,	  /**< dimension of transform*/
    msize,  /**< number of entries in affine matrices*/
    order,      /**< order of derivatives available*/
    order_max,  /**< maximum order of derivatives available (used in allocation)*/
    nvars ;     /**< number of variables*/
  double t ;    /**< kept locally for evaluation purposes*/
  double *matrix_values ; /**< transform matrix and derivatives: numerical values*/
  char **matrix_strings ;  /**< transform matrix and derivatives: expressions*/
  te_expr **matrix_expr   ; /**< transform matrix and derivatives: tinyexpr expr*/
  te_variable vars[AMC_TRANSFORM_VARIABLE_NUMBER] ; /**< variables used by transform */
  int tag, compiled ;
} ;

#define amc_transform_order(_T)           ((_T)->order)
#define amc_transform_order_max(_T)       ((_T)->order_max)
#define amc_transform_variable_number(_T) ((_T)->nvars)
#define amc_transform_dimension(_T)       (*((int *)(_T)))
#define amc_transform_matrix_size(_T)     ((_T)->msize)
#define amc_transform_matrix(_T,_i)		\
  (&((_T)->matrix_values[((_T)->msize)*(_i)]))
#define amc_transform_value(_T,_o,_i,_j)	\
  ((_T)->matrix_values[((_T)->msize)*(_o)+((_T)->dim+1)*(_i)+(_j)])
#define amc_transform_string(_T,_o,_i,_j)	\
  ((_T)->matrix_strings[((_T)->msize)*(_o)+((_T)->dim+1)*(_i)+(_j)])
#define amc_transform_expr(_T,_o,_i,_j)	\
  ((_T)->matrix_expr[((_T)->msize)*(_o)+((_T)->dim+1)*(_i)+(_j)])
#define amc_transform_variable_name(_T,_i) ((_T)->vars[(_i)].name)
#define amc_transform_variable_address(_T,_i) \
  ((double *)((_T)->vars[(_i)].address))
#define amc_transform_tag(_T) ((_T)->tag)
#define amc_transform_compiled(_T) ((_T)->compiled)
#endif /*DOXYGEN*/

#ifdef DOXYGEN
/** 
 * @typedef amc_transform_chain_t
 * 
 * Chain of transforms applied successively to build up general motion of
 * points
 */
typedef amc_transform_chain_t ;

/**
 * @brief Number of transforms in chain \a C
 */
#define amc_transform_chain_transform_number(C)     
/**
 * @brief Maximum number of transforms in chain \a C
 */
#define amc_transform_chain_transform_number_max(C) 
/**
 * @brief Pointer to \f$i\f$th transform in chain \a C
 */
#define amc_transform_chain_transform(C,i)         
/**
 * @brief TRUE if \a C is tagged
 */
#define amc_transform_chain_tag(C)

#else /*DOXYGEN*/
typedef struct _amc_transform_chain_t amc_transform_chain_t ;
struct _amc_transform_chain_t
{
  int
    ntrans,      /**< number of transforms in chain*/
    ntrans_max, /**< maximum number of transforms in chain*/
    tag ;       
  amc_transform_t **T ;		/**< transforms in chain */
} ;

#define amc_transform_chain_transform_number(_C)     ((_C)->ntrans)
#define amc_transform_chain_transform_number_max(_C) ((_C)->ntrans_max)
#define amc_transform_chain_transform(_C,_i)         ((_C)->T[(_i)])
#define amc_transform_chain_tag(_C)                  ((_C)->tag)

#define amc_transform_check(_T)						\
  do {									\
  if ( amc_transform_dimension((_T)) != 2 &&				\
       amc_transform_dimension((_T)) != 3 ) {				\
  fprintf(stderr, "%s: transform dimension (%d) not equal to 2 or 3\n", \
	  __FUNCTION__, amc_transform_dimension((_T))) ;		\
  return 1 ;								\
  }									\
  } while (0)
#endif /*DOXYGEN*/

amc_transform_t *amc_transform_alloc(int dim, int order) ;
int amc_transform_entry_set(amc_transform_t *T, int order,
			    int i, int j, double val, char *str) ;
int amc_transform_matrix_write(FILE *f, amc_transform_t *T, int order) ;
int amc_transform_expressions_compile(amc_transform_t *T) ;
int amc_transform_matrix_apply(amc_transform_t *T, int order,
			       double *xin, double *xout) ;
int amc_transform_matrix_apply_vector(amc_transform_t *T, int order,
				      double *xin, double *xout) ;
int amc_transform_matrices_evaluate(amc_transform_t *T, double t) ;
int amc_transform_variable_add(amc_transform_t *T, char *str, double *adr) ;
int amc_transform_variables_write(FILE *f, amc_transform_t *T) ;
int amc_transform_matrix_identity(amc_transform_t *T, int order) ;
int amc_transform_matrix_zero(amc_transform_t *T, int order) ;
amc_transform_chain_t *amc_transform_chain_alloc(int ntrans) ;
int amc_transform_chain_transform_add(amc_transform_chain_t *C,
				      amc_transform_t *T) ;
int amc_transform_chain_evaluate(amc_transform_chain_t *C, int order,
				 amc_transform_t *T) ;
int amc_transform_chain_derivative(amc_transform_chain_t *C, int order,
				   amc_transform_t *S, double t,
				   double *x, double *dx) ;
int amc_transform_derivatives_evaluate(amc_transform_t *T, int order) ;

int amc_transform_translation(amc_transform_t *T,
			      double dx, char *xstr,
			      double dy, char *ystr,
			      double dz, char *zstr, int order) ;
int amc_transform_rotation(amc_transform_t *T, double th, char *str,
			   int order) ;
int amc_transform_rotation_x(amc_transform_t *T, double th, char *str,
			     int order) ;
int amc_transform_rotation_y(amc_transform_t *T, double th, char *str,
			     int order) ;
int amc_transform_rotation_z(amc_transform_t *T, double th, char *str,
			     int order) ;

int amc_transform_derivative_check(amc_transform_t *T, int order,
				   double t, double dt,
				   int *imax, double *emax) ;
char *amc_sin_derivative(int i) ;
char *amc_cos_derivative(int i) ;
int amc_transform_chain_matrices_evaluate(amc_transform_chain_t *C, double t,
					  amc_transform_t *T, int order) ;

#endif /*__AMC_H_INCLUDED__*/
