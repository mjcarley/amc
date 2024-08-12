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

typedef enum {
  AMC_TRANSFORM_DEFINITION_UNKNOWN = 0,
  AMC_TRANSFORM_DEFINITION_MATRIX = 1,
  AMC_TRANSFORM_DEFINITION_ROTATION_X = 2,
  AMC_TRANSFORM_DEFINITION_ROTATION_Y = 3,
  AMC_TRANSFORM_DEFINITION_ROTATION_Z = 4
} amc_transform_definition_t ;

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
  gint dim,	  /**< dimension of transform*/
    msize,  /**< number of entries in affine matrices*/
    order,      /**< order of derivatives available*/
    order_max,  /**< maximum order of derivatives available (used in allocation)*/
    nvars ;     /**< number of variables*/
  amc_transform_definition_t def ;
  gdouble t ;    /**< kept locally for evaluation purposes*/
  gdouble *matrix_values ; /**< transform matrix and derivatives: numerical values*/
  char **matrix_strings ;  /**< transform matrix and derivatives: expressions*/
  gpointer *matrix_expr   ; /**< transform matrix and derivatives: tinyexpr expr*/
  gpointer vars ; /**< variables used by transform */
  gint tag, compiled ;
} ;

#define amc_transform_dimension(_T)       (*((gint *)(_T)))
#define amc_transform_definition(_T)      ((_T)->def)
#define amc_transform_order(_T)           ((_T)->order)
#define amc_transform_order_max(_T)       ((_T)->order_max)
#define amc_transform_variable_number(_T) ((_T)->nvars)
#define amc_transform_matrix_size(_T)     ((_T)->msize)
#define amc_transform_matrix(_T,_i)		\
  (&((_T)->matrix_values[((_T)->msize)*(_i)]))
#define amc_transform_value(_T,_o,_i,_j)	\
  ((_T)->matrix_values[((_T)->msize)*(_o)+((_T)->dim+1)*(_i)+(_j)])
#define amc_transform_string(_T,_o,_i,_j)	\
  ((_T)->matrix_strings[((_T)->msize)*(_o)+((_T)->dim+1)*(_i)+(_j)])
#define amc_transform_expr(_T,_o,_i,_j)	\
  ((_T)->matrix_expr[((_T)->msize)*(_o)+((_T)->dim+1)*(_i)+(_j)])
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
  gint
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

typedef struct _amc_frame_t amc_frame_t ;
struct _amc_frame_t {
  amc_transform_t       **T ;
  amc_transform_chain_t **C ;
  gint nt, ntmax, nchain, nc, ncmax, nv, nvmax, dim, order ;
  char **transforms, **chains ;
  gpointer vars ; /**< variables used by transforms */
  gdouble *vals ; /**< values of variables */
} ;

#define amc_frame_transform(_f,_i)               ((_f)->T[(_i)])
#define amc_frame_transform_number(_f)           ((_f)->nt)
#define amc_frame_transform_number_max(_f)       ((_f)->ntmax)
#define amc_frame_transform_chain(_f,_i)         ((_f)->C[(_i)])
#define amc_frame_transform_chain_number(_f)     ((_f)->nc)
#define amc_frame_transform_chain_number_max(_f) ((_f)->ncmax)
#define amc_frame_variable_number(_f)            ((_f)->nv)
#define amc_frame_variable_number_max(_f)        ((_f)->nvmax)
#define amc_frame_transform_name(_f,_i)          ((_f)->transforms[(_i)])
#define amc_frame_transform_chain_name(_f,_i)    ((_f)->chains[(_i)])
#define amc_frame_dimension(_f)                  ((_f)->dim)
#define amc_frame_order(_f)                      ((_f)->order)

amc_transform_t *amc_transform_alloc(gint dim, gint order) ;
gint amc_transform_entry_set(amc_transform_t *T, gint order,
			    gint i, gint j, gdouble val, char *str) ;
gint amc_transform_matrix_write(FILE *f, amc_transform_t *T, gint order) ;
gint amc_transform_expressions_compile(amc_transform_t *T) ;
gint amc_transform_matrix_apply(amc_transform_t *T, gint order,
			       gdouble *xin, gdouble *xout) ;
gint amc_transform_matrix_apply_vector(amc_transform_t *T, gint order,
				      gdouble *xin, gdouble *xout) ;
gint amc_transform_matrices_evaluate(amc_transform_t *T, gdouble t) ;
gint amc_transform_variable_add(amc_transform_t *T, char *str, gdouble *adr) ;
gint amc_transform_variables_write(FILE *f, amc_transform_t *T) ;
gint amc_transform_matrix_identity(amc_transform_t *T, gint order) ;
gint amc_transform_matrix_zero(amc_transform_t *T, gint order) ;
amc_transform_chain_t *amc_transform_chain_alloc(gint ntrans) ;
gint amc_transform_chain_transform_add(amc_transform_chain_t *C,
				      amc_transform_t *T) ;
gint amc_transform_chain_evaluate(amc_transform_chain_t *C, gint order,
				 amc_transform_t *T) ;
gint amc_transform_chain_derivative(amc_transform_chain_t *C, gint order,
				   amc_transform_t *S, gdouble t,
				   gdouble *x, gdouble *dx) ;
gint amc_transform_derivatives_evaluate(amc_transform_t *T, gint order) ;
amc_transform_definition_t amc_transform_definition_parse(char *str) ;
gint amc_transform_parse(amc_transform_t *T, char *str, gdouble val,
			 gint order) ;

gint amc_transform_translation(amc_transform_t *T,
			      gdouble dx, char *xstr,
			      gdouble dy, char *ystr,
			      gdouble dz, char *zstr, gint order) ;
gint amc_transform_rotation(amc_transform_t *T, gdouble th, char *str,
			   gint order) ;
gint amc_transform_rotation_x(amc_transform_t *T, gdouble th, char *str,
			     gint order) ;
gint amc_transform_rotation_y(amc_transform_t *T, gdouble th, char *str,
			     gint order) ;
gint amc_transform_rotation_z(amc_transform_t *T, gdouble th, char *str,
			     gint order) ;

gint amc_transform_derivative_check(amc_transform_t *T, gint order,
				   gdouble t, gdouble dt,
				   gint *imax, gdouble *emax) ;
char *amc_sin_derivative(gint i) ;
char *amc_cos_derivative(gint i) ;
gint amc_transform_chain_matrices_evaluate(amc_transform_chain_t *C, gdouble t,
					  amc_transform_t *T, gint order) ;

amc_frame_t *amc_frame_alloc(gint nt, gint nc, gint nv) ;
gint amc_frame_transform_add(amc_frame_t *f, amc_transform_t *T, char *name) ;
gint amc_frame_transform_chain_add(amc_frame_t *f,
				   amc_transform_chain_t *C, char *name) ;
gint amc_frame_transform_find(amc_frame_t *f, char *name) ;
gint amc_frame_transform_chain_find(amc_frame_t *f, char *name) ;
gint amc_frame_initialise(amc_frame_t *f) ;
gint amc_frame_evaluate(amc_frame_t *f, gdouble t) ;
gint amc_frame_write(FILE *f, amc_frame_t *F) ;
gint amc_frame_read(amc_frame_t *F, char *file) ;

#endif /*__AMC_H_INCLUDED__*/
