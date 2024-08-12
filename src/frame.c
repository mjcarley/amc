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
#include <fcntl.h>
#include <unistd.h>

#include <glib.h>

#include "tinyexpr.h"

#include "amc.h"

#define PARSER_DATA_SIZE    8
#define PARSER_DATA_FRAME   0

#define AMC_TOKEN_ERROR        -1
#define AMC_TOKEN_UNIDENTIFIED  0
#define AMC_TOKEN_VARIABLE      1
#define AMC_TOKEN_TRANSFORM     2
#define AMC_TOKEN_CHAIN         3
#define AMC_TOKEN_DIMENSION     4

static gint variable_find(te_variable *vars, gint nvars, te_variable v)

{
  gint i ;

  for ( i = 0 ; i < nvars ; i ++ ) {
    if ( strcmp(vars[i].name, v.name) == 0 ) return i ;
  }
  
  return -1 ;
}

static gint pointer_find(gpointer *ptrs, gint np, gpointer p)

{
  gint i ;

  for ( i = 0 ; i < np ; i ++ ) {
    if ( ptrs[i] == p ) return i ;
  }    
  
  return -1 ;
}

static gint token_parse(char *t)

{
  gint i ;
  char *tokens[] = {"variable", "transform", "chain", "dimension", NULL} ;

  for ( i = 0 ; tokens[i] != NULL ; i ++ ) {
    if ( strcmp(t, tokens[i]) == 0 ) {
      return i + 1 ;
    }
  }    
  
  return AMC_TOKEN_ERROR ;
}  

amc_frame_t *amc_frame_alloc(gint nt, gint nc, gint nv)

{
  amc_frame_t *F ;

  F = (amc_frame_t *)g_malloc0(sizeof(amc_frame_t)) ;

  F->T = (amc_transform_t **)g_malloc0(nt*sizeof(amc_transform_t *)) ;
  F->C = (amc_transform_chain_t **)
    g_malloc0(nt*sizeof(amc_transform_chain_t *)) ;

  F->transforms = (char **)g_malloc0((nt+nc)*sizeof(char *)) ;
  F->chains = &(F->transforms[nt]) ;
  F->vars = g_malloc0(nv*sizeof(te_variable)) ;
  F->vals = (gdouble *)g_malloc0(nv*sizeof(gdouble)) ;

  amc_frame_transform_number(F) = 0 ;
  amc_frame_transform_number_max(F) = nt ;
  amc_frame_transform_chain_number(F) = 0 ;
  amc_frame_transform_chain_number_max(F) = nc ;
  amc_frame_variable_number(F) = 0 ;
  amc_frame_variable_number_max(F) = nv ;
  
  return F ;
}

gint amc_frame_transform_add(amc_frame_t *F, amc_transform_t *T, char *name)

{
  gint i, j ;
  te_variable *Tvars, *vars ;

  if ( amc_frame_transform_number(F) >= amc_frame_transform_number_max(F) ) {
    g_error("%s: not enough space for %d transforms",
	    __FUNCTION__, amc_frame_transform_number(F)+1) ;
  }

  /*check name is not a duplicate*/
  i = amc_frame_transform_find(F, name) ;
  if ( i != -1 ) {
    g_error("%s: transform \"%s\" is already in frame", __FUNCTION__, name) ;
  }

  for ( i = 0 ; i < amc_frame_transform_number(F) ; i ++ ) {
    if ( amc_frame_transform(F,i) == T ) {
      g_error("%s: transform already in frame as \"%s\"",
	      __FUNCTION__, amc_frame_transform_name(F,i)) ;
    }
  }

  i = amc_frame_transform_number(F) ;
  amc_frame_transform(F,i)      = T ;
  amc_frame_transform_name(F,i) = g_strdup(name) ;
  amc_frame_transform_number(F) ++ ;

  /*import any variables from T and check for conflicts*/
  Tvars = (te_variable *)(T->vars) ;
  vars  = (te_variable *)(F->vars) ;

  /*start from 1 because each transform has its own local time t*/
  for ( i = 1 ; i < amc_transform_variable_number(T) ; i ++ ) {
    j = variable_find(vars, amc_frame_variable_number(F), Tvars[i]) ;
    if ( j == -1 ) {
      /*add variable main list in f*/
      if ( amc_frame_variable_number(F) >= amc_frame_variable_number_max(F) ) {
	g_error("%s: not enough space in frame for %d variables",
		__FUNCTION__, amc_frame_variable_number(F)+1) ;
      }
      vars[amc_frame_variable_number(F)].name = g_strdup(Tvars[i].name) ;
      vars[amc_frame_variable_number(F)].address = Tvars[i].address ;
      amc_frame_variable_number(F) ++ ;
    } else {
      /*check the variables really are the same*/
      if ( Tvars[i].address != vars[j].address ) {
	g_error("%s: conflicting addresses for variable \"%s\"",
		__FUNCTION__, vars[j].name) ;
      }
    }
  }

  return 0 ;
}

gint amc_frame_transform_find(amc_frame_t *F, char *name)

{
  gint i ;

  for ( i = 0 ; i < amc_frame_transform_number(F) ; i ++ ) {
    if ( strcmp(name, amc_frame_transform_name(F,i)) == 0 ) return i ;
  }
  
  return -1 ;
}

gint amc_frame_transform_chain_add(amc_frame_t *F,
				   amc_transform_chain_t *C, char *name)

{
  gint i, j ;

  if ( amc_frame_transform_chain_number(F) >=
       amc_frame_transform_chain_number_max(F) ) {
    g_error("%s: not enough space for %d chains",
	    __FUNCTION__, amc_frame_transform_chain_number(F)+1) ;
  }

  /*check name is not a duplicate*/
  i = amc_frame_transform_chain_find(F, name) ;
  if ( i != -1 ) {
    g_error("%s: chain \"%s\" is already in frame", __FUNCTION__, name) ;
  }

  for ( i = 0 ; i < amc_frame_transform_chain_number(F) ; i ++ ) {
    if ( amc_frame_transform_chain(F,i) == C ) {
      g_error("%s: transform chain already in frame as \"%s\"",
	      __FUNCTION__, amc_frame_transform_chain_name(F,i)) ;
    }
  }
  
  i = amc_frame_transform_chain_number(F) ;
  amc_frame_transform_chain(F,i)      = C ;
  amc_frame_transform_chain_name(F,i) = g_strdup(name) ;
  amc_frame_transform_chain_number(F) ++ ;

  /*check all the transforms in C are in f*/
  for ( i = 0 ; i < amc_transform_chain_transform_number(C) ; i ++ ) {
    j = pointer_find((gpointer *)(F->T), amc_frame_transform_number(F),
		     amc_transform_chain_transform(C,i)) ;
    if ( j == -1 ) {
      g_error("%s: transform %d of chain not found in frame",
	      __FUNCTION__, i) ;
    }
  }
  
  return 0 ;
}

gint amc_frame_transform_chain_find(amc_frame_t *F, char *name)

{
  gint i ;

  for ( i = 0 ; i < amc_frame_transform_chain_number(F) ; i ++ ) {
    if ( strcmp(name, amc_frame_transform_chain_name(F,i)) == 0 ) return i ;
  }
  
  return -1 ;
}

gint amc_frame_initialise(amc_frame_t *F)

{
  gint i ;
  amc_transform_t *T ;
  
  for ( i = 0 ; i < amc_frame_transform_number(F) ; i ++ ) {
    T = amc_frame_transform(F,i) ;
    amc_transform_derivatives_evaluate(T, amc_transform_order(T)) ;
    amc_transform_expressions_compile(T) ;
  }
  
  return 0 ;
}

gint amc_frame_evaluate(amc_frame_t *F, gdouble t)

{
  gint i ;
  amc_transform_t *T ;
  
  for ( i = 0 ; i < amc_frame_transform_number(F) ; i ++ ) {
    T = amc_frame_transform(F,i) ;
    amc_transform_matrices_evaluate(T, t) ;
  }
  
  return 0 ;
}

gint amc_frame_write(FILE *f, amc_frame_t *F)

{
  gint i, j, k ;
  amc_transform_t *T ;
  amc_transform_chain_t *C ;
  te_variable *vars ;
  
  vars  = (te_variable *)(F->vars) ;

  for ( i = 0 ; i < amc_frame_variable_number(F) ; i ++ ) {
    fprintf(f, "variable::%s = %lg\n",
	    vars[i].name, *((gdouble *)(vars[i].address))) ;
  }
  
  for ( i = 0 ; i < amc_frame_transform_number(F) ; i ++ ) {
    T = amc_frame_transform(F,i) ;
    fprintf(f, "transform::%s = [\n", amc_frame_transform_name(F,i)) ;
    amc_transform_matrix_write(f, T, 0) ;
    fprintf(f, "]\n") ;
  }

  for ( i = 0 ; i < amc_frame_transform_chain_number(F) ; i ++ ) {
    C = amc_frame_transform_chain(F,i) ;
    fprintf(stderr, "chain::%s = ", amc_frame_transform_chain_name(F,i)) ;
    for ( j = amc_transform_chain_transform_number(C) - 1 ; j >= 0 ; j -- ) {
      k = pointer_find((gpointer *)(F->T), amc_frame_transform_number(F),
		       amc_transform_chain_transform(C,j)) ;
      if ( k == -1 ) {
	g_error("%s: transform %d in chain %d not in frame "
		"(this should not happen)",
		__FUNCTION__, j, i) ;
      }
      fprintf(stderr, "[%s]", amc_frame_transform_name(F,k)) ;
      if ( j > 0 ) fprintf(f, "*") ;
    }
    fprintf(stderr, "\n") ;
  }

  return 0 ;
}

static gboolean token_found(GTokenType *t, gint nt, GTokenType token)

{
  gint i ;

  for ( i = 0 ; i < nt ; i ++ ) {
    if ( token == t[i] ) return TRUE ;
  }
  
  return FALSE ;
}

static GTokenType token_search(GScanner *scanner, GTokenType *t, gint nt)

{
  GTokenType token ;

  do {
    token = g_scanner_get_next_token(scanner) ;
  } while ( !token_found(t, nt, token) ) ;

  return token ;
}

static gint variable_read(GScanner *scanner, te_variable *v,
			  gdouble *val)

{
  GTokenType token, t[4] ;
  gint nt ;
  
  t[0] = G_TOKEN_EOF ; t[1] = G_TOKEN_IDENTIFIER ; nt = 2  ;
  if ( (token = token_search(scanner, t, nt)) != G_TOKEN_IDENTIFIER )
    return -1 ;
  
  /*find the variable name*/
  v->name    = g_strdup(scanner->value.v_identifier) ;
  v->address = val ;
  
 t[1] = G_TOKEN_EQUAL_SIGN ;
 if ( (token = token_search(scanner, t, nt)) != G_TOKEN_EQUAL_SIGN )
   return -1 ;

 t[1] = G_TOKEN_FLOAT ; t[2] = '-' ; nt = 3 ;
 token = token_search(scanner, t, nt) ;
 if ( token != G_TOKEN_FLOAT && token != '-' ) return -1 ;
 if ( token == '-' ) {
   nt = 2 ; token = token_search(scanner, t, nt) ;
   if ( token != G_TOKEN_FLOAT ) return -1 ;
   *val = -(scanner->value.v_float) ;
 } else {
  /*set the variable value*/
  *val = scanner->value.v_float ;
 }
 
  return 0 ;
}

static gint dimension_read(GScanner *scanner, gint *dim)

{
  GTokenType token, t[4] ;
  gint nt ;
  
  t[0] = G_TOKEN_EOF ; t[1] = G_TOKEN_EQUAL_SIGN ; nt = 2  ;
  if ( (token = token_search(scanner, t, nt)) !=
       G_TOKEN_EQUAL_SIGN ) return -1 ;
  
  t[1] = G_TOKEN_FLOAT ;
  if ( (token = token_search(scanner, t, 2)) != G_TOKEN_FLOAT )
    return -1 ;
  
  *dim = (gint)(scanner->value.v_float) ;
  
  return 0 ;
}

static gint transform_read(GScanner *scanner, amc_frame_t *F)

{
  GTokenType token, t[8] ;
  gint i, j, nt, ni, nj ;
  amc_transform_t *T ;
  gdouble sgn ;
  
  t[0] = G_TOKEN_EOF ; t[1] = G_TOKEN_IDENTIFIER ; nt = 2  ;
  if ( (token = token_search(scanner, t, nt)) != G_TOKEN_IDENTIFIER )
    return -1 ;
  
  T = amc_transform_alloc(F->dim, F->order) ; T->order = 0 ;
  amc_frame_transform_add(F, T, scanner->value.v_identifier) ;
  
  t[1] = G_TOKEN_EQUAL_SIGN ;
  if ( (token = token_search(scanner, t, nt)) != G_TOKEN_EQUAL_SIGN )
    return -1 ;

  t[1] = G_TOKEN_IDENTIFIER ; t[2] = '[' ; nt = 3 ;
  if ( (token = token_search(scanner, t, nt)) == G_TOKEN_EOF )
    return -1 ;

  if ( F->dim == 2 ) { ni = nj = 3 ; } else { ni = nj = 4 ; }
  /*token might be an opening square bracket or a function name*/
  if ( token == '[' ) {
    amc_transform_definition(T) = AMC_TRANSFORM_DEFINITION_MATRIX ;
    t[1] = G_TOKEN_STRING ; t[2] = G_TOKEN_FLOAT ; t[3] = '-' ; nt = 4 ;
    /*matrix given explicitly*/
    for ( i = 0 ; i < ni ; i ++ ) {
      for ( j = 0 ; j < nj ; j ++ ) {
	sgn = 1 ;
	token = token_search(scanner, t, nt) ;
	if ( token == '-' ) {
	  sgn = -1 ;
	  token = token_search(scanner, t, nt-1) ;
	}
	switch ( token ) {
	case G_TOKEN_FLOAT:
	  amc_transform_entry_set(T, 0, i, j,
				  sgn*(scanner->value.v_float), NULL) ;
	  break ;
	case G_TOKEN_STRING:
	  amc_transform_entry_set(T, 0, i, j, 0, scanner->value.v_string) ;
	  break ;
	default: return -1 ; break ;
	}
      }
    }
    /*should be closed with square bracket*/
    t[1] = ']' ; nt = 2 ;
    token = token_search(scanner, t, nt) ;
    if ( token != ']' ) return -1 ;
  } else {
    /*parsing functional specifiers*/
    amc_transform_definition(T) =
      amc_transform_definition_parse(scanner->value.v_identifier) ;
    /* fprintf(stderr, "function %s ", scanner->value.v_identifier) ; */
    token = g_scanner_get_next_token(scanner) ;
    if ( token != G_TOKEN_LEFT_PAREN ) return -1 ;
    t[1] = G_TOKEN_STRING ; t[2] = G_TOKEN_FLOAT ; t[3] = '-' ; nt = 4 ;
    token = token_search(scanner, t, nt) ;
    if ( token == '-' ) {
      sgn = -1 ;
      token = token_search(scanner, t, nt-1) ;
    }
    switch ( token ) {
    case G_TOKEN_FLOAT:
      /* fprintf(stderr, "(%lg", sgn*(scanner->value.v_float)) ; */
      amc_transform_parse(T, NULL, sgn*(scanner->value.v_float), T->order) ;
      break ;
    case G_TOKEN_STRING:
      /* fprintf(stderr, "(\"%s\"", scanner->value.v_string) ; */
      amc_transform_parse(T, scanner->value.v_string, 0, T->order) ;
      break ;
    default: return -1 ; break ;
    }
    token = g_scanner_get_next_token(scanner) ;
    if ( token != G_TOKEN_RIGHT_PAREN ) return -1 ;
    /* fprintf(stderr, ")\n") ;     */
  }
  
  return 0 ;
}

static gint chain_add_transform(GScanner *scanner, amc_frame_t *F,
				gint *trans, gint *nt)

{
  GTokenType tok[4], token ;
  gint ntok, i ;

  token = g_scanner_get_next_token(scanner) ;
  if ( token != '[' ) return -1 ;
  /*square brackets containing a transform name*/
  tok[0] = G_TOKEN_EOF ; tok[1] = G_TOKEN_IDENTIFIER ; ntok = 2 ;
  if ( (token = token_search(scanner, tok, ntok)) == G_TOKEN_EOF )
    return -1 ;

  i = amc_frame_transform_find(F, scanner->value.v_identifier) ;
  if ( i == -1 ) return -1 ;

  trans[(*nt)] = i ; (*nt) ++ ;
  
  token = g_scanner_get_next_token(scanner) ;
  if ( token != ']' ) return -1 ;

  /*look ahead to check for a succeeding transform*/
  token = g_scanner_peek_next_token(scanner) ;
  if ( token == '*' ) {
    g_scanner_get_next_token(scanner) ;
    return 0 ;
  }
  
  return 1 ;
}

static gint chain_read(GScanner *scanner, amc_frame_t *F)

{
  GTokenType token, t[4] ;
  gint i, trans[32], ntrans, nt ;
  amc_transform_chain_t *C ;
  char name[256] ;
  gboolean chain_read ;
  
  t[0] = G_TOKEN_EOF ; t[1] = G_TOKEN_IDENTIFIER ; nt = 2  ;
  if ( (token = token_search(scanner, t, nt)) != G_TOKEN_IDENTIFIER )
    return -1 ;
  
  /*name of transform chain*/
  strcpy(name, scanner->value.v_identifier) ;
  
  t[1] = G_TOKEN_EQUAL_SIGN ;
  if ( (token = token_search(scanner, t, nt)) != G_TOKEN_EQUAL_SIGN )
    return -1 ;

  chain_read = FALSE ; ntrans = 0 ;
  do {
    i = chain_add_transform(scanner, F, trans, &ntrans) ;
    if ( i ==  1 ) chain_read = TRUE ;
    if ( i == -1 ) return -1 ;
  } while ( !chain_read ) ;

  C = amc_transform_chain_alloc(ntrans) ;

  for ( i = ntrans - 1 ; i >= 0 ; i -- ) {
    amc_transform_chain_transform_add(C, amc_frame_transform(F, trans[i])) ;
  }

  amc_frame_transform_chain_add(F, C, name) ;
  
  return 0 ;
}

gint amc_frame_read(amc_frame_t *F, char *file)

{
  GScanner *scanner ;
  GTokenType token ;
  gint fd, amc_id, nvars, dim, i, j ;
  te_variable *vars ;
  
  vars  = (te_variable *)(F->vars) ;
  amc_frame_dimension(F) = 0 ;
  amc_frame_order(F) = 1 ;
  
  scanner = g_scanner_new(NULL) ;
  scanner->config->int_2_float = TRUE ;
  scanner->config->scan_identifier_1char = TRUE ;
  
  fd = open(file, O_RDONLY) ;
  if ( fd < 0 ) g_error("%s: cannot open file %s", __FUNCTION__, file) ;

  g_scanner_input_file(scanner, fd) ;

  amc_id = AMC_TOKEN_UNIDENTIFIED ; nvars = 0 ;

  while ( ( token = g_scanner_get_next_token(scanner) ) != G_TOKEN_EOF ) {
    if ( token == G_TOKEN_IDENTIFIER ) {
      if ( amc_id != AMC_TOKEN_UNIDENTIFIED ) {
	g_error("%s: syntax error on line %u",
		__FUNCTION__, g_scanner_cur_line(scanner)) ;
      }

      amc_id = token_parse(scanner->value.v_identifier) ;
      switch ( amc_id ) {
      default:
	g_error("%s: syntax error on line %u",
		__FUNCTION__, g_scanner_cur_line(scanner)) ;
	break ;
      case AMC_TOKEN_DIMENSION:
	if ( dimension_read(scanner, &dim) != 0 ) {
	  g_error("%s: syntax error on line %u",
		  __FUNCTION__, g_scanner_cur_line(scanner)) ;
	}
	if ( dim != 2 && dim != 3 ) {
	  g_error("%s: invalid dimension %d", __FUNCTION__, dim) ;
	}
	amc_frame_dimension(F) = dim ;
	amc_id = AMC_TOKEN_UNIDENTIFIED ;
	break ;
      case AMC_TOKEN_VARIABLE:
	if ( variable_read(scanner, &(vars[nvars]), &(F->vals[nvars])) != 0 ) {
	  g_error("%s: syntax error on line %u",
		  __FUNCTION__, g_scanner_cur_line(scanner)) ;
	}
	nvars ++ ;
	amc_id = AMC_TOKEN_UNIDENTIFIED ;
	break ;
      case AMC_TOKEN_TRANSFORM:
	if ( transform_read(scanner, F) != 0 ) {
	  g_error("%s: syntax error on line %u",
		  __FUNCTION__, g_scanner_cur_line(scanner)) ;
	}
	amc_id = AMC_TOKEN_UNIDENTIFIED ;
	break ;
      case AMC_TOKEN_CHAIN:
	if ( chain_read(scanner, F) != 0 ) {
	  g_error("%s: syntax error on line %u",
		  __FUNCTION__, g_scanner_cur_line(scanner)) ;
	}
	amc_id = AMC_TOKEN_UNIDENTIFIED ;
	break ;
      }      
    }
  }

  amc_frame_variable_number(F) = nvars ;
  
  close(fd) ;
  
  g_scanner_destroy(scanner) ;

  /*propagate variables into the transforms*/
  for ( i = 0 ; i < amc_frame_transform_number(F) ; i ++ ) {
    for ( j = 0 ; j < nvars ; j ++ ) {
      amc_transform_variable_add(amc_frame_transform(F,i),
				 (char *)(vars[j].name),
				 (gdouble *)(vars[j].address)) ;
    }
  }
  
  amc_frame_write(stderr, F) ;
  
  return 0 ;
}
