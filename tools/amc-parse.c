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

static void read_point_file(char *file, gdouble **pts, gint *npts)

{
  FILE *f ;
  gint i ;
  
  f = fopen(file, "r") ;
  if ( f == NULL ) {
    fprintf(stderr, "%s: cannot open point file \"%s\"", progname, file) ;
    exit(1) ;
  }

  fscanf(f, "%d", npts) ;
  *pts = (gdouble *)g_malloc0(3*(*npts)*sizeof(gdouble)) ;

  for ( i = 0 ; i < 3*(*npts) ; i ++ ) fscanf(f, "%lg", &((*pts)[i])) ;
    
  fclose(f) ;
  
  return ;
}

gint main(gint argc, char **argv)

{
  char ch, *frame, *ipfiles[32] = {NULL}, *opfiles[32] = {NULL} ;
  gint nt, nc, nv, nipfiles, nopfiles, npts[32], i, j, k, n, dim, order ;
  gdouble *pts[32], t, t0, t1, x[3] ;
  amc_frame_t *F ;
  amc_transform_t *T ;
  
  FILE *opf[32] ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  nt = 16 ; nc = 8 ; nv = 128 ;
  t0 = 0 ; t1 = 1 ; n = 256 ;
  dim = 3 ; order = 0 ;
  
  frame = NULL ; nipfiles = nopfiles = 0 ;
  while ( (ch = getopt(argc, argv, "f:i:o:t:")) != EOF ) {
    switch (ch ) {
    default: g_assert_not_reached() ; break ;
    case 'f': frame  = g_strdup(optarg) ; break ;
    case 'i': ipfiles[nipfiles] = g_strdup(optarg) ; nipfiles ++ ; break ;
    case 'o': opfiles[nopfiles] = g_strdup(optarg) ; nopfiles ++ ; break ;
    }
  }

  if ( frame == NULL ) {
    fprintf(stderr, "%s: no frame file specified\n", progname) ;
    return 1 ;
  } else {
    F = amc_frame_alloc(nt, nc, nv) ;
    amc_frame_read(F, frame) ;
    amc_frame_initialise(F) ;    
  }

  if ( nipfiles > amc_frame_transform_chain_number(F) ) {
    fprintf(stderr, "%s: more input files than transform chains\n",
	    progname) ;
  }
  
  if ( nipfiles == 0 ) {
    fprintf(stderr, "%s: no input files specified\n", progname) ;
    return 1 ;
  } else {
    for ( i = 0 ; i < nipfiles ; i ++ ) {
      read_point_file(ipfiles[i], &(pts[i]), &(npts[i])) ;
    }
  }

  if ( nopfiles < nipfiles ) {
    fprintf(stderr, "%s: not enough outputs (%d) for number of inputs (%d)\n",
	    progname, nopfiles, nipfiles) ;
    return 1 ;
  }
  
  if ( nopfiles == 0 ) {
    fprintf(stderr, "%s: no output files specified\n", progname) ;
    return 1 ;
  } else {
    for ( i = 0 ; i < nopfiles ; i ++ ) {
      opf[i] = fopen(opfiles[i], "w") ;
      if ( opf[i] == NULL ) {
	fprintf(stderr, "%s: cannot open output file \"%s\"",
		progname, opfiles[i]) ;
	return 1 ;
      }
    }
  }

  T = amc_transform_alloc(dim, order) ; T->order = order ;
  
  for ( i = 0 ; i < n ; i ++ ) {
    t = t0 + (t1 - t0)*i/n ;
    amc_frame_evaluate(F, t) ;
    for ( j = 0 ; j < nopfiles ; j ++ ) {
      amc_transform_chain_evaluate(amc_frame_transform_chain(F,j), 0, T) ;
      for ( k = 0 ; k < npts[j] ; k ++ ) {
	amc_transform_matrix_apply(T, 0, &((pts[j])[3*k]), x) ;
	fprintf(opf[j], "%e %e %e ", x[0], x[1], x[2]) ;
      }
      fprintf(opf[j], "\n") ;
    }
  }
  
  return 0 ;
}
