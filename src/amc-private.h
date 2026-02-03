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

#ifndef __AMC_PRIVATE_INCLUDED__
#define __AMC_PRIVATE_INCLUDED__


/*
 * y := A*x (safe to perform in place with x == y)
 */

#define amc2d_matrix_vector_mul(_y,_A,_x)		\
  do {							\
    gdouble _tmp ;					\
    _tmp  = (_A)[0]*(_x)[0] + (_A)[1]*(_x)[1] + (_A)[2] ;	\
    (_y)[1] = (_A)[3]*(_x)[0] + (_A)[4]*(_x)[1] + (_A)[5] ;	\
    (_y)[0] = _tmp ;					\
  } while (0)

/*
 * C := bt*C + al*A*B (safe to perform in place with C == A or B)
 */

#define amc2d_matrix_matrix_mul(_al,_A,_B,_bt,_C)			\
  do {									\
  gdouble _Ctmp[9] ;							\
  gint _i ;								\
  _Ctmp[0] = (_A)[0]*(_B)[0] + (_A)[1]*(_B)[3] + (_A)[2]*(_B)[6] ;	\
  _Ctmp[1] = (_A)[0]*(_B)[1] + (_A)[1]*(_B)[4] + (_A)[2]*(_B)[7] ;	\
  _Ctmp[2] = (_A)[0]*(_B)[2] + (_A)[1]*(_B)[5] + (_A)[2]*(_B)[8] ;	\
  _Ctmp[3] = (_A)[3]*(_B)[0] + (_A)[4]*(_B)[3] + (_A)[5]*(_B)[6] ;	\
  _Ctmp[4] = (_A)[3]*(_B)[1] + (_A)[4]*(_B)[4] + (_A)[5]*(_B)[7] ;	\
  _Ctmp[5] = (_A)[3]*(_B)[2] + (_A)[4]*(_B)[5] + (_A)[5]*(_B)[8] ;	\
  _Ctmp[6] = (_A)[6]*(_B)[0] + (_A)[7]*(_B)[3] + (_A)[8]*(_B)[6] ;	\
  _Ctmp[7] = (_A)[6]*(_B)[1] + (_A)[7]*(_B)[4] + (_A)[8]*(_B)[7] ;	\
  _Ctmp[8] = (_A)[6]*(_B)[2] + (_A)[7]*(_B)[5] + (_A)[8]*(_B)[8] ;	\
									\
  for ( _i = 0 ; _i < 9 ; _i ++ ) (_C)[_i] =				\
				    (_bt)*(_C)[_i] + (_al)*_Ctmp[_i] ;	\
  } while (0)

/*
 * y := A*x (safe to perform in place with x == y)
 */

#define amc3d_matrix_vector_mul(_y,_A,_x)				\
  do {									\
    gdouble _tmp[2] ;							\
  _tmp[0] = (_A)[ 0]*(_x)[0]+(_A)[ 1]*(_x)[1]+(_A)[ 2]*(_x)[2]+(_A)[ 3] ; \
  _tmp[1] = (_A)[ 4]*(_x)[0]+(_A)[ 5]*(_x)[1]+(_A)[ 6]*(_x)[2]+(_A)[ 7] ; \
  (_y)[2] = (_A)[ 8]*(_x)[0]+(_A)[ 9]*(_x)[1]+(_A)[10]*(_x)[2]+(_A)[11] ; \
  (_y)[0] = _tmp[0] ; (_y)[1] = _tmp[1] ;				\
  } while (0)

/*
 * C := bt*C + al*A*B (safe to perform in place with C == A or B)
 */

#define amc3d_matrix_matrix_mul(_al,_A,_B,_bt,_C)			\
  do {									\
  gdouble _Ctmp[16] ;							\
  gint _i ;								\
  _Ctmp[ 0] = (_A)[ 0]*(_B)[ 0]+(_A)[ 1]*(_B)[ 4]+(_A)[ 2]*(_B)[ 8]+(_A)[ 3]*(_B)[12] ; \
  _Ctmp[ 1] = (_A)[ 0]*(_B)[ 1]+(_A)[ 1]*(_B)[ 5]+(_A)[ 2]*(_B)[ 9]+(_A)[ 3]*(_B)[13] ; \
  _Ctmp[ 2] = (_A)[ 0]*(_B)[ 2]+(_A)[ 1]*(_B)[ 6]+(_A)[ 2]*(_B)[10]+(_A)[ 3]*(_B)[14] ; \
  _Ctmp[ 3] = (_A)[ 0]*(_B)[ 3]+(_A)[ 1]*(_B)[ 7]+(_A)[ 2]*(_B)[11]+(_A)[ 3]*(_B)[15] ; \
  _Ctmp[ 4] = (_A)[ 4]*(_B)[ 0]+(_A)[ 5]*(_B)[ 4]+(_A)[ 6]*(_B)[ 8]+(_A)[ 7]*(_B)[12] ; \
  _Ctmp[ 5] = (_A)[ 4]*(_B)[ 1]+(_A)[ 5]*(_B)[ 5]+(_A)[ 6]*(_B)[ 9]+(_A)[ 7]*(_B)[13] ; \
  _Ctmp[ 6] = (_A)[ 4]*(_B)[ 2]+(_A)[ 5]*(_B)[ 6]+(_A)[ 6]*(_B)[10]+(_A)[ 7]*(_B)[14] ; \
  _Ctmp[ 7] = (_A)[ 4]*(_B)[ 3]+(_A)[ 5]*(_B)[ 7]+(_A)[ 6]*(_B)[11]+(_A)[ 7]*(_B)[15] ; \
  _Ctmp[ 8] = (_A)[ 8]*(_B)[ 0]+(_A)[ 9]*(_B)[ 4]+(_A)[10]*(_B)[ 8]+(_A)[11]*(_B)[12] ; \
  _Ctmp[ 9] = (_A)[ 8]*(_B)[ 1]+(_A)[ 9]*(_B)[ 5]+(_A)[10]*(_B)[ 9]+(_A)[11]*(_B)[13] ; \
  _Ctmp[10] = (_A)[ 8]*(_B)[ 2]+(_A)[ 9]*(_B)[ 6]+(_A)[10]*(_B)[10]+(_A)[11]*(_B)[14] ; \
  _Ctmp[11] = (_A)[ 8]*(_B)[ 3]+(_A)[ 9]*(_B)[ 7]+(_A)[10]*(_B)[11]+(_A)[11]*(_B)[15] ; \
  _Ctmp[12] = (_A)[12]*(_B)[ 0]+(_A)[13]*(_B)[ 4]+(_A)[14]*(_B)[ 8]+(_A)[15]*(_B)[12] ; \
  _Ctmp[13] = (_A)[12]*(_B)[ 1]+(_A)[13]*(_B)[ 5]+(_A)[14]*(_B)[ 9]+(_A)[15]*(_B)[13] ; \
  _Ctmp[14] = (_A)[12]*(_B)[ 2]+(_A)[13]*(_B)[ 6]+(_A)[14]*(_B)[10]+(_A)[15]*(_B)[14] ; \
  _Ctmp[15] = (_A)[12]*(_B)[ 3]+(_A)[13]*(_B)[ 7]+(_A)[14]*(_B)[11]+(_A)[15]*(_B)[15] ; \
for ( _i = 0 ; _i < 16 ; _i ++ ) (_C)[_i] = (_bt)*(_C)[_i]+(_al)*_Ctmp[_i] ; \
 } while (0)

#define amc2d_matrix_identity(_A)		\
  do {						\
    (_A)[0] = 1 ; (_A)[1] = 0 ; (_A)[2] = 0 ;	\
    (_A)[3] = 0 ; (_A)[4] = 1 ; (_A)[5] = 0 ;	\
    (_A)[6] = 0 ; (_A)[7] = 0 ; (_A)[8] = 1 ;	\
  } while (0)

#define amc3d_matrix_identity(_A)				\
  do {								\
    (_A)[ 0] = 1 ; (_A)[ 1] = 0 ; (_A)[ 2] = 0 ; (_A)[ 3] = 0 ;	\
    (_A)[ 4] = 0 ; (_A)[ 5] = 1 ; (_A)[ 6] = 0 ; (_A)[ 7] = 0 ;	\
    (_A)[ 8] = 0 ; (_A)[ 9] = 0 ; (_A)[10] = 1 ; (_A)[11] = 0 ;	\
    (_A)[12] = 0 ; (_A)[13] = 0 ; (_A)[14] = 0 ; (_A)[15] = 1 ;	\
  } while (0)

#endif /*__AMC_PRIVATE_INCLUDED__*/

/* G_TOKEN_EOF 	 */

/* The end of the file. */
/* G_TOKEN_LEFT_PAREN 	 */

/* A ‘(‘ character. */
/* G_TOKEN_RIGHT_PAREN 	 */

/* A ‘)’ character. */
/* G_TOKEN_LEFT_CURLY 	 */

/* A ‘{‘ character. */
/* G_TOKEN_RIGHT_CURLY 	 */

/* A ‘}’ character. */
/* G_TOKEN_LEFT_BRACE 	 */

/* A ‘[‘ character. */
/* G_TOKEN_RIGHT_BRACE 	 */

/* A ‘]’ character. */
/* G_TOKEN_EQUAL_SIGN 	 */

/* A ‘=’ character. */
/* G_TOKEN_COMMA 	 */

/* A ‘,’ character. */
/* G_TOKEN_NONE 	 */

/* Not a token. */
/* G_TOKEN_ERROR 	 */

/* An error occurred. */
/* G_TOKEN_CHAR 	 */

/* A character. */
/* G_TOKEN_BINARY 	 */

/* A binary integer. */
/* G_TOKEN_OCTAL 	 */

/* An octal integer. */
/* G_TOKEN_INT 	 */

/* An integer. */
/* G_TOKEN_HEX 	 */

/* A hex integer. */
/* G_TOKEN_FLOAT 	 */

/* A floating point number. */
/* G_TOKEN_STRING 	 */

/* A string. */
/* G_TOKEN_SYMBOL 	 */

/* A symbol. */
/* G_TOKEN_IDENTIFIER 	 */

/* An identifier. */
/* G_TOKEN_IDENTIFIER_NULL 	 */

/* A null identifier. */
/* G_TOKEN_COMMENT_SINGLE 	 */

/* One line comment. */
/* G_TOKEN_COMMENT_MULTI */
