/* glplpf.h (LP basis factorization) */

/*----------------------------------------------------------------------
-- This code is part of GNU Linear Programming Kit (GLPK).
--
-- Copyright (C) 2000, 01, 02, 03, 04, 05, 06, 07 Andrew Makhorin,
-- Department for Applied Informatics, Moscow Aviation Institute,
-- Moscow, Russia. All rights reserved. E-mail: <mao@mai2.rcnet.ru>.
--
-- GLPK is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 2, or (at your option)
-- any later version.
--
-- GLPK is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
-- or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
-- License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with GLPK; see the file COPYING. If not, write to the Free
-- Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
-- 02110-1301, USA.
----------------------------------------------------------------------*/

#ifndef _GLPLPF_H
#define _GLPLPF_H

#ifdef _GLPLPF_PRIVATE
#include "glpscf.h"
#include "glpluf.h"

/***********************************************************************
*  The structure LPF defines the factorization of the basis mxm matrix
*  B, where m is the number of rows in corresponding problem instance.
*
*  This factorization is the following septet:
*
*     [B] = (L0, U0, R, S, C, P, Q),                                 (1)
*
*  and is based on the following main equality:
*
*     ( B  F^)     ( B0 F )       ( L0 0 ) ( U0 R )
*     (      ) = P (      ) Q = P (      ) (      ) Q,               (2)
*     ( G^ H^)     ( G  H )       ( S  I ) ( 0  C )
*
*  where:
*
*  B is the current basis matrix (not stored);
*
*  F^, G^, H^ are some additional matrices (not stored);
*
*  B0 is some initial basis matrix (not stored);
*
*  F, G, H are some additional matrices (not stored);
*
*  P, Q are permutation matrices (stored in both row- and column-like
*  formats);
*
*  L0, U0 are some matrices that defines a factorization of the initial
*  basis matrix B0 = L0 * U0 (stored in an invertable form);
*
*  R is a matrix defined from L0 * R = F, so R = inv(L0) * F (stored in
*  a column-wise sparse format);
*
*  S is a matrix defined from S * U0 = G, so S = G * inv(U0) (stored in
*  a row-wise sparse format);
*
*  C is the Schur complement for matrix (B0 F G H). It is defined from
*  S * R + C = H, so C = H - S * R = H - G * inv(U0) * inv(L0) * F =
*  = H - G * inv(B0) * F. Matrix C is stored in an invertable form.
*
*  REFERENCES
*
*  1. M.A.Saunders, "LUSOL: A basis package for constrained optimiza-
*     tion," SCCM, Stanford University, 2006.
*
*  2. M.A.Saunders, "Notes 5: Basis Updates," CME 318, Stanford Univer-
*     sity, Spring 2006.
*
*  3. M.A.Saunders, "Notes 6: LUSOL---a Basis Factorization Package,"
*     ibid. */

typedef struct LPF LPF;

struct LPF
{     /* LP basis factorization */
      /*--------------------------------------------------------------*/
      /* initial basis matrix B0 */
      int m0;
      /* the order of B0 */
      LUF *luf;
      /* LU-factorization of B0 */
      /*--------------------------------------------------------------*/
      /* current basis matrix B */
      int m;
      /* the order of B */
      double *B; /* double B[1+m*m]; */
      /* B in dense format stored by rows and used only for debugging;
         normally this array is not allocated */
      /*--------------------------------------------------------------*/
      /* augmented matrix (B0 F G H) of the order m0+n */
      int n_max;
      /* maximal number of additional rows and columns */
      int n;
      /* current number of additional rows and columns */
      /*--------------------------------------------------------------*/
      /* m0xn matrix R in column-wise format */
      int *R_ptr; /* int R_ptr[1+n_max]; */
      /* R_ptr[0] is not used;
         R_ptr[j], 1 <= j <= n, is a pointer to j-th column */
      int *R_len; /* int R_len[1+n_max]; */
      /* R_len[0] is not used;
         R_len[j], 1 <= j <= n, is the length of j-th column */
      /*--------------------------------------------------------------*/
      /* nxm0 matrix S in row-wise format */
      int *S_ptr; /* int S_ptr[1+n_max]; */
      /* S_ptr[0] is not used;
         S_ptr[i], 1 <= i <= n, is a pointer to i-th row */
      int *S_len; /* int S_len[1+n_max]; */
      /* S_len[0] is not used;
         S_len[i], 1 <= i <= n, is the length of i-th row */
      /*--------------------------------------------------------------*/
      /* Schur complement C of the order n */
      SCF *scf; /* SCF scf[1:n_max]; */
      /* factorization of the Schur complement */
      /*--------------------------------------------------------------*/
      /* matrix P of the order m0+n */
      int *P_row; /* int P_row[1+m0+n_max]; */
      /* P_row[0] is not used; P_row[i] = j means that P[i,j] = 1 */
      int *P_col; /* int P_col[1+m0+n_max]; */
      /* P_col[0] is not used; P_col[j] = i means that P[i,j] = 1 */
      /*--------------------------------------------------------------*/
      /* matrix Q of the order m0+n */
      int *Q_row; /* int Q_row[1+m0+n_max]; */
      /* Q_row[0] is not used; Q_row[i] = j means that Q[i,j] = 1 */
      int *Q_col; /* int Q_col[1+m0+n_max]; */
      /* Q_col[0] is not used; Q_col[j] = i means that Q[i,j] = 1 */
      /*--------------------------------------------------------------*/
      /* Sparse Vector Area (SVA) is a set of locations intended to
         store sparse vectors which represent columns of matrix R and
         rows of matrix S; each location is a doublet (ind, val), where
         ind is an index, val is a numerical value of a sparse vector
         element; in the whole each sparse vector is a set of adjacent
         locations defined by a pointer to its first element and its
         length, i.e. the number of its elements */
      int v_size;
      /* the SVA size, in locations; locations are numbered by integers
         1, 2, ..., v_size, and location 0 is not used */
      int v_ptr;
      /* pointer to the first available location */
      int *v_ind; /* int v_ind[1+v_size]; */
      /* v_ind[0] is not used;
         v_ind[k], 1 <= k <= v_size, is the index field of location k */
      double *v_val; /* double v_val[1+v_size]; */
      /* v_val[0] is not used;
         v_val[k], 1 <= k <= v_size, is the value field of location k */
      /*--------------------------------------------------------------*/
      int saved;
      /* if this flag is set, new column is provided */
      double *newcol; /* double newcol[1+m0+n_max]; */
      /* new column chosen to enter the basis */
      double *work1; /* double work1[1+m0+n_max]; */
      /* working array */
      double *work2; /* double work2[1+m0+n_max]; */
      /* working array */
};
#else
typedef struct { double _lpf; } LPF;
#endif

/* return codes: */
#define LPF_ELIMIT 1 /* limit reached */
#define LPF_ESING  2 /* singular basis */

#define lpf_create_it _glp_lpf_create_it
LPF *lpf_create_it(void);
/* create LP basis factorization */

#define lpf_factorize _glp_lpf_factorize
int lpf_factorize(LPF *lpf, int m, int (*col)(void *info, int j,
      int ind[], double val[]), void *info);
/* compute LP basis factorization */

#define lpf_do_ftran _glp_lpf_do_ftran
void lpf_do_ftran(LPF *lpf, double x[], int save);
/* perform forward transformation (FTRAN) */

#define lpf_do_btran _glp_lpf_do_btran
void lpf_do_btran(LPF *lpf, double x[]);
/* perform backward transformation (BTRAN) */

#define lpf_update_it _glp_lpf_update_it
int lpf_update_it(LPF *lpf, int j);
/* update LP basis factorization */

#define lpf_delete_it _glp_lpf_delete_it
void lpf_delete_it(LPF *lpf);
/* delete LP basis factorization */

#endif

/* eof */
