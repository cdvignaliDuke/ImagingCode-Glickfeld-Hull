/* glpbfi.h (basis factorization interface) */

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

#ifndef _GLPBFI_H
#define _GLPBFI_H

#ifdef _GLPBFI_PRIVATE
#include "glpinv.h"
#include "glplpf.h"

typedef struct BFI BFI;

struct BFI
{     /* basis factorization interface */
      int valid;
      /* this flag is set iff the basis factorization is valid */
      int type;
      /* type of factorization used:
         1 - INV
         2 - LPF + Bartels-Golub update for C
         3 - LPF + Givens rotation update for C */
      int m_max;
      /* maximal order of the basis (enlarged automatically) */
      INV *inv;
      /* invertable form of the basis matrix */
      LPF *lpf;
      /* LP basis factorization */
};
#else
typedef struct { double _binv; } BFI;
#endif

#define bfi_create_binv _glp_bfi_create_binv
BFI *bfi_create_binv(void);
/* create factorization of the basis matrix */

#define bfi_factorize _glp_bfi_factorize
int bfi_factorize(BFI *binv, int m, int (*col)(void *info, int j,
      int ind[], double val[]), void *info);
/* compute factorization of the basis matrix */

#define bfi_ftran _glp_bfi_ftran
void bfi_ftran(BFI *binv, double x[], int save);
/* perform forward transformation (FTRAN) */

#define bfi_btran _glp_bfi_btran
void bfi_btran(BFI *binv, double x[]);
/* perform backward transformation (BTRAN) */

#define bfi_update_binv _glp_bfi_update_binv
int bfi_update_binv(BFI *binv, int j);
/* update factorization of the basis matrix */

#define bfi_delete_binv _glp_bfi_delete_binv
void bfi_delete_binv(BFI *binv);
/* delete factorization of the basis matrix */

#endif

/* eof */
