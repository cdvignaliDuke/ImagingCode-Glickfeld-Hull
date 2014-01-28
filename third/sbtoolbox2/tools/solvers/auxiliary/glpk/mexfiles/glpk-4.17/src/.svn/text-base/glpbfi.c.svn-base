/* glpbfi.c (basis factorization interface) */

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

#define _GLPBFI_PRIVATE
#define _GLPLPF_PRIVATE
#define _GLPSCF_PRIVATE
#include "glpbfi.h"
#include "glplib.h"

BFI *bfi_create_binv(void)
{     /* create factorization of the basis matrix */
      BFI *binv;
      binv = xmalloc(sizeof(BFI));
      binv->valid = 0;
      binv->type = 1; /* default */
      binv->m_max = 0;
      binv->inv = NULL;
      binv->lpf = NULL;
      return binv;
}

int bfi_factorize(BFI *binv, int m, int (*col)(void *info, int j,
      int ind[], double val[]), void *info)
{     /* compute factorization of the basis matrix */
      int ret;
      xassert(m > 0);
      switch (binv->type)
      {  case 1:
            /* INV */
            if (binv->lpf != NULL) lpf_delete_it(binv->lpf);
            binv->lpf = NULL;
            if (binv->m_max < m)
            {  if (binv->inv != NULL) inv_delete(binv->inv);
               binv->m_max = m + 100;
               binv->inv = inv_create(binv->m_max, 50);
            }
            binv->inv->m = m;
            binv->inv->luf->n = m;
            ret = inv_decomp(binv->inv, info, col);
            break;
         case 2:
         case 3:
            /* LPF */
            binv->m_max = 0;
            if (binv->inv != NULL) inv_delete(binv->inv);
            binv->inv = NULL;
            if (binv->lpf == NULL)
            {  binv->lpf = lpf_create_it();
               if (binv->type == 2)
                  binv->lpf->scf->t_opt = SCF_TBG;
               else
                  binv->lpf->scf->t_opt = SCF_TGR;
            }
            ret = lpf_factorize(binv->lpf, m, col, info);
            break;
         default:
            xassert(binv != binv);
      }
      binv->valid = (ret == 0);
      return ret;
}

void bfi_ftran(BFI *binv, double x[], int save)
{     /* perform forward transformation (FTRAN) */
      xassert(binv->valid);
      switch (binv->type)
      {  case 1:
            inv_ftran(binv->inv, x, save);
            break;
         case 2:
         case 3:
            lpf_do_ftran(binv->lpf, x, save);
            break;
         default:
            xassert(binv != binv);
      }
      return;
}

void bfi_btran(BFI *binv, double x[])
{     /* perform backward transformation (BTRAN) */
      xassert(binv->valid);
      switch (binv->type)
      {  case 1:
            inv_btran(binv->inv, x);
            break;
         case 2:
         case 3:
            lpf_do_btran(binv->lpf, x);
            break;
         default:
            xassert(binv != binv);
      }
      return;
}

int bfi_update_binv(BFI *binv, int j)
{     /* update factorization of the basis matrix */
      int ret;
      xassert(binv->valid);
      switch (binv->type)
      {  case 1:
            ret = inv_update(binv->inv, j);
            break;
         case 2:
         case 3:
            ret = lpf_update_it(binv->lpf, j);
            break;
         default:
            xassert(binv != binv);
      }
      binv->valid = (ret == 0);
      return ret;
}

void bfi_delete_binv(BFI *binv)
{     /* delete factorization of the basis matrix */
      if (binv->inv != NULL) inv_delete(binv->inv);
      if (binv->lpf != NULL) lpf_delete_it(binv->lpf);
      xfree(binv);
      return;
}

/* eof */
