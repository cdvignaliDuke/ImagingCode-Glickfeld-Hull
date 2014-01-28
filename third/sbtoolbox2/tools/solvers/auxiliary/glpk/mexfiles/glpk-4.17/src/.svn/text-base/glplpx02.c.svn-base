/* glplpx02.c (problem retrieving routines) */

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

#define _GLPAPI_PRIVATE
#define _GLPBFI_PRIVATE
#include "glpapi.h"
#define fault xfault1

/*----------------------------------------------------------------------
-- glp_get_num_int - retrieve number of integer columns.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int glp_get_num_int(glp_prob *lp);
--
-- *Returns*
--
-- The routine glp_get_num_int returns the current number of columns,
-- which are marked as integer. */

int glp_get_num_int(glp_prob *lp)
{     GLPCOL *col;
      int j, count;
      count = 0;
      for (j = 1; j <= lp->n; j++)
      {  col = lp->col[j];
         if (col->kind == GLP_IV) count++;
      }
      return count;
}

/*----------------------------------------------------------------------
-- glp_get_num_bin - retrieve number of binary columns.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int glp_get_num_bin(glp_prob *lp);
--
-- *Returns*
--
-- The routine glp_get_num_bin returns the current number of columns,
-- which are marked as binary. */

int glp_get_num_bin(glp_prob *lp)
{     GLPCOL *col;
      int j, count;
      count = 0;
      for (j = 1; j <= lp->n; j++)
      {  col = lp->col[j];
         if (col->kind == GLP_IV &&
            (col->type == GLP_DB && col->lb == 0.0 && col->ub == 1.0))
            count++;
      }
      return count;
}

/*----------------------------------------------------------------------
-- glp_get_col_kind - retrieve column kind.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int glp_get_col_kind(glp_prob *lp, int j);
--
-- *Returns*
--
-- The routine glp_get_col_kind returns the kind of j-th column, i.e.
-- the kind of corresponding structural variable, as follows:
--
-- GLP_CV - continuous variable;
-- GLP_IV - integer variable;
-- GLP_BV - binary variable */

int glp_get_col_kind(glp_prob *lp, int j)
{     GLPCOL *col;
      int kind;
      if (!(1 <= j && j <= lp->n))
         fault("glp_get_col_kind: j = %d; column number out of range",
            j);
      col = lp->col[j];
      if (col->kind == GLP_CV)
         kind = GLP_CV;
else if (!(col->type == GLP_DB && col->lb == 0.0 && col->ub == 1.0))
         kind = GLP_IV;
      else
         kind = GLP_BV;
      return kind;
}

/*----------------------------------------------------------------------
-- lpx_get_rii - retrieve row scale factor.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double lpx_get_rii(glp_prob *lp, int i);
--
-- *Returns*
--
-- The routine lpx_get_rii returns current scale factor r[i,i] for i-th
-- row of the specified problem object. */

double lpx_get_rii(glp_prob *lp, int i)
{     if (!(1 <= i && i <= lp->m))
         fault("lpx_get_rii: i = %d; row number out of range", i);
      return lp->row[i]->rii;
}

/*----------------------------------------------------------------------
-- lpx_get_sjj - retrieve column scale factor.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double lpx_get_sjj(glp_prob *lp, int j);
--
-- *Returns*
--
-- The routine lpx_get_sjj returns current scale factor s[j,j] for j-th
-- column of the specified problem object. */

double lpx_get_sjj(glp_prob *lp, int j)
{     if (!(1 <= j && j <= lp->n))
         fault("lpx_get_sjj: j = %d; column number out of range", j);
      return lp->col[j]->sjj;
}

/*----------------------------------------------------------------------
-- lpx_is_b_avail - check if LP basis is available.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_is_b_avail(glp_prob *lp);
--
-- *Returns*
--
-- If the LP basis associated with the specified problem object exists
-- and therefore available for computations, the routine lpx_is_b_avail
-- returns non-zero. Otherwise, if the LP basis is not available, the
-- routine returns zero. */

int lpx_is_b_avail(glp_prob *lp)
{
      return lp->valid;
}

/*----------------------------------------------------------------------
-- lpx_get_b_info - retrieve LP basis information.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_get_b_info(glp_prob *lp, int i);
--
-- *Returns*
--
-- The routine lpx_get_b_info returns the ordinal number k of auxiliary
-- (1 <= k <= m) or structural (m+1 <= k <= m+n) variable, which is
-- basic variable xB[i], 1 <= i <= m, in the current basis associated
-- with the specified problem object, where m is the number of rows and
-- n is the number of columns. */

int lpx_get_b_info(glp_prob *lp, int i)
{     if (!lpx_is_b_avail(lp))
         fault("lpx_get_b_info: LP basis is not available");
      if (!(1 <= i && i <= lp->m))
         fault("lpx_get_b_info: i = %d; index out of range", i);
      return lp->basis[i];
}

/*----------------------------------------------------------------------
-- lpx_get_row_b_ind - retrieve row index in LP basis.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_get_row_b_ind(glp_prob *lp, int i);
--
-- *Returns*
--
-- The routine lpx_get_row_b_ind returns the index k of basic variable
-- xB[k], 1 <= k <= m, which is the auxiliary variable associated with
-- i-th row in the current basis for the specified problem object.
-- However, if the auxiliary variable is non-basic, the routine returns
-- zero. */

int lpx_get_row_b_ind(glp_prob *lp, int i)
{     if (!lpx_is_b_avail(lp))
         fault("lpx_get_row_b_ind: LP basis is not available");
      if (!(1 <= i && i <= lp->m))
         fault("lpx_get_row_b_ind: i = %d; row number out of range", i);
      return lp->row[i]->bind;
}

/*----------------------------------------------------------------------
-- lpx_get_col_b_ind - retrieve column index in LP basis.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_get_col_b_ind(glp_prob *lp, int j);
--
-- *Returns*
--
-- The routine lpx_get_col_b_ind returns the index k of basic variable
-- xB[k], 1 <= k <= m, which is the structural variable associated with
-- j-th column in the current basis for the specified problem object.
-- However, if the structural variable is non-basic, the routine returns
-- zero. */

int lpx_get_col_b_ind(glp_prob *lp, int j)
{     if (!lpx_is_b_avail(lp))
         fault("lpx_get_col_b_ind: LP basis is not available");
      if (!(1 <= j && j <= lp->n))
         fault("lpx_get_col_b_ind: j = %d; column number out of range",
            j);
      return lp->col[j]->bind;
}

/*----------------------------------------------------------------------
-- lpx_access_inv - access factorization of basis matrix.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- BFI *lpx_access_inv(glp_prob *lp);
--
-- *Returns*
--
-- The routine lpx_access_inv returns a pointer to a factorization of
-- the basis matrix for the specified problem object.
--
-- NOTE: This routine is intended for internal use only. */

BFI *lpx_access_inv(glp_prob *lp)
{     if (lp->binv == NULL)
      {  lp->binv = bfi_create_binv();
         lp->binv->type = lpx_get_int_parm(lp, LPX_K_BFTYPE);
      }
      return lp->binv;
}

/***********************************************************************
*  NAME
*
*  glp_get_status - retrieve generic status of basic solution
*
*  SYNOPSIS
*
*  int glp_get_status(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_status reports the generic status of the basic
*  solution for the specified problem object as follows:
*
*  GLP_OPT    - solution is optimal;
*  GLP_FEAS   - solution is feasible;
*  GLP_INFEAS - solution is infeasible;
*  GLP_NOFEAS - problem has no feasible solution;
*  GLP_UNBND  - problem has unbounded solution;
*  GLP_UNDEF  - solution is undefined. */

int glp_get_status(glp_prob *lp)
{     int status;
      status = glp_get_prim_stat(lp);
      switch (status)
      {  case GLP_FEAS:
            switch (glp_get_dual_stat(lp))
            {  case GLP_FEAS:
                  status = GLP_OPT;
                  break;
               case GLP_NOFEAS:
                  status = GLP_UNBND;
                  break;
               case GLP_UNDEF:
               case GLP_INFEAS:
                  status = status;
                  break;
               default:
                  xassert(lp != lp);
            }
            break;
         case GLP_UNDEF:
         case GLP_INFEAS:
         case GLP_NOFEAS:
            status = status;
            break;
         default:
            xassert(lp != lp);
      }
      return status;
}

/***********************************************************************
*  NAME
*
*  glp_get_prim_stat - retrieve status of primal basic solution
*
*  SYNOPSIS
*
*  int glp_get_prim_stat(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_prim_stat reports the status of the primal basic
*  solution for the specified problem object as follows:
*
*  GLP_UNDEF  - primal solution is undefined;
*  GLP_FEAS   - primal solution is feasible;
*  GLP_INFEAS - primal solution is infeasible;
*  GLP_NOFEAS - no primal feasible solution exists. */

int glp_get_prim_stat(glp_prob *lp)
{     int pbs_stat = lp->pbs_stat;
      return pbs_stat;
}

/***********************************************************************
*  NAME
*
*  glp_get_dual_stat - retrieve status of dual basic solution
*
*  SYNOPSIS
*
*  int glp_get_dual_stat(glp_prob *lp);
*
*  RETURNS
*
*  The routine glp_get_dual_stat reports the status of the dual basic
*  solution for the specified problem object as follows:
*
*  GLP_UNDEF  - dual solution is undefined;
*  GLP_FEAS   - dual solution is feasible;
*  GLP_INFEAS - dual solution is infeasible;
*  GLP_NOFEAS - no dual feasible solution exists. */

int glp_get_dual_stat(glp_prob *lp)
{     int dbs_stat = lp->dbs_stat;
      return dbs_stat;
}

/*----------------------------------------------------------------------
-- glp_get_obj_val - retrieve objective value (basic solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_get_obj_val(glp_prob *lp);
--
-- *Returns*
--
-- The routine glp_get_obj_val returns value of the objective function
-- for basic solution. */

double glp_get_obj_val(glp_prob *lp)
{     struct LPXCPS *cps = lp->cps;
      double z;
      z = lp->obj_val;
      if (cps->round && fabs(z) < 1e-9) z = 0.0;
      return z;
}

/*----------------------------------------------------------------------
-- glp_get_row_stat - retrieve row status (basic solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int glp_get_row_stat(glp_prob *lp, int i);
--
-- *Returns*
--
-- The routine glp_get_row_stat returns current status assigned to the
-- auxiliary variable associated with i-th row as follows:
--
-- GLP_BS - basic variable;
-- GLP_NL - non-basic variable on its lower bound;
-- GLP_NU - non-basic variable on its upper bound;
-- GLP_NF - non-basic free (unbounded) variable;
-- GLP_NS - non-basic fixed variable. */

int glp_get_row_stat(glp_prob *lp, int i)
{     if (!(1 <= i && i <= lp->m))
         fault("glp_get_row_stat: i = %d; row number out of range", i);
      return lp->row[i]->stat;
}

/*----------------------------------------------------------------------
-- glp_get_row_prim - retrieve row primal value (basic solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_get_row_prim(glp_prob *lp, int i);
--
-- *Returns*
--
-- The routine glp_get_row_prim returns primal value of the auxiliary
-- variable associated with i-th row. */

double glp_get_row_prim(glp_prob *lp, int i)
{     struct LPXCPS *cps = lp->cps;
      double prim;
      if (!(1 <= i && i <= lp->m))
         fault("glp_get_row_prim: i = %d; row number out of range", i);
      prim = lp->row[i]->prim;
      if (cps->round && fabs(prim) < 1e-9) prim = 0.0;
      return prim;
}

/*----------------------------------------------------------------------
-- glp_get_row_dual - retrieve row dual value (basic solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_get_row_dual(glp_prob *lp, int i);
--
-- *Returns*
--
-- The routine glp_get_row_dual returns dual value (i.e. reduced cost)
-- of the auxiliary variable associated with i-th row. */

double glp_get_row_dual(glp_prob *lp, int i)
{     struct LPXCPS *cps = lp->cps;
      double dual;
      if (!(1 <= i && i <= lp->m))
         fault("glp_get_row_dual: i = %d; row number out of range", i);
      dual = lp->row[i]->dual;
      if (cps->round && fabs(dual) < 1e-9) dual = 0.0;
      return dual;
}

/*----------------------------------------------------------------------
-- glp_get_col_stat - retrieve column status (basic solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int glp_get_col_stat(glp_prob *lp, int j);
--
-- *Returns*
--
-- The routine glp_get_col_stat returns current status assigned to the
-- structural variable associated with j-th column as follows:
--
-- GLP_BS - basic variable;
-- GLP_NL - non-basic variable on its lower bound;
-- GLP_NU - non-basic variable on its upper bound;
-- GLP_NF - non-basic free (unbounded) variable;
-- GLP_NS - non-basic fixed variable. */

int glp_get_col_stat(glp_prob *lp, int j)
{     if (!(1 <= j && j <= lp->n))
         fault("glp_get_col_stat: j = %d; column number out of range",
            j);
      return lp->col[j]->stat;
}

/*----------------------------------------------------------------------
-- glp_get_col_prim - retrieve column primal value (basic solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_get_col_prim(glp_prob *lp, int j);
--
-- *Returns*
--
-- The routine glp_get_col_prim returns primal value of the structural
-- variable associated with j-th column. */

double glp_get_col_prim(glp_prob *lp, int j)
{     struct LPXCPS *cps = lp->cps;
      double prim;
      if (!(1 <= j && j <= lp->n))
         fault("glp_get_col_prim: j = %d; column number out of range",
            j);
      prim = lp->col[j]->prim;
      if (cps->round && fabs(prim) < 1e-9) prim = 0.0;
      return prim;
}

/*----------------------------------------------------------------------
-- glp_get_col_dual - retrieve column dual value (basic solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_get_col_dual(glp_prob *lp, int j);
--
-- *Returns*
--
-- The routine glp_get_col_dual returns dual value (i.e. reduced cost)
-- of the structural variable associated with j-th column. */

double glp_get_col_dual(glp_prob *lp, int j)
{     struct LPXCPS *cps = lp->cps;
      double dual;
      if (!(1 <= j && j <= lp->n))
         fault("glp_get_col_dual: j = %d; column number out of range",
            j);
      dual = lp->col[j]->dual;
      if (cps->round && fabs(dual) < 1e-9) dual = 0.0;
      return dual;
}

/*----------------------------------------------------------------------
-- lpx_get_ray_info - retrieve row/column which causes unboundness.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_get_ray_info(glp_prob *lp);
--
-- *Returns*
--
-- The routine lpx_get_ray_info returns the number k of some non-basic
-- variable x[k], which causes primal unboundness. If such a variable
-- cannot be identified, the routine returns zero.
--
-- If 1 <= k <= m, x[k] is the auxiliary variable associated with k-th
-- row, if m+1 <= k <= m+n, x[k] is the structural variable associated
-- with (k-m)-th column, where m is the number of rows, n is the number
-- of columns in the LP problem object. */

int lpx_get_ray_info(glp_prob *lp)
{     int k;
      k = lp->some;
      return k;
}

/*----------------------------------------------------------------------
-- lpx_ipt_status - retrieve status of interior-point solution.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_ipt_status(glp_prob *lp);
--
-- *Returns*
--
-- The routine lpx_ipt_status reports the status of solution found by
-- interior-point solver as follows:
--
-- LPX_T_UNDEF - interior-point solution is undefined;
-- LPX_T_OPT   - interior-point solution is optimal. */

int lpx_ipt_status(glp_prob *lp)
{     int t_stat;
      switch (lp->ipt_stat)
      {  case GLP_UNDEF: t_stat = LPX_T_UNDEF; break;
         case GLP_OPT:   t_stat = LPX_T_OPT;   break;
         default:        xassert(lp != lp);
      }
      return t_stat;
}

/*----------------------------------------------------------------------
-- glp_ipt_obj_val - retrieve objective value (interior point).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_ipt_obj_val(glp_prob *lp);
--
-- *Returns*
--
-- The routine glp_ipt_obj_val returns value of the objective function
-- for interior-point solution. */

double glp_ipt_obj_val(glp_prob *lp)
{     struct LPXCPS *cps = lp->cps;
      double z;
      z = lp->ipt_obj;
      if (cps->round && fabs(z) < 1e-9) z = 0.0;
      return z;
}

/*----------------------------------------------------------------------
-- glp_ipt_row_prim - retrieve row primal value (interior point).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_ipt_row_prim(glp_prob *lp, int i);
--
-- *Returns*
--
-- The routine glp_ipt_row_prim returns primal value of the auxiliary
-- variable associated with i-th row. */

double glp_ipt_row_prim(glp_prob *lp, int i)
{     struct LPXCPS *cps = lp->cps;
      double pval;
      if (!(1 <= i && i <= lp->m))
         fault("glp_ipt_row_prim: i = %d; row number out of range", i);
      pval = lp->row[i]->pval;
      if (cps->round && fabs(pval) < 1e-9) pval = 0.0;
      return pval;
}

/*----------------------------------------------------------------------
-- glp_ipt_row_dual - retrieve row dual value (interior point).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_ipt_row_dual(glp_prob *lp, int i);
--
-- *Returns*
--
-- The routine glp_ipt_row_dual returns dual value (i.e. reduced cost)
-- of the auxiliary variable associated with i-th row. */

double glp_ipt_row_dual(glp_prob *lp, int i)
{     struct LPXCPS *cps = lp->cps;
      double dval;
      if (!(1 <= i && i <= lp->m))
         fault("glp_ipt_row_dual: i = %d; row number out of range", i);
      dval = lp->row[i]->dval;
      if (cps->round && fabs(dval) < 1e-9) dval = 0.0;
      return dval;
}

/*----------------------------------------------------------------------
-- glp_ipt_col_prim - retrieve column primal value (interior point).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_ipt_col_prim(glp_prob *lp, int j);
--
-- *Returns*
--
-- The routine glp_ipt_col_prim returns primal value of the structural
-- variable associated with j-th column. */

double glp_ipt_col_prim(glp_prob *lp, int j)
{     struct LPXCPS *cps = lp->cps;
      double pval;
      if (!(1 <= j && j <= lp->n))
         fault("glp_ipt_col_prim: j = %d; column number out of range",
            j);
      pval = lp->col[j]->pval;
      if (cps->round && fabs(pval) < 1e-9) pval = 0.0;
      return pval;
}

/*----------------------------------------------------------------------
-- glp_ipt_col_dual - retrieve column dual value (interior point).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_ipt_col_dual(glp_prob *lp, int j);
--
-- *Returns*
--
-- The routine glp_ipt_col_dual returns dual value (i.e. reduced cost)
-- of the structural variable associated with j-th column. */

double glp_ipt_col_dual(glp_prob *lp, int j)
{     struct LPXCPS *cps = lp->cps;
      double dval;
      if (!(1 <= j && j <= lp->n))
         fault("glp_ipt_col_dual: j = %d; column number out of range",
            j);
      dval = lp->col[j]->dval;
      if (cps->round && fabs(dval) < 1e-9) dval = 0.0;
      return dval;
}

/*----------------------------------------------------------------------
-- lpx_mip_status - retrieve status of MIP solution.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_mip_status(glp_prob *lp);
--
-- *Returns*
--
-- The routine lpx_mip_status reports the status of MIP solution found
-- by branch-and-bound solver as follows:
--
-- LPX_I_UNDEF  - MIP solution is undefined;
-- LPX_I_OPT    - MIP solution is integer optimal;
-- LPX_I_FEAS   - MIP solution is integer feasible but its optimality
--                (or non-optimality) has not been proven, perhaps due
--                to premature termination of the search;
-- LPX_I_NOFEAS - problem has no integer feasible solution (proven by
--                the solver). */

int lpx_mip_status(glp_prob *lp)
{     int i_stat;
#if 0
      if (lp->klass != LPX_MIP)
         fault("lpx_mip_status: not a MIP problem");
#endif
      switch (lp->mip_stat)
      {  case GLP_UNDEF:   i_stat = LPX_I_UNDEF; break;
         case GLP_OPT:     i_stat = LPX_I_OPT;  break;
         case GLP_FEAS:    i_stat = LPX_I_FEAS; break;
         case GLP_NOFEAS:  i_stat = LPX_I_NOFEAS; break;
         default:
            xassert(lp != lp);
      }
      return i_stat;
}

/*----------------------------------------------------------------------
-- glp_mip_obj_val - retrieve objective value (MIP solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_mip_obj_val(glp_prob *lp);
--
-- *Returns*
--
-- The routine glp_mip_obj_val returns value of the objective function
-- for MIP solution. */

double glp_mip_obj_val(glp_prob *lp)
{     struct LPXCPS *cps = lp->cps;
      double z;
      z = lp->mip_obj;
      if (cps->round && fabs(z) < 1e-9) z = 0.0;
      return z;
}

/*----------------------------------------------------------------------
-- glp_mip_row_val - retrieve row value (MIP solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_mip_row_val(glp_prob *lp, int i);
--
-- *Returns*
--
-- The routine glp_mip_row_val returns value of the auxiliary variable
-- associated with i-th row. */

double glp_mip_row_val(glp_prob *lp, int i)
{     struct LPXCPS *cps = lp->cps;
      double mipx;
      if (!(1 <= i && i <= lp->m))
         fault("glp_mip_row_val: i = %d; row number out of range", i);
      mipx = lp->row[i]->mipx;
      if (cps->round && fabs(mipx) < 1e-9) mipx = 0.0;
      return mipx;
}

/*----------------------------------------------------------------------
-- glp_mip_col_val - retrieve column value (MIP solution).
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double glp_mip_col_val(glp_prob *lp, int j);
--
-- *Returns*
--
-- The routine glp_mip_col_val returns value of the structural variable
-- associated with j-th column. */

double glp_mip_col_val(glp_prob *lp, int j)
{     struct LPXCPS *cps = lp->cps;
      double mipx;
      if (!(1 <= j && j <= lp->n))
         fault("glp_mip_col_val: j = %d; column number out of range",
            j);
      mipx = lp->col[j]->mipx;
      if (cps->round && fabs(mipx) < 1e-9) mipx = 0.0;
      return mipx;
}

/* eof */
