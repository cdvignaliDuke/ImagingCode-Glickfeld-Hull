/* glplpx03.c (control parameters and statistics routines) */

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
#include "glpapi.h"
#define fault xfault1

/*----------------------------------------------------------------------
-- lpx_reset_parms - reset control parameters to default values.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- void lpx_reset_parms(LPX *lp);
--
-- *Description*
--
-- The routine lpx_reset_parms resets all control parameters associated
-- with an LP problem object, which the parameter lp points to, to their
-- default values. */

void lpx_reset_parms(LPX *lp)
{     struct LPXCPS *cps = lp->cps;
      cps->msg_lev  = 3;
      cps->scale    = 1;
      cps->dual     = 0;
      cps->price    = 1;
      cps->relax    = 0.07;
      cps->tol_bnd  = 1e-7;
      cps->tol_dj   = 1e-7;
      cps->tol_piv  = 1e-9;
      cps->round    = 0;
      cps->obj_ll   = -DBL_MAX;
      cps->obj_ul   = +DBL_MAX;
      cps->it_lim   = -1;
      cps->it_cnt   = 0;
      cps->tm_lim   = -1.0;
      cps->out_frq  = 200;
      cps->out_dly  = 0.0;
      cps->branch   = 2;
      cps->btrack   = 2;
      cps->tol_int  = 1e-5;
      cps->tol_obj  = 1e-7;
      cps->mps_info = 1;
      cps->mps_obj  = 2;
      cps->mps_orig = 0;
      cps->mps_wide = 1;
      cps->mps_free = 0;
      cps->mps_skip = 0;
      cps->lpt_orig = 0;
      cps->presol = 0;
      cps->binarize = 0;
      cps->use_cuts = 0;
      cps->bf_type = 1;
      return;
}

/*----------------------------------------------------------------------
-- lpx_set_int_parm - set (change) integer control parameter.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- void lpx_set_int_parm(LPX *lp, int parm, int val);
--
-- *Description*
--
-- The routine lpx_set_int_parm sets (changes) the current value of an
-- integer control parameter parm. The parameter val specifies a new
-- value of the control parameter. */

void lpx_set_int_parm(LPX *lp, int parm, int val)
{     struct LPXCPS *cps = lp->cps;
      switch (parm)
      {  case LPX_K_MSGLEV:
            if (!(0 <= val && val <= 3))
               fault("lpx_set_int_parm: MSGLEV = %d; invalid value",
                  val);
            cps->msg_lev = val;
            break;
         case LPX_K_SCALE:
            if (!(0 <= val && val <= 3))
               fault("lpx_set_int_parm: SCALE = %d; invalid value",
                  val);
            cps->scale = val;
            break;
         case LPX_K_DUAL:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: DUAL = %d; invalid value",
                  val);
            cps->dual = val;
            break;
         case LPX_K_PRICE:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: PRICE = %d; invalid value",
                  val);
            cps->price = val;
            break;
         case LPX_K_ROUND:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: ROUND = %d; invalid value",
                  val);
            cps->round = val;
            break;
         case LPX_K_ITLIM:
            cps->it_lim = val;
            break;
         case LPX_K_ITCNT:
            cps->it_cnt = val;
            break;
         case LPX_K_OUTFRQ:
            if (!(val > 0))
               fault("lpx_set_int_parm: OUTFRQ = %d; invalid value",
                  val);
            cps->out_frq = val;
            break;
         case LPX_K_BRANCH:
            if (!(val == 0 || val == 1 || val == 2 || val == 3))
               fault("lpx_set_int_parm: BRANCH = %d; invalid value",
                  val);
            cps->branch = val;
            break;
         case LPX_K_BTRACK:
            if (!(val == 0 || val == 1 || val == 2 || val == 3))
               fault("lpx_set_int_parm: BTRACK = %d; invalid value",
                  val);
            cps->btrack = val;
            break;
         case LPX_K_MPSINFO:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: MPSINFO = %d; invalid value",
                  val);
            cps->mps_info = val;
            break;
         case LPX_K_MPSOBJ:
            if (!(val == 0 || val == 1 || val == 2))
               fault("lpx_set_int_parm: MPSOBJ = %d; invalid value",
                  val);
            cps->mps_obj = val;
            break;
         case LPX_K_MPSORIG:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: MPSORIG = %d; invalid value",
                  val);
            cps->mps_orig = val;
            break;
         case LPX_K_MPSWIDE:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: MPSWIDE = %d; invalid value",
                  val);
            cps->mps_wide = val;
            break;
         case LPX_K_MPSFREE:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: MPSFREE = %d; invalid value",
                  val);
            cps->mps_free = val;
            break;
         case LPX_K_MPSSKIP:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: MPSSKIP = %d; invalid value",
                  val);
            cps->mps_skip = val;
            break;
         case LPX_K_LPTORIG:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: LPTORIG = %d; invalid value",
                  val);
            cps->lpt_orig = val;
            break;
         case LPX_K_PRESOL:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: PRESOL = %d; invalid value",
                  val);
            cps->presol = val;
            break;
         case LPX_K_BINARIZE:
            if (!(val == 0 || val == 1))
               fault("lpx_set_int_parm: BINARIZE = %d; invalid value",
                  val);
            cps->binarize = val;
            break;
         case LPX_K_USECUTS:
            if (val & ~LPX_C_ALL)
               fault("lpx_set_int_parm: USECUTS = 0x%X; invalid value",
                  val);
            cps->use_cuts = val;
            break;
         case LPX_K_BFTYPE:
            if (!(1 <= val && val <= 3))
               fault("lpx_set_int_parm: BFTYPE = %d; invalid value",
                  val);
            cps->bf_type = val;
            break;
         default:
            fault("lpx_set_int_parm: parm = %d; invalid parameter",
               parm);
      }
      return;
}

/*----------------------------------------------------------------------
-- lpx_get_int_parm - query integer control parameter.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_get_int_parm(LPX *lp, int parm);
--
-- *Returns*
--
-- The routine lpx_get_int_parm returns the current value of an integer
-- control parameter parm. */

int lpx_get_int_parm(LPX *lp, int parm)
{     struct LPXCPS *cps = lp->cps;
      int val = 0;
      switch (parm)
      {  case LPX_K_MSGLEV:
            val = cps->msg_lev; break;
         case LPX_K_SCALE:
            val = cps->scale; break;
         case LPX_K_DUAL:
            val = cps->dual; break;
         case LPX_K_PRICE:
            val = cps->price; break;
         case LPX_K_ROUND:
            val = cps->round; break;
         case LPX_K_ITLIM:
            val = cps->it_lim; break;
         case LPX_K_ITCNT:
            val = cps->it_cnt; break;
         case LPX_K_OUTFRQ:
            val = cps->out_frq; break;
         case LPX_K_BRANCH:
            val = cps->branch; break;
         case LPX_K_BTRACK:
            val = cps->btrack; break;
         case LPX_K_MPSINFO:
            val = cps->mps_info; break;
         case LPX_K_MPSOBJ:
            val = cps->mps_obj; break;
         case LPX_K_MPSORIG:
            val = cps->mps_orig; break;
         case LPX_K_MPSWIDE:
            val = cps->mps_wide; break;
         case LPX_K_MPSFREE:
            val = cps->mps_free; break;
         case LPX_K_MPSSKIP:
            val = cps->mps_skip; break;
         case LPX_K_LPTORIG:
            val = cps->lpt_orig; break;
         case LPX_K_PRESOL:
            val = cps->presol; break;
         case LPX_K_BINARIZE:
            val = cps->binarize; break;
         case LPX_K_USECUTS:
            val = cps->use_cuts; break;
         case LPX_K_BFTYPE:
            val = cps->bf_type; break;
         default:
            fault("lpx_get_int_parm: parm = %d; invalid parameter",
               parm);
      }
      return val;
}

/*----------------------------------------------------------------------
-- lpx_set_real_parm - set (change) real control parameter.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- void lpx_set_real_parm(LPX *lp, int parm, double val);
--
-- *Description*
--
-- The routine lpx_set_real_parm sets (changes) the current value of
-- a real (floating point) control parameter parm. The parameter val
-- specifies a new value of the control parameter. */

void lpx_set_real_parm(LPX *lp, int parm, double val)
{     struct LPXCPS *cps = lp->cps;
      switch (parm)
      {  case LPX_K_RELAX:
            if (!(0.0 <= val && val <= 1.0))
               fault("lpx_set_real_parm: RELAX = %g; invalid value",
                  val);
            cps->relax = val;
            break;
         case LPX_K_TOLBND:
            if (!(DBL_EPSILON <= val && val <= 0.001))
               fault("lpx_set_real_parm: TOLBND = %g; invalid value",
                  val);
#if 0
            if (cps->tol_bnd > val)
            {  /* invalidate the basic solution */
               lp->p_stat = LPX_P_UNDEF;
               lp->d_stat = LPX_D_UNDEF;
            }
#endif
            cps->tol_bnd = val;
            break;
         case LPX_K_TOLDJ:
            if (!(DBL_EPSILON <= val && val <= 0.001))
               fault("lpx_set_real_parm: TOLDJ = %g; invalid value",
                  val);
#if 0
            if (cps->tol_dj > val)
            {  /* invalidate the basic solution */
               lp->p_stat = LPX_P_UNDEF;
               lp->d_stat = LPX_D_UNDEF;
            }
#endif
            cps->tol_dj = val;
            break;
         case LPX_K_TOLPIV:
            if (!(DBL_EPSILON <= val && val <= 0.001))
               fault("lpx_set_real_parm: TOLPIV = %g; invalid value",
                  val);
            cps->tol_piv = val;
            break;
         case LPX_K_OBJLL:
            cps->obj_ll = val;
            break;
         case LPX_K_OBJUL:
            cps->obj_ul = val;
            break;
         case LPX_K_TMLIM:
            cps->tm_lim = val;
            break;
         case LPX_K_OUTDLY:
            cps->out_dly = val;
            break;
         case LPX_K_TOLINT:
            if (!(DBL_EPSILON <= val && val <= 0.001))
               fault("lpx_set_real_parm: TOLINT = %g; invalid value",
                  val);
            cps->tol_int = val;
            break;
         case LPX_K_TOLOBJ:
            if (!(DBL_EPSILON <= val && val <= 0.001))
               fault("lpx_set_real_parm: TOLOBJ = %g; invalid value",
                  val);
            cps->tol_obj = val;
            break;
         default:
            fault("lpx_set_real_parm: parm = %d; invalid parameter",
               parm);
      }
      return;
}

/*----------------------------------------------------------------------
-- lpx_get_real_parm - query real control parameter.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- double lpx_get_real_parm(LPX *lp, int parm);
--
-- *Returns*
--
-- The routine lpx_get_real_parm returns the current value of a real
-- (floating point) control parameter parm. */

double lpx_get_real_parm(LPX *lp, int parm)
{     struct LPXCPS *cps = lp->cps;
      double val = 0.0;
      switch (parm)
      {  case LPX_K_RELAX:
            val = cps->relax;
            break;
         case LPX_K_TOLBND:
            val = cps->tol_bnd;
            break;
         case LPX_K_TOLDJ:
            val = cps->tol_dj;
            break;
         case LPX_K_TOLPIV:
            val = cps->tol_piv;
            break;
         case LPX_K_OBJLL:
            val = cps->obj_ll;
            break;
         case LPX_K_OBJUL:
            val = cps->obj_ul;
            break;
         case LPX_K_TMLIM:
            val = cps->tm_lim;
            break;
         case LPX_K_OUTDLY:
            val = cps->out_dly;
            break;
         case LPX_K_TOLINT:
            val = cps->tol_int;
            break;
         case LPX_K_TOLOBJ:
            val = cps->tol_obj;
            break;
         default:
            fault("lpx_get_real_parm: parm = %d; invalid parameter",
               parm);
      }
      return val;
}

/* eof */
