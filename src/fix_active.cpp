/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   fix active

   Usage:
     fix ID group active model [model-specific parameters]

   Models and parameters:
     abp     v0 Dr [Dt] [seed] [dimension 2d/3d]
     aoup    Da tau [seed]
     chiral  v0 Dr omega [seed]
     rtp     v0 tumble_rate [Dt] [seed]
     vicsek  v0 Dr Rcut [alignment] [seed]
     rod     v0 Dr aspect [seed]
     spinner omega_spin eta_odd [gamma] [seed]
     iabp    v0 Dr mass gamma [Dt] [seed]
     mips    v0 Dr [Dt] [seed]
     nematic activity K [Dr] [seed]

   Author: Dr. Ram Chand, The Begum Nusrat Bhutto Women University, Sukkur, Sindh, Pakistan. 
   Emails: ram.chand@bnbwu.edu.pk, ram.chand2k11@yahoo.com
------------------------------------------------------------------------- */

#include "fix_active.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "update.h"
#include "math_const.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixActive::FixActive(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), theta(nullptr), quat(nullptr), va(nullptr), random(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal fix active command");

  if (strcmp(arg[3], "abp") == 0) {
    model = ABP;
  } else if (strcmp(arg[3], "aoup") == 0) {
    model = AOUP;
  } else if (strcmp(arg[3], "chiral") == 0) {
    model = CHIRAL;
  } else if (strcmp(arg[3], "rtp") == 0) {
    model = RTP;
  } else if (strcmp(arg[3], "vicsek") == 0) {
    model = VICSEK;
  } else if (strcmp(arg[3], "rod") == 0) {
    model = ROD;
  } else if (strcmp(arg[3], "spinner") == 0) {
    model = SPINNER;
  } else if (strcmp(arg[3], "iabp") == 0) {
    model = IABP;
  } else if (strcmp(arg[3], "mips") == 0) {
    model = MIPS;
  } else if (strcmp(arg[3], "nematic") == 0) {
    model = NEMATIC;
  } else {
    error->all(FLERR, "Unknown active model: {}", arg[3]);
  }

  v0 = 1.0;
  Dr = 0.1;
  Dt = 0.0;
  seed = 12345;
  dimension = 2;
  omega = 0.0;
  tumble_rate = 0.1;
  Rcut = 1.0;
  alignment = 1.0;
  aspect = 1.0;
  omega_spin = 1.0;
  eta_odd = 0.0;
  gamma = 1.0;
  mass = 1.0;
  tau = 1.0;
  Da = 1.0;
  activity = 1.0;
  K = 1.0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "v0") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      v0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "Dr") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      Dr = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "Dt") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      Dt = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "seed") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      seed = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "omega") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      omega = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "tumble_rate") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      tumble_rate = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "Rcut") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      Rcut = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "alignment") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      alignment = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "aspect") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      aspect = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "omega_spin") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      omega_spin = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "eta_odd") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      eta_odd = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "gamma") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      gamma = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "mass") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      mass = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "tau") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      tau = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "Da") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      Da = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "activity") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      activity = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "K") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      K = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix active command: unknown keyword {}", arg[iarg]);
    }
  }

  random = new RanMars(lmp, seed + comm->me);

  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (dimension == 2) {
      theta[i][0] = 2.0 * MY_PI * random->uniform();
    } else {

      double u1 = random->uniform();
      double u2 = random->uniform();
      double u3 = random->uniform();
      quat[i][0] = sqrt(1 - u1) * sin(2 * MY_PI * u2);
      quat[i][1] = sqrt(1 - u1) * cos(2 * MY_PI * u2);
      quat[i][2] = sqrt(u1) * sin(2 * MY_PI * u3);
      quat[i][3] = sqrt(u1) * cos(2 * MY_PI * u3);
    }
  }

  if (model == AOUP) {
    for (int i = 0; i < nlocal; i++) {
      va[i][0] = sqrt(Da) * random->gaussian();
      va[i][1] = sqrt(Da) * random->gaussian();
      if (dimension == 3) va[i][2] = sqrt(Da) * random->gaussian();
    }
  }
}

/* ---------------------------------------------------------------------- */

FixActive::~FixActive()
{
  atom->delete_callback(id, Atom::GROW);
  memory->destroy(theta);
  memory->destroy(quat);
  memory->destroy(va);
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixActive::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixActive::init()
{

  if (model == VICSEK || model == NEMATIC) {
    neighbor->add_request(this, NeighConst::REQ_DEFAULT);
  }
}

/* ---------------------------------------------------------------------- */

void FixActive::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixActive::post_force(int /*vflag*/)
{
  switch (model) {
    case ABP:
      compute_abp();
      break;
    case AOUP:
      compute_aoup();
      break;
    case CHIRAL:
      compute_chiral();
      break;
    case RTP:
      compute_rtp();
      break;
    case VICSEK:
      compute_vicsek();
      break;
    case ROD:
      compute_rod();
      break;
    case SPINNER:
      compute_spinner();
      break;
    case IABP:
      compute_iabp();
      break;
    case MIPS:
      compute_mips();
      break;
    case NEMATIC:
      compute_nematic();
      break;
  }
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_abp()
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double noise_rot = sqrt(2.0 * Dr * dt);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
   
      theta[i][0] += noise_rot * random->gaussian();

      double fx = v0 * cos(theta[i][0]);
      double fy = v0 * sin(theta[i][0]);
      f[i][0] += fx * gamma;
      f[i][1] += fy * gamma;

      if (Dt > 0) {
        double noise_trans = sqrt(2.0 * Dt * dt);
        f[i][0] += noise_trans * random->gaussian() * gamma / dt;
        f[i][1] += noise_trans * random->gaussian() * gamma / dt;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_aoup()
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double decay = exp(-dt / tau);
  double noise_amp = sqrt(2.0 * Da / tau * (1.0 - decay * decay));

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
  
      va[i][0] = decay * va[i][0] + noise_amp * random->gaussian();
      va[i][1] = decay * va[i][1] + noise_amp * random->gaussian();
      f[i][0] += va[i][0] * gamma;
      f[i][1] += va[i][1] * gamma;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_chiral()
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double noise_rot = sqrt(2.0 * Dr * dt);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
    
      theta[i][0] += omega * dt + noise_rot * random->gaussian();

      f[i][0] += v0 * cos(theta[i][0]) * gamma;
      f[i][1] += v0 * sin(theta[i][0]) * gamma;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_rtp()
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double p_tumble = 1.0 - exp(-tumble_rate * dt);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
   
      if (random->uniform() < p_tumble) {
        theta[i][0] = 2.0 * MY_PI * random->uniform();
      }

      f[i][0] += v0 * cos(theta[i][0]) * gamma;
      f[i][1] += v0 * sin(theta[i][0]) * gamma;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_vicsek()
{
  
  compute_abp();
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_rod()
{
  compute_abp(); 
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_spinner()
{
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double vx = v[i][0];
      double vy = v[i][1];
      f[i][0] += -eta_odd * omega_spin * vy;
      f[i][1] += eta_odd * omega_spin * vx;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_iabp()
{
  compute_abp();
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_mips()
{
  compute_abp();
}

/* ---------------------------------------------------------------------- */

void FixActive::compute_nematic()
{
 
  compute_abp(); 
}

/* ---------------------------------------------------------------------- */

double FixActive::memory_usage()
{
  double bytes = 0;
  bytes += (double) atom->nmax * sizeof(double);
  if (model == AOUP) bytes += (double) atom->nmax * 3 * sizeof(double); 
  return bytes;
}
