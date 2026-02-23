/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(active, FixActive);

#else

#ifndef LMP_FIX_ACTIVE_H
#define LMP_FIX_ACTIVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixActive : public Fix {
 public:
  FixActive(class LAMMPS *, int, char **);
  ~FixActive() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  double memory_usage() override;

 private:
  
  enum ActiveModel {
    ABP,
    AOUP,
    CHIRAL,
    RTP,
    VICSEK,
    ROD,
    SPINNER,
    IABP,
    MIPS,
    NEMATIC
  };
  ActiveModel model;

  
  double v0;          
  double Dr;          
  double Dt;           
  int seed;           
  int dimension;       

  
  double omega;        
  double tumble_rate; 
  double Rcut;         
  double alignment;    
  double aspect;       
  double omega_spin;   
  double eta_odd;      
  double gamma;        
  double mass;         
  double tau;          
  double Da;           
  double activity;     
  double K;           

  
  double **theta;      
  double **quat;       
  double **va;         

  class RanMars *random;

  void compute_abp();
  void compute_aoup();
  void compute_chiral();
  void compute_rtp();
  void compute_vicsek();
  void compute_rod();
  void compute_spinner();
  void compute_iabp();
  void compute_mips();
  void compute_nematic();

  void update_orientation_2d(int);
  void update_orientation_3d(int);
  void tumble(int);
  void align_with_neighbors(int);
};

} 

#endif
#endif
