/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nufeb/growth/methanogen,FixGrowthMethanogen)

#else

#ifndef LMP_FIX_GROWTH_METHANOGEN_H
#define LMP_FIX_GROWTH_METHANOGEN_H

#include "fix_growth.h"

namespace LAMMPS_NS {

    class FixGrowthMethanogen: public FixGrowth {
    public:
        FixGrowthMethanogen(class LAMMPS *, int, char **);
        virtual ~FixGrowthMethanogen() {}

        virtual void update_atoms();
        virtual void update_cells();

    protected:
        int ih2;
        int ico2;
        int ich4;

        double h2_affinity;
        double co2_affinity;

        double growth;
        double yield;
        double maintain;
        double decay;
        double eps_yield;
        double eps_dens;
    };

}

#endif
#endif

/* ERROR/WARNING messages:
*/
