/*
 * energy_file_reader.h
 *
 *  Created on: 16 Feb 2023
 *      Author: Bowen Li
 */

#ifndef LMP_ENERGY_FILE_READER_H
#define LMP_ENERGY_FILE_READER_H

#include "pointers.h"

namespace LAMMPS_NS {
  class EnergyFileReader : public Pointers {

  public:
      EnergyFileReader(class LAMMPS *lmp, char *filename);
      virtual ~EnergyFileReader();
      virtual void read_file(char *group);

      double uptake;         // substrate uptake rate
      double decay;	         // decay rate
      double maintain;       // maintenance rate
      double dissipation;    // dissipation energy
      double biomass_gibbs;  // biomass Gibbs energy
      double max_yield;      // calculated yield
      int e_donor;           // electron donor
      double *sub_gibbs;     // substrate Gibbs energy
      double *ks_coeff;      // ks coeffs
      double *cata_coeff;    // catabolic coeffs
      double *anab_coeff;    // anabolic coeffs
      double *decay_coeff;   // decay coeffs

  private:
      int me, compressed;
      char *line, *keyword, *buffer, *style;
      FILE *fp;
      char *filename;
      char **arg;
      int narg, maxarg;

      int ngroups;    // # of groups and substrates in data file

      void open(char *);
      void header(int);
      void parse_keyword(int);
      void skip_lines(bigint);
      void parse_coeffs(char *);

      void uptake_rate(char *);
      int set_uptake(char *, char *);

      void calc_yield(char *);
      int set_yield(char *, char *);

      void electron_donor(char *);
      int set_edonor(char **, char *);

      void decay_rate(char *);
      int set_decay(char *, char *);

      void maintain_rate(char *);
      int set_maintain(char *, char *);

      void ks_coeffs(char *);
      int set_ks(int, char **, char *);

      void cata_coeffs(char *);
      int set_cata(int, char **, char *);

      void anab_coeffs(char *);
      int set_anab(int, char **, char *);

      void decay_coeffs(char *);
      int set_decay_ceoffs(int, char **, char *);

      void biomass_energy(char *);
      int set_biomass_energy(char *, char *);

      void substrate_energy();
      void set_substrate_energy(char *, int);

      void dissipation_energy(char *);
      int set_dissipation(char *, char *);
  };

} // namespace LAMMPS_NS



#endif
