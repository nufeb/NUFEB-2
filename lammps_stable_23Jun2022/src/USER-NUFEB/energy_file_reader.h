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
      virtual void read_file(char *group, char *fix_name);

      double uptake;         // substrate uptake rate
      double decay;	         // decay rate
      double dissipation;    // dissipation energy
      double biomass_gibbs;  // biomass Gibbs energy
      double *sub_gibbs;     // substrate Gibbs energy
      double *ks_coeff;      // ks coeffs
      double *cata_coeff;    // catabolic coeffs
      double *anab_coeff;    // anabolic coeffs
      double *decay_coeff;   // decay coeffs

      double **sstc_gibbs;   // Gibbs energy for the substrate in 5 forms
                             // 0:Non-Hydrated, 1:Hydrated, 2:1st Deprotonation
                             // 3:2nd Deprotonation, 4:3rd Deprotonation
      int **ncharges;         // charge numbers for substrate in 5 forms
      int *form_id;         // substrate form for microbial utilisation

      // not used
      double maintain;       // maintenance rate
      double max_yield;      // calculated yield
      int e_donor;           // electron donor

  private:
      int me, compressed;
      char *line, *keyword, *buffer, *style;
      FILE *fp;
      char *filename;
      char **arg;
      int narg, maxarg;

      int ngroups;    // # of groups in data file

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

      void substance_energy();
      void set_substance_energy(int, char **);

      void charge_numbers();
      void set_charge_numbers(int, char **);

      void uptake_form();
      void set_uptake_form(int, char **);
  };

} // namespace LAMMPS_NS



#endif
