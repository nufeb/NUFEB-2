/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include "fix_T6SS_contact.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_park.h"
#include "pair.h"
#include "update.h"
#include "grid.h"
#include "grid_masks.h"
#include <math.h>
#include <algorithm>
#include <string.h>
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixT6SSContact::FixT6SSContact(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg),effector_groups(nullptr){
  t6ss_types = nullptr;
  attacker_effectors = nullptr;
  harpoon_lens = nullptr;
  cooldowns = nullptr;
  susceptible_types = nullptr;
  susceptible_effectors = nullptr;
  effector_success_probs = nullptr;
  effector_effects = nullptr;

  if (narg < 4)
    error->all(FLERR, "Illegal fix nufeb/T6SS/contact command");

  // Initialize random number generator, same for all procs
  int seed = utils::inumeric(FLERR, arg[3], true, lmp);
  random = new RanPark(lmp, seed);

  n_t6ss_types = utils::inumeric(FLERR, arg[4], true, lmp);
  if (n_t6ss_types < 0){
      error->all(FLERR, "Specified fewer than 0 T6SS bugs in call to fix/T6SS/contact");
  }
  else {
      printf("%d types of bugs are T6SS capable. Parsing.\n",n_t6ss_types);
  }

  t6ss_types = memory->create(t6ss_types,n_t6ss_types,"fix nufeb/T6SS/contact: t6ss_types");
  attacker_effectors = memory->create(attacker_effectors,n_t6ss_types,"fix nufeb/T6SS/contact: attacker_effectors");
  harpoon_lens = memory->create(harpoon_lens,n_t6ss_types,"fix nufeb/T6SS/contact: harpoon_lens");
  cooldowns = memory->create(cooldowns,n_t6ss_types,"fix nufeb/T6SS/contact: cooldowns");

  int arg_n = 5;
  for(int i = 0; i < n_t6ss_types; i++){
    t6ss_types[i] = utils::inumeric(FLERR, arg[arg_n++], true, lmp);
    
    attacker_lookup[t6ss_types[i]]=i;
    // TODO validate bug type exists
    printf("Bug %d type code: %d\n",i,t6ss_types[i]);

    //effector id
    attacker_effectors[i] = utils::inumeric(FLERR, arg[arg_n++], true, lmp);
    printf("\t Effector ID: %d\n",attacker_effectors[i]);

    //harpoon length
    harpoon_lens[i] = utils::numeric(FLERR, arg[arg_n++], true, lmp);
    printf("\t harpoonLen: %7g\n",harpoon_lens[i]);

    //cooldown
    cooldowns[i] = utils::inumeric(FLERR, arg[arg_n++], true, lmp);
    printf("\t Cooldown ID: %d\n",cooldowns[i]);
  }

  n_susceptible_types = utils::inumeric(FLERR, arg[arg_n++], true, lmp);
  printf("%d types of bugs are susceptible to T6SS effectors. Parsing:\n", n_susceptible_types);

  susceptible_types = memory->create(susceptible_types,n_t6ss_types,"fix nufeb/T6SS/contact: susceptible_types");
  susceptible_effectors = memory->create(susceptible_effectors,n_t6ss_types,"fix nufeb/T6SS/contact: susceptible_effectors");
  effector_success_probs = memory->create(effector_success_probs,n_t6ss_types,"fix nufeb/T6SS/contact: effector_success_probs");
  effector_groups = memory->create(effector_groups,n_t6ss_types,"fix nufeb/T6SS/contact: effector_groups");
  effector_effects = memory->create(effector_effects,n_t6ss_types,"fix nufeb/T6SS/contact: effector_effects");

  for(int i = 0; i < n_susceptible_types; i++){
    susceptible_types[i] = utils::inumeric(FLERR, arg[arg_n++], true, lmp);
    printf("Bug %d type code: %d\n",i,susceptible_types[i]);
    
    //effector id
    susceptible_effectors[i] = utils::inumeric(FLERR, arg[arg_n++], true, lmp);
    printf("\t Effector ID: %d\n",susceptible_effectors[i]);

    //effect probability
    effector_success_probs[i] = utils::numeric(FLERR, arg[arg_n++], true, lmp);
    printf("\t Effect probability: %2g\n",effector_success_probs[i]);

    //effect group
    int toxgroup = -1;
    toxgroup = group->find(arg[arg_n++]);
    toxgroup = 1 | group->bitmask[toxgroup];
    if(toxgroup < 0){
        error->all(FLERR, "Can't find intoxication group.");
    }
    effector_groups[i] = toxgroup;

    //effect type
    effector_effects[i] = utils::inumeric(FLERR, arg[arg_n++], true, lmp);
    printf("\t Effect type: %d\n",effector_effects[i]);

  }
  //TODO sanity check for other args beyond this, could indicate situation
  //where a line defining a bug was added but the n_bugs style param
  //was not updated
}

/* ---------------------------------------------------------------------- */

FixT6SSContact::~FixT6SSContact() {
    delete random;
    memory->destroy(t6ss_types);
    memory->destroy(attacker_effectors);
    memory->destroy(harpoon_lens);
    memory->destroy(cooldowns);
    memory->destroy(susceptible_types);
    memory->destroy(susceptible_effectors);
    memory->destroy(effector_success_probs);
    memory->destroy(effector_groups);
    memory->destroy(effector_effects);
}

/* ---------------------------------------------------------------------- */
void FixT6SSContact::biology_nufeb()
{
    if (update -> ntimestep % nevery) return;
    do_contact();
}

/* ---------------------------------------------------------------------- */

int FixT6SSContact::setmask() {
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}
void FixT6SSContact::do_contact(){
    int i,j, ii, jj, inum, jnum, itype, jtype;
    int *ilist, *jlist, *numneigh, **firstneigh;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    std::vector<int> intoxications; 
    
    //there's going to be a series of neighbor lists we iterate thru
    //each neighbor list has one key bug and any associated neighbors
    //conventions: i is the key bug, j is the current neighbor to i
    for(ii = 0; ii < inum; ii++){
        i = ilist[ii];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        //printf("Bug %d has %d neighbors: ",i,jnum);
        for(jj =0; jj < jnum; jj++){
            j = jlist[jj];
            j &= NEIGHMASK;
            //printf("%d ",j);
            
            // nothing to do if there's no neighbors
            if(jnum <= 0) break;

            int attack_potentials = is_attack_possible(i,j);
             // if i can attack j, do attack, record intoxication 
            if(I_CANHIT_J & attack_potentials){
                intoxications.push_back(j);
            }
            // if j can attack i, do attack, record intoxication
            if(J_CANHIT_I & attack_potentials){
                intoxications.push_back(i);
            }
        }
        //printf("\n");
    }
    
    int *mask = atom->mask;
    int *type = atom->type;
    //apply intoxication transitions
    std::vector<int>::const_iterator intox_i;
    for(intox_i = intoxications.begin(); intox_i != intoxications.end(); intox_i++){
        //lookup the correct index value for intoxicated bug type
        int effect_lookup = -1;
        for(int st = 0; st < n_susceptible_types; st++){
            //printf("st %d, susceptible_type: %d\n",st, susceptible_types[st]);
            //when true, st is an index into susceptibility data associated 
            //with the intoxicated bug's type
            if(type[*intox_i] == susceptible_types[st]){
               effect_lookup = st; 
            }
        }
        if(effect_lookup != -1){
        //printf("Intoxication application on Bug:%d of type %d. Lookup %d  mask %d  group %d\n",*intox_i, type[*intox_i], effect_lookup, effector_groups[effect_lookup],effector_effects[effect_lookup]);
        mask[*intox_i] = effector_groups[effect_lookup];
        type[*intox_i] = effector_effects[effect_lookup];
        }
    } 
/*    double **conc = grid->conc;
    double **reac = grid->reac;
    double **dens = grid->dens;
    for(int i = 0; i < grid->ncells; i++){
        if(grid->mask[i] & GRID_MASK){
            printf("Proc grid %d dens %7g\n",i,dens[igroup][i]);
        }
    }
 */   
}


//All T6SS attack specifications are associated with a bug type
//this returns the appropriate index used to retrieve attack specifications
//based on the type associated with the individual bug id
//returns -1 if the bug type does not have a T6SS specification
// TODO rethink this so that we can multiple effectors per bug type
// TODO arrange code so this isn't called multiple times
int FixT6SSContact::get_attack_index(int bug_id){
    int *type = atom->type;
    
    int i = -1;
    //printf("Looking up attack index of bug %d of type %d\n",bug_id,type[bug_id]);
    auto lookup = attacker_lookup.find(type[bug_id]);
    bool is_attacker = lookup != attacker_lookup.end();
    if(is_attacker){
        i = lookup->second;
    }
    return i;
}

// for bug i to be susceptible to bug j, bug i's type  must have a listed 
// susceptibility to an effector which is listed with an attack associated with
// but j's type
// TODO, could probably return a list - save work when applying intoxication
bool FixT6SSContact::is_susceptible(int i, int j){
    int *type = atom->type;
    
    //get the effector associated with j
    int j_index = get_attack_index(j);
    if(j_index == -1) return false;
    int j_effector = attacker_effectors[j_index ];
   
    //check each associated susceptibility associated with i against effector
    for(int k = 0; k < n_susceptible_types; k++){
        //when true, k is an index into susceptibility data associated 
        //with bug i's type
        if(type[i] == susceptible_types[k]){
            //printf("Bug %d type %d is susceptible to effector %d. Bug %d is of type %d and has effector %d\n",i,type[i],susceptible_effectors[k],j,type[j],j_effector);
            if(susceptible_effectors[k] == j_effector) return true;
        }
    }
    return false;
}

// give the distance between SURFACES of two atoms  along dimension k, 
// taking into account radius, domain size and periodicity in that direction
double FixT6SSContact::periodic_displacement(int i, int j, int k){
    double **x = atom->x;
    double *radius = atom->radius;
    double irad = radius[i];
    double jrad = radius[j];

    //centerpoint distance
    double d = fabs(x[i][k] - x[j][k]);
    //printf("%d vs %d displacement %d is %7g\n",i,j,k,d);
    if (d > domain->prd_half[k]) d = domain->prd[k]-d;
    //printf("%d vs %d perioidc displacement %d is %7g\n",i,j,k,d);
    return d; 
}

// For a pair of bugs, determine if an attack is possible from i to j or
// j to i.  Encodes outcome in integer intended to be masked with:
// I_CANHIT_J and J_CANHIT_I
int FixT6SSContact::is_attack_possible(int i, int j){
    int *type = atom->type;
    int *mask = atom->mask;

    int possible_attacks = NO_ATTACKS;
  
    //at times, we'll have eliminated one attack route early
    //these are used to avoid unnecessary checks once we determine
    //an attack route is impossible
    bool keep_checking_i_attack_j = true;
    bool keep_checking_j_attack_i = true;

    // If either bug isn't part of a T6SS contact fix, then there's no attack
    bool is_i_T6SS = (mask[i] & groupbit);
    bool is_j_T6SS = (mask[j] & groupbit);
    //if(is_i_T6SS) printf("%d type %d is T6SS participant\n",i,type[i]); 
    //if(is_j_T6SS) printf("%d type %d is T6SS participant\n",j,type[j]); 
    if(!(is_i_T6SS || is_j_T6SS)) return NO_ATTACKS;

    // Bugs of the same type cannot intoxicate each other
    if(type[i] == type[j]) return NO_ATTACKS;

    //if((i==1 && j==4) || (j==4 && i==1))
        //printf("Got to rough dist check\n");
    //Based on quick x, y, z displacement calculations and harpoon lengths,
    //determine if the bugs have any chance of reaching each other
    double dx = periodic_displacement(i,j,0);
    double dy = periodic_displacement(i,j,1);
    double dz = periodic_displacement(i,j,2);
    double max_disp = std::max({dx,dy,dz});
    
   
    int i_ai = get_attack_index(i);
    int j_ai = get_attack_index(j); 

    //keep_checking_i_attack_j = (i_ai != -1) && (harpoon_lens[i_ai] > max_disp);
    //keep_checking_j_attack_i = (j_ai != -1) && (harpoon_lens[j_ai] > max_disp);
    keep_checking_i_attack_j = (harpoon_lens[i_ai] > max_disp);
    keep_checking_j_attack_i = (harpoon_lens[j_ai] > max_disp);
  

    //are there any effector-susceptibility pairs?
    if(keep_checking_i_attack_j){
        keep_checking_i_attack_j = is_susceptible(j, i);
    } 
     if(keep_checking_j_attack_i){
        keep_checking_j_attack_i = is_susceptible(i, j);
    }  
    
    // do the full distance calculation
    //avoiding sqrt call and compariing against distance squared
    //double dist_i_j = sqrt(dx*dx+dy*dy+dz*dz);
    double dist_i_jsq = dx*dx+dy*dy+dz*dz;
    double *radius = atom->radius;
    double radius_offset = radius[i]+radius[j];

    double r2 = radius_offset*radius_offset;
    double harplen = 0;
    if(keep_checking_i_attack_j){
    //    keep_checking_i_attack_j = (harpoon_lens[i_ai]>(dist_i_j-radius_offset));
        harplen=harpoon_lens[i_ai];
        keep_checking_i_attack_j = (harplen*harplen+2*harplen*radius_offset + r2)  >  dist_i_jsq;
    }
    if(keep_checking_j_attack_i){
     //   keep_checking_j_attack_i = (harpoon_lens[j_ai]>(dist_i_j-radius_offset));
        harplen=harpoon_lens[j_ai];
        keep_checking_j_attack_i = (harplen*harplen+2*harplen*radius_offset + r2)  >  dist_i_jsq;
}

    //TODO be better about taking advantage of boolean for masking
    if(keep_checking_i_attack_j) possible_attacks = possible_attacks | I_CANHIT_J; 
    if(keep_checking_j_attack_i) possible_attacks = possible_attacks | J_CANHIT_I; 
    
    return possible_attacks;
}


void FixT6SSContact::init_list(int id, NeighList *ptr){
    list = ptr;
}

void FixT6SSContact::init(){
    neighbor->add_request(this,NeighConst::REQ_DEFAULT);
}
