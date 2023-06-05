#ifndef LMP_GIT_VERSION_H
#define LMP_GIT_VERSION_H
bool LAMMPS_NS::LAMMPS::has_git_info() { return false; }
const char *LAMMPS_NS::LAMMPS::git_commit() { return "(unknown)"; }
const char *LAMMPS_NS::LAMMPS::git_branch() { return "(unknown)"; }
const char *LAMMPS_NS::LAMMPS::git_descriptor() { return "(unknown)"; }
#endif
