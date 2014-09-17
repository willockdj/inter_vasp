#include <stdio.h>
#include <string.h>
#include "maxima.h"
#include "structures.h"

int find_atom_label(atom *p_molecule, int num_atoms, char *p_ref_lab,
                    int *p_mol_index)
{
int index, iatom;

index= -1;
for (iatom= 0; iatom <= num_atoms; iatom++)
  {
    if (strcmp(p_molecule->label, p_ref_lab) == 0)
     {
       index = iatom;
       *p_mol_index = p_molecule->mol;
       break;
     }
    p_molecule++;
  }

return index;
}
