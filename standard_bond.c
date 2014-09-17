/*************************************************/
/*                                               */
/*   look up standard bond length                */
/*                                               */
/*************************************************/

#include <stdio.h>
#include <string.h>
#include "maxima.h"
#include "data.h"

int compare_strings( char *p_ichar1, char *p_ichar2 );

double standard_bond( char *atom1, char *atom2 )
{
 int iloop;

 for (iloop = 0; iloop < num_bond_list; iloop++)
  {
     if (   strlen(bond_table[iloop].label1) == strlen(atom1)
         && strlen(bond_table[iloop].label2) == strlen(atom2) )
          {
             if (    compare_strings( &(bond_table[iloop].label1[0]), atom1 )
                  && compare_strings( &(bond_table[iloop].label2[0]), atom2 ) )
                         return bond_table[iloop].bond_length;
           }

       
     if (   strlen(bond_table[iloop].label2) == strlen(atom1)
         && strlen(bond_table[iloop].label1) == strlen(atom2) )
          {
            if (    compare_strings( &(bond_table[iloop].label1[0]), atom2 )
                 && compare_strings( &(bond_table[iloop].label2[0]), atom1 ) )
               return bond_table[iloop].bond_length;
          }
  }
 return -1;
}

