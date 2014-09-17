#
# development makefile 
#
# for SGI
#
#CC = cc
#CFLAGS = -g -mips2 -xansi -wlint -fullwarn -prototypes 
#LFLAGS = -lm
#
BIN = /home/sacdjw/bin
#
# for LINUX
#
CC = gcc
#
#
CFLAGS = -Wall
#
LFLAGS = -lm
#
COBJS = main.o read_input.o read_car.o locate_string.o read_line.o \
        read_atom_data.o \
        get_doub.o put_string.o copy_int.o int_to_string.o next_space.o \
        next_none_space.o get_int.o generate_neighbours.o standard_bond.o \
        atom_separation_squared.o compare_strings.o min_image.o open_file.o \
        cart_to_fract.o fract_to_cart.o cart_latt_vecs.o vec_cross.o vec_dot.o \
        find_chunk.o move_molecule.o rotate_with_flags.o write_car.o \
        write_atom_data.o unit_vector.o size_vector.o \
        move_molecule_with_flags.o  \
        string_to_int.o real_random.o tok_get.o get_double.o get_integer.o \
        find_kind.o read_gulp.o read_atom_data_glp.o write_gulp.o \
        write_atom_data_gulp.o \
        read_poscar.o read_potcar.o read_atom_data_vasp.o \
        write_atom_data_poscar.o \
        write_poscar.o sort_by_elem.o read_outcar.o read_incar.o \
        cut_to_uscore.o cut_to_dot.o latt_vecs_from_cart.o \
        find_line.o find_line_with_stopper.o find_mol.o atomic_mass_list.o centre_of_mass.o \
        find_atom_label.o centre_of_mass_flagged.o flag_chosen_atoms.o \
        move_selected_atoms.o write_pdb.o read_force_block.o \
        react_coord.o mol_dot.o mol_dist.o read_doscar.o write_doscsv.o \
        smear_dos.o read_xdatcar.o read_fdf.o read_siesta_vectors.o \
        count_doscar.o reorientate_cell.o rotate_vector.o read_star_format.o

#


main.o             : main.c read_input.c read_outcar.c read_incar.c read_car.c \
                     cart_latt_vecs.c \
                     generate_neighbours.c \
                     compare_strings.c unit_vector.c move_molecule_with_flags.c \
                     size_vector.c string_to_int.c open_file.c read_poscar.c \
                     read_potcar.c  read_fdf.c\
                     structures.h header.h data.h maxima.h 
read_input.o       : read_input.c tok_get.c find_kind.c get_double.c get_integer.c \
                     structures.h maxima.h reader.h global_values.h
tok_get.o          : tok_get.c find_kind.c maxima.h global_values.h reader.h 
find_kind.o        : find_kind.c maxima.h reader.h
get_double.o       : get_double.c tok_get.c global_values.h
get_integer.o      : get_integer.c tok_get.c global_values.h
read_car.o         : read_car.c locate_string.c read_atom_data.c get_doub.c \
                     put_string.c structures.h global_values.h
locate_string.o    : locate_string.c 
read_atom_data.o   : read_atom_data.c locate_string.c next_none_space.c \
                     next_space.c get_int.c get_doub.c copy_int.c put_string.c \
                     int_to_string.c structures.h
read_atom_data_vasp.o : read_atom_data_vasp.c
get_doub.o         : get_doub.c get_int.c global_values.h
put_string.o       : put_string.c global_values.h
copy_int.o         : copy_int.c 
next_none_space.o  : next_none_space.c 
int_to_string.o    : int_to_string.c 
next_space.o       : next_space.c 
get_int.o          : get_int.c global_values.h
generate_neighbours.o : generate_neighbours.c standard_bond.c  atom_separation_squared.c \
                        structures.h maxima.h data.h
standard_bond.o    : standard_bond.c compare_strings.c maxima.h data.h
atom_separation_squared.o : atom_separation_squared.c min_image.c structures.h global_values.h
compare_strings.o  : compare_strings.c 
min_image.o        : min_image.c cart_to_fract.c fract_to_cart.c maxima.h global_values.h structures.h
cart_to_fract.o    : cart_to_fract.c
fract_to_cart.o    : fract_to_cart.c
find_chunk.o       : find_chunk.c maxima.h global_values.h structures.h
move_molecule.o    : move_molecule.c structures.h
rotate_with_flags.o : rotate_with_flags.c move_molecule.c maxima.h global_values.h structures.h
write_car.o        : write_car.c write_atom_data.c structures.h
write_atom_data.o  : write_atom_data.c structures.h
unit_vector.o      : unit_vector.c size_vector.c
size_vector.o      : size_vector.c
move_molecule_with_flags.o : move_molecule_with_flags.c structures.h
string_to_int.o    : string_to_int.c
real_random.o      : real_random.c 
open_file.o        : open_file.c
read_gulp.o        : read_gulp.c read_line.c locate_string.c read_atom_data_glp.c \
                     get_doub.c get_int.c string_to_int.c next_none_space.c \
                     next_space.c int_to_string.c \
                     structures.h maxima.h global_values.h
read_poscar.o      : read_poscar.c read_line.c locate_string.c \
                     read_atom_data_vasp.c get_doub.c latt_vecs_from_cart.c \
                     put_string.c structures.h global_values.h
read_fdf.o         : read_fdf.c
read_potcar.o      : read_potcar.c tok_get.c structures.h
read_outcar.o      : read_outcar.c tok_get.c latt_vecs_from_cart.c structures.h \
                     find_line.c read_force_block.c
read_incar.o       : read_incar.c tok_get.c structures.h
write_gulp.o       : write_gulp.c cart_to_fract.c write_atom_data_gulp.c put_string.c \
                     structures.h global_values.h
write_atom_data_gulp.o : write_atom_data_gulp.c structures.h
write_atom_data_poscar.o : write_atom_data_poscar.c structures.h
write_poscar.o : write_poscar.c write_atom_data_poscar.c structures.h
sort_by_elem.o : sort_by_elem.c structures.h
cut_to_uscore.o : cut_to_uscore.c 
cut_to_dot.o : cut_to_dot.c 
latt_vecs_from_cart.o : latt_vecs_from_cart.c vec_cross.c vec_dot.c constants.h
find_line.o       : find_line.c tok_get.c structures.h maxima.h
find_line_with_stopper.o       : find_line_with_stopper.c tok_get.c structures.h maxima.h
find_mol.o        : find_mol.c structures.h
atomic_mass_list.o : atomic_mass_list.c maxima.h data.h
centre_of_mass.o   : centre_of_mass.c
centre_of_mass_flagged.o   : centre_of_mass_flagged.c structures.h
find_atom_label.o : find_atom_label.c structures.h maxima.h
flag_chosen_atoms.o : flag_chosen_atoms.c structures.h maxima.h
move_selected_atoms.o : move_selected_atoms.c structures.h maxima.h
write_pdb.o  : write_pdb.c maxima.h structures.h constants.h
read_force_block.o : read_force_block.c structures.h
react_coord.o      : react_coord.c structures.h
mol_dot.o          : mol_dot.c structures.h
mol_dist.o          : mol_dist.c structures.h
read_doscar.o      : read_doscar.c structures.h
count_doscar.o     : count_doscar.c structures.h
write_doscsv.o     : write_doscsv.c structures.h
smear_dos.o        : smear_dos.c structures.h
read_xdatcar.o     : read_xdatcar.c get_integer.c structures.h
read_siesta_vectors.o : read_siesta_vectors.c
reorientate_cell.o : reorientate_cell.c
rotate_vector.o : rotate_vector.c
read_star_format.o : read_star_format.c

inter_vasp : $(COBJS) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/inter_vasp $(COBJS) $(LFLAGS) 

inter_vasp_dos : $(COBJS) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/inter_vasp_dos $(COBJS) $(LFLAGS) 

inter_vasp_test : $(COBJS) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/inter_vasp_test $(COBJS) $(LFLAGS) 

trial      : $(COBJS) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/trial $(COBJS) $(LFLAGS) 

