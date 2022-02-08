#!/usr/bin/python3
""" 
    Converts AMBER NetCDF files into Gaussian inputs. 
    Written by Mariano Curti, Institut Català d'Investigació Química.
    Funding from the EU through the I2: ICIQ Impulsion programme is gratefully acknowledged.

"""
# To test:
# - Save files to current directory instead of that of the prmtop
# is sanitize_name necessary?

# To do:
# all paths should be given from CLI
# Check extra spaces in prints
# Export pdb & xyz
# Apply type conversion before
# "\n" in amber2gaussian.log
# - Improve water selection algorithm
# - Seg fault when file order wrong
# - Improve distance to capping
# - Check input files actually exist
# - Take all information (including route section) from input file
# - "C": "H4" problem. Read "capping_type"
# - Check whether parm corresponds to coordinates
# - Option to save parameters in external file [do it for first iteration]
# - Add error catching
# Sanitize input so that ints or strings are type-checked
# even spacing in parameters section
# User should be warned when link_data is not available
# Notice of the ~125 frames bug in parmed
# What happens when capping atom is connected to more than one frag?
# bottom.txt / TIP3P
# >>> sort parameters section (to ensure diff works)

import fileinput
import sys
import argparse
import datetime
import parmed as pmd
from elementnames import ElementName

VERSION = "v0.1.7"


def parse_args():
    """ Parses command line arguments """
    parser = argparse.ArgumentParser(description='Converts AMBER topology/\
                                                  coordinates pairs into Gaussian inputs')
    parser.add_argument('topfile', help='AMBER topology file')
    parser.add_argument('coordfile', help='AMBER NetCDF trajectory or restart file')
    parser.add_argument('inputfile', help='Input file')
    parser.add_argument('-o', '--output', help='Output Gaussian job filename. \
                         It is saved as amber2gaussian.<jobtype>.<frame>.gjf. \
                         Default: <topfile>.<jobtype>.<frame>.gjf')
    parser.add_argument('-p', '--param_file', help='External file containing the MM parameters.')
    parser.add_argument('-f', '--frames', type=int, help='Maximum number of frames to process.')
    return parser.parse_args()


def parse_input_file(file):
    """ Function to parse the input file with user options."""
    with open(file) as inp_file:
        data = inp_file.readlines()
        opts = {}

        for inp_line in data:
            # parse input, assign values to variables
            inp_line = inp_line.split("=")

            if len(inp_line) == 2:
                key, value = inp_line[0], inp_line[1]
                key = key.strip().lower()

                if "mask" not in key:
                    # Sanitize non-mask inputs
                    value = value.strip().lower()

                opts[key] = value

    # Default values for non-set parameters
    default_options = {"mask": "@*",
                       "high_level_mask": "@*",
                       "job_type": "oniom",
                       "freeze_mask": "!@*",
                       "print_connectivity": True,
                       "freeze_capping": False,
                       "use_charges": False}

    for key, default in default_options.items():
        if key not in opts:
            opts[key] = default

    opts["print_connectivity"] = opts["print_connectivity"] in ['true', '1', True]
    opts["freeze_capping"] = opts["freeze_capping"] in ['true', '1', True]
    opts["use_charges"] = opts["use_charges"] in ['true', '1', True]

    assert "job_top" in opts, "A filename must be given for the top file."

    assert opts["job_type"] in ["oniom", "eet", "regular"], "Job type must be either \
                                                             ONIOM, EET, or regular."

    if opts["job_type"] == "eet":
        assert "fragments" in opts, "An EET calculation has been requested, \
            but the number of fragments (\"fragments\") has not been defined.\
            Program stopped."

        try:
            opts["fragments"] = int(opts["fragments"])
        except TypeError:
            print("The number of \
            fragments (\"fragments\") is not valid. Program stopped.")
            sys.exit()

        assert opts["fragments"] > 0 and opts["fragments"] <= 30, "The\
            number of fragments (\"fragments\") must be between 1 and \
            30. Program stopped."

        for frag in range(1, opts["fragments"] + 1):
            if f"fragment{frag}_mask" not in opts:
                print(f'An EET calculation has been requested, but the mask \
                      defining fragment {frag} ("fragment{frag}_mask") has not been defined. \
                      Program stopped.')
                sys.exit()

            if frag == 1:
                opts["high_level_mask"] = opts[f"fragment{frag}_mask"]
            else:
                opts["high_level_mask"] += "|" + opts[f"fragment{frag}_mask"]

    #opts["job_bottom"] = "bottom_short.txt" # !!!

    return opts


def get_n_waters():
    """Expands a distance-based mask until options['n_waters'] water molecules
        are selected"""
    step = 1.0  # Step, in A, in which the mask will be expanded (variable)
    max_dist = 20.0  # Largest distance that will be tested

    distance = 0.0

    while distance < max_dist:
        mask = f"(({options['waters_from']}<:{str(distance)})&:WAT)"

        selected_atoms = pmd.amber.mask.AmberMask(parm, mask).Selection()

        n_selected_waters = sum(selected_atoms) / 3

        if n_selected_waters == int(options["n_waters"]):
            print(f"    Distance {distance} worked. {int(n_selected_waters)} water molecules selected.")
            return f"|{mask}"
        else:
            print(f"    Distance {distance} didn't work. Target water molecules: {options['n_waters']}. Selected: {int(n_selected_waters)}.")

            if n_selected_waters > int(options["n_waters"]):
                # If we went over the target number, step back and decrease
                # the step size
                distance = round(distance - step, 5)
                step = step / 2
                distance = round(distance + step, 5)
            else:
                distance = round(distance + step, 5)

    # If options['n_waters'] was not reached before the limit, just return an
    # empty string
    return ""


def format_coord(coord):
    """Give a nice format to x,y,z coordinates and transform to string"""
    return ("%.6f" % coord).rjust(9, ' ')


def sanitize_name(name):
    """Remove - and + symbols from atom or residue names (e.g. from Na+) to
    comply with Gaussian requirements"""
    remove = "-+"
    for char in remove:
        name = name.replace(char, "")
    return name


def link_data(atm):
    """Returns a string specifying link atom information (if any) for ONIOM"""
    # It could be desirable to give this information as input
    capping_type = {"C": "H4", "N": "H"}
    if hasattr(atm, "isCapping"):
        if atm.QMtype in capping_type:
            return f" H-{capping_type[atm.QMtype]} {atm.QMbond}   0.0000"
        else:
            return ' ! This is a capping atom, but there is no link type.'
    else:
        return ""


def freeze_code(atm, atm_index):
    """ Returns the freeze code for a given atom """
    if frozen_list[atm_index] == 0:
        return " 0"

    if hasattr(atm, "isCapping") and not options["freeze_capping"]:
        return " 0"

    return "-1"


def replace_nonamber_names(atm_types):
    """ Replaces non-AMBER FF names (recogized by small caps in the first \
        character) by made-up ones """
    conversion = {}
    suffix_index = 0
    suffixes = ["E", "F", "D", "G", "I", "J", "K", "L", "M", "Q", "S", "U",
                "Y", "Z", "6", "7", "9", "0"]

    for atm_type in sorted(atm_types):
        if atm_type[0].isupper() is False and atm_type[0].isdigit() is False:
            conversion[atm_type] = atm_type[0].upper() + suffixes[suffix_index]
            suffix_index += 1
        elif atm_type == "Na+":
            conversion[atm_type] = "QN"  # To avoid Gaussian confusion w/NA
        elif atm_type == "Cl-":
            conversion[atm_type] = "QC"

        else:
            conversion[atm_type] = atm_type
    return conversion


def sort_dihedral_atoms(dih):
    """ Sort (alphabetically) the atoms in a dihedral.
        This prevents printing the same dihedral twice.
    """
    orig_atom1 = type_conversion[dih.atom1.type]
    orig_atom2 = type_conversion[dih.atom2.type]
    orig_atom3 = type_conversion[dih.atom3.type]
    orig_atom4 = type_conversion[dih.atom4.type]

    if orig_atom1 != orig_atom4:
        new_atom1 = sorted([orig_atom1, orig_atom4])[0]
        new_atom2 = orig_atom2
        new_atom3 = orig_atom3
        new_atom4 = sorted([orig_atom1, orig_atom4])[1]

        # if order was flipped, do the same for the inner atoms
        if new_atom1 != orig_atom1:
            new_atom2 = orig_atom3
            new_atom3 = orig_atom2

    else:  # If atom1 and atom4 are the same, sort the inner atoms
        new_atom1 = orig_atom1
        new_atom2 = sorted([orig_atom2, orig_atom3])[0]
        new_atom3 = sorted([orig_atom2, orig_atom3])[1]
        new_atom4 = orig_atom4

    return new_atom1, new_atom2, new_atom3, new_atom4


if __name__ == '__main__':
    print(f"amber2gaussian {VERSION}\r")
    args = parse_args()

    if args.output is None:
        # Default output file uses topfile name (but current path)
        args.output = "amber2gaussian"

    options = parse_input_file(args.inputfile)

    # Load coordinates file
    try:
        trajectory = pmd.load_file(args.coordfile)
    except FileNotFoundError as error:  # except FileNotFoundError as error
        print("A problem occured while loading the trajectory. \
              Program stopped.")
        print(error)
        sys.exit()

    sys.stdout.write(f"{len(trajectory.coordinates)} frame(s) detected.\n")
    if args.frames is not None:
        sys.stdout.write(f"{args.frames} frame(s) will be processed.\n")

    # Loop over all frames in the trajectory
    for frame_id, frame in enumerate(trajectory.coordinates):
        # If requested, stop at the appropiate frame
        if frame_id == args.frames:
            break

        parm = pmd.load_file(args.topfile, xyz=frame)  # This is the main bottleneck

        sys.stdout.write(f"Processing frame {frame_id + 1}...\r")
        sys.stdout.flush()

        output_file = f"{args.output}.{options['job_type']}.{frame_id}.gjf"

        geometry = []                       # Geometry section
        connect_table = []                  # Used in ONIOM
        vdw_parameters_table = set()        # MM parameters for ONIOM
        bonds_parameters_table = set()      # MM parameters for ONIOM
        angles_parameters_table = set()     # MM parameters for ONIOM
        impropers_parameters_table = set()  # MM parameters for ONIOM
        dihedrals_parameters_table = set()  # MM parameters for ONIOM
        charges_table = []                  # Used with Gaussian CHARGE keyword
        capping_atoms = []
        atom_types = set()

        # This avoids problems when the masked atoms change through frames
        high_level_mask = options["high_level_mask"]

        if "n_waters" in options and "waters_from" in options:
            high_level_mask += get_n_waters()

        # Select all atoms matching mask
        sel_atoms = pmd.amber.mask.AmberMask(parm, options["mask"]).Selected()

        # A list is needed to iterate more than once
        atoms_list = list(sel_atoms)

        # List of atoms in QM region
        qm_list = pmd.amber.mask.AmberMask(parm, high_level_mask).Selection()

        # Lists of atoms in each fragment for eet calculations
        if options["job_type"] == "eet":
            fragments_list = []

            for f in range(1, options["fragments"] + 1):
                mask = options[f"fragment{f}_mask"]
                fragments_list.append(pmd.amber.mask.AmberMask(parm, mask).Selection())

        # List of frozen atoms
        frozen_list = pmd.amber.mask.AmberMask(parm, options["freeze_mask"]).Selection()

        # Loop over all selected atoms to set QM, fragment, and capping status,
        # and save unique atom types; and store total charge
        total_charge = 0.0
        qm_charge = 0.0

        # original_index corresponds to all atoms in the parm file, while
        # new_index corresponds to the list of selected atoms
        for new_index, original_index in enumerate(atoms_list):
            atom = parm.atoms[original_index]
            total_charge += atom.charge

            if options["job_type"] == "oniom":
                # save (unique) atom types
                atom_types.add(atom.type)

                # save new index to use in connectivity table
                atom.new_index = new_index

            if qm_list[original_index] == 1:
                atom.isQM = True
                qm_charge += atom.charge

                # Assign to fragment
                if options["job_type"] == "eet":
                    for i, fragment in enumerate(fragments_list):
                        if fragment[original_index] == 1:
                            atom.fragment = i + 1

                # Check if atoms connected to the QM region are included in it
                # If not, mark them as "capping"
                for partner in atom.bond_partners:
                    if qm_list[partner.idx] == 0:
                        partner.isCapping = True
                        # and remember to which QM atom it is bonded
                        partner.QMbond = new_index + 1
                        # and the type (to find the appropiate H type)
                        partner.QMtype = atom.type
                        if options["job_type"] == "eet":
                            # The capping atom belongs to only one fragment. \
                            # Could give problems if it is connected to more \
                            # than one fragment.
                            partner.fragment = atom.fragment
            else:
                atom.isQM = False

        # This dictionary allows to replace nonamber types by made up ones, without repeating
        type_conversion = replace_nonamber_names(atom_types)

        # Loop (again) over all selected atoms to obtain geometry,
        # connectivity and charges information
        for index in atoms_list:
            atom = parm.atoms[index]
            atom_id = index + 1

            if options["job_type"] == "oniom":
                # Store geometry
                first = (f"{ElementName[atom.element]}-"
                         f"{type_conversion[atom.type]}-{atom.charge:.5f}"
                         f"(PDBName={sanitize_name(atom.name)},"
                         f"ResName={sanitize_name(atom.residue.name)},"
                         f"ResNum={atom.residue.number + 1})")
                first = first.ljust(60, ' ')

                second = (f"   {freeze_code(atom, index)}   "
                          f"{format_coord(atom.xx)}  "
                          f"{format_coord(atom.xy)}  "
                          f"{format_coord(atom.xz)}")

                Level = {True: " H", False: " L"}  # ONIOM layer specification

                third = Level[atom.isQM] + link_data(atom)

                geometry.append(first + second + third + "\n")

                if options["print_connectivity"]:
                    # Store connectivity table
                    connect_table.append(str(atom.new_index + 1))

                    for partner in atom.bond_partners:
                        # Avoid repeating connectivity info. The second check
                        # makes sure the atom is selected.
                        if partner.idx > atom_id - 1 and hasattr(partner, "new_index"):
                            connect_table.append(
                                f" {partner.new_index + 1} 1.0")

                    connect_table.append("\n")

            elif options["job_type"] == "eet" or options["job_type"] == "regular":
                if atom.isQM:
                    # Store geometry of QM atoms
                    fragment_data = ""
                    if options["job_type"] == "eet":
                        fragment_data = f"(Fragment={atom.fragment})"

                    first = f"{ElementName[atom.element]}{fragment_data}"
                    first = first.ljust(16, ' ')
                    first += f"  {freeze_code(atom, index)}   "

                    second = (f"{format_coord(atom.xx)}  "
                              f"{format_coord(atom.xy)}  "
                              f"{format_coord(atom.xz)}")

                    geometry.append(first + second + "\n")
                else:
                    # If it is a capping atom, store it as such.
                    # If not, store it as a point charge
                    if hasattr(atom, "isCapping"):
                        capping_atoms.append(atom)
                    else:
                        first = (f"{format_coord(atom.xx)}  "
                                 f"{format_coord(atom.xy)}  "
                                 f"{format_coord(atom.xz)}   ")

                        second = ("%.5f" % atom.charge).rjust(10, ' ')

                        charges_table.append(first + second + "\n")

        # Add capping atoms as H atoms to geometry
        if options["job_type"] == "eet" or options["job_type"] == "regular":
            for capping_atom in capping_atoms:
                fragment_data = ""
                if options["job_type"] == "eet":
                    fragment_data = f"(Fragment={capping_atom.fragment})"

                fragment_data = fragment_data.ljust(16, ' ')

                coords = (f"{freeze_code(capping_atom, capping_atom.idx)}   "
                          f"{format_coord(capping_atom.xx)}  "
                          f"{format_coord(capping_atom.xy)}  "
                          f"{format_coord(capping_atom.xz)}")

                geometry.append(f"H{fragment_data} {coords}")

                cap_data = (f" ! Capping H added in place of "
                            f"{ElementName[capping_atom.element]}@{capping_atom.idx}\n")

                geometry.append(cap_data)

        # Write Gaussian route section
        with open(output_file, 'w') as fout, fileinput.input(options["job_top"]) as fin:
            for line in fin:
                line = line.replace('{total_charge}', f'{round(total_charge)}')
                line = line.replace('{qm_charge}', f'{round(qm_charge)}')

                fout.write(line)

        # Print geometry (all job types)
        geometry.append("\n")
        open(output_file, 'a').write("".join(geometry))

        if options["job_type"] == "oniom":
            if options["print_connectivity"]:
                # Print connectivity table
                open(output_file, 'a').write("".join(connect_table))

            if args.param_file is None:
                # Extract parameters information from prmtop
                # This could actually be done once, instead of repeating
                # for every frame
                for at_type in sorted(atom_types):
                    vdw_data = (f'VDW "{type_conversion[at_type]}" '
                                f'{parm.LJ_radius[parm.LJ_types[at_type]-1]:.4f}  '
                                f'{parm.LJ_depth[parm.LJ_types[at_type]-1]:.4f}\n')

                    vdw_parameters_table.add(vdw_data)

                # instead of looping over all bonds (or angles, or dihedrals)
                # it could be more efficient to loop over parm.bond_types
                for bond in parm.bonds:
                    if hasattr(bond.atom1, "isQM") and hasattr(bond.atom2, "isQM"): # check all these
                        # Sorting atom1 and atom2 prevents printing the same bond twice
                        atom1 = sorted(
                            [type_conversion[bond.atom1.type], type_conversion[bond.atom2.type]])[0]
                        atom2 = sorted(
                            [type_conversion[bond.atom1.type], type_conversion[bond.atom2.type]])[1]
                        bonds_parameters_table.add(
                            f'HrmStr1 "{atom1}" "{atom2}" {bond.type.k:.2f} {bond.type.req:.4f}\n')

                for angle in parm.angles:
                    if hasattr(angle.atom1, "isQM") and hasattr(angle.atom2, "isQM") and hasattr(angle.atom3, "isQM"):
                        atom1 = sorted(
                            [type_conversion[angle.atom1.type], type_conversion[angle.atom3.type]])[0]
                        atom2 = type_conversion[angle.atom2.type]
                        atom3 = sorted(
                            [type_conversion[angle.atom1.type], type_conversion[angle.atom3.type]])[1]
                        angles_parameters_table.add(f'HrmBnd1 "{atom1}" "{atom2}" "{atom3}" '
                                                    f'{angle.type.k:.1f} {angle.type.theteq:.1f}\n')

                for dihedral in parm.dihedrals:
                    if dihedral.improper and hasattr(dihedral.atom1, "isQM") and hasattr(dihedral.atom2, "isQM") and hasattr(dihedral.atom3, "isQM") and hasattr(dihedral.atom4, "isQM"):
                        atom1 = type_conversion[dihedral.atom1.type]
                        atom2 = type_conversion[dihedral.atom2.type]
                        atom3 = type_conversion[dihedral.atom3.type]
                        atom4 = type_conversion[dihedral.atom4.type]

                        imp = (f'ImpTrs "{atom1}" "{atom2}" "{atom3}" "{atom4}"'
                               f'{dihedral.type.phi_k:.1f} {dihedral.type.phase:.1f} '
                               f'{dihedral.type.per:.1f} \n')

                        impropers_parameters_table.add(imp)

                phases = [0, 0, 0, 0]
                force_constants = [0, 0, 0, 0]

                for index, dihedral in enumerate(parm.dihedrals):
                    if not dihedral.improper and hasattr(dihedral.atom1, "isQM") and hasattr(dihedral.atom2, "isQM") and hasattr(dihedral.atom3, "isQM") and hasattr(dihedral.atom4, "isQM"):
                        term = dihedral.type.per - 1

                        phases[term] = dihedral.type.phase
                        force_constants[term] = dihedral.type.phi_k

                        atom1, atom2, atom3, atom4 = sort_dihedral_atoms(dihedral)

                        next_dihedral = parm.dihedrals[index + 1]
                        
                        if hasattr(next_dihedral.atom1, "isQM") and hasattr(next_dihedral.atom2, "isQM") and hasattr(next_dihedral.atom3, "isQM") and hasattr(next_dihedral.atom4, "isQM"):
                        
                            next_1, next_2, next_3, next_4 = sort_dihedral_atoms(next_dihedral)

                            if (atom1, atom2, atom3, atom4) != (next_1, next_2, next_3, next_4):
                                # This check makes sure multi-term dihedrals are
                                # only printed when the last term is reached (i.e.
                                # when the next dihedral contains different atoms).

                                phases_data = (f"{phases[0]:.0f} {phases[1]:.0f} "
                                               f"{phases[2]:.0f} {phases[3]:.0f}")
                                force_constants_data = (f"{force_constants[0]:.3f} "
                                                        f"{force_constants[1]:.3f} "
                                                        f"{force_constants[2]:.3f} "
                                                        f"{force_constants[3]:.3f}")
                                dihedrals_parameters_table.add(f'AmbTrs "{atom1}" "{atom2}" '
                                                               f'"{atom3}" "{atom4}"    {phases_data}'
                                                               f'    {force_constants_data}    1.0\n')

                                # Clear data to make room for next dihedral
                                phases = [0, 0, 0, 0]
                                force_constants = [0, 0, 0, 0]

                open(output_file, 'a').write("NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2\n\n")
                
                # We convert all sets to lists to sort them
                vdw_parameters_table = list(vdw_parameters_table)
                bonds_parameters_table = list(bonds_parameters_table)
                angles_parameters_table = list(angles_parameters_table)
                impropers_parameters_table = list(impropers_parameters_table)
                dihedrals_parameters_table = list(dihedrals_parameters_table)
                
                vdw_parameters_table.sort()
                bonds_parameters_table.sort()
                angles_parameters_table.sort()
                impropers_parameters_table.sort()
                dihedrals_parameters_table.sort()
                
                open(output_file, 'a').write("".join(vdw_parameters_table))
                open(output_file, 'a').write("".join(bonds_parameters_table))
                open(output_file, 'a').write("".join(angles_parameters_table))
                open(output_file, 'a').write("".join(impropers_parameters_table))
                open(output_file, 'a').write("".join(dihedrals_parameters_table))

                # Print extra MM parameters from external file
                with open(output_file, 'a') as fout, fileinput.input(options["job_bottom"]) as fin:
                    for line in fin:
                        fout.write(line)
                    fout.write("\n\n\n\n\n")

            elif args.param_file != "":
                # Print link to external file containing MM parameters, if
                # given in the input
                open(output_file, 'a').write(f"@{args.param_file}/N\n\n")

        elif (options["job_type"] == "eet" or options["job_type"] == "regular") \
                and options["use_charges"]:
            open(output_file, 'a').write("".join(charges_table))

    sys.stdout.write("Conversion complete.   \n")

    # Save main input information in a log file
    with open("amber2gaussian.log", "a") as out:
        out.write(f'# amber2gaussian {VERSION} - '
                  f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
        out.write(f'# Command: {" ".join(sys.argv)}\n')
        out.write('# Input arguments:\n')
        out.write(f'{options}\n\n')
