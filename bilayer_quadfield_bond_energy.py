import transform_quadraticfield_cif_v2
import pymatgen.io.cif
import numpy as np
import math

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('Top_cif_filename', type=str, help='The cif file name for top surface.')
    parser.add_argument('m1', type=int, help='m1 for m1*I + m2*Omega quadfield transformation for top surface')
    parser.add_argument('m2', type=int, help='m2 for m1*I + m2*Omega quadfield transformation for top surface')
    parser.add_argument('Bottom_cif_filename', type=str, help='The cif file name for bottom surface.')
    parser.add_argument('n1', type=int, help='n1 for n1*I + n2*Omega quadfield transformation for top surface')
    parser.add_argument('n2', type=int, help='n2 for n1*I + n2*Omega quadfield transformation for top surface')
    parser.add_argument('d', type=int, help='d for quadfield')
    parser.add_argument('-w', '--workdir', type=str, default='', help='working directory to export POSCAR and cif files')
    parser.add_argument('-i', '--init', type=float, default='0.0', help='initial layer distance')
    parser.add_argument('-f', '--final', type=float, default='0.0', help='final layer distance')
    parser.add_argument('-s', '--step', type=float, default='0.0', help='step for layer distance')

    args = parser.parse_args()

    print(f'{args.Top_cif_filename} is transformed by {args.m1}*I + {args.m2}*Omega')
    print(f'{args.Bottom_cif_filename} is transformed by {args.n1}*I + {args.n2}*Omega')

    top_structure = pymatgen.io.cif.CifParser(args.Top_cif_filename).get_structures()[0]
    bottom_structure = pymatgen.io.cif.CifParser(args.Bottom_cif_filename).get_structures()[0]

    top_supercell_structure = transform_quadraticfield_cif_v2.get_transformed_quadfield_structure(top_structure, args.m1, args.m2, args.d)
    bottom_supercell_structure = transform_quadraticfield_cif_v2.get_transformed_quadfield_structure(bottom_structure, args.n1, args.n2, args.d)

    top_strain_structure, bottom_strain_structure = transform_quadraticfield_cif_v2.match_planar_coordinates_both(top_supercell_structure, bottom_supercell_structure)
    top_component_structure, bottom_compoment_structure = transform_quadraticfield_cif_v2.get_out_of_plain_coords_in_heterostructure(top_strain_structure, bottom_strain_structure)

    concatonated_structure = transform_quadraticfield_cif_v2.concat_structures(top_component_structure, bottom_compoment_structure)

    transform_quadraticfield_cif_v2.export_structure('top_nostrain', top_supercell_structure, output_dir=args.workdir, export_type='vdir')
    transform_quadraticfield_cif_v2.export_structure('bottom_nostrain', bottom_supercell_structure, output_dir=args.workdir, export_type='vdir')
    transform_quadraticfield_cif_v2.export_structure('top_c', top_component_structure, output_dir=args.workdir, export_type='vdir')
    transform_quadraticfield_cif_v2.export_structure('bottom_c', bottom_compoment_structure, output_dir=args.workdir, export_type='vdir')
    transform_quadraticfield_cif_v2.export_structure('concat_c', concatonated_structure, output_dir=args.workdir, export_type='vdir')

    if args.step != 0.0:
        interlayer_distance_list = np.linspace(args.init, args.final, int(np.floor((args.final - args.init)/args.step)) + 1)
    else:
        interlayer_distance_list = [args.init]
    for interlayer_distance in interlayer_distance_list:
        top_component_structure, bottom_compoment_structure = transform_quadraticfield_cif_v2.get_out_of_plain_coords_in_heterostructure(top_strain_structure, bottom_strain_structure, distance=interlayer_distance)
        distance_structure = transform_quadraticfield_cif_v2.concat_structures(top_component_structure, bottom_compoment_structure)
        transform_quadraticfield_cif_v2.export_structure(f'distance_{interlayer_distance:.2f}', distance_structure, export_type='cif_vdir')
    
if __name__ == "__main__":
    main()