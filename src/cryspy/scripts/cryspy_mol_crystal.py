#!/usr/bin/env python3
'''
分子のxyz座標から、その分子で構成される最適化分子性結晶のcifファイルを生成するスクリプト
Script to generate optimized molecular crystal CIF files from molecular xyz coordinates
'''

import argparse
from logging import getLogger
import os
import tempfile

import numpy as np
from pymatgen.core import Molecule, Structure
from pymatgen.io.cif import CifWriter
from pyxtal import pyxtal
from pyxtal.tolerance import Tol_matrix
from ase import Atoms
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from ase.constraints import FixSymmetry
from ase.calculators.emt import EMT
from pymatgen.io.ase import AseAtomsAdaptor

from cryspy.util.struc_util import get_mol_data, out_cif, set_mindist
from cryspy.util.utility import set_logger
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


logger = getLogger('cryspy')


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='分子のxyz座標から最適化分子性結晶のCIFファイルを生成'
    )
    parser.add_argument(
        'xyz_files', 
        nargs='+', 
        help='分子のxyzファイル（複数指定可能）'
    )
    parser.add_argument(
        '-o', '--output',
        default='optimized_crystal.cif',
        help='出力CIFファイル名 (default: optimized_crystal.cif)'
    )
    parser.add_argument(
        '-n', '--nstruct',
        type=int,
        default=1,
        help='生成する結晶構造の数 (default: 1)'
    )
    parser.add_argument(
        '--nmol',
        type=int,
        nargs='+',
        help='各分子の数（分子種ごとに指定、未指定時は自動設定）'
    )
    parser.add_argument(
        '--spgnum',
        type=int,
        nargs='+',
        default=[1, 2, 14, 15, 19, 61, 62, 63, 64, 65, 92, 96, 142, 143, 144, 145],
        help='空間群番号のリスト (default: 分子性結晶でよく使われる空間群)'
    )
    parser.add_argument(
        '--vol-factor',
        type=float,
        default=1.1,
        help='体積因子 (default: 1.1)'
    )
    parser.add_argument(
        '--fmax',
        type=float,
        default=0.05,
        help='構造最適化の力の閾値 (eV/Å) (default: 0.05)'
    )
    parser.add_argument(
        '--steps',
        type=int,
        default=1000,
        help='構造最適化の最大ステップ数 (default: 1000)'
    )
    parser.add_argument(
        '--calculator',
        choices=['EMT'],
        default='EMT',
        help='構造最適化に使用する計算機 (default: EMT)'
    )
    parser.add_argument(
        '--mindist-factor',
        type=float,
        default=0.8,
        help='最小距離の因子 (default: 0.8)'
    )
    parser.add_argument(
        '--no-optimization',
        action='store_true',
        help='構造最適化を行わない'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='デバッグモード'
    )
    parser.add_argument(
        '--optimize-density',
        action='store_true',
        help='全ての分子性結晶空間群を試して最も密度が高い構造を選択（構造最適化は無効化）'
    )
    
    return parser.parse_args()


def load_molecules(xyz_files):
    """Load molecules from xyz files."""
    molecules = []
    atom_types = set()
    
    for xyz_file in xyz_files:
        if not os.path.exists(xyz_file):
            raise FileNotFoundError(f"XYZファイルが見つかりません: {xyz_file}")
        
        mol = Molecule.from_file(xyz_file)
        molecules.append(mol)
        
        # 原子種を収集
        for species in mol.species:
            atom_types.add(species.symbol)
    
    atype = tuple(sorted(atom_types))
    logger.info(f"読み込んだ分子数: {len(molecules)}")
    logger.info(f"原子種: {atype}")
    
    return molecules, atype


def generate_molecular_crystal_all_symmetries(molecules, atype, nmol, spgnum, vol_factor, mindist_factor, max_attempts_per_spg=10):
    """Generate molecular crystal structures for all space groups and return the one with highest density."""
    
    from contextlib import redirect_stdout, redirect_stderr
    from io import StringIO
    
    # Set minimum distances
    if len(atype) == 1:
        mindist = ((2.0,),)  # 単原子種の場合
    else:
        # 複数原子種の場合、簡単な距離マトリックスを作成
        n_types = len(atype)
        mindist = []
        for i in range(n_types):
            row = []
            for j in range(n_types):
                # 典型的な共有結合半径に基づく最小距離
                if atype[i] == 'H' or atype[j] == 'H':
                    dist = 1.5
                elif atype[i] in ['C', 'N', 'O'] and atype[j] in ['C', 'N', 'O']:
                    dist = 2.5
                else:
                    dist = 3.0
                row.append(dist * mindist_factor)
            mindist.append(tuple(row))
        mindist = tuple(mindist)
    
    logger.info(f"最小距離設定: {mindist}")
    
    # Set tolerance matrix
    tolmat = Tol_matrix(prototype="molecular")
    
    successful_structures = []
    
    # Try all space groups
    for spg in spgnum:
        logger.info(f"空間群 {spg} で結晶生成を試行中...")
        
        # Use space-group-specific volume factors to allow natural density differences
        # Different space groups have different packing efficiencies, so we need to
        # let them settle to their natural volumes rather than forcing optimal packing
        base_factors = [0.8, 1.0, 1.2, 1.5, 2.0]  # Wider range, less biased toward high density
        
        # Add space-group-specific scaling to account for different symmetry efficiencies
        # Higher symmetry groups often pack more efficiently, so they get slightly higher starting volumes
        spg_multipliers = {
            range(1, 3): 1.0,      # P1, P-1 (triclinic) - low symmetry, keep as reference
            range(3, 16): 1.1,     # monoclinic 
            range(16, 75): 1.2,    # orthorhombic
            range(75, 143): 1.3,   # tetragonal  
            range(143, 168): 1.4,  # trigonal
            range(168, 195): 1.5,  # hexagonal
            range(195, 231): 1.6,  # cubic - highest symmetry, often most efficient packing
        }
        
        # Find the multiplier for this space group
        spg_mult = 1.0
        for spg_range, mult in spg_multipliers.items():
            if spg in spg_range:
                spg_mult = mult
                break
        
        volume_factors = [vol_factor * base_factor * spg_mult for base_factor in base_factors]
        
        # Instead of finding the "best" (highest density), take the first successful structure
        # This allows each space group to settle to its natural density
        structure_found = False
        
        for vol_f in volume_factors:
            if structure_found:
                break
                
            attempt = 0
            pyxtal_fail_count = 0  # Track PyXtal internal failures
            while attempt < max_attempts_per_spg:
                try:
                    # Create pyxtal structure
                    crystal = pyxtal(molecular=True)
                    
                    # Capture stdout/stderr to detect PyXtal internal failures
                    f = StringIO()
                    with redirect_stdout(f):
                        with redirect_stderr(f):
                            crystal.from_random(
                                dim=3,
                                group=spg,
                                species=molecules,
                                numIons=nmol,
                                factor=vol_f,
                                conventional=False,
                                tm=tolmat
                            )
                    
                    output = f.getvalue()
                    if output and "Cannot generate crystal after max attempts" in output:
                        pyxtal_fail_count += 1
                        logger.debug(f"PyXtal内部で最大試行回数に達しました (空間群 {spg}, vol_factor={vol_f:.2f}, 試行 {attempt + 1})")
                        # If PyXtal itself failed 3 times for this vol_factor, skip to next vol_factor
                        if pyxtal_fail_count >= 3:
                            logger.debug(f"空間群 {spg} (vol_factor={vol_f:.2f}) をスキップ - PyXtal内部失敗が多すぎます")
                            break  # Break to try next volume factor
                    
                    if crystal.valid:
                        # Convert to pymatgen structure
                        structure = crystal.to_pymatgen()
                        density = structure.density
                        
                        # Take this structure (first successful one, not optimizing for highest density)
                        successful_structures.append({
                            'structure': structure,
                            'space_group': spg,
                            'density': density,
                            'volume': structure.volume,
                            'vol_factor_used': vol_f
                        })
                        logger.info(f"空間群 {spg} で結晶構造生成成功: "
                                   f"密度 = {density:.3f} g/cm³, vol_factor={vol_f:.2f}")
                        structure_found = True
                        break  # Success, move to next space group
                    else:
                        logger.debug(f"無効な結晶構造: 空間群 {spg}, vol_factor={vol_f:.1f}, 試行 {attempt + 1}")
                        
                except Exception as e:
                    logger.debug(f"結晶生成失敗 (空間群 {spg}, vol_factor={vol_f:.1f}, 試行 {attempt + 1}): {e}")
                
                attempt += 1
        
        if not structure_found:
            logger.warning(f"空間群 {spg} で全ての体積因子で結晶構造生成に失敗")
    
    if not successful_structures:
        raise RuntimeError("全ての空間群で結晶構造生成に失敗しました")
    
    # Sort by density (highest first)
    successful_structures.sort(key=lambda x: x['density'], reverse=True)
    
    # Log all successful structures
    logger.info("成功した結晶構造:")
    for i, struct_info in enumerate(successful_structures):
        logger.info(f"  {i+1}. 空間群 {struct_info['space_group']}: "
                   f"密度 = {struct_info['density']:.3f} g/cm³, "
                   f"体積 = {struct_info['volume']:.2f} Å³, "
                   f"vol_factor = {struct_info['vol_factor_used']:.2f}")
    
    # Return structure with highest density
    best_structure_info = successful_structures[0]
    logger.info(f"最高密度構造を選択: 空間群 {best_structure_info['space_group']}, "
               f"密度 = {best_structure_info['density']:.3f} g/cm³")
    
    return best_structure_info['structure']


def generate_molecular_crystal(molecules, atype, nmol, spgnum, vol_factor, mindist_factor, max_attempts=100):
    """Generate molecular crystal structure using pyxtal.
    
    Args:
        molecules: list of pymatgen Molecule objects
        atype: tuple of atom types
        nmol: tuple of number of molecules per type
        spgnum: list of space group numbers to try
        vol_factor: volume factor for structure generation
        mindist_factor: minimum distance factor
        max_attempts: maximum number of attempts (default: 100)
    
    Returns:
        pymatgen Structure object
    """
    
    from contextlib import redirect_stdout, redirect_stderr
    from io import StringIO
    
    # Set minimum distances
    if len(atype) == 1:
        mindist = ((2.0,),)  # 単原子種の場合
    else:
        # 複数原子種の場合、簡単な距離マトリックスを作成
        n_types = len(atype)
        mindist = []
        for i in range(n_types):
            row = []
            for j in range(n_types):
                # 典型的な共有結合半径に基づく最小距離
                if atype[i] == 'H' or atype[j] == 'H':
                    dist = 1.5
                elif atype[i] in ['C', 'N', 'O'] and atype[j] in ['C', 'N', 'O']:
                    dist = 2.5
                else:
                    dist = 3.0
                row.append(dist * mindist_factor)
            mindist.append(tuple(row))
        mindist = tuple(mindist)
    
    logger.info(f"最小距離設定: {mindist}")
    
    # Set tolerance matrix
    tolmat = Tol_matrix(prototype="molecular")
    
    # Generate crystal structure
    attempt = 0
    failed_spgs = {}  # Track failed space groups to avoid repeating
    last_log_attempt = 0
    log_interval = 10  # Log every 10 attempts to reduce output
    
    while attempt < max_attempts:
        try:
            # Choose random space group
            spg = np.random.choice(spgnum)
            
            # Track failures per space group
            if spg not in failed_spgs:
                failed_spgs[spg] = 0
            
            # Skip space groups that failed too many times
            if failed_spgs[spg] >= 5:
                continue
            
            # Log progress less frequently
            if attempt - last_log_attempt >= log_interval:
                logger.info(f"試行 {attempt + 1}/{max_attempts}: 結晶生成を継続中...")
                last_log_attempt = attempt
            elif attempt == 0:
                logger.info(f"試行 {attempt + 1}: 空間群 {spg} で結晶生成中...")
            
            # Create pyxtal structure
            crystal = pyxtal(molecular=True)
            
            # Capture stdout/stderr to detect PyXtal internal failures
            f = StringIO()
            with redirect_stdout(f):
                with redirect_stderr(f):
                    crystal.from_random(
                        dim=3,
                        group=spg,
                        species=molecules,
                        numIons=nmol,
                        factor=vol_factor,
                        conventional=False,
                        tm=tolmat
                    )
            
            output = f.getvalue()
            if output and "Cannot generate crystal after max attempts" in output:
                failed_spgs[spg] += 1
                logger.debug(f"PyXtal内部で最大試行回数に達しました (空間群 {spg}, 試行 {attempt + 1})")
            
            if crystal.valid:
                # Convert to pymatgen structure
                structure = crystal.to_pymatgen()
                logger.info(f"結晶構造生成成功: 空間群 {spg} (試行回数: {attempt + 1})")
                return structure
            else:
                failed_spgs[spg] += 1
                logger.debug(f"無効な結晶構造: 空間群 {spg}")
                
        except Exception as e:
            failed_spgs[spg] = failed_spgs.get(spg, 0) + 1
            logger.debug(f"結晶生成失敗 (空間群 {spg}): {e}")
        
        attempt += 1
    
    # Log summary of failures
    logger.warning(f"{max_attempts}回の試行後も結晶構造生成に失敗しました")
    logger.warning(f"失敗した空間群の統計: {failed_spgs}")
    raise RuntimeError(f"{max_attempts}回の試行後も結晶構造生成に失敗しました")


def optimize_structure(structure, fmax=0.05, steps=1000, calculator='EMT'):
    """Optimize crystal structure using ASE."""
    logger.info("構造最適化を開始...")
    
    # Convert to ASE atoms
    atoms = AseAtomsAdaptor.get_atoms(structure)
    
    # Set calculator
    if calculator == 'EMT':
        atoms.calc = EMT()
    else:
        raise ValueError(f"未対応の計算機: {calculator}")
    
    # Set constraints to maintain symmetry
    atoms.set_constraint([FixSymmetry(atoms)])
    
    # Apply cell filter for cell optimization
    cell_filter = FrechetCellFilter(atoms)
    
    # Optimize
    optimizer = BFGS(cell_filter)
    
    try:
        converged = optimizer.run(fmax=fmax, steps=steps)
        
        # Get optimized structure
        lattice = cell_filter.atoms.cell[:]
        species = cell_filter.atoms.get_chemical_symbols()
        coords = cell_filter.atoms.get_scaled_positions()
        opt_structure = Structure(lattice=lattice, species=species, coords=coords)
        
        energy = cell_filter.atoms.get_total_energy()
        
        logger.info(f"構造最適化完了: 収束={converged}, エネルギー={energy:.4f} eV")
        return opt_structure, energy, converged
        
    except Exception as e:
        logger.error(f"構造最適化失敗: {e}")
        return structure, np.nan, False


def write_cif(structure, output_file, structure_id=1):
    """Write structure to CIF file."""
    try:
        # Use pymatgen CifWriter
        cif_writer = CifWriter(structure)
        cif_string = str(cif_writer)
        
        # Modify the title for identification and add density info
        lines = cif_string.split('\n')
        density = structure.density
        for i, line in enumerate(lines):
            if line.startswith('_chemical_formula_sum'):
                lines[i] = f"_chemical_formula_sum   'Structure_{structure_id}_density_{density:.3f}_g_cm3'"
                break
        
        # Add density as a comment at the beginning
        density_comment = f"# Density: {density:.3f} g/cm³\n# Volume: {structure.volume:.2f} Å³\n"
        cif_string = density_comment + '\n'.join(lines)
        
        # Write to file
        mode = 'w' if structure_id == 1 else 'a'
        with open(output_file, mode) as f:
            f.write(cif_string)
            f.write('\n')
        
        logger.info(f"CIFファイルに書き込み完了: {output_file} (密度: {density:.3f} g/cm³)")
        
    except Exception as e:
        logger.error(f"CIF書き込み失敗: {e}")
        raise


def write_structure_info(structure, nmol, output_file, structure_id=1, symprec=0.01):
    """Write crystal structure information to a text file.
    
    Args:
        structure: pymatgen Structure object
        nmol: tuple of number of molecules per type
        output_file: output file name
        structure_id: structure ID number
        symprec: symmetry precision for space group determination
    """
    try:
        # Get space group information
        analyzer = SpacegroupAnalyzer(structure, symprec=symprec)
        spg_symbol = analyzer.get_space_group_symbol()
        spg_number = analyzer.get_space_group_number()
        
        # Get density
        density = structure.density
        
        # Total number of molecules per unit cell
        total_nmol = sum(nmol)
        
        # Write information
        mode = 'w' if structure_id == 1 else 'a'
        with open(output_file, mode) as f:
            f.write(f"# Structure {structure_id}\n")
            f.write(f"Number of molecules per unit cell: {total_nmol}\n")
            f.write(f"  Molecules per type: {nmol}\n")
            f.write(f"Space group number: {spg_number}\n")
            f.write(f"Space group symbol: {spg_symbol}\n")
            f.write(f"Density: {density:.4f} g/cm³\n")
            f.write(f"Volume: {structure.volume:.4f} Å³\n")
            f.write(f"{'='*60}\n\n")
        
        logger.info(f"構造情報ファイルに書き込み完了: {output_file}")
        
    except Exception as e:
        logger.error(f"構造情報書き込み失敗: {e}")
        raise


def main():
    """Main function."""
    args = parse_args()
    
    # Set up logging
    set_logger(debug=args.debug, logfile='log_mol_crystal')
    logger = getLogger('cryspy')
    
    logger.info("分子性結晶CIF生成スクリプトを開始")
    logger.info(f"入力XYZファイル: {args.xyz_files}")
    logger.info(f"出力CIFファイル: {args.output}")
    
    try:
        # Load molecules
        molecules, atype = load_molecules(args.xyz_files)
        
        # Set number of molecules if not specified
        if args.nmol is None:
            # Default: 4 molecules per type for reasonable packing
            nmol = tuple([4] * len(molecules))
        else:
            if len(args.nmol) != len(molecules):
                raise ValueError(f"--nmolの数({len(args.nmol)})が分子ファイル数({len(molecules)})と一致しません")
            nmol = tuple(args.nmol)
        
        logger.info(f"分子数設定: {nmol}")
        
        if args.optimize_density:
            logger.info("密度最適化モード: 全ての対称性を試して最も密度が高い構造を選択")
            logger.info("構造最適化は無効化されます（分子構造が壊れるのを防ぐため）")
        
        # Determine structure info output file name
        info_output = args.output.replace('.cif', '_info.txt')
        if info_output == args.output:
            info_output = args.output + '_info.txt'
        
        # Generate structures
        for i in range(args.nstruct):
            logger.info(f"構造 {i+1}/{args.nstruct} を生成中...")
            
            # Generate initial crystal structure
            if args.optimize_density:
                structure = generate_molecular_crystal_all_symmetries(
                    molecules, atype, nmol, args.spgnum, 
                    args.vol_factor, args.mindist_factor
                )
            else:
                structure = generate_molecular_crystal(
                    molecules, atype, nmol, args.spgnum, 
                    args.vol_factor, args.mindist_factor
                )
            
            # Optimize structure if requested (但し、--optimize-densityが指定されている場合は無効化)
            if not args.no_optimization and not args.optimize_density:
                structure, energy, converged = optimize_structure(
                    structure, args.fmax, args.steps, args.calculator
                )
                if not converged:
                    logger.warning("構造最適化が収束しませんでした")
            elif args.optimize_density:
                logger.info("密度最適化モードのため構造最適化をスキップ")
            
            # Write to CIF file
            write_cif(structure, args.output, i + 1)
            
            # Write structure information to text file
            write_structure_info(structure, nmol, info_output, i + 1)
        
        logger.info(f"全ての構造生成が完了しました: {args.output}")
        logger.info(f"構造情報ファイル: {info_output}")
        
    except Exception as e:
        logger.error(f"エラーが発生しました: {e}")
        raise


if __name__ == '__main__':
    main()