#!/usr/bin/env python
"""
The :mod:`hotspots.calculation` handles the main Fragment Hotspot Maps algorithm. In addition, an alternative pocket burial method, Ghecom, is provided.

The main classes of the :mod:`hotspots.calculation` module are:

- :class:`hotspots.calculation.Buriedness`
- :class:`hotspots.calculation.Runner`

More information about the Fragment Hotspot Maps method is available from:
    - Radoux, C.J. et. al., Identifying the Interactions that Determine Fragment Binding at Protein Hotspots J. Med. Chem. 2016, 59 (9), 4314-4325 [dx.doi.org/10.1021/acs.jmedchem.5b01980]

More information about the Ghecom method is available from:
    - Kawabata T, Go N. Detection of pockets on protein surfaces using small and large probe spheres to find putative ligand binding sites. Proteins 2007; 68: 516-529
"""
from __future__ import print_function, division

import multiprocessing
import operator
import random
import shutil
import sys
import os
import tempfile
import time
from functools import reduce
from os import system, environ
from os.path import join, dirname

ghecom =join(dirname(dirname(__file__)),"win_ghecom/Ghecom.exe")
if os.path.exists(ghecom) and not 'GHECOM_EXE' in environ:
    environ['GHECOM_EXE'] = ghecom

import numpy as np
import pkg_resources

from scipy import ndimage
from skimage.morphology import ball
from skimage.transform import resize
# from hotspots.protoss import Protoss
from tqdm import tqdm
from time import perf_counter

# Development note. We import the CCDC modules after 3rd part modules where possible as
# we generally find that the ccdc package is less fussy about the underlying compiled
# libraries used.
#  

from ccdc.cavity import Cavity
from ccdc.io import MoleculeWriter, MoleculeReader
from ccdc.molecule import Molecule, Coordinates
from ccdc.protein import Protein
from ccdc.utilities import PushDir

from hotspots.atomic_hotspot_calculation import _AtomicHotspot
from hotspots.grid_extension import Grid
from hotspots.hs_utilities import Helper
from hotspots.wrapper_pdb import PDBResult
from hotspots.wrapper_ghecom import Ghecom
from hotspots.result import Results


class ExpBuriedness(object):

    def __init__(self, prot, out_grid, max_probe_radius=11):
        self.prot = prot
        self.out_grid = out_grid
        self.max_probe = max_probe_radius
        self.probe_selem_dict = self._generate_probe_selem()
        self.protein_grid = self._from_molecule(prot)
        if self.out_grid is None:
            self.out_grid = self.protein_grid.copy_and_clear()
        print("buriedness init")


    def _generate_probe_selem(self):
        """
        Generates a set of structuring elements (i.e. probe spheres) starting from a 3x3 cross (0.5 A radius) probe up
        to the max probe radius
        :return:
        """
        probe_selem_dict = {}

        # Small probe

        probe_selem_dict[2] = ball(2)
        # r= probe radius
        for r in range(3,self.max_probe +1):
            probe_selem_dict[r] = ball(r)
        print("generating probes")
        return probe_selem_dict

    def _from_molecule(self, mol, scaling=1):
        """
        generate a molecule mask where gp within the vdw radius of the molecule heavy atoms are set to 1.0
        :param mol: `ccdc.molecule.Molecule`
        :param padding: int
        :param scaling: float
        :return: `hotspots.grid_extension.Grid`
        """

        coords = [a.coordinates for a in mol.atoms]
        g = Grid.initalise_grid(coords=coords, padding= 10, spacing=1)

        for probe in sorted(self.probe_selem_dict.keys(), reverse=True):
            for a in mol.heavy_atoms:
                g.set_sphere(point=a.coordinates,
                             radius=probe * scaling,
                             value=probe,
                             mode='replace',
                             scaling='None')

        for a in mol.heavy_atoms:
            g.set_sphere(point=a.coordinates,
                         radius=a.vdw_radius,
                         value=100,
                         mode='replace',
                         scaling='None')

        try:
            out_bound_box = self.out_grid.bounding_box
            print(out_bound_box)
            origin_indices = g.point_to_indices(out_bound_box[0])
            far_indices = g.point_to_indices(out_bound_box[1])
            region = origin_indices + far_indices
            print(region)
            return g.sub_grid(region)
        except AttributeError:
            return g

    def _close_grid(self, g, probe):

        g_array = g.get_array().astype(int)
        closed_array = ndimage.binary_erosion(g_array, structure=self.probe_selem_dict[probe])
        return Grid.array_to_grid(closed_array.astype(int), g)

    def _multiscale_closing(self,g):


        #probe_sizes.reverse()
        all_g = None
        prot_mask = g > 90

        for probe in sorted(self.probe_selem_dict.keys())[1:]:
            if all_g is None:
                all_g = prot_mask.copy()

            probe_mask = ((g < probe+0.1) & (g>0)) | (prot_mask)


            probe_mask = self._close_grid(probe_mask,probe)
            novel = (-all_g) * (g>0)
            all_g += (probe_mask * novel * (self.max_probe - probe))

        return all_g * (prot_mask < 1)

    def _open_grid(self, g, probe):

        g_array = g.get_array().astype(int)
        opened_array = ndimage.binary_opening(g_array, structure=self.probe_selem_dict[probe])
        return Grid.array_to_grid(opened_array.astype(int), g)

    def buriedness_grid(self):

        closed_g = self._multiscale_closing(self.protein_grid)
        print("close_g")
        out_g = self._open_grid(closed_g,2) * closed_g
        print("out_g")
        out_array = out_g.get_array()
        print("out array")
        scaled_g = Grid.initalise_grid(self.out_grid.bounding_box, padding=0, spacing=0.5)
        print("scaled_g")
        scaled_array = resize(out_array, scaled_g.nsteps ,anti_aliasing=False)

        # Future tweaking here
        final_array = scaled_array
        print("array to grid")
        return Grid.array_to_grid(final_array.astype(int), scaled_g)


class _WeightedResult(object):
    """
    private class

    A class to handle the weighted grid results

    :param str identifier: identifier, default is the probe identifier assigned at the Atomic Hotspot Calculation stage
    :param `ccdc.utilities.Grid` grid: result grid
    """

    def __init__(self, identifier, grid):
        self.identifier = identifier
        self.grid = grid


class _SampleGrid(object):
    """
    private class

    class to handle sampled grids

    :param name: str, name of probe (donor, acceptor, apolar, positive, negative)
    :param grid: a :class: `ccdc.utilities.Grid` instance
    :param atom_predicate: atom_predicate will be used to select atoms of a molecule for sampling
    """

    def __init__(self, name, grid, atom_predicate):
        self.name = name
        self.grid = grid
        self.atom_predicate = atom_predicate
        self.mol = None

    @staticmethod
    def add_coordinates(coord, trans):
        """
        provides a coordinate list of atoms to be scored be scored in the SampleGrid
        :param coord: tup, (float(x), float(y), float(z)), set of atomic coordinates for "active" coordinates
        :param trans: tup, (float(x), float(y), float(z)), set of translations to translate probes to points
        above threshold
        :return: list of tup
        """

        return [coord[i] + trans[i] for i in range(0, len(trans))]

    def sample(self, coordinate_list, trans):
        """
        score Molecule in grids for which it has active atoms
        :param coordinate_list: list, set of coordinates of translated and rotated probe atoms to be scored
        :param trans:
        :return:
        """

        try:
            return [self.grid._grid.value(self.add_coordinates(c, trans)) for c in coordinate_list]
        except RuntimeError:
            return [0]

    def set_molecule(self, mol, polar_contribution):
        """
        set which atoms match the grid
        probes with polar atoms contribute do not contribute to apolar maps as this leads to artefacts

        :param mol:
        :param polar_contribution:
        :return:
        """
        self.mol = mol
        if self.name == 'apolar' and not polar_contribution and len([a for a in mol.atoms
                                                                     if a.is_donor or a.is_acceptor]) > 0:
            self._active_atoms = []
        elif not hasattr(self, '_active_atoms'):
            self._active_atoms = [a for a in self.mol.atoms if self.atom_predicate(a)]

    @staticmethod
    def is_donor(a):
        """
        returns true if a given atom is a donor

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom classification is "donor"
        """

        if a.is_donor and a.formal_charge == 0:
            return True
        else:
            return False

    @staticmethod
    def is_acceptor(a):
        """
        returns true if a given atom is a acceptor

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom classification is "acceptor"
        """
        if a.is_acceptor and a.formal_charge == 0:
            return True
        else:
            return False

    @staticmethod
    def is_apolar(a):
        """
        returns true if a given atom is a apolar

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom classification is "apolar"
        """
        if a.is_donor or a.is_acceptor or a.formal_charge != 0 or a.atomic_symbol == "Xe":
            return False
        else:
            return True

    @staticmethod
    def is_positive(a):
        """
        returns true if a given atom is a positively charged

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom is positively charged
        """
        return a.formal_charge > 0

    @staticmethod
    def is_negative(a):
        """
        returns true if a given atom is a negatively charged

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom is negatively charged
        """
        return a.formal_charge < 0

    @staticmethod
    def is_aromatic(a):
        """
        returns true if a given atom is aromatic

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom is aromatic
        """
        return a.atomic_symbol == 'C' and a.is_cyclic and any(b.atom_type == 'aromatic' for b in a.bonds)


class Runner(object):
    """
    A class for running the Fragment Hotspot Map calculation
    """

    class Settings(object):
        """
        adjusts the default settings for the calculation

        :param int nrotations: number of rotations (keep it below 10**6)
        :param float apolar_translation_threshold: translate probe to grid points above this threshold. Give lower values for greater sampling. Default 15
        :param float polar_translation_threshold: translate probe to grid points above this threshold. Give lower values for greater sampling. Default 15
        :param bool polar_contributions: allow carbon atoms of probes with polar atoms to contribute to the apolar output map.
        :param bool return_probes: Generate a sorted list of molecule objects, corresponding to probe poses
        :param bool sphere_maps: When setting the probe score on the output maps, set it for a sphere (radius 1.5) instead of a single point.
        :param bool fix_seeds: use fixed random seeds when sampling map
        """

        def __init__(self, nrotations=3000, apolar_translation_threshold=15, polar_translation_threshold=15,
                     polar_contributions=False, return_probes=False, sphere_maps=False, fix_seeds=False):
            self.nrotations = nrotations
            self.apolar_translation_threshold = apolar_translation_threshold
            self.polar_translation_threshold = polar_translation_threshold
            self.polar_contributions = polar_contributions
            self.return_probes = return_probes
            self.sphere_maps = sphere_maps
            self.fix_seeds = fix_seeds

        @property
        def _num_gp(self):
            """
            number of grid point for a given volume
            :return:
            """
            return int(float(500) / self.nrotations ** 3)


    class _Sampler(object):
        """
        Samples one or more grids with a probe molecule
        """

        def __init__(self, *grids, **kw):
            """
            Settings used to run fragment-hotspot-maps script

            :param grids: list, list of :class:`ccdc.utilities.Grid` instances
            :param kw: 'settings'
            """

            self.grids = grids
            self.grid_dic = {}
            for g in grids:
                self.grid_dic[g.name] = g
            settings = kw.get('settings')
            if settings is None:
                settings = self.Settings()

            self.settings = settings
            self.probe_grids = [_SampleGrid(g.name, g.grid.copy_and_clear(), g.atom_predicate) for g in self.grids]

        def get_priority_atom(self, molecule):
            """
            Select priority atom. Select polar atom. If multiple polar atoms, select the one furthest from the centre of
            geometry. If no polar atoms, select atom furthest from centre of geometry

            :param molecule: a :class: `ccdc.molecule.Molecule` instance
            :return: a :class: `ccdc.molecule.Molecule` instance, str atom type
            """
            c = molecule.centre_of_geometry()
            polar_atoms = [a for a in molecule.atoms if a.is_donor or a.is_acceptor]
            atom_by_distance = {}
            if len(polar_atoms) > 0:
                for a in polar_atoms:
                    d = Helper.get_distance(c, a.coordinates)
                    atom_by_distance[d] = a
            else:
                for a in molecule.atoms:
                    d = Helper.get_distance(c, a.coordinates)
                    atom_by_distance[d] = a

            greatest_distance = sorted(atom_by_distance.keys())[0]
            priority_atom = atom_by_distance[greatest_distance]

            pa_type = None
            if priority_atom.formal_charge != 0:
                if priority_atom.formal_charge < 0:
                    pa_type = "acceptor"
                elif priority_atom.formal_charge > 0:
                    pa_type = "donor"
            else:
                if priority_atom.is_acceptor:
                    pa_type = "acceptor"
                elif priority_atom.is_donor:
                    pa_type = "donor"
                else:
                    pa_type = "apolar"

            return priority_atom, pa_type

        def get_translation_points(self, priority_atom_type):
            """
            returns a list of coordinates that are greater than the threshold, that the probe will be translated to

            :param priority_atom_type: str, atomic interaction type
            :return: list, list of :class: `ccdc.molecule.Molecule` instances
            """
            translate_probe = []
            wg = self.grid_dic[priority_atom_type]

            if priority_atom_type == 'apolar':
                translation_threshold = self.settings.apolar_translation_threshold
            else:
                translation_threshold = self.settings.polar_translation_threshold
            hs = wg.grid.islands(translation_threshold)
            for g in hs:
                nx, ny, nz = g.nsteps

                maxima = [g.indices_to_point(i, j, k)
                          for i in range(nx) for j in range(ny) for k in range(nz)
                          if g.value(i, j, k) >= translation_threshold]

                translate_probe = translate_probe + maxima
            return translate_probe

        def generate_rand_quaternions(self):
            """
            Returns a list of random quaternions. Length matches settings.nrotations

            :return: tup, (a,b,c,d)
            """
            quaternions = []

            # Based on averages from several runs
            seed = 1741620750.3308594
            delta = 0.000005472091

            i = 1
            if self.settings.nrotations > 1:
                while i <= self.settings.nrotations:
                    if self.settings.fix_seeds:
                        seed += delta
                        random.seed(seed)
                    r1 = random.uniform(-1, 1)

                    if self.settings.fix_seeds:                  
                        seed += delta
                        random.seed(seed)
                    r2 = random.uniform(-1, 1)

                    s1 = r1 * r1 + r2 * r2
                    if s1 < 1:
                        if self.settings.fix_seeds:
                            seed += delta
                            random.seed(seed)
                        r3 = random.uniform(-1, 1)

                        if self.settings.fix_seeds:
                            seed += delta
                            random.seed(seed)
                        r4 = random.uniform(-1, 1)

                        s2 = r3 * r3 + r4 * r4
                        if s2 < 1:
                            q = (r1, r2, r3 * (np.sqrt((1 - s1) / s2)), r4 * (np.sqrt((1 - s1) / s2)))
                            quaternions.append(q)

                            i += 1

            return quaternions

        @staticmethod
        def score(values):
            """
            Calculate geometric mean of scores

            :param values: float, scores of atoms in probe
            :return: float, geometric mean of probe atom scores
            """

            return reduce(operator.__mul__, values, 1.0) ** (1. / len(values))
            # return sum(values)

        def sample_pose(self, trans, active_atoms_dic, probe):
            """
            Return a pose score (as defined by the score(self,dic) function) and a dictionary of atom:scores

            :param trans: list of translations
            :param active_atoms_dic: dict {"interaction type": "atoms to be scored"}
            :param probe: str, interaction_type
            :return:
            """
            if probe == "negative" or probe == "positive":
                atoms = [g.mol.atoms for g in self.grids if g.name == "positive"][0]
                weight = int(6 / len([a for a in atoms if str(a.atomic_symbol) == "C"]))

                apolar_values = [g.sample(active_atoms_dic[g.name], trans) for g in self.grids if
                                 len(active_atoms_dic[g.name]) > 0 and g.name == "apolar"] * weight
                charged_values = [g.sample(active_atoms_dic[g.name], trans) for g in self.grids if
                                  len(active_atoms_dic[g.name]) > 0 and g.name != "apolar"]
                values = apolar_values + charged_values

            else:
                values = [g.sample(active_atoms_dic[g.name], trans) for g in self.grids if
                          len(active_atoms_dic[g.name]) > 0]
            scores = [item for sublist in values for item in sublist]
            return self.score(scores)

        def update_out_grids(self, score, active_coordinates_dic, trans):
            """
            For active atoms for a given grid, set closest grid point value to score, unless already set to a higher
            value

            :param score: float, score of a given probe
            :param active_coordinates_dic:
            :param trans: list of translations
            :return:
            """

            for pg in self.probe_grids:
                actives = active_coordinates_dic[pg.name]
                if len(actives) == 0:
                    continue
                for a in actives:
                    coords = pg.add_coordinates(a, trans)
                    i, j, k = pg.grid.point_to_indices(coords)
                    orig_value = pg.grid.value(i, j, k)
                    #
                    if self.settings.sphere_maps:
                        # pg.grid.set_sphere(coords, 2, score, mode='max')
                        if score > orig_value:
                            pg.grid.set_sphere(coords, 1.5, score, mode='max', scaling='None')
                    else:
                        pg.grid.set_value(i, j, k, max(score, orig_value))

        def get_active_coordinates(self):
            """
            Returns a dictionary of {grid_name:[Coordinates]}

            :return:
            """
            active_coords_dic = {
                g.name: [a.coordinates for a in g._active_atoms]
                for g in self.grids
            }

            return active_coords_dic

        def sample(self, molecule, probe, update_grids=True):
            """
            Sample the grids according to the settings

            :param molecule:
            :param probe: str, interaction type, (donor, acceptor, negative, positive, apolar)
            :return:
            """
            priority_atom, priority_atom_type = self.get_priority_atom(molecule)
            translate_points = self.get_translation_points(priority_atom_type)
            molecule.remove_hydrogens()
            quaternions = self.generate_rand_quaternions()
            high_scoring_probes = {}
            print("\n    nRotations:", len(quaternions), "nTranslations:", len(translate_points), "probename:", probe)

            for g in self.grids:
                g.set_molecule(molecule, True)

            for g in self.probe_grids:
                g.set_molecule(molecule, self.settings.polar_contributions)

            for q in tqdm(quaternions):
                molecule.apply_quaternion(q)
                priority_atom_coordinates = priority_atom.coordinates
                active_coordinates_dic = self.get_active_coordinates()

                for priority_atom_point in translate_points:

                    translation = [priority_atom_point[i] - priority_atom_coordinates[i]
                                   for i in range(0, len(priority_atom_coordinates))]

                    score = self.sample_pose(translation, active_coordinates_dic, probe)

                    try:
                        if score < 1:
                            continue
                    except TypeError:
                        continue
                    if update_grids:
                        self.update_out_grids(score, active_coordinates_dic, translation)

                    if self.settings.return_probes is True:
                        if score > 10:
                            m = molecule.copy()
                            m.translate(translation)
                            m.identifier = "{}".format(score)

                            try:
                                high_scoring_probes[score].append(m)
                            except KeyError:
                                high_scoring_probes[score] = [m]

            if self.settings.return_probes is True:
                sampled_probes = []
                for key in sorted(high_scoring_probes.keys(), reverse=True)[:100]:
                    sampled_probes.extend(high_scoring_probes[key])
                print('Returned probes = ', len(sampled_probes))
                return sampled_probes

                # if len(sampled_probes) > 100:
                #     return sampled_probes[:100]
                # else:
                #     return sampled_probes

    def __init__(self, settings=None, output_path=None, pdb_id=None, save_cavity=False):
        self.out_grids = {}
        self.super_grids = {}
        self.buriedness = None
        self.sampled_probes = {}

        if settings is None:
            self.sampler_settings = self.Settings()
        else:
            self.sampler_settings = settings
            
        self.output_path = output_path
        self.pdb_id = pdb_id
        self.save_cavity = save_cavity

    @property
    def protein(self):
        """
        protein submitted for calculation
        :return:
        """
        return self._protein

    @protein.setter
    def protein(self, prot):
        """
        protein submitted for calculation
        :param `ccdc.protein.Protein` prot: protein
        :return:
        """
        if isinstance(prot, Protein):
            self._protein = prot
        else:
            raise TypeError("`ccdc.protein.Protein` must be supplied. Hint: Use Protein.from_file()")

    @property
    def charged_probes(self):
        """
        optional settings, True if charged features are to be calculated
        :return:
        """
        return self._charged_probes

    @charged_probes.setter
    def charged_probes(self, option):
        """
        optional settings, True if charged features are to be calculated
        :param bool option: (default = True)
        :return:
        """
        if type(option) is bool:
            self._charged_probes = option
        else:
            raise TypeError("Expecting a bool, got {} instead".format(type(option)))

    @property
    def probe_size(self):
        """
        optional settings, change probe size. NB: feature has not been tested.
        :return:
        """
        return self._probe_size

    @probe_size.setter
    def probe_size(self, size):
        """
        optional settings, change probe size. NB: feature has not been tested.
        :param int size: number of heavy atoms contain within the molecular probes (integer between 3 and 8)
        :return:
        """
        if size in range(3, 8):
            self._probe_size = size
        else:
            raise ValueError("Probe size must be an integer between 3-7")

    @property
    def buriedness_method(self):
        """
        optional settings, pocket detection method. (default = "ghecom")
        :return:
        """
        return self._buriedness_method

    @buriedness_method.setter
    def buriedness_method(self, method):
        """
        optional settings, pocket detection method. (default = "ghecom")
        :param str method: either 'ghecom' or 'ligsite'
        :return:
        """
        method = method.lower()
        if method == 'ghecom':

            if 'GHECOM_EXE' in environ:
                self._buriedness_method = method
            else:
                raise EnvironmentError("Must set Ghecom environment variable to use Ghecom")

        elif method == 'ghecom_internal':
            self._buriedness_method = method

        elif method == 'ligsite':
            if sys.platform == 'linux' or sys.platform == 'linux2':
                print("RECOMMENDATION: you have chosen LIGSITE as buriedness method, ghecom is recommended")
            self._buriedness_method = method

        else:
            raise ValueError("Buriedness method must be 'ghecom' (default) or 'ligsite")

    @property
    def cavities(self):
        """
        optional settings, cavities supplied to the calculation
        :return:
        """
        return self._cavities

    @cavities.setter
    def cavities(self, obj):
        """
        optional settings, cavities supplied to the calculation
        :param list or Coordinate or `ccdc.molecule.Molecule` or `ccdc.cavity.Cavity`: cavity information provided
        :return:
        """
        if obj is not None:
            if isinstance(obj, list) or isinstance(obj, tuple):
                if isinstance(obj, Coordinates):
                    try:
                        print(obj.x)
                        self._cavities = [obj]
                    except AttributeError:
                        self._cavities = obj

                    self._cavities = [obj]
                elif isinstance(obj[0], Molecule):
                    self._cavities = [m.centre_of_geometry() for m in obj]
                elif isinstance(obj[0], Cavity):
                    self._cavities = [Helper.cavity_centroid(c) for c in obj]
                else:
                    try:
                        self._cavities = [Coordinates(x=obj[0], y=obj[1], z=obj[2])]
                    except RuntimeError:
                        print("WARNING! Failed to detected cavity, Atomic Hotspot detection to run on whole protein")
                        self._cavities = None

            elif isinstance(obj, Molecule):
                self._cavities = [obj.centre_of_geometry()]
            elif isinstance(obj, Cavity):
                self._cavities = [Helper.cavity_centroid(obj)]

            else:
                print("WARNING! Failed to detected cavity, Atomic Hotspot detection to run on whole protein")
                self._cavities = None
        else:
            self._cavities = None

    @property
    def nprocesses(self):
        """
        number of processes used in parallelisation
        :return:
        """
        return self._nprocesses

    @nprocesses.setter
    def nprocesses(self, num):
        """
        number of processes used in parallelisation
        :param int num: number of processor
        :return:
        """
        num = int(num)
        if num in range(0, int(multiprocessing.cpu_count())+1):
            self._nprocesses = num
        else:
            raise OSError("CPU count = {}".format(multiprocessing.cpu_count()))

    @property
    def sampler_settings(self):
        """
        sampler settings
        :return:
        """
        return self._sampler_settings

    @sampler_settings.setter
    def sampler_settings(self, settings):
        """
        sampler settings
        :param `hotspots.calculation.Runner.Settings` settings: adjust run settings
        :return:
        """
        if isinstance(settings, self.Settings):
            self._sampler_settings = settings
        else:
            self._sampler_settings = None

    def _get_weighted_maps(self):
        """
        private method

        weight superstar output by burriedness
        :return: a list of :class:`WeightedResult` instances
        """
        results = []
        for s in self.superstar_grids:
            g, b = Grid.common_grid([s.grid, self.buriedness], padding=1)
            weighted_grid = g * b
            results.append(_WeightedResult(s.identifier, weighted_grid))

        return results

    def _dock_probe(self,mol, grid_dict):

        donor_grid = _SampleGrid('donor', grid_dict['donor'], _SampleGrid.is_donor)
        acceptor_grid = _SampleGrid('acceptor', grid_dict['acceptor'], _SampleGrid.is_acceptor)
        apolar_grid = _SampleGrid('apolar', grid_dict['apolar'], _SampleGrid.is_apolar)

        kw = {'settings': self.sampler_settings}
        self.sampler = self._Sampler(apolar_grid, donor_grid, acceptor_grid, **kw)

        probes = self.sampler.sample(mol, probe='probe', update_grids=False)
        return probes

    def _get_out_maps(self, probe, grid_dict, return_probes=False):
        """
        private method

        organises the sampling of weighted superstar maps by molecular probes
        :param str probe: probe identifier set in the Atomic Hotspot calculation
        :param dict grid_dict: dictionary with key = probe identifier and value = `hotspots.grid_extension.Grid`
        :param bool return_probes: optional, bool indicating if probe molecules should be returned
        :return:
        """
        donor_grid = _SampleGrid('donor', grid_dict['donor'], _SampleGrid.is_donor)
        acceptor_grid = _SampleGrid('acceptor', grid_dict['acceptor'], _SampleGrid.is_acceptor)
        apolar_grid = _SampleGrid('apolar', grid_dict['apolar'], _SampleGrid.is_apolar)

        if self.charged_probes:
            negative_grid = _SampleGrid('negative', grid_dict['negative'], _SampleGrid.is_negative)
            positive_grid = _SampleGrid('positive', grid_dict['positive'], _SampleGrid.is_positive)

        kw = {'settings': self.sampler_settings}
        if self.charged_probes:
            self.sampler = self._Sampler(apolar_grid, donor_grid, acceptor_grid, negative_grid, positive_grid, **kw)
        else:
            self.sampler = self._Sampler(apolar_grid, donor_grid, acceptor_grid, **kw)

        probe_path = pkg_resources.resource_filename('hotspots', 'probes/')

        if self.charged_probes:
            if probe == "negative" or probe == "positive":
                mol = MoleculeReader(join(probe_path, "rotate-{}_{}_flat.mol2".format(probe, "test")))[0]
            else:
                mol = MoleculeReader(join(probe_path, "rotate-{}_{}_flat.mol2".format(probe, self.probe_size)))[0]
        else:
            mol = MoleculeReader(join(probe_path, "rotate-{}_{}_flat.mol2".format(probe, self.probe_size)))[0]

        probes = self.sampler.sample(mol, probe=probe)

        for pg in self.sampler.probe_grids:
            if pg.name.lower() == probe:
                try:
                    self.out_grids[pg.name].append(pg.grid)
                except KeyError:
                    self.out_grids[pg.name] = [pg.grid]

        if return_probes is True:
            return probes

    def _calc_hotspots(self, return_probes=False):
        """
        handles the organisation of the hotspot calculation
        :param return_probes: optional, bool indicating if probe molecules should be returned
        :return:
        """
        print("Start atomic hotspot detection\n        Processors: {}".format(self.nprocesses))
        t_ah1 = perf_counter()
        a = _AtomicHotspot()
        a.settings.atomic_probes = {"apolar": "AROMATIC CH CARBON",
                                    "donor": "UNCHARGED NH NITROGEN",
                                    "acceptor": "CARBONYL OXYGEN"}
        if self.charged_probes:
            a.settings.atomic_probes = {"negative": "CARBOXYLATE OXYGEN", "positive": "CHARGED NH NITROGEN"}

        probe_types = a.settings.atomic_probes.keys()
        self.superstar_grids = a.calculate(protein=self.protein,
                                           nthreads=self.nprocesses,
                                           cavity_origins=self.cavities)

        if self.clear_tmp:
            shutil.rmtree(a.settings.temp_dir)

        t_ah2 = perf_counter()
        print(f"Atomic hotspot detection complete in {t_ah2-t_ah1:.2f} seconds\n")

        print("Start buriedness calculation")
        t_bur = perf_counter()
        if self.buriedness_method.lower() == 'ghecom' and self.buriedness is None:
            print("    method: Ghecom")
            out_grid = self.superstar_grids[0].buriedness.copy_and_clear()
            
            cavity_file = "ghecom_out.pdb"
            if self.pdb_id:
                cavity_file = f"{self.pdb_id}_{cavity_file}"
            if self.output_path:
                cavity_file = os.path.join(self.output_path, cavity_file)

            ghecom = Ghecom()
            
            path = cavity_file if self.save_cavity and os.path.isfile(cavity_file) else None
            if not path:
                path = ghecom.run(self.protein)
                if self.save_cavity:
                    shutil.copy2(path, cavity_file)
            else:
                print(f"    reusing cavity file: {path}")


            t_p2g = perf_counter()
            self.buriedness = ghecom.pdb_to_grid(path, out_grid)
            print(f"Time to convert pdb to grid = {perf_counter()-t_p2g:.2f} seconds")

            shutil.rmtree(ghecom.temp) # probs not needed

        elif self.buriedness_method.lower() == 'ghecom_internal' and self.buriedness is None:
            print("    method: Internal version Ghecom")
            out_grid = self.superstar_grids[0].buriedness.copy_and_clear()
            b = ExpBuriedness(prot=self.protein, out_grid=None)
            self.buriedness = b.buriedness_grid()

        elif self.buriedness_method.lower() == 'ligsite' and self.buriedness is None:
            print("    method: LIGSITE")
            self.buriedness = Grid.get_single_grid(grd_dict={s.identifier: s.buriedness for s in self.superstar_grids},
                                                   mask=False)

        self.weighted_grids = self._get_weighted_maps()

        print(f"Buriedness calculation complete in {perf_counter() - t_bur:.2f} seconds\n")

        print("Start sampling")
        t_samp = perf_counter()
        grid_dict = {w.identifier: w.grid for w in self.weighted_grids}

        for probe in probe_types:
            t_probe = perf_counter()
            if return_probes is True:
                ps = self._get_out_maps(probe, grid_dict, return_probes=True)
                print(len(ps))
                self.sampled_probes.update({probe: ps})

            else:
                self._get_out_maps(probe, grid_dict)
            print(f"Sampling with {probe} complete in {perf_counter()-t_probe:.1f} seconds")

        print(f"All sampling complete in {perf_counter()-t_samp:.1f} seconds\n")

    def _prepare_protein(self, protoss=False):
        """
        default protein preparation settings on the protein
        :return:
        """
        self.protein.remove_all_waters()
        for lig in self.protein.ligands:
            self.protein.remove_ligand(lig.identifier)
        self.protein.remove_all_metals()

        if protoss is False:
            self.protein.add_hydrogens()

    def from_superstar(self, protein, superstar_grids, buriedness, charged_probes=False, probe_size=7,
                       settings=None, clear_tmp=True):
        """
        calculate hotspot maps from precalculated superstar maps. This enables more effective parallelisation and reuse
        of object such as the Buriedness grids

        :param protein: a :class:`ccdc.protein.Protein` instance
        :param superstar_grids: a :class:`hotspots.atomic_hotspot_calculation._AtomicHotspotResult` instance
        :param buriedness: a :class:`hotspots.grid_extension.Grid` instance
        :param bool charged_probes: If True, include positive and negative probes
        :param int probe_size: Size of probe in number of heavy atoms (3-8 atoms)
        :param settings: `hotspots.calculation.Runner.Settings` settings: holds the sampler settings
        :param bool clear_tmp: If True, clear the temporary directory
        :return:
        """
        start = time.time()
        self.super_grids = {}
        self.superstar_grids = superstar_grids
        self.probe_types = [p.identifier for p in self.superstar_grids]
        self.buriedness = buriedness

        self.protein = protein
        self.charged_probes = charged_probes
        self.probe_size = probe_size
        self.clear_tmp = clear_tmp

        if settings is None:
            self.sampler_settings = self.Settings()
        else:
            self.sampler_settings = settings

        self.weighted_grids = self._get_weighted_maps()

        print("Start sampling")
        grid_dict = {w.identifier: w.grid for w in self.weighted_grids}

        for probe in self.probe_types:
            self._get_out_maps(probe, grid_dict)

        self.super_grids = {p: g[0] for p, g in self.out_grids.items()}

        print("Sampling complete\n")
        print("Runtime = {}seconds".format(time.time() - start))

        return Results(super_grids=self.super_grids,
                       protein=self.protein,
                       buriedness=self.buriedness)


    def from_protein(self, protein, charged_probes=False, probe_size=7, buriedness_method='ghecom',
                     cavities=None, nprocesses=1, settings=None, buriedness_grid=None, clear_tmp=True):
        """
        generates a result from a protein

        :param protein: a :class:`ccdc.protein.Protein` instance
        :param bool charged_probes: If True include positive and negative probes
        :param int probe_size: Size of probe in number of heavy atoms (3-8 atoms)
        :param str buriedness_method: Either 'ghecom' or 'ligsite'
        :param cavities: Coordinate or `ccdc.cavity.Cavity` or `ccdc.molecule.Molecule` or list specifying the cavity or cavities on which the calculation should be run
        :param int nprocesses: number of CPU's used
        :param `hotspots.calculation.Runner.Settings` settings: holds the sampler settings
        :param `ccdc.utilities.Grid` buriedness_grid: pre-calculated buriedness grid
        :return: a :class:`hotspots.result.Results` instance


        >>> from ccdc.protein import Protein
        >>> from hotspots.calculation import Runner

        >>> protein = Protein.from_file(<path_to_protein>)

        >>> runner = Runner()
        >>> settings = Runner.Settings()
        >>> settings.nrotations = 1000  # fewer rotations increase speed at the expense of accuracy
        >>> runner.from_protein(protein, nprocesses=3, settings=settings)
        Result()

        """
        start = time.time()
        self.super_grids = {}
        self.buriedness = buriedness_grid
        self.protein = protein
        self.charged_probes = charged_probes
        self.probe_size = probe_size
        self.buriedness_method = buriedness_method
        self.cavities = cavities
        self.clear_tmp = clear_tmp

        self.nprocesses = nprocesses
        if settings is None:
            self.sampler_settings = self.Settings()
        else:
            self.sampler_settings = settings
        self._calc_hotspots()  # return probes = False by default
        self.super_grids = {p: g[0] for p, g in self.out_grids.items()}
        print(f"Runtime for 'from_protein' function = {time.time() - start:.2f} seconds")

        return Results(super_grids=self.super_grids,
                       protein=self.protein,
                       buriedness=self.buriedness,
                       superstar={x.identifier: x.grid for x in self.superstar_grids},
                       weighted_superstar={x.identifier: x.grid for x in self.weighted_grids})

    def from_pdb(self, pdb_code, charged_probes=False, probe_size=7, buriedness_method='ghecom', nprocesses=3,
                 cavities=False, settings=None, clear_tmp=True):
        """
        generates a result from a pdb code

        :param str pdb_code: PDB code
        :param bool charged_probes: If True include positive and negative probes
        :param int probe_size: Size of probe in number of heavy atoms (3-8 atoms)
        :param str buriedness_method: Either 'ghecom' or 'ligsite'
        :param int nprocesses: number of CPU's used
        :param `hotspots.calculation.Runner.Settings` settings: holds the calculation settings
        :return: a :class:`hotspots.result.Result` instance


        >>> from hotspots.calculation import Runner

        >>> runner = Runner()
        >>> runner.from_pdb("1hcl")
        Result()

        """
        protoss = False

        tmp = tempfile.mkdtemp()
        if clear_tmp is False:
            print("Detailed output written to",tmp)
        # if  protoss is True:
        #     protoss = Protoss(out_dir=tmp)
        #     self.protein = protoss.add_hydrogens(pdb_code).protein
        #
        # else:
        PDBResult(identifier=pdb_code).download(out_dir=tmp)
        fname = join(tmp, "{}.pdb".format(pdb_code))
        self.protein = Protein.from_file(fname)
        self._prepare_protein(protoss)
        self.charged_probes = charged_probes
        self.probe_size = probe_size
        self.buriedness_method = buriedness_method
        self.clear_tmp = clear_tmp
        self.cavities = None
        if cavities is True:
            self.cavities = Cavity.from_pdb_file(fname)
        self.nprocesses = nprocesses

        if settings is None:
            self.sampler_settings = self.Settings()
        else:
            self.sampler_settings = settings

        if self.sampler_settings.return_probes is True:
            self._calc_hotspots(return_probes=True)

        else:
            self._calc_hotspots()

        self.super_grids = {p: g[0] for p, g in self.out_grids.items()}

        if clear_tmp == True:
            shutil.rmtree(tmp)

        return Results(super_grids=self.super_grids,
                       protein=self.protein,
                       buriedness=self.buriedness,
                       superstar={x.identifier: x.grid for x in self.superstar_grids},
                       weighted_superstar={x.identifier: x.grid for x in self.weighted_grids})

    def global_dock(self, hr:Results, mol, settings = None):
        self.super_grids = hr.super_grids
        self.buriedness = hr.buriedness
        self.protein = hr.protein
        self.charged_probes = False
        self.probe_size = 7
        self.buriedness_method = 'ligsite'
        self.cavities = None
        self.clear_tmp = True

        if settings is None:
            self.sampler_settings = self.Settings()

        self.sampler_settings.return_probes = True
        self.sampler_settings.nrotations = 2000
        self.apolar_translation_threshold = 17
        self.polar_translation_threshold = 17

        print(hr.weighted_superstar)

        grid_dict = hr.weighted_superstar


        ps = self._dock_probe(mol, grid_dict)
        return ps

