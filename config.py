# -----------------------------------------------------------------------------
# KMCLib version 1.1-b1
# Distributed under the GPLv3 license
# Copyright (C)  2012-2014  Mikael Leetmaa
# Developed by Mikael Leetmaa <leetmaa@kth.se>
#
# This program is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE and README files, and the source code, for details.
#
# You should have received a copy of the GNU General Public License version 3
# (GPLv3) along with this program. If not, see <http://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------

from KMCLib import *

# -----------------------------------------------------------------------------
# Unit cell

cell_vectors = [[   2.870000e+00,   0.000000e+00,   0.000000e+00],
                [   0.000000e+00,   2.870000e+00,   0.000000e+00],
                [   0.000000e+00,   0.000000e+00,   2.870000e+00]]

basis_points = [[   0.000000e+00,   0.000000e+00,   0.000000e+00],
                [   5.000000e-01,   5.000000e-01,   5.000000e-01]]

unit_cell = KMCUnitCell(
    cell_vectors=cell_vectors,
    basis_points=basis_points)

# -----------------------------------------------------------------------------
# Lattice

lattice = KMCLattice(
    unit_cell=unit_cell,
    repetitions=(5,5,5),
    periodic=(True, True, True))

# -----------------------------------------------------------------------------
# Configuration

types = ['Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','V','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','V','Fe',
         'Fe','Fe','Fe','V','Fe','Fe','Fe','V','Fe','V','Fe',
         'Fe','Fe','V','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','V','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','V',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','V','Fe','Fe','Fe','V','Fe','Fe','V','Fe',
         'Fe','Fe','V','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','Fe','V','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','V','Fe','Fe','Fe','Fe','V','Fe','Fe','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','V','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','V','Fe','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','V','Fe',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe',
         'Fe','Fe','Fe','V','Fe','Fe','Fe','Fe','Fe','Fe','V',
         'Fe','Fe','Fe','Fe','Fe','Fe','Fe','Fe']

possible_types = ['Fe','V']

configuration = KMCConfiguration(
    lattice=lattice,
    types=types,
    possible_types=possible_types)
