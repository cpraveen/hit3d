# hit3d #

This code is a pseudo-spectral DNS code for simulation of homogeneous isotropic incompressible turbulence in three dimensional space. It is aspiring to be a standard code for DNS of isotropic homogeneous turbulence in triple-periodic box.

The features include:
  * MPI framework
    * currently has been tested for OpenMPI and MVAPICH on variety of NSF and DOE clusters
  * FFTW3 is used for Fourier transforms. See http://www.fftw.org for details.
  * Lagrangian particles (tracers) to gather Lagrangian statistics.
  * Ability to perform Large-Eddy Simulation (LES), several models are implemented.

The code is released under the GNU Public License.

With any questions, contact the project owner or visit http://schumakov.info

**When using HIT3D for your research, please cite the following papers and provide the URL:**
  * SG Chumakov, A priori study of subgrid-scale flux of a passive scalar in turbulence, Phys. Rev. E, 78 15563 (2008)
  * SG Chumakov, Scaling properties of subgrid-scale energy dissipation, Phys. Fluids, 19 058104 (2007)
  * http://hit3d.googlecode.com
