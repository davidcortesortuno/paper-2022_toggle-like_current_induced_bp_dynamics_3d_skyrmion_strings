Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4384569.svg)](https://doi.org/10.5281/zenodo.4384569)

# Data Set

This repository contains both the experimental data and the simulation scripts
to reproduce the results of [1]: *Toggle-like current-induced Bloch point
dynamics of 3D skyrmion strings in a room-temperature nanowire* by M. T. Birch,
D. Cortés-Ortuño, K. Litzius, S. Wintz, F. Schulz, M. Weigand, A. Štefančič, D.
Mayoh, G. Balakrishnan, P.D. Hatton, G. Schütz.

![](images/sk_transition.jpg)

# Experimental data

Please refer to the README in the [Experimental Data](Experimental_Data)
directory.

Cloning this repository will not load the experimental data files
automatically. You need to install `git lfs` in the cloned repository and then
`git pull`.

# Simulations

Simulations are based on the MuMax3 [2] code and can be found in the
`sims/mumax3` directory. Simulations are written directly in Go (`go` files) or
in MuMax's scripting language (`mx3` files). Simulations written in Go require
to be compiled with a working installation of MuMax3 by running, for example,

```shell
go build skyrmion_lattice_Co8Zn9Mn3/sk_lattice_D_0d6e-3_A_5d7295e-12_LD_120nm_ROT.go
```

There are different simulations in the [sims/mumax3/](sims/mumax3) directory

- **skyrmion_lattice_Co8Zn9Mn3** : Simulation of a skyrmion lattice. In the
  manuscript it is used the rotated lattice given by the script ending in `ROT`

- **skyrmion_memory_device** : Sequence of skyrmion strings of different length

- **two_tubes_center_Co8Zn9Mn3** : 2 skyrmion strings at the centre of the
  sample using different magnitudes of the exchange and DMI constant such that
  the helical length of the system is kept at 120 nm. This is to find the
  optimal values for the parameters to make the skyrmion string size comparable
  to experimental data.

- **two_tubes_field-sweep_Co8Zn9Mn3_anisotropy** : Starting from two tubes at
  the centre of the sample, separated by a small distance of 5% of the length
  of the film, at a field of By=-100 mT, these simulations do a field sweep
  decreasing the field magnitude to 0 mT, for different values of a uniaxial
  anisotropy with hard axis in the x-direction. This is done in order to check
  a threshold value where helices propagating in the x-direction are formed.

- **two_tubes_separated_field-sweep_Co8Zn9Mn3_anisotropy** : Similarly than
  before but using a specific magnitude of anisotropy and using different
  separation of the initial two skyrmion strings. The anisotropy value was
  chosen from the threshold value of the previous simulations.

- **two_tubes_field-sweep_Co8Zn9Mn3_NEW** : Field sweep starting from two
  skyrmion strings at the centre of the sample, at a field of By=-100mT, for
  two different initial separation of the skyrmions. The field sweeps are
  performed both increasing and decreasing the field magnitude.

- **two_tubes_separation_Co8Zn9Mn3** : Energy minimisation of two skyrmion
  strings solutions initially separated by different distances, at different
  applied field strengths. Field is applied along the film width, in the
  y-direction. This folder includes simulations of the conical state.


# Jupyter Notebooks

Notebooks where the energy of the configurations is analysed and where
visualisations of the configurations are shown, can be found in the
[data_analysis](data_analysis) directory. Data analysis of the MuMax3 [2]
simulation outputs is performed using the OOMMFPy library [5] Plots of the
simulations are performed using Matplotlib [3]. Data analysis is done via Numpy
arrays [4]. 

In the notebooks the main method is to compute the mean value of the out of
plane component of the magnetization across the film thickness, which is in the
z-direction. This is done to compare the results with the experimental data.
The average is computed by loading the magnetization data from `ovf` files
using OOMMFPy and then using Numpy array operations.

- **two_tubes_DMI_Exchange_vis.ipynb** : Data analysis of the simulation of two
  skyrmion strings using different values of exchange and DMI. Additionally, it
  is included visualisation of the skyrmion strings using PyVista. This
  notebook is well documented and can be used as a reference for the methods
  used in the other notebooks.

- **conical_state.ipynb** : Analysis of the conical state for a field applied
  in the direction of the film width. Results with and without anisotropy.

- **sk_lattice_vs_field.ipynb** : Simulation of skyrmion lattice in the
  xy-plane.

- **skyrmion_device.ipynb** : Visualisation of skyrmion strings of different
  length.

- **two_tubes_separation_vs_field.ipynb** : Compute images of the field sweep
  applied to the film with two skyrmion strings, starting from skyrmions
  separated by increasing distances.

- **two_tubes_separation_vs_field-CONICAL.ipynb** : Like before but looking
  at simulations where a conical background was specified in the initial state.

- **two_tubes_vs_anisotropies.ipynb** : Field sweep on a film with two skyrmion
  strings starting at a field of By=-100 mT, using different values of uniaxial
  anisotropy with hard axis in the x-direction.


# Cite

If you find this material useful please cite us (you might need the LaTeX's
`url` package)

    @Misc{Birch2021,
      author       = {M. T. Birch and D. Cort\'es-Ortu\~no},
      title        = {{Data set for: Toggle-like current-induced Bloch point dynamics of 3D skyrmion strings in a room-temperature nanowire}},
      howpublished = {Zenodo \url{doi:10.5281/zenodo.XXXX}. Github: \url{https://github.com/davidcortesortuno/paper-2022_toggle-like_current_induced_bp_dynamics_3d_skyrmion_strings}},
      year         = {2022},
      doi          = {10.5281/zenodo.XXX},
      url          = {https://doi.org/10.5281/zenodo.XXXX},
    }

# References

[1]  Birch M, Cortés-Ortuño D, Litzius K, et al. Toggle-like current-induced
Bloch point dynamics of 3D skyrmion strings in a room-temperature nanowire.
Research Square; 2022. DOI: 10.21203/rs.3.rs-1235546/v1.
[https://europepmc.org/article/ppr/ppr440740]

[2] Vansteenkiste, A. et al. The design and verification of MuMax3. AIP
Advances 4, 107133 (2014).

[3] J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science
& Engineering, vol. 9, no. 3, pp. 90-95, 2007.

[4] Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming
with NumPy. Nature 585, 357–362 (2020). DOI: 0.1038/s41586-020-2649-2.

[5] Cortés-Ortuño, D. OOMMFPy Python module. doi: 10.5281/zenodo.2611194 (2019).
[https://github.com/davidcortesortuno/oommfpy]
