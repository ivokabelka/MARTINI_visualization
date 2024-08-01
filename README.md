# MARTINI v2 visualization
**martiniviz.py**

This script is designed for visualizing systems using the MARTINI v2 force field. It reads a GROMACS *.top topology file to generate the bonds between beads and outputs a PSF file compatible with VMD for visualization. The script is written in python3 and has no dependencies.

Running the script:
```
python3 martiniviz.py -i topol.top -o system.psf
```

Using the output:
```
vmd system.gro traj.xtc system.psf
```
