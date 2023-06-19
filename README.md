# CAM-Nor Physics

This repository holds source code which have been modified from the version in ESCOMP/CAM.

There are 3 directories in this repository:

- fv: These are files which modify the finite-volume dycore and shadow src/dynamics/fv
- control: This file shadows src/control/camsrfexch.F90 and modifies the surface exchange physics
- physics: This code modifies various aspects of the physics squence including the ZM scheme. It shadows src/physics/cam
