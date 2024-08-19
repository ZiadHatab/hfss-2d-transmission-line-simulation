
# Python-controlled 2D HFSS Simulation of Transmission Lines

This is something I have wanted to do for a while... still a work in progress. The goal is to automatically simulate typical transmission lines (maybe even waveguides) and get the 2D parameters, i.e. characteristic impedance and propagation constant, from which RLGC parameters could then be calculated.

I also plan to include Jacobian calculation, since HFSS supports derivation during simulation.

Surface roughness and platting are on the TODO list.

BTW, if anyone has already done something similar, please get in touch!
## Installation

There is no installtion. Just run the code provided in repository. You need to install [`pyaedt`](https://github.com/ansys/pyaedt). Of course, you also need Ansys EDT installed, e.g., <https://www.ansys.com/academic/students/ansys-electronics-desktop-student>

[`pyaedt`](https://github.com/ansys/pyaedt) package is constantly being updated, and I wouldn't be surprised if some of the commands get deprecated. So, it is what it is 😒

```cmd
  python -m pip install pyaedt -U
```

## License

[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)

