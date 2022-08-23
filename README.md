# SPADExp
SPADExp is the abbreviation of "Simulator of photoemission angular distribution for experiment."
SPADExp calculates the photoemission angular distribution (PAD), which is the momentum dependence of spectrum intensity in angle-resolved photoemission spectroscopy (ARPES).
The software can directly load the output of the first-principles software package OpenMX, so users do not need to construct tight-binding models as previous studies did.
As a result, we can calculate the PADs of large systems such as quasicrystals and slab systems.


# Requirements
## C, C++
- make
- HDF5 (compiled with the ```--enable-cxx``` option)
- OpenMP
- BLAS

## Python3
- PyQt5
- pyqtgraph
- h5py
- numpy
- scipy

# [Docs](https://github.com/Hiroaki-Tanaka-0606/SPADExp/tree/main/Docs)
- [Docs (ja)](https://github.com/Hiroaki-Tanaka-0606/SPADExp/raw/main/Docs/SPADExp_docs_ja.pdf)

We will serve docs both in Japanese and English.

# [Examples](https://github.com/Hiroaki-Tanaka-0606/SPADExp/tree/main/Examples)
See also docs to learn what keywords in the input file mean.

# Contact
Hiroaki Tanaka

The Institute for Solid State Physics, The University of Tokyo

E-mail: hiroaki-tanaka_at_issp.u-tokyo.ac.jp (replace \_at\_ with @)

# Links
- OpenMX [Japanese](http://www.openmx-square.org/openmx_man3.9jp/index.html) [English](http://www.openmx-square.org/openmx_man3.9/index.html)
- ADPACK [Japanese](http://www.openmx-square.org/adpack_man2.2_jp/adpack2_2_jp.html) [English](http://www.openmx-square.org/adpack_man2.2/adpack2_2.html)
- [HDF5](https://www.hdfgroup.org/)
- Paper about SPADExp (in preparation)
