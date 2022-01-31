# Description

Bindings code originally written by Lukas Mosimann (lukasm91).
This whole work is very *experimental* and must be understood as a proof of concept only!

- `src/atlas4py/_atlas4py.cpp` is the file providing the atlas binding
- `Atlas4py-Test.ipynb` is a Python notebook showing how to use the atlas bindings in the current experimental version.

If the atlas library is already installed (in a non-standard location), the path to atlas can be specified by
```
export ATLAS_INSTALL_DIR=<path-to-atlas-installation>
```
If atlas is not found, it will be installed by atlas4py.

# Future work

- We should think about which features we are going to support in this package

 
