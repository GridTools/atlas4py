# Description

Bindings code originally written by Lukas Mosimann (@lukasm91).
This whole work is very *experimental* and must be understood as a proof of concept only!

- `src/atlas4py/_atlas4py.cpp` is the file providing the atlas binding
- `Atlas4py-Test.ipynb` is a Python notebook showing how to use the atlas bindings in the current experimental version.

If the atlas library is already installed in a non-standard location, the path to atlas can be specified by:
```
export atlas_ROOT=<path-to-atlas-installation>
```
or other environment variables recognized by CMake [find_package(atlas)](https://cmake.org/cmake/help/latest/command/find_package.html)

If atlas is not found it will be automatically downloaded and compiled by atlas4py.

Experimental pre-built packages can be found in [TestPyPI](https://test.pypi.org/project/atlas4py/) and installed with `pip`:
```
pip install -i https://test.pypi.org/simple/ atlas4py
```

# Future work

- We should think about which features we are going to support in this package

