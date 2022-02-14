import os
import re
import sys
import subprocess

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

HERE = os.path.abspath(os.path.dirname(__file__))

VERSIONS = dict(
    cmake="3.14.0",
    ecbuild="3.6.3",
    eckit="1.17.1",
    atlas="0.26.0",
    pybind11="2.5.0",
    python="3.6",
)

# Package meta-data.
NAME = "atlas4py"
DESCRIPTION = "Python bindings for Atlas: a ECMWF library for parallel data-structures"
URL = "https://github.com/ecmwf/atlas.git"
EMAIL = "willem.deconinck@ecmwf.int"
AUTHOR = "Willem Deconinck"
LICENSE = "Apache License 2.0"
CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Apache Software License",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: Implementation :: CPython",
    "Topic :: Scientific/Engineering :: Atmospheric Science",
]


PYTHON_REQUIRES = f">={VERSIONS['python']}"

CMAKE_MIN_VERSION = f"{VERSIONS['cmake']}"
BUILD_JOBS = os.cpu_count()


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        cmake_version = LooseVersion(re.search(r"version\s*([\d.]+)", out.decode()).group(1))
        if cmake_version < CMAKE_MIN_VERSION:
            raise RuntimeError(f"CMake >= {CMAKE_MIN_VERSION} is required on Windows")

        # Activate ccache
        # Taken from: https://github.com/h5py/h5py/pull/1382
        # This allows ccache to recognise the files when pip builds in a temp
        # directory. It speeds up repeatedly running tests through tox with
        # ccache configured (CC="ccache gcc"). It should have no effect if
        # ccache is not in use.
        os.environ["CCACHE_BASEDIR"] = HERE
        os.environ["CCACHE_NOHASHDIR"] = "1"

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        os.makedirs(self.build_temp, exist_ok=True)
        cfg = "Debug" if self.debug else "Release"

        # Run CMake configure
        print("-" * 10, "Running CMake prepare", "-" * 40)
        RPATH = "'@loader_path'" if sys.platform == "darwin" else "'$ORIGIN'"
        cmake_args = [
            "-DCMAKE_BUILD_WITH_INSTALL_RPATH=ON",
            "-DCMAKE_INSTALL_RPATH=" + RPATH,
            "-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON",
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
            "-DCMAKE_BUILD_TYPE=" + cfg,
            "-DATLAS4PY_CMAKE_MINIMUM_REQUIRED_VERSION=" + VERSIONS["cmake"],
            "-DATLAS4PY_ECBUILD_VERSION=" + VERSIONS["ecbuild"],
            "-DATLAS4PY_ECKIT_VERSION=" + VERSIONS["eckit"],
            "-DATLAS4PY_ATLAS_VERSION=" + VERSIONS["atlas"],
            "-DATLAS4PY_PYBIND11_VERSION=" + VERSIONS["pybind11"],
        ]
        # print(f"./{self.build_temp}$ " + " ".join(["cmake", ext.sourcedir] + cmake_args))
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp)

        # Run CMake build
        print("-" * 10, "Building extensions", "-" * 40)
        build_args = ["--config", cfg, "-j", str(BUILD_JOBS)]
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)


setup(
    name=NAME,
    version=VERSIONS["atlas"],
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    long_description="",
    license=LICENSE,
    classifiers=CLASSIFIERS,
    python_requires=PYTHON_REQUIRES,
    packages=find_packages("src"),
    package_dir={"": "src"},
    ext_modules=[CMakeExtension("atlas4py._atlas4py", sourcedir="src/atlas4py")],
    cmdclass={"build_ext": CMakeBuild},
    include_package_data=True,
    zip_safe=False,
)
