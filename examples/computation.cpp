/*cppimport
<%
cfg['compiler_args'] = ['-std=c++17', '-fopenmp', '-O0']
cfg['linker_args'] = [
    '-L/home/lukas/documents/work/eckit/install/lib/',
    '-L/home/lukas/documents/work/atlas/install/lib/',
    '-Wl,-rpath,/home/lukas/documents/work/atlas/install/lib:/home/lukas/documents/work/eckit/install/lib',
    '-fopenmp'
    ]
cfg['include_dirs'] = [
    '/home/lukas/documents/work/atlas4py/build/_deps/pybind11_fetch-src/include',
    '/home/lukas/documents/work/atlas/install/include',
    '/home/lukas/documents/work/eckit/install/include',
    '/home/lukas/documents/work/dawn/dawn/prototype/',
    '/home/lukas/documents/work/dawn/dawn/src/'
    ]
cfg['libraries'] = ['eckit', 'atlas']

setup_pybind11(cfg)
%>
*/
#include <pybind11/pybind11.h>

#include "atlas_interface.hpp"
#include "generated_Diffusion.hpp"

PYBIND11_MODULE(computation, m) {
    m.def("run_computation", [](atlas::Mesh const& mesh, int ksize,
                                atlas::Field in_field, atlas::Field out_field) {
        auto in_view = atlasInterface::Field<double>(
            atlas::array::make_view<double, 2>(in_field));
        auto out_view = atlasInterface::Field<double>(
            atlas::array::make_view<double, 2>(out_field));
        dawn_generated::cxxnaiveico::generated<atlasInterface::atlasTag>(
            mesh, ksize, in_view, out_view)
            .run();
    });
}
