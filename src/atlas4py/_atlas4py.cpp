#include <functional>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h> // Required for std::string support
#include <nanobind/stl/array.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/optional.h> // required for std::optional arguments
#include <nanobind/stl/tuple.h> // required for std::optional arguments

#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildNode2CellConnectivity.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/meshgenerator.h"
#include "atlas/option/Options.h"
#include "atlas/output/Gmsh.h"
#include "atlas/library.h"

#if __has_include("pluto/pointer_info.h")
#include "pluto/pointer_info.h"
#define ATLAS4PY_HAS_PLUTO
#else
#warning "pluto/pointer_info.h not found, it is only available starting from atlas version 0.42.0. Assuming all arrays are host accessible and not pinned, device, or managed"
namespace pluto {
    inline bool is_pinned(const void*) { return false; }
    inline bool is_host(const void*) { return true; }
    inline bool is_device(const void*) { return false; }
    inline bool is_managed(const void*) { return false; }
    inline bool is_host_accessible(const void*) { return true; }
    inline bool is_device_accessible(const void*) { return false; }
}
#endif

#include "eckit/value/Value.h"
#include "eckit/config/Configuration.h"

namespace nb = ::nanobind;

using namespace atlas;
using namespace nb::literals;

namespace nanobind {
namespace detail {
template <>
struct type_caster<atlas::array::ArrayStrides>
    : public type_caster<std::vector<atlas::array::ArrayStrides::value_type>> {};
template <>
struct type_caster<atlas::array::ArrayShape> : public type_caster<std::vector<atlas::array::ArrayShape::value_type>> {};

}  // namespace detail
}  // namespace nanobind


namespace atlas4py::dtype {
std::string to_python_name( DataType const& dt ) {
    switch ( dt.kind() ) {
        case DataType::KIND_INT32:
            return "int32";
        case DataType::KIND_INT64:
            return "int64";
        case DataType::KIND_REAL32:
            return "float32";
        case DataType::KIND_REAL64:
            return "float64";
        case DataType::KIND_UINT64:
            return "uint64";
        default:
            return "";
    }
}

::nb::dlpack::dtype to_nb_dtype( DataType const& dt ) {
    switch ( dt.kind() ) {
        case DataType::KIND_INT32:
            return nb::dtype<int32_t>();
        case DataType::KIND_INT64:
            return nb::dtype<int64_t>();
        case DataType::KIND_REAL32:
            return nb::dtype<float>();
        case DataType::KIND_REAL64:
            return nb::dtype<double>();
        case DataType::KIND_UINT64:
            return nb::dtype<uint64_t>();
        default:
            return nb::dlpack::dtype();
    }
}

DataType from_python_name( std::string const& dtype_name ) {
    if (dtype_name == "int32") {
        return DataType::KIND_INT32;
    }
    else if (dtype_name == "int64") {
        return DataType::KIND_INT64;
    }
    else if (dtype_name == "float32" || dtype_name == "real32") {
        return DataType::KIND_REAL32;
    }
    else if (dtype_name == "float64" || dtype_name == "real64") {
        return DataType::KIND_REAL64;
    }
    else if (dtype_name == "uint64") {
        return DataType::KIND_UINT64;
    }
    else {
        throw std::out_of_range( "unsupported dtype: " + dtype_name );
    }
}

DataType from_python_object( nb::object const& dtype ) {
    if (nb::hasattr(dtype, "__name__")) {
        return from_python_name( nb::cast<std::string>(dtype.attr("__name__")) );
    }
    else if (nb::hasattr(dtype, "name")) {
        return from_python_name( nb::cast<std::string>(dtype.attr("name")) );
    }
    else if ( nb::isinstance<nb::str>( dtype ) ) {
        return from_python_name( nb::cast<std::string>(dtype) );
    }
    else {
        throw std::out_of_range( "unsupported dtype" );
    }
}
}


namespace {

void atlasInitialise() {
    static bool already_initialised = false;
    if (already_initialised)
        return;
    already_initialised = true;

    nb::module_ sys = nb::module_::import_("sys");
    nb::list sys_argv = sys.attr("argv");
    int argc = sys_argv.size();
    char** argv = new char*[argc];
    for (int i = 0; i < argc; ++i) {
        argv[i] = (char*) PyUnicode_AsUTF8(sys_argv[i].ptr());
    }

    atlas::initialise(argc, argv);
    // Warning: Don't delete argv, because this is being referenced within
    // the atlas library and deleting it may cause issues.
}

nb::object toPyObject( eckit::Configuration const& v );
nb::object toPyObject( eckit::Configuration const& v, std::string const& key );

nb::object toPyObject(bool v) {
    return nb::bool_(v);
}
nb::object toPyObject(long v) {
    return nb::int_(v);
}
nb::object toPyObject(double v) {
    return nb::float_(v);
}
nb::object toPyObject(std::string const& v) {
    return nb::str(v.c_str());
}
template <typename T>
nb::object toPyObject( std::vector<T> const& v ) {
    nb::list ret;
    for ( auto const& val : v ) {
        ret.append( toPyObject( val ) );
    }
    return ret;
}

nb::object toPyObject( eckit::Configuration const& v, std::string const& key ) {
    if ( v.isSubConfiguration ( key ) ) {
        return toPyObject( v.getSubConfiguration( key ) );
    }
    else if (v.isBoolean( key )) {
        return toPyObject( v.getBool( key ) );
    }
    else if (v.isIntegral( key )) {
        return toPyObject( v.getLong( key ) );
    }
    else if (v.isFloatingPoint( key )) {
        return toPyObject( v.getDouble( key ) );
    }
    else if (v.isString( key )) {
        return toPyObject( v.getString( key ) );
    }
    else if (v.isSubConfigurationList( key )) {
        std::vector<eckit::LocalConfiguration> subconfigs = v.getSubConfigurations( key );
        return toPyObject( subconfigs );
    }
    else if (v.isIntegralList( key )) {
        std::vector<long> values = v.getLongVector( key );
        return toPyObject( values );
    }
    else if (v.isFloatingPointList( key )) {
        std::vector<double> values = v.getDoubleVector( key );
        return toPyObject( values );
    }
    else if (v.isStringList( key )) {
        std::vector<std::string> values = v.getStringVector( key );
        return toPyObject( values );
    }
    else if (v.isBooleanList( key )) {
        throw std::out_of_range( "boolean lists not supported for key " + key );
    }
    else {
        throw std::out_of_range( "type of value unsupported for key " + key );
    }
}

nb::object toPyObject( eckit::Configuration const& v ) {
    nb::dict ret;
    for ( auto const& key : v.keys()) {
        ret[ key.c_str() ] = toPyObject( v, key );
    }
    return ret;
}

int get_nb_device_type( const void* ptr ) {
    if (pluto::is_pinned(ptr)) {
        return nb::device::cuda_host::value;
    }
    else if (pluto::is_host(ptr)) {
        return nb::device::cpu::value;
    }
    else if (pluto::is_device(ptr)) {
        return nb::device::cuda::value;
    }
    else if (pluto::is_managed(ptr)) {
        return nb::device::cuda_managed::value;
    }
    else {
        return nb::device::none::value;
    }
}

enum class MemorySpace {
    host,
    device
};
template<typename Framework>
auto make_ndarray(atlas::array::Array& array, MemorySpace memory_space) {
    void* data_ptr;
    if (memory_space == MemorySpace::host) {
        data_ptr = array.host_data<void>();
        if (not pluto::is_host_accessible(data_ptr)) {
            throw std::runtime_error( "array data is not host accessible" );
        }
    }
    else {
        data_ptr = array.device_data<void>();
        if (not pluto::is_device_accessible(data_ptr)) {
            throw std::runtime_error( "array data is not device accessible" );
        }
    }
    constexpr int max_ndim = 8;
    const auto ndim = array.rank();
    if (ndim > max_ndim) {
        throw std::runtime_error( "array rank exceeds maximum supported by atlas4py" );
    }
    std::array<size_t,max_ndim> shape;
    std::array<int64_t,max_ndim> strides;
    for( int i=0; i<ndim; ++i ) {
        shape[i]   = array.shape()[i];
        strides[i] = array.strides()[i];
    }
    return nb::ndarray<Framework>(
        data_ptr, // pointer to data
        ndim,  // ndim
        shape.data(), // shape
        nb::handle{}, // owner
        strides.data(), // strides
        atlas4py::dtype::to_nb_dtype( array.datatype()), // dtype
        get_nb_device_type(data_ptr), // device_type
        0, // device_id
        'C' // order
    );
}

}  // namespace

#define STRINGIFY(s) STRINGIFY_HELPER(s)
#define STRINGIFY_HELPER(s) #s

NB_MODULE( _atlas4py, m ) {
    m.def("_initialise", atlasInitialise)
     .def("_finalise",   atlas::finalise);
    m.attr("__version__") = STRINGIFY(ATLAS4PY_VERSION_STRING);

    nb::class_<PointLonLat>( m, "PointLonLat" )
        .def( nb::init<double,double>(), "lon"_a, "lat"_a )
        .def_prop_ro( "lon", nb::overload_cast<>( &PointLonLat::lon, nb::const_ ) )
        .def_prop_ro( "lat", nb::overload_cast<>( &PointLonLat::lat, nb::const_ ) )
        .def( "__repr__", []( PointLonLat const& p ) {
            return "_atlas4py.PointLonLat(lon=" + std::to_string( p.lon() ) + ", lat=" + std::to_string( p.lat() ) + ")";
        } );
    nb::class_<PointXY>( m, "PointXY" )
        .def( nb::init<double,double>(), "x"_a, "y"_a )
        .def_prop_ro( "x", nb::overload_cast<>( &PointXY::x, nb::const_ ) )
        .def_prop_ro( "y", nb::overload_cast<>( &PointXY::y, nb::const_ ) )
        .def( "__repr__", []( PointXY const& p ) {
            return "_atlas4py.PointXY(x=" + std::to_string( p.x() ) + ", y=" + std::to_string( p.y() ) + ")";
        } );

    nb::class_<Projection>( m, "Projection" )
        .def( "__repr__", []( Projection const& p ) {
        return "_atlas4py.Projection("_s + nb::str( toPyObject( p.spec() ) ) + ")"_s;
        } );

    nb::class_<Domain>( m, "Domain" )
        .def_prop_ro( "type", &Domain::type )
        .def_prop_ro( "is_global", &Domain::global )
        .def_prop_ro( "units", &Domain::units )
        .def( "__repr__", []( Domain const& d ) {
            if (d) {
                return nb::str("_atlas4py.Domain("_s + nb::str( toPyObject( d.spec() ) ) + ")"_s);
            }
            return nb::str("_atlas4py.Domain()"_s);
        } );
    nb::class_<RectangularDomain, Domain>( m, "RectangularDomain" )
        .def( nb::init<const atlas::RectangularDomain::Interval&, const atlas::RectangularDomain::Interval&, const std::string&>(), "x_interval"_a, "y_interval"_a, "units"_a = "degrees" );

    nb::class_<Grid>( m, "Grid" )
        .def( nb::init<const std::string&>(), "name"_a  )
        .def_prop_ro( "name", &Grid::name )
        .def_prop_ro( "uid", &Grid::uid )
        .def_prop_ro( "size", &Grid::size )
        .def_prop_ro( "projection", &Grid::projection )
        .def_prop_ro( "domain", &Grid::domain )
        .def( "__repr__",
              []( Grid const& g ) { return "_atlas4py.Grid("_s + nb::str( toPyObject( g.spec() ) ) + ")"_s; } );

    nb::class_<grid::Spacing>( m, "Spacing" )
        .def( "__len__", &grid::Spacing::size )
        .def( "__getitem__", &grid::Spacing::operator[])
        .def( "__repr__", []( grid::Spacing const& spacing ) {
            return "_atlas4py.Spacing("_s + nb::str( toPyObject( spacing.spec() ) ) + ")"_s;
        } );
    nb::class_<grid::LinearSpacing, grid::Spacing>( m, "LinearSpacing" )
        .def( nb::init<double,double,long,bool>(), "start"_a, "stop"_a, "N"_a, "endpoint_included"_a = true );


    nb::class_<grid::GaussianSpacing, grid::Spacing>( m, "GaussianSpacing" )
        .def( nb::init<long>(), "N"_a );

    nb::class_<StructuredGrid, Grid>( m, "StructuredGrid" )
        .def( nb::init<const Grid&, Domain>(), "grid"_a, "domain"_a = Domain() )
        .def( nb::init<const std::string&, Domain>(), "name"_a, "domain"_a = Domain() )
        .def( nb::init<const grid::LinearSpacing&, const grid::Spacing&>(), "x_spacing"_a, "y_spacing"_a )
        .def("__init__", [](StructuredGrid *g, std::vector<grid::LinearSpacing> xLinearSpacings, grid::Spacing ySpacing, Domain const& d) {
            std::vector<grid::Spacing> xSpacings;
            std::copy( xLinearSpacings.begin(), xLinearSpacings.end(), std::back_inserter( xSpacings ) );
            new (g) StructuredGrid( xSpacings, ySpacing, Projection(), d );
            }, "x_spacings"_a, "y_spacing"_a, "domain"_a = Domain() )
        .def( "__enter__", []( StructuredGrid& self ) { return &self; } )
        .def( "__exit__", []( StructuredGrid& self, nb::object exc_type, nb::object exc_value, nb::object traceback ) { self.reset(nullptr); }, "exc_type"_a = nb::none(), "exc_value"_a = nb::none(), "traceback"_a = nb::none() )
        .def("__bool__", &StructuredGrid::valid )
        .def_prop_ro( "valid", &StructuredGrid::valid )
        .def_prop_ro( "ny", &StructuredGrid::ny )
        .def_prop_ro( "nx", nb::overload_cast<>( &StructuredGrid::nx, nb::const_ ) )
        .def_prop_ro( "nxmax", &StructuredGrid::nxmax )
        .def_prop_ro( "y", nb::overload_cast<>( &StructuredGrid::y, nb::const_ ) )
        .def_prop_ro( "x", &StructuredGrid::x )
        .def( "xy", nb::overload_cast<idx_t, idx_t>( &StructuredGrid::xy, nb::const_ ), "i"_a, "j"_a )
        .def( "lonlat", nb::overload_cast<idx_t, idx_t>( &StructuredGrid::lonlat, nb::const_ ), "i"_a, "j"_a )
        .def_prop_ro( "reduced", &StructuredGrid::reduced )
        .def_prop_ro( "regular", &StructuredGrid::regular )
        .def_prop_ro( "periodic", &StructuredGrid::periodic );

    nb::class_<eckit::Configuration>( m, "eckit.Configuration" );

    nb::class_<eckit::LocalConfiguration, eckit::Configuration>( m, "eckit.LocalConfiguration" );

    // TODO This is a duplicate of metadata below (because same base class)
    nb::class_<util::Config, eckit::LocalConfiguration>( m, "Config" )
        .def( nb::init() )
        .def( "__setitem__",
              []( util::Config& config, std::string const& key, nb::object value ) {
                  if ( nb::isinstance<nb::bool_>( value ) )
                      config.set( key, nb::cast<bool>( value ) );
                  else if ( nb::isinstance<nb::int_>( value ) )
                      config.set( key, nb::cast<long long>( value ) );
                  else if ( nb::isinstance<nb::float_>( value ) )
                      config.set( key, nb::cast<double>( value ) );
                  else if ( nb::isinstance<nb::str>( value ) )
                      config.set( key, nb::cast<std::string>( value ) );
                  else
                      throw std::out_of_range( "type of value unsupported" );
              } )
        .def( "__getitem__",
              []( util::Config& config, std::string const& key ) -> nb::object {
                  if ( !config.has( key ) )
                      throw std::out_of_range( "key <" + key + "> could not be found" );

                  // TODO: We have to query metadata.get() even though this should
                  // not be done (see comment in Config::get). We cannot
                  // avoid this right now because otherwise we cannot query
                  // the type of the underlying data.
                  return toPyObject( config, key );
              } )
        .def( "__repr__", []( util::Config const& config ) {
            return "_atlas4py.Config("_s + nb::str( toPyObject( config ) ) + ")"_s;
        } );

    nb::class_<Field>( m, "Field" )
        .def_prop_ro( "name", &Field::name )
        .def_prop_ro( "strides", &Field::strides )
        .def_prop_ro( "shape", nb::overload_cast<>( &Field::shape, nb::const_ ) )
        .def_prop_ro( "size", &Field::size )
        .def_prop_ro( "rank", &Field::rank )
        .def_prop_ro( "dtype", []( Field& f ) { return atlas4py::dtype::to_python_name(f.datatype()); } )
        .def_prop_ro( "metadata", nb::overload_cast<>( &Field::metadata, nb::const_ ), nb::rv_policy::reference_internal )
        .def("host_array", [](Field& self) {
            return make_ndarray<nb::array_api>(self, MemorySpace::host);
            }, nb::rv_policy::reference_internal)
        .def("device_array", [](Field& self) {
            return make_ndarray<nb::array_api>(self, MemorySpace::device);
            }, nb::rv_policy::reference_internal)
        // Numpy array interface, see https://numpy.org/doc/stable/reference/arrays.interface.html
        .def("__array__", [](Field &self, nb::handle dtype, nb::handle copy) {
            return make_ndarray<nb::numpy>(self, MemorySpace::host);
            }, "dtype"_a = nb::none(), "copy"_a = nb::none(),
            nb::rv_policy::reference_internal)
        // CuPy array interface, see https://docs.cupy.dev/en/stable/reference/cupy.ndarray.html#cupy.ndarray.__array_interface__
        .def("__cuda_array_interface__", [](Field& self) {
            // WARNING: not tested
            return make_ndarray<nb::cupy>(self, MemorySpace::device);
            }, nb::rv_policy::reference_internal)
        // DLPack interface, see https://data-apis.org/array-api/latest/API_specification/generated/array_api.array.__dlpack__.html
        .def("__dlpack__", [](nb::pointer_and_handle<Field> self, nb::kwargs kwargs) {
            nb::object aa = nb::cast( make_ndarray<nb::array_api>(*self.p, MemorySpace::host), nb::rv_policy::reference_internal, self.h);
            return aa.attr("__dlpack__")(**kwargs);
            })
        .def("__dlpack_device__", [](nb::handle /*self*/) {
            return std::make_pair(nb::device::cpu::value, 0);
        });

    nb::class_<Mesh>( m, "Mesh" )
        .def_prop_ro( "grid", &Mesh::grid )
        .def_prop_ro( "projection", &Mesh::projection )
        .def_prop_ro( "nodes", nb::overload_cast<>( &Mesh::nodes, nb::const_ ))
        .def_prop_ro( "edges", nb::overload_cast<>( &Mesh::edges, nb::const_ ))
        .def_prop_ro( "cells", nb::overload_cast<>( &Mesh::cells, nb::const_ ));

    nb::class_<StructuredMeshGenerator>( m, "StructuredMeshGenerator" )
        // TODO in FunctionSpace below we expose config options, not the whole config object
        .def( nb::init<const util::Config&>(), "config"_a )
        .def( nb::init() )
        .def( "generate", nb::overload_cast<Grid const&>( &StructuredMeshGenerator::generate, nb::const_ ) );

    m.def( "build_edges", []( Mesh& mesh, const eckit::Configuration& config ) {
        mesh::actions::build_edges( mesh, config);
    }, "mesh"_a, "config"_a = util::Config() );
    m.def( "build_node_to_edge_connectivity",
           nb::overload_cast<Mesh&>( &mesh::actions::build_node_to_edge_connectivity ) );
    m.def( "build_element_to_edge_connectivity",
           nb::overload_cast<Mesh&>( &mesh::actions::build_element_to_edge_connectivity ) );
    m.def( "build_node_to_cell_connectivity",
           nb::overload_cast<Mesh&>( &mesh::actions::build_node_to_cell_connectivity ) );
    m.def( "build_median_dual_mesh", nb::overload_cast<Mesh&>( &mesh::actions::build_median_dual_mesh ) );
    m.def( "build_periodic_boundaries", nb::overload_cast<Mesh&>( &mesh::actions::build_periodic_boundaries ) );
    m.def( "build_halo", nb::overload_cast<Mesh&, int>( &mesh::actions::build_halo ) );
    m.def( "build_parallel_fields", nb::overload_cast<Mesh&>( &mesh::actions::build_parallel_fields ) );

    nb::class_<mesh::IrregularConnectivity>( m, "IrregularConnectivity" )
        .def( "__getitem__",
              []( mesh::IrregularConnectivity const& c, std::tuple<idx_t, idx_t> const& pos ) {
                  auto const& [row, col] = pos;
                  return c( row, col );
              } )
        .def_prop_ro( "rows", &mesh::IrregularConnectivity::rows )
        .def( "cols", &mesh::IrregularConnectivity::cols, "row_idx"_a )
        .def_prop_ro( "maxcols", &mesh::IrregularConnectivity::maxcols )
        .def_prop_ro( "mincols", &mesh::IrregularConnectivity::mincols );

    nb::class_<mesh::BlockConnectivity>( m, "BlockConnectivity" )
        .def( "__getitem__",
              []( mesh::BlockConnectivity const& c, std::tuple<idx_t, idx_t> const& pos ) {
                  auto const& [row, col] = pos;
                  return c( row, col );
              } )
        .def_prop_ro( "rows", &mesh::BlockConnectivity::rows )
        .def_prop_ro( "cols", &mesh::BlockConnectivity::cols );

    nb::class_<mesh::MultiBlockConnectivity>( m, "MultiBlockConnectivity" )
        .def( "__getitem__",
              []( mesh::MultiBlockConnectivity const& c, std::tuple<idx_t, idx_t> const& pos ) {
                  auto const& [row, col] = pos;
                  return c( row, col );
              } )
        .def( "__getitem__",
              []( mesh::MultiBlockConnectivity const& c, std::tuple<idx_t, idx_t, idx_t> const& pos ) {
                  auto const& [block, row, col] = pos;
                  return c( block, row, col );
              } )
        .def_prop_ro( "rows", &mesh::MultiBlockConnectivity::rows )
        .def( "cols", &mesh::MultiBlockConnectivity::cols, "row_idx"_a )
        .def_prop_ro( "maxcols", &mesh::MultiBlockConnectivity::maxcols )
        .def_prop_ro( "mincols", &mesh::MultiBlockConnectivity::mincols )
        .def_prop_ro( "blocks", &mesh::MultiBlockConnectivity::blocks )
        .def( "block", nb::overload_cast<idx_t>( &mesh::MultiBlockConnectivity::block, nb::const_ ), nb::rv_policy::reference_internal );

    nb::class_<mesh::Nodes>( m, "Nodes" )
        .def_prop_ro( "size", &mesh::Nodes::size )
        .def_prop_ro( "edge_connectivity", nb::overload_cast<>( &mesh::Nodes::edge_connectivity, nb::const_ ) )
        .def_prop_ro( "cell_connectivity", nb::overload_cast<>( &mesh::Nodes::cell_connectivity, nb::const_ ) )
        .def_prop_ro( "lonlat", nb::overload_cast<>( &Mesh::Nodes::lonlat, nb::const_ ) )
        .def("field", []( mesh::Nodes const& n, std::string const& name ) { return n.field( name ); }, "name"_a, nb::rv_policy::reference_internal )
        .def( "flags", []( mesh::Nodes const& n ) { return n.flags(); }, nb::rv_policy::reference_internal);

    nb::class_<mesh::HybridElements>( m, "HybridElements" )
        .def_prop_ro( "size", &mesh::HybridElements::size )
        .def( "nb_nodes", &mesh::HybridElements::nb_nodes )
        .def( "nb_edges", &mesh::HybridElements::nb_edges )
        .def_prop_ro( "node_connectivity", nb::overload_cast<>( &mesh::HybridElements::node_connectivity, nb::const_ ) )
        .def_prop_ro( "edge_connectivity", nb::overload_cast<>( &mesh::HybridElements::edge_connectivity, nb::const_ ) )
        .def_prop_ro( "cell_connectivity", nb::overload_cast<>( &mesh::HybridElements::cell_connectivity, nb::const_ ) )
        .def( "field", []( mesh::HybridElements const& he, std::string const& name ) { return he.field( name ); }, "name"_a, nb::rv_policy::reference_internal )
        .def( "flags", []( mesh::HybridElements const& he ) { return he.flags(); }, nb::rv_policy::reference_internal );

    auto m_fs = m.def_submodule( "functionspace" );
    nb::class_<FunctionSpace>( m_fs, "FunctionSpace" )
        .def_prop_ro( "size", &FunctionSpace::size )
        .def_prop_ro( "type", &FunctionSpace::type )
        .def(
            "create_field",
            []( FunctionSpace const& fs
                , nb::object dtype
                , std::optional<std::string> const& name
                , std::optional<int> levels
                , std::optional<int> variables
            ) {
                util::Config config;
                if ( name )
                    config = config | option::name( *name );
                // TODO what does it mean in atlas if levels is not set?
                if ( levels )
                    config = config | option::levels( *levels );
                if ( variables )
                    config = config | option::variables( *variables );
                config = config | option::datatype( atlas4py::dtype::from_python_object( dtype ) );
                return fs.createField( config );
            }, "dtype"_a,  "name"_a = std::nullopt, "levels"_a = std::nullopt, "variables"_a = std::nullopt );

    nb::class_<functionspace::EdgeColumns, FunctionSpace>( m_fs, "EdgeColumns" )
        .def("__init__", [](functionspace::EdgeColumns *t, const Mesh&m, int halo) { new (t) functionspace::EdgeColumns(m, util::Config()( "halo", halo )); }, "mesh"_a, "halo"_a = 0 )
        .def_prop_ro( "nb_edges", &functionspace::EdgeColumns::nb_edges )
        .def_prop_ro( "mesh", &functionspace::EdgeColumns::mesh )
        .def_prop_ro( "edges", &functionspace::EdgeColumns::edges )
        .def_prop_ro( "valid", &functionspace::EdgeColumns::valid );

    nb::class_<functionspace::NodeColumns, FunctionSpace>( m_fs, "NodeColumns" )
        .def("__init__", [](functionspace::NodeColumns *t, const Mesh&m, int halo) { new (t) functionspace::NodeColumns(m, util::Config()( "halo", halo )); }, "mesh"_a, "halo"_a = 0 )
        .def_prop_ro( "nb_nodes", &functionspace::NodeColumns::nb_nodes )
        .def_prop_ro( "mesh", &functionspace::NodeColumns::mesh )
        .def_prop_ro( "nodes", &functionspace::NodeColumns::nodes )
        .def_prop_ro( "valid", &functionspace::NodeColumns::valid );
    nb::class_<functionspace::CellColumns, FunctionSpace>( m_fs, "CellColumns" )
        .def("__init__", [](functionspace::CellColumns *t, const Mesh&m, int halo) { new (t) functionspace::CellColumns(m, util::Config()( "halo", halo )); }, "mesh"_a, "halo"_a = 0 )
        .def_prop_ro( "nb_cells", &functionspace::CellColumns::nb_cells )
        .def_prop_ro( "mesh", &functionspace::CellColumns::mesh )
        .def_prop_ro( "cells", &functionspace::CellColumns::cells )
        .def_prop_ro( "valid", &functionspace::CellColumns::valid );

    nb::class_<util::Metadata>( m, "Metadata" )
        .def_prop_ro( "keys", &util::Metadata::keys )
        .def( "__setitem__",
              []( util::Metadata& metadata, std::string const& key, nb::object value ) {
                  if ( nb::isinstance<nb::bool_>( value ) )
                      metadata.set( key, nb::cast<bool>(value) );
                  else if ( nb::isinstance<nb::int_>( value ) )
                      metadata.set( key, nb::cast<long long>(value) );
                  else if ( nb::isinstance<nb::float_>( value ) )
                      metadata.set( key, nb::cast<double>(value) );
                  else if ( nb::isinstance<nb::str>( value ) )
                      metadata.set( key, nb::cast<std::string>(value) );
                  else
                      throw std::out_of_range( "type of value unsupported" );
              } )
        .def( "__getitem__",
              []( util::Metadata& metadata, std::string const& key ) -> nb::object {
                  if ( !metadata.has( key ) )
                      throw std::out_of_range( "key <" + key + "> could not be found" );

                  // TODO: We have to query metadata.get() even though this should
                  // not be done (see comment in Config::get). We cannot
                  // avoid this right now because otherwise we cannot query
                  // the type of the underlying data.
                  return toPyObject( metadata, key );
              } )
        .def( "__repr__", []( util::Metadata const& metadata ) {
            return "_atlas4py.Metadata("_s + nb::str( toPyObject( metadata ) ) + ")"_s;
        } );

    nb::class_<mesh::Nodes::Topology> topology( m, "Topology" );
    topology.attr( "NONE" )     = nb::cast( int( mesh::Nodes::Topology::NONE ) );
    topology.attr( "GHOST" )    = nb::cast( int( mesh::Nodes::Topology::GHOST ) );
    topology.attr( "PERIODIC" ) = nb::cast( int( mesh::Nodes::Topology::PERIODIC ) );
    topology.attr( "BC" )       = nb::cast( int( mesh::Nodes::Topology::BC ) );
    topology.attr( "WEST" )     = nb::cast( int( mesh::Nodes::Topology::WEST ) );
    topology.attr( "EAST" )     = nb::cast( int( mesh::Nodes::Topology::EAST ) );
    topology.attr( "NORTH" )    = nb::cast( int( mesh::Nodes::Topology::NORTH ) );
    topology.attr( "SOUTH" )    = nb::cast( int( mesh::Nodes::Topology::SOUTH ) );
    topology.attr( "PATCH" )    = nb::cast( int( mesh::Nodes::Topology::PATCH ) );
    topology.attr( "POLE" )     = nb::cast( int( mesh::Nodes::Topology::POLE ) );
    topology.def_static( "reset", &mesh::Nodes::Topology::reset );
    topology.def_static( "set", &mesh::Nodes::Topology::set );
    topology.def_static( "unset", &mesh::Nodes::Topology::unset );
    topology.def_static( "toggle", &mesh::Nodes::Topology::toggle );
    topology.def_static( "check", &mesh::Nodes::Topology::check );
    topology.def_static( "check_all", &mesh::Nodes::Topology::check_all );
    topology.def_static( "check_any", &mesh::Nodes::Topology::check_any );

    nb::class_<output::Gmsh>( m, "Gmsh" )
        .def( nb::init<std::string const&>(), "path"_a )
        .def( "__enter__", []( output::Gmsh& self ) { return &self; } )
        .def( "__exit__", []( output::Gmsh& self, nb::object exc_type, nb::object exc_value, nb::object traceback ) { self.reset(nullptr); }, "exc_type"_a = nb::none(), "exc_value"_a = nb::none(), "traceback"_a = nb::none() )
        .def( "write", []( output::Gmsh& gmsh, Mesh const& mesh ) { gmsh.write( mesh ); }, "mesh"_a )
        .def( "write", []( output::Gmsh& gmsh, Field const& field ) { gmsh.write( field ); }, "field"_a )
        .def( "write", []( output::Gmsh& gmsh, Field const& field, FunctionSpace const& fs ) { gmsh.write( field, fs ); }, "field"_a, "functionspace"_a )
        .def("__repr__", []( output::Gmsh const& gmsh ) { return "_atlas4py.output.Gmsh()"; } );
}
