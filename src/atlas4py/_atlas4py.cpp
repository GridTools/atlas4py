#include <functional>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
#include "atlas/trans/Trans.h"
#include "atlas/interpolation.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/util/function/SphericalHarmonic.h"

#include "eckit/value/Value.h"
#include "eckit/config/Configuration.h"


namespace py = ::pybind11;
using namespace atlas;
using namespace pybind11::literals;

namespace pybind11 {
namespace detail {
template <>
struct type_caster<atlas::array::ArrayStrides>
    : public type_caster<std::vector<atlas::array::ArrayStrides::value_type>> {};
template <>
struct type_caster<atlas::array::ArrayShape> : public type_caster<std::vector<atlas::array::ArrayShape::value_type>> {};

}  // namespace detail
}  // namespace pybind11

namespace {

struct PySys {
    int argc;
    char** argv;
    static const PySys& instance() {
        static PySys _instance;
        return _instance;
    }
private:
    PySys() {
        py::module sys = py::module::import("sys");
        py::list sys_argv = sys.attr("argv");
        argc = (int)sys_argv.size();
        argv = (char**)malloc(argc * sizeof(char*));
        for (int i = 0; i < argc; ++i) {
            argv[i] = (char*)PyUnicode_AsUTF8(sys_argv[i].ptr());
        }
    }
};

void initialise_sys_argv() {
    atlas::initialise(PySys::instance().argc,PySys::instance().argv);
};

void config_set( util::Config& config, const std::string& key, py::handle value ) {
        if ( py::isinstance<py::bool_>( value ) ) {
            config.set(key,value.cast<bool>());
        }
        else if ( py::isinstance<py::int_>( value ) ) {
            config.set(key,value.cast<long long>());
        } 
        else if ( py::isinstance<py::float_>( value ) ) {
            config.set(key,value.cast<double>());
        }
        else if ( py::isinstance<py::str>( value ) ) {
            config.set(key, value.cast<std::string>());
        }
        else {
            throw std::out_of_range( "type of value unsupported" );
        }
}

util::Config to_config( py::kwargs kwargs ) {
    util::Config config;
    for( const auto& pair : kwargs ) {
        const auto key = pair.first.cast<std::string>();
        const auto& value = pair.second;
        config_set(config, key, value);
    }
    return config;
}

py::object toPyObject( eckit::Value const& v ) {
    if ( v.isBool() )
        return py::bool_( v.as<bool>() );
    else if ( v.isNumber() )
        return py::int_( v.as<long long>() );
    else if ( v.isDouble() )
        return py::float_( v.as<double>() );
    else if ( v.isMap() ) {
        py::dict ret;
        auto const& map = v.as<eckit::ValueMap>();
        for ( auto const& [k, v] : map ) {
            ret[k.as<std::string>().c_str()] = toPyObject( v );
        }
        return ret;
    }
    else if ( v.isList() ) {
        py::list ret;
        auto const& list = v.as<eckit::ValueList>();
        for ( auto const& v : list )
            ret.append( toPyObject( v ) );
        return ret;
    }
    else if ( v.isString() )
        return py::str( v.as<std::string>() );
    else
        throw std::out_of_range( "type of value unsupported (" + v.typeName() + ")" );
}
std::string atlasToPybind( array::DataType const& dt ) {
    switch ( dt.kind() ) {
        case array::DataType::KIND_INT32:
            return py::format_descriptor<int32_t>::format();
        case array::DataType::KIND_INT64:
            return py::format_descriptor<int64_t>::format();
        case array::DataType::KIND_REAL32:
            return py::format_descriptor<float>::format();
        case array::DataType::KIND_REAL64:
            return py::format_descriptor<double>::format();
        case array::DataType::KIND_UINT64:
            return py::format_descriptor<uint64_t>::format();
        default:
            return "";
    }
}
array::DataType pybindToAtlas( py::dtype const& dtype ) {
    if ( dtype.is( py::dtype::of<int32_t>() ) )
        return array::DataType::KIND_INT32;
    else if ( dtype.is( py::dtype::of<int64_t>() ) )
        return array::DataType::KIND_INT64;
    else if ( dtype.is( py::dtype::of<float>() ) )
        return array::DataType::KIND_REAL32;
    else if ( dtype.is( py::dtype::of<double>() ) )
        return array::DataType::KIND_REAL64;
    else if ( dtype.is( py::dtype::of<uint64_t>() ) )
        return array::DataType::KIND_UINT64;
    else
        return { 0 };
}


}  // namespace

PYBIND11_MODULE( _atlas4py, m ) {
    auto m_library = m.def_submodule( "library" );
    m_library.def("initialize", []() { atlas::initialise(PySys::instance().argc, PySys::instance().argv);})
             .def("initialise", []() { atlas::initialise(PySys::instance().argc, PySys::instance().argv);})
             .def("finalize",   []() { atlas::finalize(); })
             .def("finalise",   []() { atlas::finalize(); });
    m_library.attr("version") = atlas::Library::instance().version();

    m.def("initialize", []() { atlas::initialise(PySys::instance().argc, PySys::instance().argv);})
     .def("initialise", []() { atlas::initialise(PySys::instance().argc, PySys::instance().argv);})
     .def("finalize",   []() { atlas::finalize(); })
     .def("finalise",   []() { atlas::finalize(); });

    py::class_<PointLonLat>( m, "PointLonLat" )
        .def( py::init( []( double lon, double lat ) {
                  return PointLonLat( { lon, lat } );
              } ),
              "lon"_a, "lat"_a )
        .def_property_readonly( "lon", py::overload_cast<>( &PointLonLat::lon, py::const_ ) )
        .def_property_readonly( "lat", py::overload_cast<>( &PointLonLat::lat, py::const_ ) )
        .def( "__repr__", []( PointLonLat const& p ) {
            return "_atlas4py.PointLonLat(lon=" + std::to_string( p.lon() ) + ", lat=" + std::to_string( p.lat() ) +
                   ")";
        } );
    py::class_<PointXY>( m, "PointXY" )
        .def( py::init( []( double x, double y ) {
                  return PointXY( { x, y } );
              } ),
              "x"_a, "y"_a )
        .def_property_readonly( "x", py::overload_cast<>( &PointXY::x, py::const_ ) )
        .def_property_readonly( "y", py::overload_cast<>( &PointXY::y, py::const_ ) )
        .def( "__repr__", []( PointXY const& p ) {
            return "_atlas4py.PointXY(x=" + std::to_string( p.x() ) + ", y=" + std::to_string( p.y() ) + ")";
        } );

    py::class_<Projection>( m, "Projection" ).def( "__repr__", []( Projection const& p ) {
        return "_atlas4py.Projection("_s + py::str( toPyObject( p.spec().get() ) ) + ")"_s;
    } );
    py::class_<Domain>( m, "Domain" )
        .def_property_readonly( "type", &Domain::type )
        .def_property_readonly( "global", &Domain::global )
        .def_property_readonly( "units", &Domain::units )
        .def( "__repr__", []( Domain const& d ) {
            return "_atlas4py.Domain("_s + ( d ? py::str( toPyObject( d.spec().get() ) ) : "" ) + ")"_s;
        } );
    py::class_<RectangularDomain, Domain>( m, "RectangularDomain" )
        .def( py::init( []( std::tuple<double, double> xInterval, std::tuple<double, double> yInterval ) {
                  auto [xFrom, xTo] = xInterval;
                  auto [yFrom, yTo] = xInterval;
                  return RectangularDomain( { xFrom, xTo }, { yFrom, yTo } );
              } ),
              "x_interval"_a, "y_interval"_a );

    py::class_<Grid>( m, "Grid" )
        .def( py::init<const std::string&>() )
        .def_property_readonly( "name", &Grid::name )
        .def_property_readonly( "uid", &Grid::uid )
        .def_property_readonly( "size", &Grid::size )
        .def_property_readonly( "projection", &Grid::projection )
        .def_property_readonly( "domain", &Grid::domain )
        .def( "lonlat", [](Grid const& self) {
                 return self.lonlat().begin();
             })
        .def( "__repr__",
              []( Grid const& g ) { return "_atlas4py.Grid("_s + py::str( toPyObject( g.spec().get() ) ) + ")"_s; } );

    py::class_<grid::IteratorLonLat>(m, "IteratorLonLat")
        .def("__iter__", [](grid::IteratorLonLat& it) { return it; })
        .def("__next__", [](grid::IteratorLonLat& it) {
            PointLonLat p;
            if( !it.next(p) ) {
                throw py::stop_iteration();
            }
            return p;});

    py::class_<grid::Spacing>( m, "Spacing" )
        .def( "__len__", &grid::Spacing::size )
        .def( "__getitem__", &grid::Spacing::operator[])
        .def( "__repr__", []( grid::Spacing const& spacing ) {
            return "_atlas4py.Spacing("_s + py::str( toPyObject( spacing.spec().get() ) ) + ")"_s;
        } );
    py::class_<grid::LinearSpacing, grid::Spacing>( m, "LinearSpacing" )
        .def( py::init( []( double start, double stop, long N, bool endpoint ) {
                  return grid::LinearSpacing{ start, stop, N, endpoint };
              } ),
              "start"_a, "stop"_a, "N"_a, "endpoint_included"_a = true );
    py::class_<grid::GaussianSpacing, grid::Spacing>( m, "GaussianSpacing" )
        .def( py::init( []( long N ) { return grid::GaussianSpacing{ N }; } ), "N"_a );

    py::class_<StructuredGrid, Grid>( m, "StructuredGrid" )
        .def( py::init( []( std::string const& s, Domain const& d ) {
                  return StructuredGrid{ s, d };
              } ),
              "gridname"_a, "domain"_a = Domain() )
        .def( py::init( []( grid::LinearSpacing xSpacing, grid::Spacing ySpacing ) {
                  return StructuredGrid{ xSpacing, ySpacing };
              } ),
              "x_spacing"_a, "y_spacing"_a )
        .def(
            py::init( []( std::vector<grid::LinearSpacing> xLinearSpacings, grid::Spacing ySpacing, Domain const& d ) {
                std::vector<grid::Spacing> xSpacings;
                std::copy( xLinearSpacings.begin(), xLinearSpacings.end(), std::back_inserter( xSpacings ) );
                return StructuredGrid{ xSpacings, ySpacing, Projection(), d };
            } ),
            "x_spacings"_a, "y_spacing"_a, "domain"_a = Domain() )
        .def_property_readonly( "valid", &StructuredGrid::valid )
        .def_property_readonly( "ny", &StructuredGrid::ny )
        .def_property_readonly( "nx", py::overload_cast<>( &StructuredGrid::nx, py::const_ ) )
        .def_property_readonly( "nxmax", &StructuredGrid::nxmax )
        .def_property_readonly( "y", py::overload_cast<>( &StructuredGrid::y, py::const_ ) )
        .def_property_readonly( "x", &StructuredGrid::x )
        .def( "xy", py::overload_cast<idx_t, idx_t>( &StructuredGrid::xy, py::const_ ), "i"_a, "j"_a )
        .def( "lonlat", py::overload_cast<idx_t, idx_t>( &StructuredGrid::lonlat, py::const_ ), "i"_a, "j"_a )
        .def_property_readonly( "reduced", &StructuredGrid::reduced )
        .def_property_readonly( "regular", &StructuredGrid::regular )
        .def_property_readonly( "periodic", &StructuredGrid::periodic );

    py::class_<GaussianGrid, StructuredGrid>( m, "GaussianGrid" )
        .def( py::init( []( Grid const& g ) {
                  return GaussianGrid{ g };
              } ),
              "grid"_a )
        .def_property_readonly( "valid", &GaussianGrid::valid )
        .def("__bool__", &GaussianGrid::valid);

    py::class_<eckit::Configuration>( m, "eckit.Configuration" );
    py::class_<eckit::LocalConfiguration, eckit::Configuration>( m, "eckit.LocalConfiguration" );

    // TODO This is a duplicate of metadata below (because same base class)
    py::class_<util::Config, eckit::LocalConfiguration>( m, "Config" )
        .def( py::init() )
        .def( py::init( []( py::kwargs kwargs) {
            return to_config(kwargs);
        } ) )
        .def( "__setitem__",
              []( util::Config& config, std::string const& key, py::object value ) {
                  config_set(config,key,value);
              } )
        .def( "__getitem__",
              []( util::Config& config, std::string const& key ) -> py::object {
                  if ( !config.has( key ) )
                      throw std::out_of_range( "key <" + key + "> could not be found" );
                  return toPyObject( config.get().element( key ) );
              } )
        .def( "__repr__", []( util::Config const& config ) {
            return "_atlas4py.Config("_s + py::str( toPyObject( config.get() ) ) + ")"_s;
        } );

    py::class_<grid::Partitioner>( m, "Partitioner" )
        .def( py::init( []( py::kwargs kwargs ) { return grid::Partitioner( to_config(kwargs)); } ) )
        .def( py::init( []( util::Config const& config ) { return grid::Partitioner( config ); } ) )
        .def( py::init( []( const std::string& type ) { return grid::Partitioner(type); } ) ) 
        .def( py::init( [](){ return grid::Partitioner(); } ) );

    py::class_<grid::MatchingPartitioner,grid::Partitioner>( m, "MatchingPartitioner" )
        .def( py::init( []( Mesh const& mesh, py::kwargs kwargs ) { return grid::MatchingPartitioner(mesh, to_config(kwargs)); } ) )
        .def( py::init( []( FunctionSpace const& functionspace, py::kwargs kwargs ) { return grid::MatchingPartitioner(functionspace, to_config(kwargs)); } ) );

    py::class_<MeshGenerator>( m, "MeshGenerator" )
        .def( py::init( []( py::kwargs kwargs ) { return MeshGenerator( to_config(kwargs)); } ) )
        .def( py::init( []( util::Config const& config ) { return MeshGenerator( config ); } ) )
        .def( py::init( []( const std::string& type ) { return MeshGenerator(type); } ) ) 
        .def( py::init( [](){ return MeshGenerator(); } ) )
        .def( "generate", [](MeshGenerator const& self, Grid const& grid){ return self.generate(grid); } )
        .def( "generate", [](MeshGenerator const& self, Grid const& grid, grid::Partitioner const& partitioner ){ return self.generate(grid, partitioner); } );

    py::class_<StructuredMeshGenerator>( m, "StructuredMeshGenerator" )
        // TODO in FunctionSpace below we expose config options, not the whole config object
        .def( py::init( []( util::Config const& config ) { return StructuredMeshGenerator( config ); } ) )
        .def( py::init() )
        .def( "generate", py::overload_cast<Grid const&>( &StructuredMeshGenerator::generate, py::const_ ) );

    py::class_<Mesh>( m, "Mesh" )
        .def( py::init( []( const Grid& grid ) { return Mesh(grid); } ) )
        .def( py::init( []( const Grid& grid, const grid::Partitioner& partitioner ) { return Mesh(grid,partitioner); } ) )
        .def_property_readonly( "grid", &Mesh::grid )
        .def_property_readonly( "projection", &Mesh::projection )
        .def_property( "nodes", py::overload_cast<>( &Mesh::nodes, py::const_ ), py::overload_cast<>( &Mesh::nodes ) )
        .def_property( "edges", py::overload_cast<>( &Mesh::edges, py::const_ ), py::overload_cast<>( &Mesh::edges ) )
        .def_property( "cells", py::overload_cast<>( &Mesh::cells, py::const_ ), py::overload_cast<>( &Mesh::cells ) );
    m.def( "build_edges", []( Mesh& mesh, std::optional<std::reference_wrapper<eckit::Configuration>> const& config ) {
        if(config)
            mesh::actions::build_edges( mesh, config.value().get());
        else
            mesh::actions::build_edges( mesh);
    }, "mesh"_a, "config"_a = std::nullopt );
    m.def( "build_node_to_edge_connectivity",
           py::overload_cast<Mesh&>( &mesh::actions::build_node_to_edge_connectivity ) );
    m.def( "build_element_to_edge_connectivity",
           py::overload_cast<Mesh&>( &mesh::actions::build_element_to_edge_connectivity ) );
    m.def( "build_node_to_cell_connectivity",
           py::overload_cast<Mesh&>( &mesh::actions:: build_node_to_cell_connectivity ) );
    m.def( "build_median_dual_mesh", py::overload_cast<Mesh&>( &mesh::actions::build_median_dual_mesh ) );
    m.def( "build_periodic_boundaries", py::overload_cast<Mesh&>( &mesh::actions::build_periodic_boundaries ) );
    m.def( "build_halo", py::overload_cast<Mesh&, int>( &mesh::actions::build_halo ) );
    m.def( "build_parallel_fields", py::overload_cast<Mesh&>( &mesh::actions::build_parallel_fields ) );

    py::class_<mesh::IrregularConnectivity>( m, "IrregularConnectivity" )
        .def( "__getitem__",
              []( mesh::IrregularConnectivity const& c, std::tuple<idx_t, idx_t> const& pos ) {
                  auto const& [row, col] = pos;
                  return c( row, col );
              } )
        .def_property_readonly( "rows", &mesh::IrregularConnectivity::rows )
        .def( "cols", &mesh::IrregularConnectivity::cols, "row_idx"_a )
        .def_property_readonly( "maxcols", &mesh::IrregularConnectivity::maxcols )
        .def_property_readonly( "mincols", &mesh::IrregularConnectivity::mincols );
    py::class_<mesh::MultiBlockConnectivity>( m, "MultiBlockConnectivity" )
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
        .def_property_readonly( "blocks", &mesh::MultiBlockConnectivity::blocks )
        .def( "block", py::overload_cast<idx_t>( &mesh::MultiBlockConnectivity::block, py::const_ ), py::return_value_policy::reference_internal)
        .def_property_readonly( "rows", &mesh::MultiBlockConnectivity::rows )
        .def( "cols", &mesh::MultiBlockConnectivity::cols, "row_idx"_a )
        .def_property_readonly( "maxcols", &mesh::MultiBlockConnectivity::maxcols )
        .def_property_readonly( "mincols", &mesh::MultiBlockConnectivity::mincols );
    py::class_<mesh::BlockConnectivity>( m, "BlockConnectivity" )
        .def( "__getitem__",
              []( mesh::BlockConnectivity const& c, std::tuple<idx_t, idx_t> const& pos ) {
                  auto const& [row, col] = pos;
                  return c( row, col );
              } )
        .def_property_readonly( "rows", &mesh::BlockConnectivity::rows )
        .def_property_readonly( "cols", &mesh::BlockConnectivity::cols );

    py::class_<mesh::Nodes>( m, "Nodes" )
        .def_property_readonly( "size", &mesh::Nodes::size )
        .def_property_readonly( "edge_connectivity",
                                py::overload_cast<>( &mesh::Nodes::edge_connectivity, py::const_ ) )
        .def_property_readonly( "cell_connectivity",
                                py::overload_cast<>( &mesh::Nodes::cell_connectivity, py::const_ ) )
        .def_property_readonly( "lonlat", py::overload_cast<>( &Mesh::Nodes::lonlat, py::const_ ) )
        .def_property_readonly( "xy", py::overload_cast<>( &Mesh::Nodes::xy, py::const_ ) )
        .def("field", []( mesh::Nodes const& n, std::string const& name ) { return n.field( name ); },
             "name"_a, py::return_value_policy::reference_internal )
        .def( "flags", []( mesh::Nodes const& n ) { return n.flags(); },
              py::return_value_policy::reference_internal);

    py::class_<mesh::HybridElements>( m, "HybridElements" )
        .def_property_readonly( "size", &mesh::HybridElements::size )
        .def( "nb_nodes", &mesh::HybridElements::nb_nodes )
        .def( "nb_edges", &mesh::HybridElements::nb_edges )
        .def_property_readonly( "node_connectivity",
                                py::overload_cast<>( &mesh::HybridElements::node_connectivity, py::const_ ) )
        .def_property_readonly( "edge_connectivity",
                                py::overload_cast<>( &mesh::HybridElements::edge_connectivity, py::const_ ) )
        .def_property_readonly( "cell_connectivity",
                                py::overload_cast<>( &mesh::HybridElements::cell_connectivity, py::const_ ) )

        .def("field", []( mesh::HybridElements const& he, std::string const& name ) { return he.field( name ); },
            "name"_a, py::return_value_policy::reference_internal )
        .def( "flags", []( mesh::HybridElements const& he ) { return he.flags(); },
              py::return_value_policy::reference_internal );

    auto m_fs = m.def_submodule( "functionspace" );
    py::class_<FunctionSpace>( m_fs, "FunctionSpace" )
        .def_property_readonly( "size", &FunctionSpace::size )
        .def_property_readonly( "type", &FunctionSpace::type )
        .def(
            "create_field",
            []( FunctionSpace const& fs, std::optional<py::object> dtype, py::kwargs kwargs ) {
                util::Config config = to_config(kwargs);
                if ( dtype )
                    config.set( option::datatype( pybindToAtlas( py::dtype::from_args( *dtype ) ) ));
                else
                    config.set( option::datatypeT<double>() );
                return fs.createField( config );
            },
            "dtype"_a = std::nullopt)
        .def_property_readonly("lonlat", &FunctionSpace::lonlat );
    py::class_<functionspace::EdgeColumns, FunctionSpace>( m_fs, "EdgeColumns" )
        .def( py::init( []( Mesh const& m, int halo, int levels ) {
                  return functionspace::EdgeColumns( m, util::Config()( "halo", halo )("levels", levels) );
              } ),
              "mesh"_a, "halo"_a = 0, "levels"_a = 0 )
        .def_property_readonly( "nb_edges", &functionspace::EdgeColumns::nb_edges )
        .def_property_readonly( "mesh", &functionspace::EdgeColumns::mesh )
        .def_property_readonly( "edges", &functionspace::EdgeColumns::edges )
        .def_property_readonly( "valid", &functionspace::EdgeColumns::valid );
    py::class_<functionspace::NodeColumns, FunctionSpace>( m_fs, "NodeColumns" )
        .def( py::init( []( FunctionSpace fs ) { return functionspace::NodeColumns{fs}; } ) )
        .def( py::init( []( Mesh const& m, int halo, int levels ) {
                  return functionspace::NodeColumns( m, util::Config()( "halo", halo )("levels",levels) );
              } ),
              "mesh"_a, "halo"_a = 0, "levels"_a = 0 )
        .def_property_readonly( "nb_nodes", &functionspace::NodeColumns::nb_nodes )
        .def_property_readonly( "mesh", &functionspace::NodeColumns::mesh )
        .def_property_readonly( "nodes", &functionspace::NodeColumns::nodes )
        .def_property_readonly( "valid", &functionspace::NodeColumns::valid )
        .def( "__bool__", [](const functionspace::NodeColumns& self) { return self.valid(); } );

    py::class_<functionspace::CellColumns, FunctionSpace>( m_fs, "CellColumns" )
        .def( py::init( []( Mesh const& m, int halo, int levels ) {
                  return functionspace::CellColumns( m, util::Config()( "halo", halo )("levels",levels) );
              } ),
              "mesh"_a, "halo"_a = 0, "levels"_a = 0 )
        .def_property_readonly( "nb_cells", &functionspace::CellColumns::nb_cells )
        .def_property_readonly( "mesh", &functionspace::CellColumns::mesh )
        .def_property_readonly( "cells", &functionspace::CellColumns::cells )
        .def_property_readonly( "valid", &functionspace::CellColumns::valid );

    py::class_<functionspace::Spectral, FunctionSpace>( m_fs, "Spectral" )
        .def( py::init( []( FunctionSpace fs ) { return functionspace::Spectral{fs}; } ) )
        .def( py::init( []( int truncation ) { return functionspace::Spectral( truncation ); } ), "truncation"_a )
        .def_property_readonly( "nb_spectral_coefficients", &functionspace::Spectral::nb_spectral_coefficients )
        .def_property_readonly( "nb_spectral_coefficients_global", &functionspace::Spectral::nb_spectral_coefficients_global )
        .def_property_readonly( "truncation", &functionspace::Spectral::truncation )
        .def_property_readonly( "valid", &functionspace::Spectral::valid )
        .def( "__bool__", [](const functionspace::Spectral& self) { return self.valid(); } )
        .def("parallel_for", []( const functionspace::Spectral& self, const py::function &f) {
            self.parallel_for<std::function<void(idx_t,idx_t,int,int)>>(f);});


    py::class_<functionspace::StructuredColumns, FunctionSpace>( m_fs, "StructuredColumns" )
        .def( py::init( []( Grid const& g, grid::Partitioner const& p, py::kwargs kwargs ) {
                  return functionspace::StructuredColumns( g, p, to_config(kwargs) ); } ),
              "grid"_a, "partitioner"_a )
        .def( py::init( []( Grid const& g, py::kwargs kwargs ) { return functionspace::StructuredColumns( g, to_config(kwargs) ); } ), "grid"_a )
        .def_property_readonly( "grid", &functionspace::StructuredColumns::grid )
        .def_property_readonly( "valid", &functionspace::StructuredColumns::valid );

    py::class_<util::Metadata>( m, "Metadata" )
        .def_property_readonly( "keys", &util::Metadata::keys )
        .def( "__setitem__",
              []( util::Metadata& metadata, std::string const& key, py::object value ) {
                  if ( py::isinstance<py::bool_>( value ) )
                      metadata.set( key, value.cast<bool>() );
                  else if ( py::isinstance<py::int_>( value ) )
                      metadata.set( key, value.cast<long long>() );
                  else if ( py::isinstance<py::float_>( value ) )
                      metadata.set( key, value.cast<double>() );
                  else if ( py::isinstance<py::str>( value ) )
                      metadata.set( key, value.cast<std::string>() );
                  else
                      throw std::out_of_range( "type of value unsupported" );
              } )
        .def( "__getitem__",
              []( util::Metadata& metadata, std::string const& key ) -> py::object {
                  if ( !metadata.has( key ) )
                      throw std::out_of_range( "key <" + key + "> could not be found" );
                  return toPyObject( metadata.get().element( key ) );
              } )
        .def( "__repr__", []( util::Metadata const& metadata ) {
            return "_atlas4py.Metadata("_s + py::str( toPyObject( metadata.get() ) ) + ")"_s;
        } );

    py::class_<Field>( m, "Field", py::buffer_protocol() )
        .def_property_readonly( "name", &Field::name )
        .def_property_readonly( "strides", &Field::strides )
        .def_property_readonly( "shape", py::overload_cast<>( &Field::shape, py::const_ ) )
        .def_property_readonly( "size", &Field::size )
        .def_property_readonly( "rank", &Field::rank )
        .def_property_readonly( "levels", &Field::levels )
        .def_property_readonly( "datatype", []( Field& f ) { return atlasToPybind( f.datatype() ); } )
        .def_property( "metadata", py::overload_cast<>( &Field::metadata, py::const_ ),
                       py::overload_cast<>( &Field::metadata ) )
        .def_property_readonly( "functionspace", py::overload_cast<>( &Field::functionspace, py::const_ ) )
        .def( "halo_exchange", []( Field& f) { f.haloExchange(); })
        .def_property_readonly( "halo_dirty", []( Field& f) { f.haloExchange(); })
        .def_property( "halo_dirty", &Field::dirty, &Field::set_dirty, py::return_value_policy::copy)
        .def_buffer( []( Field& f ) {
            auto strides = f.strides();
            std::transform( strides.begin(), strides.end(), strides.begin(),
                            [&]( auto const& stride ) { return stride * f.datatype().size(); } );
            return py::buffer_info( f.storage(), f.datatype().size(), atlasToPybind( f.datatype() ), f.rank(),
                                    f.shape(), strides );
        } );

    py::class_<mesh::Nodes::Topology> topology( m, "Topology" );
    topology.attr( "NONE" )     = py::cast( int( mesh::Nodes::Topology::NONE ) );
    topology.attr( "GHOST" )    = py::cast( int( mesh::Nodes::Topology::GHOST ) );
    topology.attr( "PERIODIC" ) = py::cast( int( mesh::Nodes::Topology::PERIODIC ) );
    topology.attr( "BC" )       = py::cast( int( mesh::Nodes::Topology::BC ) );
    topology.attr( "WEST" )     = py::cast( int( mesh::Nodes::Topology::WEST ) );
    topology.attr( "EAST" )     = py::cast( int( mesh::Nodes::Topology::EAST ) );
    topology.attr( "NORTH" )    = py::cast( int( mesh::Nodes::Topology::NORTH ) );
    topology.attr( "SOUTH" )    = py::cast( int( mesh::Nodes::Topology::SOUTH ) );
    topology.attr( "PATCH" )    = py::cast( int( mesh::Nodes::Topology::PATCH ) );
    topology.attr( "POLE" )     = py::cast( int( mesh::Nodes::Topology::POLE ) );
    topology.def_static( "reset", &mesh::Nodes::Topology::reset );
    topology.def_static( "set", &mesh::Nodes::Topology::set );
    topology.def_static( "unset", &mesh::Nodes::Topology::unset );
    topology.def_static( "toggle", &mesh::Nodes::Topology::toggle );
    topology.def_static( "check", &mesh::Nodes::Topology::check );
    topology.def_static( "check_all", &mesh::Nodes::Topology::check_all );
    topology.def_static( "check_any", &mesh::Nodes::Topology::check_any );

    py::class_<output::Gmsh>( m, "Gmsh" )
        .def( py::init( []( std::string const& path ) { return output::Gmsh{ path }; } ), "path"_a )
        .def( py::init( []( std::string const& path, eckit::Configuration const& config, py::kwargs kwargs ) {
            util::Config cfg = util::Config(config);
            cfg.set(to_config(kwargs));
            return output::Gmsh(path,cfg);
          }), "path"_a, "config"_a )
        .def( py::init( []( std::string const& path, py::kwargs kwargs ) {return output::Gmsh{ path, to_config(kwargs) }; } ), "path"_a )
        .def( "__enter__", []( output::Gmsh& gmsh ) { return gmsh; } )
        .def( "__exit__", []( output::Gmsh& gmsh, py::object exc_type, py::object exc_val,
                              py::object exc_tb ) { gmsh.reset( nullptr ); } )
        .def(
            "write", []( output::Gmsh& gmsh, Mesh const& mesh ) { gmsh.write( mesh ); return gmsh; }, "mesh"_a )
        .def(
            "write", []( output::Gmsh& gmsh, Field const& field ) { gmsh.write( field );  return gmsh;}, "field"_a )
        .def(
            "write", []( output::Gmsh& gmsh, Field const& field, FunctionSpace const& fs ) { gmsh.write( field, fs );  return gmsh; },
            "field"_a, "functionspace"_a );


    py::class_<trans::Trans>( m, "Trans" )
        .def( py::init( [](const FunctionSpace& gp, const FunctionSpace& sp){ return trans::Trans(gp,sp);} ), "gp"_a, "sp"_a )
        .def( py::init( [](const Grid& grid, int truncation){ return trans::Trans(grid,truncation);} ), "grid"_a, "truncation"_a )
        .def( "dirtrans", []( trans::Trans& trans, const Field& gpfield, Field& spfield) { trans.dirtrans(gpfield,spfield);} )
        .def( "invtrans", []( trans::Trans& trans, const Field& spfield, Field& gpfield) { trans.invtrans(spfield,gpfield);} )
        .def_property_readonly( "truncation", &trans::Trans::truncation )
        .def_property_readonly( "nb_spectral_coefficients", &trans::Trans::spectralCoefficients )
        .def_static("backend", [] ( const std::string& backend ){ trans::Trans::backend(backend); } )
        .def_static("has_backend", [] ( const std::string& backend ){ return trans::Trans::hasBackend(backend); } );

    py::class_<Interpolation>( m, "Interpolation" )
        .def( py::init( [](const std::string& type, const FunctionSpace& source, const FunctionSpace& target, py::kwargs kwargs){
                auto config = to_config(kwargs);
                config.set("type",type);
                return Interpolation(config,source,target);
            } ), "type"_a, "source"_a, "target"_a )
        .def( py::init( [](const std::string& type, const Grid& source, const Grid& target, py::kwargs kwargs){
                auto config = to_config(kwargs);
                config.set("type",type);
                return Interpolation(config,source,target);
            } ), "type"_a, "source"_a, "target"_a )
        .def( "execute", []( Interpolation const& self, const Field& source, Field& target) { return self.execute(source,target);} )
        .def_property_readonly( "source", &Interpolation::source )
        .def_property_readonly( "target", &Interpolation::target );


    auto m_function = m.def_submodule( "function" );
    m_function.def("vortex_rollup", [](double lon, double lat, double t) { return util::function::vortex_rollup(lon,lat,t); } );
    m_function.def("spherical_harmonic", [](double lon, double lat,int n, int m ) { return util::function::spherical_harmonic(n,m,lon,lat); }, "lon"_a, "lat"_a, "n"_a, "m"_a );

}
