#include "_atlas4py_Config.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "eckit/config/Configuration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/util/Config.h"

namespace nb = ::nanobind;

namespace {

nb::object _toPyObject( eckit::Configuration const& v );
nb::object _toPyObject( eckit::Configuration const& v, std::string const& key );

nb::object _toPyObject(bool v) {
    return nb::bool_(v);
}
nb::object _toPyObject(long v) {
    return nb::int_(v);
}
nb::object _toPyObject(double v) {
    return nb::float_(v);
}
nb::object _toPyObject(std::string const& v) {
    return nb::str(v.c_str());
}
template <typename T>
nb::object _toPyObject( std::vector<T> const& v ) {
    nb::list ret;
    for ( auto const& val : v ) {
        ret.append( _toPyObject( val ) );
    }
    return ret;
}

template <>
nb::object _toPyObject( std::vector<bool> const& v ) {
    nb::list ret;
    for ( auto const& val : v ) {
        ret.append( _toPyObject( val ) );
    }
    return ret;
}

nb::object _toPyObject( eckit::Configuration const& v, std::string const& key ) {
    if ( v.isSubConfiguration ( key ) ) {
        return _toPyObject( v.getSubConfiguration( key ) );
    }
    else if (v.isBoolean( key )) {
        return _toPyObject( v.getBool( key ) );
    }
    else if (v.isIntegral( key )) {
        return _toPyObject( v.getLong( key ) );
    }
    else if (v.isFloatingPoint( key )) {
        return _toPyObject( v.getDouble( key ) );
    }
    else if (v.isString( key )) {
        return _toPyObject( v.getString( key ) );
    }
    else if (v.isSubConfigurationList( key )) {
        std::vector<eckit::LocalConfiguration> subconfigs = v.getSubConfigurations( key );
        return _toPyObject( subconfigs );
    }
    else if (v.isBooleanList( key )) {
        auto getBooleanVector = [&]() {
            std::vector<long> values = v.getLongVector( key );
            std::vector<bool> bools;
            bools.reserve( values.size() );
            for ( auto const& val : values ) {
                bools.push_back( static_cast<bool>( val ) );
            }
            return bools;
        };
        return _toPyObject( getBooleanVector() );
    }
    else if (v.isIntegralList( key )) {
        std::vector<long> values = v.getLongVector( key );
        return _toPyObject( values );
    }
    else if (v.isFloatingPointList( key )) {
        std::vector<double> values = v.getDoubleVector( key );
        return _toPyObject( values );
    }
    else if (v.isStringList( key )) {
        std::vector<std::string> values = v.getStringVector( key );
        return _toPyObject( values );
    }
    else {
        throw std::out_of_range( "type of value unsupported for key " + key );
    }
}

nb::object _toPyObject( eckit::Configuration const& v ) {
    nb::dict ret;
    for ( auto const& key : v.keys()) {
        ret[ key.c_str() ] = _toPyObject( v, key );
    }
    return ret;
}


void config_set( eckit::LocalConfiguration& config, const std::string& key, nb::handle value ) {
    auto py_type_str = [](nb::handle h) -> std::string {
        nb::object type_obj = nb::borrow<nb::object>(Py_TYPE(h.ptr()));
        nb::object name = type_obj.attr("__name__");
        return nb::cast<std::string>(name);
    };

    if (nb::isinstance<nb::bool_>(value)) {
        config.set(key, nb::cast<bool>(value));
    } else if (nb::isinstance<nb::int_>(value)) {
        config.set(key, nb::cast<long long>(value));
    } else if (nb::isinstance<nb::float_>(value)) {
        config.set(key, nb::cast<double>(value));
    } else if (nb::isinstance<nb::str>(value)) {
        config.set(key, nb::cast<std::string>(value));
    } else if (nb::isinstance<nb::sequence>(value)) {
        nb::object seq = nb::cast<nb::object>(value);
        size_t n = nb::len(seq);
        if (n == 0) {
            config.set(key, std::vector<long>{});
            return;
        }
        auto elem = seq[0];
        auto handle_sequence = [&](auto&& conv) {
            using VecType = std::vector<std::decay_t<decltype(conv(elem))>>;
            VecType vec;
            for (size_t i = 0; i < n; ++i) vec.push_back(conv(seq[i]));
            config.set(key, vec);
        };
        if (nb::isinstance<nb::bool_>(elem)) {
            std::vector<int> vec;
            for (size_t i = 0; i < n; ++i) vec.push_back( nb::cast<bool>(seq[i]) ? 1 : 0 );
            config.set(key, vec);
        } else if (nb::isinstance<nb::int_>(elem)) {
            handle_sequence([](nb::handle v) { return nb::cast<long long>(v); });
        } else if (nb::isinstance<nb::float_>(elem)) {
            handle_sequence([](nb::handle v) { return nb::cast<double>(v); });
        } else if (nb::isinstance<nb::str>(elem)) {
            handle_sequence([](nb::handle v) { return nb::cast<std::string>(v); });
        } else if (nb::isinstance<nb::mapping>(elem)) {
            std::vector<eckit::LocalConfiguration> vec;
            for (size_t i = 0; i < n; ++i) {
                eckit::LocalConfiguration subconfig;
                nb::object mapping = nb::cast<nb::object>(seq[i]);
                for (nb::handle item : mapping.attr("keys")()) {
                    config_set(subconfig, nb::cast<std::string>(item), mapping[item]);
                }
                vec.push_back(subconfig);
            }
            config.set(key, vec);
        } else {
            throw std::out_of_range("Unsupported sequence element type for key '" + key + "': got type '" + py_type_str(elem) + "'");
        }
    } else if (nb::isinstance<nb::mapping>(value)) {
        eckit::LocalConfiguration subconfig;
        nb::object mapping = nb::cast<nb::object>(value);
        for (nb::handle item : mapping.attr("keys")()) {
            config_set(subconfig, nb::cast<std::string>(item), mapping[item]);
        }
        config.set(key, subconfig);
    } else {
        throw std::out_of_range("type of value unsupported for key '" + key + "': got type '" + py_type_str(value) + "'");
    }
}

std::string to_lowercase(const std::string& str) {
    std::string lowercase_str = str;
    std::transform(lowercase_str.begin(), lowercase_str.end(), lowercase_str.begin(),
                    [](unsigned char c) { return std::tolower(c); });
    return lowercase_str;
}

} // namespace

nb::object atlas4py::make_object( eckit::Configuration const& v ) {
    return _toPyObject( v );
}


atlas::util::Config atlas4py::make_Config( nb::kwargs const& kwargs ) {
    atlas::util::Config config;
    for( const auto& pair : kwargs ) {
        const auto key = nb::cast<std::string>(pair.first);
        const auto& value = pair.second;
        config_set(config, key, value);
    }
    return config;
}

void atlas4py::bind_Config( nb::module_& m ) {
    using namespace nanobind::literals;

    nb::class_<eckit::Configuration>( m, "eckit.Configuration" )
        .def( "keys", &eckit::Configuration::keys )
        .def( "__contains__",
              []( eckit::Configuration const& config, std::string const& key ) {
                  return config.has( key );
              } )
        .def( "__getitem__",
              []( eckit::Configuration const& config, std::string const& key ) -> nb::object {
                  if ( !config.has( key ) )
                      throw std::out_of_range( "key <" + key + "> could not be found" );
                  return _toPyObject( config, key );
              } )
        .def( "__iter__",
              []( eckit::Configuration const& config ) {
                  return nb::iter(nb::cast(config.keys()));
              } )
        .def( "__len__",
              []( eckit::Configuration const& config ) {
                  return config.keys().size();
              } )
        .def( "__repr__", []( eckit::Configuration const& config ) {
            return "_atlas4py.eckit.Configuration("_s + nb::str( make_object( config ) ) + ")"_s;
        } );

    nb::class_<eckit::LocalConfiguration, eckit::Configuration>( m, "eckit.LocalConfiguration" )
        .def( nb::init() )
        .def( "__setitem__",
              []( eckit::LocalConfiguration& config, std::string const& key, nb::object value ) {
                  config_set(config,key,value);
              } )
        .def( "__repr__", []( eckit::LocalConfiguration const& config ) {
            return "_atlas4py.eckit.LocalConfiguration("_s + nb::str( make_object( config ) ) + ")"_s;
        } );

     nb::class_<atlas::util::Config, eckit::LocalConfiguration>( m, "Config" )
        .def( nb::init() )
        .def_static( "from_kwargs", []( nb::kwargs kwargs ) {
            return make_Config(kwargs);
        } )
        .def_static( "from_yaml", []( std::string const& yaml ) {
            return atlas::util::Config( eckit::YAMLConfiguration(yaml) );
        }, "yaml"_a )
        .def_static( "from_json", []( std::string const& json ) {
            return atlas::util::Config( eckit::YAMLConfiguration(json) );
        }, "json"_a )
        .def_static( "from_file", []( const nb::object path, std::string const& format ) {
            // path accepts string but also path-like objects (e.g. pathlib.Path)
            if ( !format.empty() ) {
                auto format_lowercase = to_lowercase(format);
                if (format_lowercase != "yaml" && format_lowercase != "json") {
                    throw std::runtime_error("Only 'yaml' or 'json' format is supported");
                }
            }
            return atlas::util::Config( eckit::PathName{nb::cast<std::string>( nb::str(path) )} );
        }, "path"_a, "format"_a = "" )
        .def( "__repr__", []( atlas::util::Config const& config ) {
            return "_atlas4py.Config("_s + nb::str( make_object( config ) ) + ")"_s;
        } );
}
