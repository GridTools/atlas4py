
#include "_mpi.hpp"

#include <sstream>
#include <map>
#include <stdexcept>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "atlas/library/config.h"
#include "atlas/parallel/mpi/mpi.h"
#include "eckit/mpi/Comm.h"
#include "_atlas4py_common.hpp"

namespace atlas4py {
namespace mpi {
namespace {

const eckit::mpi::Comm& comm_world() {
    static const eckit::mpi::Comm& world = atlas::mpi::comm("world");
    return world;
}
const eckit::mpi::Comm& comm_self() {
    static const eckit::mpi::Comm& self = atlas::mpi::comm("self");
    return self;
}
int comm_world_int() {
    static int world = comm_world().communicator();
    return world;
}
int comm_self_int() {
    static int self = comm_self().communicator();
    return self;
}

std::string comm_name_from_int( int comm_as_int ) {
    if ( comm_as_int == comm_world_int() ) {
        return "world";
    }
    if ( comm_as_int == comm_self_int() ) {
        return "self";
    }
    std::ostringstream s;
    s << "int." << comm_as_int;
    return s.str();
}

enum class RegisterIfNotFound {
    yes,
    no
};
const eckit::mpi::Comm& lookup_comm_in_eckit( int comm_as_int, RegisterIfNotFound register_if_not_found ) {
    if ( comm_as_int == comm_world_int() ) {
        return comm_world();
    }
    else if ( comm_as_int == comm_self_int() ) {
        return comm_self();
    }
    else {
        for( std::string name : eckit::mpi::listComms() ) {
            const auto& comm = atlas::mpi::comm(name);
            if ( comm.communicator() == comm_as_int ) {
                return comm;
            }
        }
    }
    if ( register_if_not_found == RegisterIfNotFound::yes ) {
        std::string name = comm_name_from_int(comm_as_int);
        eckit::mpi::addComm( name.c_str(), comm_as_int );
        return atlas::mpi::comm( name );
    }
    throw std::out_of_range( "comm with int value " + std::to_string(comm_as_int) + " not found in eckit and register_if_not_found is false" );
}

int comm_to_int( const nb::object& comm ) {
    if ( nb::isinstance<nb::int_>( comm ) ) {
        return nb::cast<int>(comm);
    }
    else if (nb::hasattr(comm, "toint")) { // This is supporting mpi4py.MPI.Comm
        return nb::cast<int>(comm.attr("toint")());
    }
    else if (nb::hasattr(comm, "to_int")) { // This is supporting atlas4py.mpi.Comm
        return nb::cast<int>(comm.attr("to_int")());
    }
    throw std::out_of_range( "unsupported type for comm" );
}

}

bool has_comm( const std::string& name ) {
    return eckit::mpi::hasComm( name.c_str() );
}

Comm register_comm( const std::string& name, int comm ) {
    eckit::mpi::addComm( name.c_str(), comm );
    return Comm(name);
}

#if ATLAS_ECKIT_VERSION_AT_LEAST(2,0,0)
void unregister_comm( const std::string& name ) {
    eckit::mpi::unregisterComm( name.c_str() );
}
#endif

void delete_comm( const std::string& name ) {
    eckit::mpi::deleteComm( name.c_str() ); // Note, this calls MPI_Comm_free!
}

void set_default_comm( const std::string& name) {
    eckit::mpi::setCommDefault( name.c_str() );
}

Comm::Comm() : comm_( atlas::mpi::comm() ) {}
Comm::Comm( const std::string& name ) : comm_( atlas::mpi::comm(name) ) {}
Comm::Comm( const eckit::mpi::Comm& comm ) : comm_( comm ) {}
Comm::Comm( int comm ) : comm_( lookup_comm_in_eckit(comm, RegisterIfNotFound::yes) ) {}

int Comm::size() const {
    if (size_ < 0) {
        size_ = comm_.size();
    }
    return size_;
}
int Comm::rank() const {
    if (rank_ < 0) {
        rank_ = comm_.rank();
    }
    return rank_;
}
void Comm::barrier() const { comm_.barrier(); }
void Comm::abort(int errorcode) const { comm_.abort(errorcode); }
const std::string& Comm::name() const {
    if (name_.empty()) {
        name_ = comm_.name();
    }
    return name_;
}
int Comm::to_int() const {
    if (comm_int_ == std::numeric_limits<int>::min()) {
        comm_int_ = comm_.communicator();
    }
    return comm_int_;
}

const Comm& Comm::WORLD() {
    static const Comm world = mpi::Comm( "world" );
    return world;
}
const Comm& Comm::SELF() {
    static const Comm self = mpi::Comm( "self" );
    return self;
}

Comm Comm::split(int color, const std::string& name) const {
    return Comm( comm_.split(color, name.c_str()) );
}

class Scope {
    public:
        Scope( const std::string& name ) {
            push(name);
        }
        Scope( const mpi::Comm& comm ) : Scope(comm.name()) {}
        Scope( const int comm ) : Scope(mpi::Comm(comm).name()) {}
        Scope( const nb::object& comm ) : Scope(mpi::Comm(mpi::comm_to_int(comm)).name()) {}
        ~Scope() {
            if (not ended_) {
                pop();
            }
            ended_ = true;
        }

        void reset( const std::string& name ) {
            pop();
            push(name);
        }
        void reset( const mpi::Comm& comm ) {
            pop();
            push(comm.name());
         }
        void reset( std::nullptr_t ) {
            pop();
            ended_ = true;
        }
        static void push(const std::string& name) {
            #if ATLAS4PY_ATLAS_VERSION_AT_LEAST(0,46,0)
                atlas::mpi::scope::push(name);
            #else
                atlas::mpi::push(name);
            #endif
        }
        static void pop() {
            #if ATLAS4PY_ATLAS_VERSION_AT_LEAST(0,46,0)
                atlas::mpi::scope::pop();
            #else
                atlas::mpi::pop();
            #endif
        }
        bool ended_ = false;
};

} // namespace mpi

namespace nb = ::nanobind;
void bind_submodule_mpi(nb::module_ &m) {
    using namespace nb::literals;
    nb::module_ m_mpi = m.def_submodule("mpi", "mpi python binding.");


    m_mpi.def("set_default_comm", [](const std::string& name) { mpi::set_default_comm(name); }, "name"_a);
    m_mpi.def("has_comm", &mpi::has_comm, "name"_a );
    m_mpi.def("list_comms", &eckit::mpi::listComms);
    m_mpi.def("register_comm", [](const std::string& name, int comm) { return mpi::register_comm(name, comm); }, "name"_a, "comm"_a );
    m_mpi.def("register_comm", [](const std::string& name, nb::object comm) { return mpi::register_comm(name, mpi::comm_to_int(comm)); }, "name"_a, "comm"_a );
#if ATLAS_ECKIT_VERSION_AT_LEAST(2,0,0)
    m_mpi.def("unregister_comm", [](const std::string& name) { mpi::unregister_comm(name); }, "name"_a );
    m_mpi.def("unregister_comm", [](mpi::Comm& comm) { mpi::unregister_comm(comm.name()); }, "comm"_a );
#endif
    m_mpi.def("delete_comm", &mpi::delete_comm, "name"_a );
    m_mpi.def("comm", []() { return mpi::Comm(); } );
    m_mpi.def("comm", [](const std::string& name) { return mpi::Comm(name); } );
    m_mpi.def("comm", [](int comm) { return mpi::Comm(comm); } );
    m_mpi.def("comm", [](nb::object comm) { return mpi::Comm(mpi::comm_to_int(comm)); } );
    m_mpi.def("size", []() { return atlas::mpi::size(); } );
    m_mpi.def("rank", []() { return atlas::mpi::rank(); } );
    m_mpi.def("barrier", []() { atlas::mpi::comm().barrier(); } );
    m_mpi.def("abort", [](int errorcode = -1) { atlas::mpi::comm().abort(errorcode); }, "errorcode"_a = -1 );
    m_mpi.def("finalize", []() { atlas::mpi::finalize(); } );
    m_mpi.def("push", [](const std::string& name) { mpi::Scope::push(name); }, "name"_a );
    m_mpi.def("pop",  []() { mpi::Scope::pop(); } );

    nb::class_<mpi::Comm>( m_mpi, "Comm" )
        .def( nb::init<>() )
        .def_prop_ro("size", &mpi::Comm::size )
        .def_prop_ro("rank", &mpi::Comm::rank )
        .def_prop_ro("name", &mpi::Comm::name )
        .def("barrier", &mpi::Comm::barrier )
        .def("abort", &mpi::Comm::abort, "errorcode"_a = -1 )
        .def("split", &mpi::Comm::split, "color"_a, "name"_a )
        .def("to_int", &mpi::Comm::to_int )
        .def("__int__", &mpi::Comm::to_int )
        .def("__eq__", [](const mpi::Comm &self, const mpi::Comm &other) {
            return self.to_int() == other.to_int();
        } )
        .def("__repr__", []( const mpi::Comm& comm ) {
            return "_atlas4py.mpi.Comm(name=" + comm.name() + ", size=" + std::to_string( comm.size() ) +
                   ", rank=" + std::to_string( comm.rank() ) + ", int=" + std::to_string( comm.to_int() ) + ")";
        } );

    nb::class_<mpi::Scope>( m_mpi, "Scope" )
        .def( nb::init<const std::string&>(), "name"_a )
        .def( nb::init<const mpi::Comm&>(), "comm"_a )
        .def( nb::init<int>(), "comm"_a )
        .def( nb::init<nb::object>(), "comm"_a )
        .def( "__enter__", []( mpi::Scope& self ) { return &self; } )
        .def( "__exit__", []( mpi::Scope& self, nb::object exc_type, nb::object exc_value, nb::object traceback ) { self.reset(nullptr); }, "exc_type"_a = nb::none(), "exc_value"_a = nb::none(), "traceback"_a = nb::none() );

    m_mpi.def("__getattr__", [m_mpi](nb::handle name) -> nb::object {
        // Check if the requested attribute is the one we want to lazily initialize
        std::string name_str = nb::cast<std::string>(name);
        if (name_str == "COMM_WORLD" || name_str == "COMM_SELF") {
            const mpi::Comm& comm = [name_str]() -> const mpi::Comm& {
                if (name_str == "COMM_WORLD") {
                    return mpi::Comm::WORLD();
                }
                else { // name_str == "COMM_SELF"
                    return mpi::Comm::SELF();
                }
            }();
            nb::object value = nb::cast(comm, nb::rv_policy::reference); // policy::reference tells nanobind NOT to manage the memory.

            // 2. Cache it in the module's __dict__ for fast subsequent access
            // This makes it act like a read-only property after the first call
            nb::getattr(nb::handle(m_mpi), "__dict__")[name] = value;
            return value;
        }
        // Raise AttributeError for other missing attributes
        throw nb::attribute_error();
    });
}


}  // namespace atlas4py
