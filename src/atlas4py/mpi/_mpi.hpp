#pragma once

/*
Module atlas4py.mpi

Interoperability with mpi4py or other MPI libraries is possible via the underlying MPI integer handle
```python
import atlas4py
from mpi4py import MPI

mpi4py_comm = MPI.Comm.fromint( int(atlas4py.mpi.comm()) )
atlas4py_comm = atlas4py.mpi.comm( mpi4py_comm.toint() )

# You can also register the communicator by name
atlas4py.mpi.add_comm( "my_comm", mpi4py_comm.toint() )

# and retrieve it later
atlas4py_my_comm = atlas4py.mpi.comm( "my_comm" )
```
*/

#include <string>

#include <nanobind/nanobind.h>

namespace eckit::mpi {
class Comm;
}
namespace atlas4py {
namespace nb = ::nanobind;
void bind_submodule_mpi(nb::module_ &m);

namespace mpi {
class Comm {
public:
    Comm();
    Comm( const std::string& name );
    Comm( const eckit::mpi::Comm& comm );
    explicit Comm( int comm );
    int size() const;
    int rank() const;
    void barrier() const;
    void abort(int errorcode = -1) const;
    Comm split(int color, const std::string& name) const;
    const std::string& name() const;
    int to_int() const;
    static const Comm& WORLD();
    static const Comm& SELF();

private:
    const eckit::mpi::Comm& comm_;
    mutable int size_ = -1;
    mutable int rank_ = -1;
    mutable int comm_int_ = std::numeric_limits<int>::min();
    mutable std::string name_;
};

} // namespace mpi

}  // namespace atlas4py
