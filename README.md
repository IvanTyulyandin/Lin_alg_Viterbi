# Lin_alg_Viterbi

Viterbi algorithm implementations using linear algebra.
The main purpose is to check if specialization can give performance improvement.

## Building

To build and run code, one have to install
[SuiteSparse:GraphBLAS](https://github.com/DrTimothyAldenDavis/GraphBLAS).
Since that library uses `OpenMP`, user also should install it.

The script `build.sh` is available and can be used as a building guide.

## Implementations

There are implementations of the Viterbi algorithm in terms of linear algebra:

1. Default: `LA_Viterbi.h` and `LA_Viterbi.cpp`.
2. Specialized: `LA_Viterbi_spec.h` and `LA_Viterbi_spec.cpp`
3. Specialized: `LA_Viterbi_assoc_spec.h` and `LA_Viterbi_assoc_spec.cpp`

The second version precalculates matrix multipications only for one observation.
The third one is an updated version of the second.
Since Min_plus semiring is associative, observation matrices can be combined and precalculated.

## Testing

Folder `tests` contains some tests. To run them,
use script `run_tests.sh`.
If flag `-v` is passed to the script,
all tests will run with `Valgrind` to detect memory-leaks.
