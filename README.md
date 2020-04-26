# SuperQuadP2

# Requirements
- OpenMP

# To Build
- `make`

Or if any changes were made to dependencies:

- `cmake CMakeLists.txt`
- `make` or `cmake --build .`


# To Run

- `./Superquad <num threads> <numnodes> [filename | -vol] [filename]
    - Specify a file to save all statistics to or use `... -vol filename` to only record the volume (saves CPU time as a 1-threaded run is not needed to calculate speedup).
- See  `azure-pipelines.yml` for examples.

