
# Installation

## Prerequisites

* an implementation of the Message Passing Interface (MPI) standard such as [mpich](https://www.mpich.org/), [open-mpi](https://www.open-mpi.org/) or [Intel® MPI Library](https://software.intel.com/en-us/mpi-library)
* [autoconf](https://www.gnu.org/software/autoconf/)
* [automake 1.15](https://www.gnu.org/software/automake/)
* [make](https://www.gnu.org/software/make/)

The MPI compiler but be available in your PATH or set with the CC environment variable.

If you do not have automake 1.15 but a former version (such as 1.13), you can edit in the `configure.ac` file and change the line `AM_INIT_AUTOMAKE([1.15 foreign -Wall])` with `AM_INIT_AUTOMAKE([1.13 foreign -Wall])`.

If `automake` and `autoconf` have been installed in custom directories, be sure their are available in your PATH:

`export PATH=path_to_automake/automake-1.15/bin:path_to_autoconf/autoconf-2.69/bin:${PATH}`

If needed, you can set your PATH according to your configuration directly in your `${HOME}/.bashrc` file.


Custom options can be used with `configure` such as `--prefix` to set the destination installation path or `CC` for the MPI compiler, for example:
`./configure CC=mpi_bin_path --prefix`

This should be only what you need to know about how to use `./configure` but if you are interested, more details are available in the [README-configure](README-configure) and on the command line `./configure --help`.

## Build from the git repository

```
git clone https://github.com/bioinfo-pf-curie/mpiBWA.git
cd mpiBWA
# Checkout the branch of the version you want to install, for example:
# git checkout version-1.0.0
aclocal
autoconf
automake --add-missing
# If not yet in your PATH, you can provide the PATH to `mpicc`
# or your favourite MPI compiler at the configure stage
# using the CC environment variable, for example:
#./configure CC=/usr/lib64/mpich/bin/mpicc
./configure --prefix=${HOME}/local/mpiBWA
make
make install
```

## Build from a tar.gz archive

Download  the source code archive of the version you want to install from [https://github.com/bioinfo-pf-curie/mpiBWA/releases](https://github.com/bioinfo-pf-curie/mpiBWA/releases).

```
tar xzf mpibwa-mpibwaidx-1.0.tar.gz
cd mpbwa-1.0
# If not yet in your PATH, you can provide the PATH to `mpicc`
# or your favourite MPI compiler at the configure stage
# using the CC environment variable, for example:
#./configure CC=/usr/lib64/mpich/bin/mpicc
./configure --prefix=${HOME}/local/mpiBWA
make
make install
```


## Package the source code into tar.gz archive


```
git clone https://github.com/bioinfo-pf-curie/mpiBWA.git
cd mpiBWA
# Checkout the branch of the version you want to package, for example:
# git checkout version-1.0.0
aclocal
autoconf
automake --add-missing
# If not yet in your PATH, you can provide the PATH to `mpicc`
# or your favourite MPI compiler at the configure stage
# using the CC environment variable, for example:
#./configure CC=/usr/lib64/mpich/bin/mpicc
./configure --prefix=${HOME}/local/mpiBWA
make dist
```



## Integration further BWA-MEM release


```
# Clone the BWA repository
git clone https://github.com/lh3/bwa

# Clone the mpiBWA repository
git clone https://github.com/bioinfo-pf-curie/mpiBWA.git
# Checkout the branch of the version you want to install, for example:
# git checkout version-1.0.0

# Copy all the .c and .h from bwa repo into the src/ of mpibwa

cp bwa/*.c bwa/*.h mpiBWA/src/

# Now compile mpiBWA normally
cd mpiBWA
aclocal
autoconf
automake --add-missing
# If not yet in your PATH, you can provide the PATH to `mpicc`
# or your favourite MPI compiler at the configure stage
# using the CC environment variable, for example:
#./configure CC=/usr/lib64/mpich/bin/mpicc
./configure --prefix=${HOME}/local/mpiBWA
make
make install

```


