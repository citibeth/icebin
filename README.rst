Introduction
============

IceBin is a coupler for GCMs and Ice Models.  It also provides
capabilities for regridding with elevation classes, and with
generalized grids.

License
=======
See the file LICENSE for details on IceBin licensing.

Requirements
============

We support building and running IceBin on Linux.  IceBin should also
work on Macintosh; however, there is not currently an easy way to
build it.

IceBin is built with Spack, which requires a functional Python 2.6 or
Python 2.7 on your system.  Other standard build tools are assumed as
well: make, git, etc.  All other prerequisites are built with Spack.


Download and Install
====================

In order to install IceBin's numerous dependencies, we recommend using
Spack to install IceBin.  Instructions for Linux are provided here.  For more information on Spack, see: http://software.llnl.gov/spack

Install Spack
--------------

1. Download::

    cd ~
    # git clone git@github.com:citibeth/spack.git -b efischer/icebin
    git clone https://github.com/citibeth/spack.git -b efischer/icebin

2. Add to your ``.bashrc`` file::

    export SPACK_ROOT=$HOME/spack
    . $SPACK_ROOT/share/spack/setup-env.sh


Set up Spack Compilers
----------------------

IceBin is known to work with GCC 4.9.3.  In theory, it works with
other C++11 compilers as well.  Enter the following, to see what
compilers Spack has found on your system::

    spack compilers

This produces a file ``~/.spack/compilers.yaml``, which looks as
follows on CentOS 7::

    compilers:
      linux-x86_64:
        gcc@4.8.5:
          cc: /usr/bin/gcc
          cxx: /usr/bin/g++
          f77: /usr/bin/gfortran
          fc: /usr/bin/gfortran

If you are happy with the compilers Spack found, you can proceed.  Otherwise, you need to build GCC 4.9.3::

    spack install gcc@4.9.3

Once that completes, add GCC 4.9.3 to the ``compilers.yaml`` file:

1. Identify the location of your new GCC with::

    spack location -i gcc@4.9.3

2. Add it to your compilers.yaml file, to look something like::

    compilers:
      linux-x86_64:
        gcc@4.8.5:
          cc: /usr/bin/gcc
          cxx: /usr/bin/g++
          f77: /usr/bin/gfortran
          fc: /usr/bin/gfortran
        gcc@4.9.3:
          cc: /home/rpfische/spack/opt/spack/linux-x86_64/gcc-4.8.5/gcc-4.9.3-layphctulnk3omsbjpzftqv6dlxpfe3d/bin/gcc
          cxx: /home/rpfische/spack/opt/spack/linux-x86_64/gcc-4.8.5/gcc-4.9.3-layphctulnk3omsbjpzftqv6dlxpfe3d/bin/g++
          f77: /home/rpfische/spack/opt/spack/linux-x86_64/gcc-4.8.5/gcc-4.9.3-layphctulnk3omsbjpzftqv6dlxpfe3d/bin/gfortran
          fc: /home/rpfische/spack/opt/spack/linux-x86_64/gcc-4.8.5/gcc-4.9.3-layphctulnk3omsbjpzftqv6dlxpfe3d/bin/gfortran


Configure Spack
---------------

Now it is time to tell Spack which compilers and package versions are preferred.  Do this by creating the file ``~/.spack/packages.yaml``.  It should look like this::

    packages:
        python:
            version: [3.5.1]
        py-cython:
            version: [0.23.4]

        netcdf-cxx4:
            version: [ecdf914]

        ibmisc:
            version: [0.1.0]

        icebin:
            version: [0.1.0]

        openssl:
            paths:
                openssl@system: /usr
            buildable: False

        all:
            compiler: [gcc@4.9.3]
            providers:
                mpi: [openmpi]
                blas: [atlas]
                lapack: [atlas]

A few things to note here:

1. The ``compiler`` section tells Spack which compilers to use, in
   preferred order.

2. The ``openssl`` section tells Spack to use the OS version of the
   OpenSSL library, rather than building one itself.  This is for
   security reasons.

   If you choose this route, Spack will later give you
   spurious warnings that look like::

        ==> Warning: This installation depends on an old version of OpenSSL,
                     which may have known security issues.
        ==> Warning: Consider updating to the latest version of this package.
        ==> Warning: More details at http://www.openssl.org

   You can safely ignore these warnings because they are false.

Install Git
-----------

You might wish to install the latest, greatest version of git.  Do this with::

    spack install git+curl+expat

Once Git is installed, make it available to Bash via::

    spack load git



Install IBMisc
-----------------

Spack can install packages automatically, or assist in building packages manually.  We will use Spack to automatically install all of IceBin's prerequisites, and then manually install IceBin and its support library IBMisc from GitHub.

1. Download the IBMisc library (support for IceBin)::

    cd ~
    git clone https://github.com/citibeth/ibmisc.git -b v0.1.0
    cd ibmisc

2. Ask Spack about the prerequisites for IBMisc::

    spack spec ibmisc@local +python +netcdf ^netcdf+mpi ^eigen~suitesparse ^py-numpy+lapack ^atlas ^python@3:

3. If this looks good, install the prerequisites (change ``spec`` to ``spconfig`` on the command line)::

    spack spconfig ibmisc@local +python +netcdf ^netcdf+mpi ^eigen~suitesparse ^py-numpy+lapack ^atlas ^python@3:

4. Now build IBMisc itself::

    mkdir build
    cd build
    ../spconfig.py ..
    make
    make install

Install IceBin
--------------

The manual install of IceBin itself is similar::

    cd ~
    git clone https://github.com/citibeth/ibmisc.git -b v0.1.0
    cd ibmisc

    spack spec icebin@local +gridgen +python ~coupler ~pism ^ibmisc@local ^netcdf+mpi ^eigen~suitesparse ^py-numpy+lapack ^atlas ^python@3:
    spack spconfig icebin@local +gridgen +python ~coupler ~pism ^ibmisc@local ^netcdf+mpi ^eigen~suitesparse ^py-numpy+lapack ^atlas ^python@3:

    mkdir build
    cd build
    ../spconfig.py ..
    make
    make install

Set Up Spack Python
-------------------

IceBin produces a Python extension.  The following Spack commands will install the Python modules necessary to run that extension::

    spack install py-cython ^python@3:
    spack activate py-cython
    spack install py-numpy+blas+lapack ^atlas ^python@3:
    spack activate py-numpy
    spack install py-scipy ^atlas ^python@3:
    spack activate py-scipy
