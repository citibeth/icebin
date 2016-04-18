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

Ubuntu 12.04
-------------

The following installation was enough to prepare a bare Ubuntu 12.04
machine to build the IceBin software stack with Spack:

    sudo apt-get install build-essential git curl m4 bison gettext \
        libssl-dev libjpeg-dev libpthread-stubs0-dev libxau-dev libqt4-dev



Download and Install
====================

In order to install IceBin's numerous dependencies, we recommend using
Spack to install IceBin.  Installing IceBin therefore involves three steps:

1. Install essential build tools (Spack, Environment Modules, GCC and git).

2. Use Spack to install IceBin and related packages.

Step 1 is general for any package one might install with Spack;
whereas step 2 is specific to IceBin.  Instructions for Linux are
provided here; this has not been tested or debugged on Macintosh.

Install Build Tools
``````````````````````

The build tools consist of Spack, which may be used to install
Environment Modules, GCC and git if necessary.  By using Spack to
bootstrap the build envrionment, one ensures that an up-to-date
version of the build tools is available, no matter what host operating
system is being used.  For more information on Spack, see:
http://software.llnl.gov/spack


Install Spack
--------------

1. Download::

    cd ~
    # git clone git@github.com:citibeth/spack.git -b efischer/icebin
    git clone https://github.com/citibeth/spack.git -b efischer/icebin

2. Add to your ``.bashrc`` file::

    export SPACK_ROOT=$HOME/spack
    . $SPACK_ROOT/share/spack/setup-env.sh

3. Remove non-system stuff from your ``PATH``, ``LD_LIBRARY_PATH`` and
   other environment variables, which can cause strange errors when
   building with Spack.

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


**NOTE**: ``spack compilers`` does not always detect all installed compiilers.  Make sure to check ``compilers.yaml`` after it has been auto-generated.  Errors like "Compiler 'gcc@4.6' does not support compiling C++ programs" indicate a problem in ``compilers.yaml``.

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

Create the file ``~/.spack/packages.yaml``.  It can look like this for now::

    packages:
        openssl:
            paths:
                openssl@system: /usr
            buildable: False

        all:
            compiler: [gcc@4.9.3]

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


Install Environment Modules
-------------------------------

In order to use Spack's generated environment modules, you must have
installed the *Environment Modules* package.  On many Linux
distributions, this can be installed from the vendor's repository.
For example: ```yum install environment-modules``
(Fedora/RHEL/CentOS).  If your Linux distribution does not have
Environment Modules, you can get it with Spack:

1. Install with::

    spack install environment-modules

2. Activate with::

    TMP=`tempfile`
    echo >$TMP
    MODULE_HOME=`spack location -i environment-modules`
     MODULE_VERSION=`ls -1 $MODULE_HOME/Modules | head -1`
     ${MODULE_HOME}/Modules/${MODULE_VERSION}/bin/add.modules <$TMP
     cp .bashrc $TMP
     echo "MODULE_VERSION=${MODULE_VERSION}" > .bashrc
     cat $TMP >>.bashrc

This adds to your ``.bashrc`` (or similar) files, enabling Environment
Modules when you log in.  Re-load your .bashrc (or log out and in
again), and then test that the ``module`` command is found with:

    module avail

Enable Spack Shell Support
--------------------------------

Spack shell support allows Spack to work with Environment Modules.
You can enable it by sourcing some files in the ``/share/spack``
directory.  Add this to your ``.bashrc`` file, after the environment
modules setup section:

.. code-block:: sh

    export SPACK_ROOT=$HOME/spack
    . $SPACK_ROOT/share/spack/setup-env.sh


You can put the above code in your ``.bashrc`` or ``.cshrc``, and
Spack's shell support will be available on the command line.

Log out and in again; you can now test this with a simple command like::

    spack find
    spack load gcc

You can also use environment modules directly, for example::

    module avail


Install Git
-----------

Older system-supplied versions of git do not provide features that are
necessary today.  You might wish to install the latest, greatest
version of git.  Do this with::

    spack install git+curl+expat

Once Git is installed, make it available to Bash via::

    spack load git

Install Curl
-------------

Older system-supplied versions of ``curl`` may not work to download
some packages, particularly those using the *https* protocol.  If this
is the case, you may use a newer Spack-installed ``curl`` instead::

    spack install curl
    spack load curl

Updated Binutils
-----------------

On Ubuntu 12.04, an updated ``binutils`` may be required to build some packages (eg, ``py-numpy``)::

    spack load binutils


Install IceBin Application Packages
````````````````````````````````````

The IceBin library has many build and run dependencies.  The
instructions below will install them all.

Configure Package Versions
-----------------------------

Now it is time to tell Spack which compilers and package versions are
preferred.  Do this by adding to ``~/.spack/packages.yaml`` so it
looks like this::

    packages:
        python:
            version: [3.5.1]
        py-cython:
            version: [0.23.4]
        py-proj:
            # Normal released version is buggy
            version: [citibeth-latlong2]

        netcdf-cxx4:
            version: [ecdf914]

        ibmisc:
            version: [0.1.0]

        icebin:
            version: [0.1.0]

        # Recommended for security reasons
        # Do not install OpenSSL as non-root user.
        openssl:
            paths:
                openssl@system: /usr
            buildable: False

        # Recommended, unless your system doesn't provide Qt4
        qt:
            paths:
                qt@system: /usr
            buildable: False

        all:
            compiler: [gcc@4.9.3]
            providers:
                mpi: [openmpi]
                blas: [openblas]
                lapack: [openblas]


Install IceBin
-----------------

.. code-block:: bash

    spack install icebin@0.1.3 +gridgen +python ~coupler ~pism \
        ^ibmisc@0.1.1 ^netcdf+mpi ^eigen~suitesparse ^py-numpy+lapack \
        ^openblas~shared ^python@3:

Additionally, download the IceBin source code for testing purposes::

    cd ~
    git clone https://github.com/citibeth/icebin.git -b develop
    cd icebin

Spack Python Stack
-------------------

IceBin produces a Python extension.  The following Spack commands will install the Python modules necessary to run that extension::

    spack install py-basemap ^py-matplotlib+gui+ipython ^py-numpy+blas+lapack ^openblas~shared ^python@3:
    spack install py-giss ^py-matplotlib+gui+ipython ^py-numpy+blas+lapack ^openblas~shared ^py-proj@citibeth-latlong2 ^python@3:
    spack activate py-numpy


Activate Stuff You Need
-----------------------


The following command will load the Spack-installed packages needed
for basic Python use of IceBin::

    module load `spack module find tcl icebin netcdf cmake@3.5.1`
    module load `spack module find --dependencies tcl py-basemap py-giss`

You can speed up shell startup by turning these into ``module load`` commands.

1. Cut-n-paste the script ``make_spackenv``::

    #!/bin/sh
    #
    # Generate commands to load the Spack environment

    SPACKENV=$HOME/spackenv.sh

    spack module find --shell tcl git icebin@local ibmisc netcdf cmake@3.5.1 >$SPACKENV
    spack module find --dependencies --shell tcl py-basemap py-giss >>$SPACKENV

2. Add the following to your ``.bashrc`` file::

    source $HOME/spackenv.sh
    # Preferentially use your checked-out Python source
    export PYTHONPATH=$HOME/icebin/pylib:$PYTHONPATH

3. Run ``sh make_spackenv`` whenever your Spack installation changes (including right now).



Test the Activation
----------------------

The loaded packages may be tested as follows::

    # These do not produce output
    python3 -c 'import cython'
    python3 -c 'import numpy'
    python3 -c 'import scipy'
    python3 -c 'import netCDF4'
    python3 -c 'import matplotlib'
    python3 -c 'from mpl_toolkits.basemap import Basemap'
    python3 -c 'import ibmisc'
    python3 -c 'import icebin'
    python3 -c 'import giss'

    # This does produce output...
    python3 -c 'import pyproj; pyproj.test()'
    which overlap

IceBin Demos
============

IceBin comes with two demonstrations of its regridding capabilities,
in the folders ``tests/test_conserv`` and ``tests/test_issm``.  Each
demo is self-contained, and is executed with the ``make`` command.
The user may inspect ``makefile`` to determine the steps being
demonstrated.

test_conserv
-------------

This tests all six regridding matrices between a ModelE (global) GCM
grid and the SeaRISE (local) ice grid.  Conservation is checked
quantitatively.  Plots are also generated showing the grids involved,
and the regridded fields.

test_issm
----------

This tests regridding bewteen an ISSM mesh (I) and the SeaRISE Grid
(A).  Plots are generated showing the grids involved and regridded
fields.  However:

1. This software is not able to directly plot functions on an ISSM
   mesh.  Only versions of the function on the SeaRISE grid are
   plotted at this time.

2. No quantiative measures of conservation are taken.

3. Construction of the regridding matrix ``AvI`` is currently only in
   Python.  It has not been ported into the C++ code.
