Build on Ubuntu
===============

* Install all needed packages
    sudo apt-get install git cmake-curses-gui libboost-all-dev libgsl0-dev libopenbabel-dev gfortran

* Clone cluster optimization project and PaGMO
    git clone https://github.com/enikulenkov/magistracy-proj.git
    git clone https://github.com/enikulenkov/PaGMO.git

* Build and install PaGMO
    cd PaGMO
    mkdir build
    cd build
    ccmake ../
    ------------------------------------------------------------------
    Press 'c'
    Choose following options:
    BUILD_EXAMPLES                   OFF
    BUILD_MAIN                       ON
    BUILD_PYGMO                      OFF
    CMAKE_BACKWARDS_COMPATIBILITY    2.4
    CMAKE_BUILD_TYPE                 Debug
    CMAKE_INSTALL_PREFIX             /usr/local
    ENABLE_GSL                       OFF
    ENABLE_GTOP_DATABASE             OFF
    ENABLE_IPOPT                     OFF
    ENABLE_MPI                       OFF
    ENABLE_NLOPT                     OFF
    ENABLE_SNOPT                     OFF
    ENABLE_TESTS                     OFF
    EXECUTABLE_OUTPUT_PATH           
    INSTALL_HEADERS                  ON
    LIBRARY_OUTPUT_PATH              
    SYSTEM_M_LIBRARY                 /usr/lib/i386-linux-gnu/libm.so
    Press 'c' again
    Press 'g'
    ------------------------------------------------------------------
    make
    sudo make install
    cd ../src/Eigen/
    find . -maxdepth 1 -type f | grep -v CmakeLists.txt | while read i ; do sudo cp $i /usr/local/include/pagmo/Eigen ; done
    cd ../../..

* Build cluster optimization program
   cd magistracy-proj
   ln -s build/Options-linux Options
   cd TGROUP/Ir_N/
   make

Running cluster optimization
============================
  cd magistracy-proj/TGROUP/Ir_N

  There are several test data files on which cluster optimization can be
  performed. They are located in test_data directory and have .mat extension.
  To run cluster optimization for cluster with 40 atoms, run:

  ./run.sh test_data/FCC_40.mat

  Script will do cluster optimization using two methods: genetic algorithm and
  molecular dynamics method. Results will be located in res/YYYY-MM-DD-HH:MM
  folder. Content of the folder:
    - c_res.mat          Coordinates of atoms in resulting cluster configuration
                         obtained by genetic algorithm
    - c_vis.cml          Visualization of cluster configuration obtained by genetic
                         algorithm
                         using genetic algorithm
    - FCC_40.mat         Initial cluster configuration
    - orig.cml           Visualization of initial cluster configuration
    - gen_alg_prefs.ini  Settings of genetic algorithm
    - out_f.cml          Visualization of cluster configuration obtained by
                         molecular dynamics method
    - summary.txt        Contains found minimum value and consumed time by
                         both algorithms
   
    Not all result files of molecular dynamics method are copied to res
    folder. Files with *.ang *.ene *.xyz *.res extensions are remained in
    Ir_N folder.
