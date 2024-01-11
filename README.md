# simpleTBEM
## Configuration and Dependencies
The library can be configured by running 
~~~
cmake CMakeLists.txt
~~~
in the base directory or
~~~
cmake ..
~~~
in a subdirectory for an out-of-source build.
This should automatically generate a target with which the Eigen library for matrix computations can be built by running
~~~
make Eigen
~~~
if it is not already available.
The [complex_bessel library](https://github.com/joeydumont/complex_bessel)
that is used for passing complex arguments to the Hankel and Bessel functions 
aswell as the [arpackpp library](https://github.com/m-reuter/arpackpp) which gives a <tt>c++</tt> interface to the [arpack library](https://github.com/opencollab/arpack-ng) are installed automatically and don't need to be built.
<tt>arpack</tt> and <tt>lapack</tt> need to be installed separately and can usually be done so with your distributions packagemanager.

For <tt>arch</tt> based distros:
~~~
sudo pacman -S arpack
sudo pacman -S lapack
~~~
For <tt>debian</tt> based distros:
~~~
sudo apt install libboost-all-dev
sudo apt install libarpack2-dev 
sudo apt install liblapack3-dev
~~~

Running CMake also configures some examples of how to use the library as <tt>make</tt> targets.
These can then be built by running 

~~~
make <target_name>
~~~

The compiled binary can be found in the <tt>bin</tt> directory.

