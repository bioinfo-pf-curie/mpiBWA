Installation
---------

git clone https://github.com/bioinfo-pf-curie/mpiBWA.git <br />
git checkout master  <br />
git pull <br />
export PATH=/PATH_TO/automake-1.15/bin:/PATH_TO/autoconf-2.69/bin:$PATH <br />
aclocal <br />
automake --add-missing <br />
autoconf <br />
 ./configure CC=/PATH_TO/mpicc <br />
make && make intall <br />

If you don't have 1.15 change in the configure.ac the line  <br />
AM_INIT_AUTOMAKE([1.15 foreign -Wall]) with AM_INIT_AUTOMAKE([1.13 foreign -Wall])  <br />
