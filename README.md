# GREMLIN

*G*et *R*ams*E*s si*M*u*L*ation *I*nformatio*N*

Quick python class and functions for handling RAMSES simulations

Point it to the run path at initialization and it will try and open info files [first and last, but if you need the aexps it will open all or specified snaps], the .nml file [will get confused if there are several - can specify file name for these cases], hydro_file_descriptor.txt, count the available snapshots and get their numbers.

Everything gets put in intuitive-ish names etc, full .nml namelist accessible from object using f90nml python library.
