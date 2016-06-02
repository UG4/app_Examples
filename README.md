# Examples #

**UG4-App** providing example scripts for various problems.

Copyright 2010-2016 Goethe Center for Scientific Computing, University Frankfurt

Please install/clone this repository through UG4's package manager
[ughub](https://github.com/UG4/ughub):

    ughub install Examples


## How to use the example scripts ##
Once you compiled UG4 and after sourcing 'ugbash'
(cf. https://github.com/UG4/ughub/wiki#compilation-of-ug4),
you may execute a script e.g. like this (from any folder you like):

    ugshell -ex Examples/laplace.lua

A list of available parameters is shown using the -help option:

    ugshell -ex Examples/laplace.lua -help


To see a demonstration of refinement projectors you may execute the following:
	
	ugshell -ex Examples/laplace.lua -grid grids/sphere_2d.ugx

When opening the output (sol_laplace_2d.vtu) in Paraview or Visit you'll notice
the smooth boundary approximation and circular inner grid lines.