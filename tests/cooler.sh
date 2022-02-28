#!/bin/sh -x

rm -f sol_cooler.pvd sol_cooler_t00*.vtu

${UGSHELL-ugshell} -numRefs 0 -ex ../cooler.lua && \
tail -n7 sol_cooler.pvd sol_cooler_t00*.vtu
