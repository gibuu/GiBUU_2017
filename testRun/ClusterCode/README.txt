############################################################################
#     This directory contains two models of nuclear fragmentation          #
############################################################################
 
(1) code/coalescence: 
Simple phenomenological Phase-Space Coalescence.
Only Evaporation is optionaly included.

(2) code/smm:
Full Statistical Multifragmentation Model (SMM).
Fermi break-up, evaporation/fision included.
(provided by A. Botvina)


Usage of (1):
(a) Run GiBUU code as long as freeze-out has been achieved.
    For X+X@SIS-energies this can take up to time_max~100-200 fm/c
    For p+X@SIS-energies           -//-      time_max~150 fm/c
    During the run it is convenient to write the (real) 
    particle vector at different times (e.g. in time sequence of 5 or 10 fm/c) 
    in the following format (N=number of ensemple, M=number of subsequent run):
    ID,charge,Mass(GeV),x,y,z (fm),px,py,pz (GeV),N,M
    Optionally also the impact parameter is needed for each subsequent run. 
    In this case one should provide also this information in the file 
    containing the particle vector.
(b) The coalescence-code reads then the file containing the particle vector 
    and perfomes phase-space coalescence in an event-by-event basis.

Usage of (2):
(a) Run GiBUU code as long as the system has arrived to equilibrium. 
    For X+X@SIS-energies this can take up to time_max~50-60 fm/c
    For p+X@SIS-energies           -//-      time_max~100-150 fm/c
    This is automatically done during GiBUU-runs by choosing in the 
    GiBUU-jobCard SMM_Flag=true. After equilibration, information 
    on source(s) properties, e.g. mass, charge, excitation energy, etc. 
    of fireball and spectators are printed out. Also, the particle 
    vector is automatically printed out at different times after the 
    onset of equilibration.
(b) The SMM-Code needs only mass, charge, excitation energy and (in some 
    cases, e.g. expanding fireball) radial flow energy, and finally number 
    of MC-events. 
