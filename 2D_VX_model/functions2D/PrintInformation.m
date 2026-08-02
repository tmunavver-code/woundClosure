fprintf(1,'Finished timestep %d of %d\n',tstep,param.Nsteps);
fprintf(1,'Free energy of the system : %3.4f\n', energymat(tstep));
fprintf(1,'Relative change in energy : %3.4E\n', relenergychange);
fprintf(1,'Number of T1 transitions occured : %d\n', param.nT1);
fprintf(1,'Total simulation time (min): %2.3f\n\n', toc(totalsimtime)/60);