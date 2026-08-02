% Plotting evolution of tissue
for tstep = 1:param.Nsteps
    % Applying stretch (does nothing when not wanted)
    applyStretchX;
    
    % Plotting every nVisualisation steps
    if mod(tstep,param.nVisualisation) == 0
        p1 = plot3Dtissue(celldata.nCells,Coordinates(:,2*tstep-1:2*tstep),Connectivity(:,tstep),param);
        rectangle('Position',[0 0 param.Lx param.Ly]);
        fprintf(1,'Plotting time step %d\n',tstep);
        pause(0.1)
        
        if tstep < param.Nsteps %% keeping the last patch
            for i = 1:length(p1)
                delete(p1(i)); % deleting patches to not overload figure
            end
        end
    end
end
