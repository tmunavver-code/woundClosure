function p = plot3Dtissue(nCells, Coords, Connec,param)

[Cnew, Vnew] = getCellDataforPlottingwithoutPeriJumps(nCells, Coords, Connec, param);

p = gobjects(nCells, 1); % Pre-allocate graphics object array

for cellID = 1:nCells
    
    % --- NEW: Fix for "white hole" bug ---
    % Only try to draw valid polygons (3+ vertices)
    if length(Cnew{cellID}) > 2
        
        if ismember(cellID,param.SelfPropellingCellIDs) || ismember(cellID,param.cellIDtoContract) || ismember(cellID,param.cellIDstoTrack)
          p(cellID) =  patch('faces',Cnew{cellID},'vertices',Vnew , 'FaceColor','red','FaceAlpha',1) ;
        else
           p(cellID) = patch('faces',Cnew{cellID},'vertices',Vnew , 'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',1) ;
        end
        
    end
    % --- End of fix ---
    
    axis image off;

end
end