function p = plot3Dtissue(nCells, Coords, Connec,param)
% Draws the cells as patches; tracked and wound-margin cells are highlighted.

[Cnew, Vnew] = getCellDataforPlottingwithoutPeriJumps(nCells, Coords, Connec, param);

p = gobjects(nCells, 1); % Pre-allocate graphics object array

for cellID = 1:nCells
    
    % Draw only valid polygons (3 or more vertices)
    if length(Cnew{cellID}) > 2
        
        if ismember(cellID,param.SelfPropellingCellIDs) || ismember(cellID,param.cellIDtoContract) || ismember(cellID,param.cellIDstoTrack)
          p(cellID) =  patch('faces',Cnew{cellID},'vertices',Vnew , 'FaceColor','red','FaceAlpha',1) ;
        else
           p(cellID) = patch('faces',Cnew{cellID},'vertices',Vnew , 'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',1) ;
        end
        
    end
    
    axis image off;

end
end