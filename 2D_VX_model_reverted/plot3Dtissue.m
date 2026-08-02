function p = plot3Dtissue(nCells, Coords, Connec,param)

[Cnew, Vnew] = getCellDataforPlottingwithoutPeriJumps(nCells, Coords, Connec, param);

for cellID = 1:nCells
    if ismember(cellID,param.SelfPropellingCellIDs) || ismember(cellID,param.cellIDtoContract) || ismember(cellID,param.cellIDstoTrack)
      p(cellID) =  patch('faces',Cnew{cellID},'vertices',Vnew , 'FaceColor','red','FaceAlpha',1) ;
    else
       p(cellID) = patch('faces',Cnew{cellID},'vertices',Vnew , 'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',1) ;
    end    
    axis image off;

end