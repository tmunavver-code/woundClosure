function [VertexSlaveToMastermap, augmentedCellMSmap] = findMasters(CoordsSlaves,CoordsMasters,Vsubset,Csubset,celldata, param)
% Given a slave vertex, Vertexmap will give the global ID of it's
%master

Lx = param.Lx;
Ly = param.Ly;
augmentedCellMSmap = cell(size(Csubset));
nMasters = length(CoordsMasters);
AllCoords = [CoordsMasters;CoordsSlaves];

vertexCoordsLxpos = [CoordsSlaves(:,1) + Lx, CoordsSlaves(:,2)];
vertexCoordsLxneg = [CoordsSlaves(:,1) - Lx, CoordsSlaves(:,2)];
vertexCoordsLypos = [CoordsSlaves(:,1), CoordsSlaves(:,2) + Ly];
vertexCoordsLyneg = [CoordsSlaves(:,1), CoordsSlaves(:,2) - Ly];

vertexCoordsTR = [CoordsSlaves(:,1) + Lx, CoordsSlaves(:,2) + Ly];
vertexCoordsBR = [CoordsSlaves(:,1) + Lx, CoordsSlaves(:,2) - Ly];
vertexCoordsTL = [CoordsSlaves(:,1) - Lx, CoordsSlaves(:,2) + Ly];
vertexCoordsBL = [CoordsSlaves(:,1) - Lx, CoordsSlaves(:,2) - Ly];

[a,shiftedRight] = ismembertol(vertexCoordsLxpos,CoordsMasters,1e-5,'ByRows',true);
[b,shiftedLeft] = ismembertol(vertexCoordsLxneg,CoordsMasters,1e-5,'ByRows',true);
[c,shiftedUp] = ismembertol(vertexCoordsLypos,CoordsMasters,1e-5,'ByRows',true);
[d,shiftedDown] = ismembertol(vertexCoordsLyneg,CoordsMasters,1e-5,'ByRows',true);

[e,shiftedTR] = ismembertol(vertexCoordsTR,CoordsMasters,1e-5,'ByRows',true);
[f,shiftedBR] = ismembertol(vertexCoordsBR,CoordsMasters,1e-5,'ByRows',true);
[g,shiftedTL] = ismembertol(vertexCoordsTL,CoordsMasters,1e-5,'ByRows',true);
[h,shiftedBL] = ismembertol(vertexCoordsBL,CoordsMasters,1e-5,'ByRows',true);



% just checkig that we don't have repitions of blanks
if find(a + b + c + d + e + f + g + h > 1 | a + b + c + d + e + f + g + h== 0) > 0
    error('Master-Slave map not done properly: either a master is missing, or there are two masters for a slave')
end

VertexSlaveToMastermap = shiftedRight + shiftedUp + shiftedDown + shiftedLeft + shiftedTR + shiftedBR + shiftedTL + shiftedBL;

% finding to which cell these slave vertices belonged to and adding the
% master ones
[~, idxOfSlaveInVsubset] = ismembertol(CoordsSlaves,Vsubset,1e-5,'ByRows',true);
[~, idxOfMasterInVsubset] = ismembertol(CoordsMasters,Vsubset,1e-5,'ByRows',true);
globalIdxOfMasters = 1:length(CoordsMasters) ;
globalIdxOfSlaves = [1:length(CoordsSlaves)] + nMasters;
for cellID = 1:length(Csubset)
    verticesmat = Csubset{cellID};
    slavebelongedTocell = ismember(idxOfSlaveInVsubset,verticesmat);
    masterbelongedTocell = ismember(idxOfMasterInVsubset,verticesmat);
    augmentedCellMSmap{cellID} = [ augmentedCellMSmap{cellID}, globalIdxOfMasters(masterbelongedTocell), globalIdxOfSlaves(slavebelongedTocell)]; 
    Barycenters = getCellBarycenter(cellID,augmentedCellMSmap,AllCoords);
    augmentedCellMSmap = rearrangeAnticlockwise(cellID,augmentedCellMSmap,AllCoords,Barycenters);
end
    

end

