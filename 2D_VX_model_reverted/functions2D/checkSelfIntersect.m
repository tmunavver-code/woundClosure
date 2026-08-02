function [ celldata ] = checkSelfIntersect( celldata)
%CHECKSELFINTERSECT Perform node swap in the case where internal elements self
%intersect

% for cellID = 1:celldata.nCells
%      Barycenters = getCellBarycenter( cellID, celldata.connec, celldata.r );
%      celldata.connec = rearrangeAnticlockwise( cellID,  celldata.connec,  celldata.r, Barycenters );
% end

