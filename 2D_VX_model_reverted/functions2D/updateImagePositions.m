function [celldata, augT1flagVec,augT1relaxstepcountVec]   = updateImagePositions(celldata,augT1flagVec,augT1relaxstepcountVec,param)
%% updateImagePositions: update image vertex positions
% The images displace as their masters

celldata.slave_r = celldata.slave_r0 + celldata.r(celldata.StoMmap,:) - celldata.r0(celldata.StoMmap,:);
celldata.allVertices = [celldata.r;celldata.slave_r];
