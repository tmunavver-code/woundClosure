function [coordinates, connectivity, edgedata,verttocell, T1flagVec,T1relaxstepcountVec, nT1]   = checkT1transitions(nCells,coordinates, connectivity, edgedata,verttocell,param,T1flagVec,T1relaxstepcountVec,imAugmented)
%% checkT1transitions: Checks for and performs a T1 transition
nT1 = param.nT1;

if imAugmented == 1
    label = 'augmented';
else
    label = 'reduced';
end

for cellID = 1:nCells
    verticesmat = connectivity{cellID};
    edgelenmat  = edgedata{cellID};
    for j = 1:length(verticesmat)
        if edgelenmat(j) < param.T1_TOL && T1flagVec(j) == 0
            vertID1 = verticesmat(j);
            if j < length(verticesmat)
                vertID2 = verticesmat(j+1);
            else % Final edge
                vertID2 = verticesmat(1);
            end
            % Find the surrounding cells to the edge to be flipped
            [cellID1, cellID2, cellID3, cellID4] = getT1cellIDs(vertID1, vertID2,  verttocell,connectivity);
            % Perform the T1 transition
            [coordinates, connectivity, edgedata,verttocell] = performT1swap(cellID1, cellID2, cellID3, cellID4, coordinates, connectivity, edgedata,verttocell, vertID1, vertID2, param,imAugmented);
            % Display status
            fprintf(1,'-------------\nT1 transition performed in %s \n',label);
            fprintf(1,'V1 = %3d, V2 = %3d\n',vertID1, vertID2);
            fprintf(1,'C1 = %3d, C2 = %3d, C3 = %3d, C4 = %3d]\n',cellID1, cellID2, cellID3, cellID4);
            fprintf(1,'-------------\n');
            nT1 = nT1 + 1;
            T1flagVec(j) = 1;
        else
            % Do nothing
        end
    end
end


for i = 1:length(T1flagVec)
    if T1flagVec(i) == 1
        T1relaxstepcountVec(i) = T1relaxstepcountVec(i) + 1;
    end
    
    if T1relaxstepcountVec(i) > param.maxT1relaxsteps
        T1flagVec(i) = 0;
        T1relaxstepcountVec(i) = 0;
    end
end