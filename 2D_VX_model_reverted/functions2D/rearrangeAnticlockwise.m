function [ Connectivity ] = rearrangeAnticlockwise( cellID, Connectivity, AllverticesCoords, Barycenter )
%REARRANGEANTICLOCKWISE Rearrange the cell vertices in an anticlockwise order

cellvertices = Connectivity{cellID};
angles = zeros(size(cellvertices));
for j = 1:length(cellvertices)
    dx = AllverticesCoords(cellvertices(j),1) - Barycenter(1);
    dy = AllverticesCoords(cellvertices(j),2) - Barycenter(2);
    angles(j) = atan2(dy, dx)/pi;
end
[~,sortind] = sort(angles);
Connectivity{cellID} = cellvertices(sortind);

end



