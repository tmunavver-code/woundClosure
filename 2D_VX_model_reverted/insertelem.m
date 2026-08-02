function vect = insertelem(insertval, insertind, vect)
% Assuming vect is a 1D vector, the function adds insertval at the
% localtion insertind in the vect. If insertind is zero, then the element
% is added at the end of the array
if insertind == 0
	insertind = length(vect) + 1;
end

vect = cat(2,  vect(1:insertind - 1), insertval, vect(insertind:length(vect)));