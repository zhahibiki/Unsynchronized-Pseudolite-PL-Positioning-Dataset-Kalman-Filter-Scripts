function [E, N, U] = geo2utm(latitude, longitude, h, i)



utmZone = findUtmZone(latitude, longitude);

ecefP = geo2cartd(latitude, longitude, h, i);


[E, N, U] = cart2utm(ecefP(1), ecefP(2), ecefP(3), utmZone);





% ans =
%           39.8998834126787
% ans =
%           116.299939906796
% navSolutions.E(currMeasNr), ...
%          navSolutions.N(currMeasNr), ...
%          navSolutions.U(currMeasNr)
% ans =
%           440036.122761453
% ans =
%           4417112.13856293
% ans =
%           931.284742046148
% navSolutions.height(currMeasNr)
% ans =
%            1039.3546040887