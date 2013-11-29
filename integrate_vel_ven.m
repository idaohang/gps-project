% integrate_vel_ven.m	(actual file name: integrate_vel_ven.m)
%
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% this function performs trapezoidal integration on the velocity
% time history in vertical/east/north coordinates to produce
% the accumulated vertical/east/north position change time
% history.
%
% input:  'thist'        Nvel-by-1 vector of times (sec) at which the
%                        vertical-east-north velocity vector components
%                        are available in the corresponding rows
%                        of velhist_VEN
%
%						   [ t; t; t; ...; t]
%
%		  'velhist_VEN'  Nvel-by-3 matrix of velocities (meters/sec)
%                        in vertical/east/north coordinates, each
%                        row is a VEN velocity vector at the
%                        time in the corresponding row of thist.
%
%						   [ velVertical velEast velNorth; ...
%                            velVertical velEast velNorth; ...
%                                .
%                                .
%                                .
%                            velVertical velEast velNorth]
%
% output: 'poshist_VEN'  Nvel-by-3 matrix of positions relative
%                        to the starting point (meters) in
%                        vertical/east/north coordinates, each
%                        row is a VEN position vector at the
%                        time in the corresponding row of thist.
%
%						   [ posVertical posEast posNorth; ...
%                            posVertical posEast posNorth; ...
%                                .
%                                .
%                                .
%                            posVertical posEast posNorth]
%
%                        Note: because these positions are measured
%                        relative to the start, the first row
%                        of poshist_VEN contains all zeros.
%          
%
function poshist_VEN = integrate_vel_ven(thist,velhist_VEN)
% determine the number of time points and one less than it.
Nvel = size(thist,1);
Nvelm1 = Nvel - 1;
% initialize the output with all zeros.
poshist_VEN = zeros(Nvel,3);
% compute the time intervals.
diff_thist = diff(thist);
% do the Vertical component trapezoidal integration
% followed by the East component and concluding with the North
% component, one per iteration of the loop.
for jj = 1:3
   ddelposhist = diff_thist.*(velhist_VEN(1:Nvelm1,jj) + ...
                              velhist_VEN(2:Nvel,jj))*0.5;
   poshist_VEN(2:Nvel,jj) = cumsum(ddelposhist);
end
return