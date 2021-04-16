function y = cdc(x,N)
% PURPOSE: adaptating boundary conditions to cope with checkboards limits.
% The checkboard becomes a torus.
% ------------------------------------------------------------
% SYNTAX:  y = cdc(x,N);
% ------------------------------------------------------------
% OUTPUT: y: 1x1 -> final position
% ------------------------------------------------------------
% INPUT: x: 1x1 -> initial position (row or column)
%        N: 1x1 -> checkboard dimension
% ------------------------------------------------------------

% written by:
%  Enrique M. Quilis
%  <equilis@gmail.com>

% Version 1.0 [August 2015]

switch x
    case 0 %Moving to the left: 0 becomes N
        y = N;
    case N+1 %Moving to the right: N+1 becomes 1
        y = 1; 
    otherwise %Preserving remaining positions
        y = x;
end
