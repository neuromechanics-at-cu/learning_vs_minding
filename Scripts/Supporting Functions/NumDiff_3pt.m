function [dX] = NumDiff_3pt(X, t)

%==========================================================================
%==========================================================================
% [dX] = NumDiff_3pt(X, t) 
%
% This uses the three-point estimation (central difference) method of 
%  numerical differentiation to calculate the derivatives of a series of
%  points.  Forward and backward differences are used for the first and 
%  last points in the set.
%
%
% Author: Ben K. Bradley
% Date: 09/27/2013
%
%
% INPUTS:         Description                                         Units
%
%  X          - points to take the derivative of                    (n x m)
%  t          - time                                                (n x 1)
%
% OUTPUT:       
%    
%  dX         - derivative of X                                     (n x m)
%
%
%
% Coupling:
%
%  none
%
% References:
% 
%  [1] Kharab, A., An Introduction to Numerical Methods: A MATLAB 
%	     Approach, 2nd ed., 2005.
%
%  [2] Fornberg, B., "Generation of Finite Difference Formulas on 
%        Arbitrarily Spaced Grid," Mathematics of Computation, Vol. 51, 
%        No. 184, pp. 699-706, October, 1988.
%
%  [3] "Finite Difference Coefficient," Wikipedia, 2013.
%
%==========================================================================
%==========================================================================


% Use Forward Divided Difference for first data entry
dX(1,:) = (-1.5*X(1,:) + 2*X(2,:) - 0.5*X(3,:)) ./ (t(2,1) - t(1,1));


% Loop through all but last data entries
sz = length(X(:,1)); %number of rows

for m = 2:(sz-1)
    
    dX(m,:) = (X(m+1,:) - X(m-1,:)) ./ (2*(t(m+1,1) - t(m,1)));

end

% Use Backward Divided Difference for last data entry 
dX(sz,:) = (-1.5*X(sz,:) + 2*X(sz-1,:) - 0.5*X(sz-2,:)) ./ (t(sz,1) - t(sz-1,1));










