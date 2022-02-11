function G = fitgain(f0,v0)

% For given data vectors v(t) and f(t) of equal length, this function fits
% an estimated gain 'g' such that f = g*v. (For curl field catch trials.)

% fh     Fit function, f(t)
% errfh  Error function (want to minimize this; can also use norm, etc.)
% g0     Initial guess
% g      Final estimated gain

options = optimset('MaxFunEvals',1e6, 'MaxIter',1e6);

n = length(f0);
f = zeros(1,n); f(1:n) = f0(1:n);
v = zeros(1,n); v(1:n) = v0(1:n);

g0 = max(f)/max(v);

fh = @(g) ( v.*g );
errfh = @(g,dat) sum( (dat - fh(g)).^2 );
G = fminsearch(errfh,g0,options,f);

