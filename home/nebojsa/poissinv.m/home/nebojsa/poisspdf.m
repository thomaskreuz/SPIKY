function x = poissinv(p,lambda)
%POISSINV Inverse of the Poisson cumulative distribution function (cdf).
%	X = POISSINV(P,LAMBDA) returns the inverse of the Poisson cdf 
%	with parameter lambda. Since the Poisson distribution is discrete,
%	POISSINV returns the smallest value of X, such that the poisson 
%	cdf evaluated, at X, equals or exceeds P.
%
%	The size of X is the common size of P and LAMBDA. A scalar input   
%	functions as a constant matrix of the same size as the other input.	 

%	B.A. Jones 1-15-93
%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:55:53 $

if nargin < 2, 
    error('Requires two input arguments.'); 
end

[errorcode p lambda] = distchck(2,p,lambda);

if errorcode > 0
    error('The arguments must be the same size or be scalars.');
end

x = zeros(size(p));

cumdist = poisspdf(0,lambda);
count = 0;

% Compare P to the poisson cdf.
k = find(lambda > 0 & p >= 0 & p <= 1);
while any(any(p(k) > cumdist(k)))
    count = count + 1;
    cumdist(k) = cumdist(k) + poisspdf(count,lambda(k));
    if count == 1
            x(k) = x(k) + 1;
    end
    index = find(cumdist(k) < p(k));
    x(k(index)) = x(k(index)) + 1;  
end

% Return NaN if the arguments are outside their respective limits.
k = find(lambda <= 0 | p < 0 | p > 1);
if any(k)
    x(k) = NaN * ones(size(k));
end
