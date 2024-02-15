function I = cauchy_integral(F, X, K,im)
% This function calculates Cauchy type integral I(k) = int(F(x)/(x-k), x=-A..+A)
% The integration is performed along the line segment [-A,+A]. X is an
% array, corresponding to this segment X = [-A,+A]
% The values of k are taken from array K.
% parameter "im" defines the imaginary part of z
% It is assumed that function F decays fast enough to be approximated with
% 0 close the ends of the integration inerval
%   Detailed explanation goes here
% 


A = X(end);

    nX = length(X);
    nK = length(K);

if abs(im)>1e-13
    F2 = [fliplr(F) F(2:end)];
    spline_F = spline(X, F2);
    func_F = @(x) ppval(spline_F,x);
    
    I = zeros(1,nK);
    
    % tic
    for j=1:nK
        I(j) = quadgk(@(x) func_F(x)./(x-(K(j)+im)),X(1),X(end),'Waypoints',X);
    end
else


    spline_F = spline(X, F);
    % func_F_approx = @(x) ppval(spline_F_approx,x);
   
    integrand = zeros(nX,nK);

    for j=1:nK
        ii = find(abs(X-K(j))<1e-13);
%         j
%         K(j)
%         ii
        integrand(:,j) = (F - F(ii))./(X-K(j));
        integrand(ii,j) = spline_F.coefs(ii,3);
    end
    
    I = zeros(1,nK);
    
    % tic
    for j=1:nK
        spline_integrand = spline(X, squeeze(integrand(:,j)));
        func_integrand = @(x) ppval(spline_integrand, x);
        I(j) = quadgk(@(x) func_integrand(x),X(1),X(end),'Waypoints',X);
    end
    
end
end

