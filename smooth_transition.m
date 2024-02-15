function F = smooth_transition(AA,BB,X,f,asympt_p,asympt_m)

syms X1 X2 X3 X4;
[k1_p, k2_p, k3_p, k4_p] = solve(...
    AA^3*X1 + AA^2*X2 + AA*X3 + X4 == 1, ...
    BB^3*X1 + BB^2*X2 + BB*X3 + X4 == 0,...
    3*AA^2*X1 + 2*AA*X2 + X3 == 0, ...
    3*BB^2*X1 + 2*BB*X2 + X3 == 0);
k1_p = double(k1_p);
k2_p = double(k2_p);
k3_p = double(k3_p);
k4_p = double(k4_p);

[k1_m, k2_m, k3_m, k4_m] = solve(...
    -AA^3*X1 + AA^2*X2 - AA*X3 + X4 == 1, ...
    -BB^3*X1 + BB^2*X2 - BB*X3 + X4 == 0,...
    3*AA^2*X1 - 2*AA*X2 + X3 == 0, ...
    3*BB^2*X1 - 2*BB*X2 + X3 == 0);
k1_m = double(k1_m);
k2_m = double(k2_m);
k3_m = double(k3_m);
k4_m = double(k4_m);

fff_p = @(x) k1_p.*x.^3 + k2_p.*x.^2 + k3_p.*x + k4_p;
fff_m = @(x) k1_m.*x.^3 + k2_m.*x.^2 + k3_m.*x + k4_m;

PSI_p = @(x) (1-heaviside(x-AA)) + fff_p(x).*heaviside(x-AA).*(1-heaviside(x-BB));
PSI_m = @(x) fff_m(x).*heaviside(x+BB).*(1-heaviside(x+AA)) + heaviside(x+AA);

f_p = f(X>=0).*PSI_p(X(X>=0)) + asympt_p.*(1-PSI_p(X(X>=0)));
f_m = f(X<=0).*PSI_m(X(X<=0)) + asympt_m.*(1-PSI_m(X(X<=0)));

F = [f_m f_p(2:end)];
end