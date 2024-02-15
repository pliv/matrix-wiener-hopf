function [DISCR, DIFF, SOL] = benchmark(phi, gamma, n_iter, far_point, END, X)

time = cputime;
%%% Here we use the general algorithm ()
%% INPUT DATA
I = 1i;

a.func = @(x) (x.^2+10)./(x.^2+1);
a.array = a.func(X);
b.func = @(x) gamma./(x.^2+1);
b.array = b.func(X);

U_minus.analytic.func = @(x) 1./(x-2i);
U_minus.analytic.array = U_minus.analytic.func(X);

V_minus.analytic.func = @(x) 1./(x-3i);
V_minus.analytic.array = V_minus.analytic.func(X);

U_plus.analytic.func = @(x) 1./(x+2i);
U_plus.analytic.array = U_plus.analytic.func(X);

V_plus.analytic.func = @(x) 1./(x+3i);
V_plus.analytic.array = V_plus.analytic.func(X);


U_mm = U_minus.analytic.array;
U_pp = U_plus.analytic.array;
V_mm = V_minus.analytic.array;
V_pp = V_plus.analytic.array;


%% RHS


f1.func = @(x) U_minus.analytic.func(x) + a.func(x).*U_plus.analytic.func(x) + ...
    b.func(x).*exp(1i.*x.*phi).*V_plus.analytic.func(x);
f2.func = @(x) V_minus.analytic.func(x) + b.func(x).*exp(-1i.*x.*phi).*U_plus.analytic.func(x) +...
    a.func(x).*V_plus.analytic.func(x);

f1.array = f1.func(X);
f2.array = f2.func(X);

g1.func = @(x) f1.func(x);
g2.func = @(x) -b.func(x).*exp(-1i.*x.*phi).*f1.func(x)./a.func(x) + f2.func(x);

g1.array = g1.func(X);
g2.array = g2.func(X);

m_1_minus.func = @(x) (x-1i*sqrt(10)) ./ (x-1i);
m_1_plus.func = @(x) (x+1i*sqrt(10)) ./ (x+1i);
m_2_minus.func = @(x)  (x-1i*sqrt(10-gamma)).*(x-1i*sqrt(10+gamma))./ (x-1i) ./ (x-1i*sqrt(10));
m_2_plus.func = @(x)  (x+1i*sqrt(10-gamma)).*(x+1i*sqrt(10+gamma))./ (x+1i) ./ (x+1i*sqrt(10));

m_1_minus.array = m_1_minus.func(X);
m_1_plus.array = m_1_plus.func(X);
m_2_minus.array = m_2_minus.func(X);
m_2_plus.array = m_2_plus.func(X);

ix_0 = find(X==0);


%% C_1 C_2

C_1.func = @(x) g1.func(x) ./ m_1_minus.func(x);
C_2.func = @(x) g2.func(x) ./ m_2_minus.func(x);

C_1.array = C_1.func(X);
C_2.array = C_2.func(X);

ext1 = END/10;

END_big = END + ext1;
step = X(end)-X(end-1);

N_big = floor(length(X(abs(X)>END-ext1))/2);
X_p1 = linspace(END+step,END_big,N_big);
X_p2 = linspace(END_big+step, END_big+ext1, N_big);

X_m1 = -fliplr(X_p1);
X_m2 = -fliplr(X_p2);

X_p3 = [X(X>=0), X_p2];
X_m3 = [X_m2, X(X<=0)];

X_big = [X_m1, X, X_p1];

C_1_plus.func = @(k) ...
    9 ./ ((sqrt(10) + 1) .* (k + I)) ...
    - 6 ./ ((sqrt(10) + 2) .* (k + 2 .* I)) ...
    + gamma.*(exp(I .* k .* phi) ./ ((2 .* sqrt(10) + 2) .* (k + I)) ...
    - exp(I .* k .* phi) ./ ((2 .* sqrt(10) + 6) .* (k + 3 .* I)) ...
    - (exp(I .* k .* phi) - exp(-sqrt(10) .* phi)) ./ ...
      ((sqrt(10) + 3) .* (sqrt(10) + 1) .* (k - I .* sqrt(10))));


C_1_plus.array = C_1_plus.func(X);
    
C_1_minus.func = @(x) ...
    ((sqrt(10) - 1) ./ (x - I .* sqrt(10)) - 1 ./ (x - 2 .* I)) ./ (sqrt(10) - 2) ...
    - gamma.*exp(-sqrt(10) .* phi) ./ ((sqrt(10) + 3) .* (sqrt(10) + 1) .* (x - I .* sqrt(10)));


C_1_minus.array = C_1_minus.func(X);


g2_plus = @(k) ...
    gamma .* exp(-sqrt(10) .* phi) .* sqrt(10) ./ ((k + I .* sqrt(10)) .* (20 .* sqrt(10) + 40)) ...
    + (1 ./ 8) .* (gamma^2 - 1) ./ (k + 3 .* I) ...
    + (1 ./ 4) .* (-gamma^2 + 81) ./ ((k + I) .* 9) ...
    - (1 ./ 20) .* gamma^2 .* sqrt(10) ./ ((9 .* (sqrt(10) - 3)) .* (k + I .* sqrt(10)));

g2_minus = @(x) ...
    -gamma .* exp(-I .* x .* phi) .* (1 ./ (6 .* x - 12 .* I) - sqrt(10) ./ ((x - I .* sqrt(10)) .* (20 .* sqrt(10) - 40))) ...
    + gamma .* (exp(-I .* x .* phi) - exp(-sqrt(10) .* phi)) .* sqrt(10) ./ ((x + I .* sqrt(10)) .* (20 .* sqrt(10) + 40)) ...
    + 1 ./ (x - 3 .* I) ...
    + (1 ./ 8) .* (gamma^2 - 81) ./ (9 .* x - 9 .* I) ...
    - (1 ./ 20) .* gamma^2 .* sqrt(10) ./ ((9 .* sqrt(10) + 27) .* (x - I .* sqrt(10)));

k1 = 1i*sqrt(10-gamma);
k2 = 1i*sqrt(10+gamma);

mu_1 = (k1-1i)*(k1-1i*sqrt(10))/(k1-1i*sqrt(10+gamma));
mu_2 = (k2-1i)*(k2-1i*sqrt(10))/(k2-1i*sqrt(10-gamma));

g_s_plus = @(k) g2_plus(k)./m_2_minus.func(k) + g2_plus(k1).*mu_1./(k1-k) + g2_plus(k2).*mu_2./(k2-k);
g_s_minus = @(k) - g2_plus(k1).*mu_1./(k1-k) - g2_plus(k2).*mu_2./(k2-k);


C_2_plus.func = @(k) g_s_plus(k);
C_2_minus.func = @(k) g_s_minus(k) + g2_minus(k)./m_2_minus.func(k);

C_2_plus.array = C_2_plus.func(X);
C_2_minus.array = C_2_minus.func(X);


%% initial solution ()
U_m = m_1_minus.array.*(C_1_minus.array);
U_p = (C_1_plus.array)./m_1_plus.array;
V_m = m_2_minus.array.*(C_2_minus.array);
V_p = (C_2_plus.array)./m_2_plus.array;

%%  check

LHS_1 = U_m ./ m_1_minus.array + U_p .* m_1_plus.array;
RHS_1 = C_1.array - b.array.*exp(1i.*X.*phi).*V_p./m_1_minus.array;

LHS_2 = V_m ./ m_2_minus.array + V_p .* m_2_plus.array;
RHS_2 = C_2.array + b.array.*exp(-1i.*X.*phi).*U_m./m_2_minus.array./a.array;

discr_1_0 = abs(LHS_1 - RHS_1);
discr_2_0 = abs(LHS_2 - RHS_2);

%% ITERATIONS (start from step 2)

U_p_big.array = zeros(length(n_iter)+1,length(X));
U_m_big.array = U_p_big.array;
V_p_big.array = U_p_big.array;
V_m_big.array = U_p_big.array;

diff_U_m = U_p_big.array;
diff_U_p = U_p_big.array;
diff_V_m = U_p_big.array;
diff_V_p = U_p_big.array;

U_p_big.array(1,:) = U_p;
U_m_big.array(1,:) = U_m;
V_p_big.array(1,:) = V_p;
V_m_big.array(1,:) = V_m;

SOL(1, 1,:) = U_m;
SOL(2, 1,:) = U_p;
SOL(3, 1,:) = V_m;
SOL(4, 1,:) = V_p;


discr_1_big = zeros(length(n_iter)+1,length(X));
discr_2_big = zeros(length(n_iter)+1,length(X));

discr_1_big(1,:) = discr_1_0;
discr_2_big(1,:) = discr_2_0;

E_p = exp(1i.*X.*phi);
E_m = exp(-1i.*X.*phi);

ITERS = 0:n_iter;

diff_norm_U_m = 0.*ITERS;
diff_norm_U_p = 0.*ITERS;
diff_norm_V_m = 0.*ITERS;
diff_norm_V_p = 0.*ITERS;

diff_norm_U_m(1) = norm(abs(U_m-U_mm),2);
diff_norm_U_p(1) = norm(abs(U_p-U_pp),2);
diff_norm_V_m(1) = norm(abs(V_m-V_mm),2);
diff_norm_V_p(1) = norm(abs(V_p-V_pp),2);

discr_norm_1=zeros(1,length(n_iter+1));
discr_norm_2=discr_norm_1;

    discr_norm_1(1) = norm(discr_1_big(1,(X~=0)),2);
    discr_norm_2(1) = norm(discr_2_big(1,(X~=0)),2);

for m_step=1:n_iter
    
    D_1.array = b.array.*E_p.*V_p_big.array(m_step,:)./m_1_minus.array;
    D_2.array = b.array.*E_m.*U_m_big.array(m_step,:)./m_2_minus.array./a.array;
    
    ff1 = b.array.*V_p_big.array(m_step,:)./m_1_minus.array;
    ff2 = b.array.*U_m_big.array(m_step,:)./m_2_minus.array./a.array;
    
    ff1_re = least_squares_inf(real(ff1(X>=far_point)), X(X>=far_point),3);
    ff1_im = least_squares_inf(imag(ff1(X>=far_point)), X(X>=far_point),4);
    ff2_re = least_squares_inf(real(ff2(X>=far_point)), X(X>=far_point),3);
    ff2_im = least_squares_inf(imag(ff2(X>=far_point)), X(X>=far_point),4);
    
    D_1_inf = (ff1_re./(X+I).^3 + 0.*(3i*ff1_re + ff1_im)./(X+I).^4).*E_p;
    D_1.star = D_1.array - D_1_inf;
    
    D_2_inf = (ff2_re./(X-I).^3 + 0.*(-3i*ff2_re + ff2_im)./(X-I).^4).*E_m;
    D_2.star = D_2.array - D_2_inf;
    
    D1_XP = [D_1.star(X>=0), zeros(1,length(X_p2))];
    D1_XM = [zeros(1,length(X_m2)), D_1.star(X<=0)];
    
    D2_XP = [D_2.star(X>=0), zeros(1,length(X_p2))];
    D2_XM = [zeros(1,length(X_m2)), D_2.star(X<=0)];
    
    D1_big_p = interp1(X_p3,D1_XP,X_big(X_big>=0),'spline');
    D1_big_m = interp1(X_m3,D1_XM,X_big(X_big<=0),'spline');
    D1_big_pre = [D1_big_m, D1_big_p(2:end)];
    
    D2_big_p = interp1(X_p3,D2_XP,X_big(X_big>=0),'spline');
    D2_big_m = interp1(X_m3,D2_XM,X_big(X_big<=0),'spline');
    D2_big_pre = [D2_big_m, D2_big_p(2:end)];
    
    D_1.big = smooth_transition(END,END_big,X_big,D1_big_pre,0,0);
    D_2.big = smooth_transition(END,END_big,X_big,D2_big_pre,0,0);
    
    tic
    
    I_D1 = cauchy_integral(D_1.big,X_big,X,0);
    I_D2 = cauchy_integral(D_2.big,X_big,X,0);
    toc
    
    Q_LN_D1 = D_1.star.*log(abs((X_big(end)-X)./(X_big(end)+X)));
    Q_LN_D2 = D_2.star.*log(abs((X_big(end)-X)./(X_big(end)+X)));
    
    integral_D1 = I_D1 + Q_LN_D1;
    integral_D2 = I_D2 + Q_LN_D2;
    
    D_1_plus.array = 0.5.*D_1.star + integral_D1./2./pi./1i + D_1_inf;
    D_1_minus.array = 0.5.*D_1.star - integral_D1./2./pi./1i;
    D_2_plus.array = 0.5.*D_2.star + integral_D2./2./pi./1i;
    D_2_minus.array = 0.5.*D_2.star - integral_D2./2./pi./1i + D_2_inf;
    
    %%%%%%%
        %%%%% Solution
    U_m_big.array(m_step+1,:) = m_1_minus.array.*(C_1_minus.array - D_1_minus.array);
    U_p_big.array(m_step+1,:) = (C_1_plus.array - D_1_plus.array)./m_1_plus.array;
    V_m_big.array(m_step+1,:) = m_2_minus.array.*(C_2_minus.array + D_2_minus.array);
    V_p_big.array(m_step+1,:) = (C_2_plus.array + D_2_plus.array)./m_2_plus.array;
    
    U_p = U_p_big.array(m_step+1,:);
    U_m = U_m_big.array(m_step+1,:);
    V_p = V_p_big.array(m_step+1,:);
    V_m = V_m_big.array(m_step+1,:);
    
    SOL(1, m_step+1,:) = U_m;
    SOL(2, m_step+1,:) = U_p;
    SOL(3, m_step+1,:) = V_m;
    SOL(4, m_step+1,:) = V_p;
    
    LHS_1 = U_m ./ m_1_minus.array + U_p .* m_1_plus.array;
    RHS_1 = C_1.array - b.array.*exp(1i.*X.*phi).*V_p./m_1_minus.array;
    
    LHS_2 = V_m ./ m_2_minus.array + V_p .* m_2_plus.array;
    RHS_2 = C_2.array + b.array.*exp(-1i.*X.*phi).*U_m./m_2_minus.array./a.array;
    
    
    discr_1_big(m_step+1,:) = abs(LHS_1 - RHS_1);
    discr_2_big(m_step+1,:) = abs(LHS_2 - RHS_2);
    
    discr_norm_1(m_step+1) = norm(discr_1_big(m_step+1,(X~=0)),2);
    discr_norm_2(m_step+1) = norm(discr_2_big(m_step+1,(X~=0)),2);
    
    discr_norm_1(m_step+1)
    discr_norm_2(m_step+1)
    
    diff_U_m(m_step+1,:) = abs(U_mm - U_m);
    diff_U_p(m_step+1,:) = abs(U_pp - U_p);
    diff_V_m(m_step+1,:) = abs(V_mm - V_m);
    diff_V_p(m_step+1,:) = abs(V_pp - V_p);
    
    diff_norm_U_m(m_step+1) = norm(diff_U_m(m_step+1,(X~=0)),2);
    diff_norm_U_p(m_step+1) = norm(diff_U_p(m_step+1,(X~=0)),2);
    diff_norm_V_m(m_step+1) = norm(diff_V_m(m_step+1,(X~=0)),2);
    diff_norm_V_p(m_step+1) = norm(diff_V_p(m_step+1,(X~=0)),2);
    diff_norm_U_m(m_step+1)
    diff_norm_U_p(m_step+1)
    diff_norm_V_m(m_step+1)
    diff_norm_V_p(m_step+1)
    
    disp(['iteration ',num2str(m_step),' done'])
end


DISCR(1,:) = discr_norm_1;
DISCR(2,:) = discr_norm_2;
DIFF(1,:) = diff_norm_U_m;
DIFF(2,:) = diff_norm_U_p;
DIFF(3,:) = diff_norm_V_m;
DIFF(4,:) = diff_norm_V_p;

end


