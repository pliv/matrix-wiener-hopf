L = 0.5;
d = 0.5;

r_a = 0.1;
r_b = 0.2;
a = L + r_a;
b = L + r_b;

a0 = 2;
  % xi/k
% f_eps = @(y) -(exp(-1i.*y.*3)) ./ sqrt(eps-1i.*y) ./ sqrt(eps+1i.*y) ;
% g = @(y) -(exp(-1i.*y.*3)-1) ./ abs(y) ;
% g_eps = @(y) -(exp(-1i.*y.*3)-1) ./ sqrt(eps-1i.*y) ./ sqrt(eps+1i.*y);
EPS = 1e-15;

P = @(t) load_P(t,r_a,r_b);%1i.*(exp(-b.*1i.*t)-exp(-a.*1i.*t))./t;  %% Fourier transform of the applied load p_j(x-L)=g_j(x)
G = @(t) P(t)-P(0).*exp(1i.*L.*t)./(1-1i.*a0.*t);
q = @(t) -d.*G(t)./abs(t) + P(0).*exp(1i.*L.*t)./(1-1i.*a0.*t);
f = @(t) -P(t).*exp(-L.*1i.*t);

load('results05a','discrepancy','R_diff','Z_diff','Dif_R','Dif_Z','coef_b_r','coef_c_r','coef_b_z','coef_c_z','ZZ', 'RR', 'Zzero_p', 'END', 'xx')
xx_big = xx;
% END = 60; %%%endpoint of the integration interval
A = END;
h = 1;
% END_big = END+h;  %%%extension of the integration interval for calculation of Cauchy-type integral (avoiding the singularity on the edge)
END_big = END+h;

% END = 300; %%%endpoint of the integration interval
END_small = 60;
h = 1;
% END_big = END+h;  %%%extension of the integration interval for calculation of Cauchy-type integral (avoiding the singularity on the edge)
END_big = END+h;


N = length(xx)
i0 = find(abs(xx)<1e-13);

A = xx_big(end);
n_stop = 15;
N_END = 15;
n = length(discrepancy);
steps = 1:n;

xp_big = END_big+0.1:0.1:END_big+50;
xm_big = -END_big-50:0.1:-END_big-0.1;

xx_ext = [xm_big xx_big xp_big];
xx_ext_p = xx_ext(xx_ext>=0);
xx_ext_m = xx_ext(xx_ext<=0);

xx_big_p = xx_big(xx_big>=0);
xx_big_m = xx_big(xx_big<=0);
xx_p = xx(xx>=0);
xx_m = xx(xx<=0);

XP = [xx(xx>0) xp_big];
XM = [xm_big xx(xx<0)];

AA = END-30;
BB = END_big;
syms X1 X2 X3 X4;
[k1_p, k2_p, k3_p, k4_p] = solve(...
    AA^3*X1 + AA^2*X2 + AA*X3 + X4 == 1, ...
    END_big^3*X1 + END_big^2*X2 + END_big*X3 + X4 == 0,...
    3*AA^2*X1 + 2*AA*X2 + X3 == 0, ...
    3*END_big^2*X1 + 2*END_big*X2 + X3 == 0);
k1_p = double(k1_p);
k2_p = double(k2_p);
k3_p = double(k3_p);
k4_p = double(k4_p);

[k1_m, k2_m, k3_m, k4_m] = solve(...
    -AA^3*X1 + AA^2*X2 - AA*X3 + X4 == 1, ...
    -END_big^3*X1 + END_big^2*X2 - END_big*X3 + X4 == 0,...
    3*AA^2*X1 - 2*AA*X2 + X3 == 0, ...
    3*END_big^2*X1 - 2*END_big*X2 + X3 == 0);
k1_m = double(k1_m);
k2_m = double(k2_m);
k3_m = double(k3_m);
k4_m = double(k4_m);

fff_p = @(x) k1_p.*x.^3 + k2_p.*x.^2 + k3_p.*x + k4_p;
fff_m = @(x) k1_m.*x.^3 + k2_m.*x.^2 + k3_m.*x + k4_m;
% dff = @(x) 3.*k1_p.*x.^2 + 2.*k2_p.*x + k3_p;
PSI_p = @(x) (1-heaviside(x-AA)) + fff_p(x).*heaviside(x-AA).*(1-heaviside(x-END_big));
PSI_m = @(x) fff_m(x).*heaviside(x+END_big).*(1-heaviside(x+AA)) + heaviside(x+AA);

FS1 = 16;
% title(['R inf ',num2str(iter)])
figure
ax1 = subplot(2,2,[1 2]);
semilogy(ax1,steps(2:N_END), discrepancy(2:N_END),'b-*', 'LineWidth', 1)
grid on
% hold on
% plot(ax1,steps(2:N_END), discrepancy(2:N_END),'b-*')
xlabel('n','interpreter','LaTeX','FontSize',FS1)
title(['$$\left\Vert F-F_n\right\Vert$$, $L=$' num2str(L),', $d=$',num2str(d), ', $A=$',num2str(A-1)],'interpreter','LaTeX','FontSize',FS1)
hold on
ax2 = subplot(2,2,3);
semilogy(ax2,steps(2:N_END), R_diff(2:N_END), 'LineWidth', 1)
grid on
hold on
plot(ax2,steps(2:N_END), R_diff(2:N_END),'b*')
xlabel('n','interpreter','LaTeX','FontSize',FS1)
title('$$\left\Vert\tilde{R}^-_{n+1}-\tilde{R}^-_{n}\right\Vert$$','interpreter','LaTeX','FontSize',FS1)
hold on
ax3 = subplot(2,2,4);
semilogy(ax3,steps(2:N_END), Z_diff(2:N_END), 'LineWidth', 1)
grid on
hold on
plot(ax3,steps(2:N_END), Z_diff(2:N_END),'b*')
xlabel('n','interpreter','LaTeX','FontSize',FS1)
title('$$\left\Vert\tilde{Z}^+_{n+1}-\tilde{Z}^+_{n}\right\Vert$$','interpreter','LaTeX','FontSize',FS1)
grid on
hold off


% figure
% bx1 = subplot(1,2,1);
% for i=2:4
%     txt = ['iteration = ',num2str(i)];
%     semilogy(bx1,xx, abs(squeeze(Dif_R(i,:))), '--','DisplayName',txt)
%     hold on
% end
% title('abs(R(n+1)-R(n))')
% xlabel('t')
% hold on
% legend show
% 
% bx2 = subplot(1,2,2);
% for i=2:4
%     txt = ['iteration = ',num2str(i)];
%     semilogy(bx2,xx, abs(squeeze(Dif_Z(i,:))), '--','DisplayName',txt)
%     hold on
% end
% title('abs(Z(n+1)-Z(n))')
% xlabel('t')
% hold off
% legend show
END_ext = 100;
Xm = -END_big-END_ext:0.1:-END_big-0.1;
Xp = END_big+0.1:0.1:END_big+END_ext;
X = [Xm xx_big Xp];
figure
cx1 = subplot(1,2,1);
for i=1:5
    d_b_r = coef_b_r(i+1)-coef_b_r(i);
    d_c_r = coef_c_r(i+1)-coef_c_r(i);
    f_R = @(x) d_b_r.*exp(-1i.*x.*L)./sqrt(1-1i.*x)./(1-1i.*x.*a0) + d_c_r./(1-1i.*x.*a0);
%     diff_R = [f_R(Xm), squeeze(Dif_R(i+1,:)), f_R(Xp)];
    
    diff_R_p = squeeze(Dif_R(i+1,i0:end));
    diff_R_m = squeeze(Dif_R(i+1,1:i0));
    
    d_R_XP = [diff_R_p(2:end) zeros(1,length(xp_big))];
    d_R_XM = [zeros(1,length(xm_big)) diff_R_m(1:end-1)];
    
    d_R_p = interp1(XP,d_R_XP,xx_big(xx_big>0),'spline');
    d_R_p = [diff_R_p(1) d_R_p];
    d_R_m = interp1(XM,d_R_XM,xx_big(xx_big<0),'spline');
    d_R_m = [d_R_m diff_R_p(end)];
    
    d_R_p_psi = d_R_p.*PSI_p(xx_big_p) + f_R(xx_big_p).*(1-PSI_p(xx_big_p));
    d_R_m_psi = d_R_m.*PSI_m(xx_big_m) + f_R(xx_big_m).*(1-PSI_m(xx_big_m));
    
    diff_R = [d_R_m_psi d_R_p_psi(2:end)];
    diff_R = [f_R(Xm) diff_R f_R(Xp)];
    
    txt = ['n=',num2str(i+1)];
    semilogy(cx1,X, abs(diff_R),'--','LineWidth',1,'DisplayName',txt)
    hold on
end
legend show
line(cx1,[-A -A],ylim,'LineStyle','--','Color','k')
line(cx1,[A A],ylim,'LineStyle','--','Color','k')
title('$$ \left| \tilde{R}^-_{n+1}-\tilde{R}^-_{n} \right|  $$','interpreter','latex','FontSize',FS1)
xlabel('t','interpreter','LaTeX','FontSize',FS1)
hold on


cx2 = subplot(1,2,2);
for i=1:5
    d_b_z = coef_b_z(i+1)-coef_b_z(i);
    d_c_z = coef_c_z(i+1)-coef_c_z(i);
    f_Z = @(x) d_b_z.*exp(1i.*x.*L)./(1-1i.*x.*a0) + d_c_z./sqrt(1-1i.*x)./(1-1i.*x.*a0);
    %     diff_Z = [f_Z(Xm), squeeze(Dif_Z(i+1,:)), f_Z(Xp)];
    
    diff_Z_p = squeeze(Dif_Z(i+1,i0:end));
    diff_Z_m = squeeze(Dif_Z(i+1,1:i0));
    
    d_Z_XP = [diff_Z_p(2:end) zeros(1,length(xp_big))];
    d_Z_XM = [zeros(1,length(xm_big)) diff_Z_m(1:end-1)];
    
    d_Z_p = interp1(XP,d_Z_XP,xx_big(xx_big>0),'spline');
    d_Z_p = [diff_Z_p(1) d_Z_p];
    d_Z_m = interp1(XM,d_Z_XM,xx_big(xx_big<0),'spline');
    d_Z_m = [d_Z_m diff_Z_p(end)];
    
    d_Z_p_psi = d_Z_p.*PSI_p(xx_big_p) + f_Z(xx_big_p).*(1-PSI_p(xx_big_p));
    d_Z_m_psi = d_Z_m.*PSI_m(xx_big_m) + f_Z(xx_big_m).*(1-PSI_m(xx_big_m));
    
    diff_Z = [d_Z_m_psi d_Z_p_psi(2:end)];
    diff_Z = [f_Z(Xm) diff_Z f_Z(Xp)];

    txt = ['n=',num2str(i+1)];
    semilogy(cx2,X, abs(diff_Z),'--','LineWidth',1,'DisplayName',txt)
    hold on
end
legend show
line(cx2,[-A -A],ylim,'LineStyle','--','Color','k')
line(cx2,[A A],ylim,'LineStyle','--','Color','k')
title('$$\left| \tilde{Z}^+_{n+1}-\tilde{Z}^+_{n}\right|$$','interpreter','latex','FontSize',FS1)
xlabel('t','interpreter','LaTeX','FontSize',FS1)
% set(get(gca,'ylabel'),'rotation',0)
hold off



norm_RR = zeros(1,length(steps));
for i=1:n-1
    d_b_r = coef_b_r(i+1)-coef_b_r(i);
    d_c_r = coef_c_r(i+1)-coef_c_r(i);
    f_R = @(x) d_b_r.*exp(-1i.*x.*L)./sqrt(1-1i.*x)./(1-1i.*x.*a0) + d_c_r./(1-1i.*x.*a0);
%     diff_R = [f_R(Xm), squeeze(Dif_R(i+1,:)), f_R(Xp)];
    
    diff_R_p = squeeze(Dif_R(i+1,i0:end));
    diff_R_m = squeeze(Dif_R(i+1,1:i0));
    
    d_R_XP = [diff_R_p(2:end) zeros(1,length(xp_big))];
    d_R_XM = [zeros(1,length(xm_big)) diff_R_m(1:end-1)];
    
    d_R_p = interp1(XP,d_R_XP,xx_big(xx_big>0),'spline');
    d_R_p = [diff_R_p(1) d_R_p];
    d_R_m = interp1(XM,d_R_XM,xx_big(xx_big<0),'spline');
    d_R_m = [d_R_m diff_R_p(end)];
    
    d_R_p_psi = d_R_p.*PSI_p(xx_big_p) + f_R(xx_big_p).*(1-PSI_p(xx_big_p));
    d_R_m_psi = d_R_m.*PSI_m(xx_big_m) + f_R(xx_big_m).*(1-PSI_m(xx_big_m));
    
    RR = [d_R_m_psi d_R_p_psi(2:end)];
%     RR = [f_R(Xm) RR f_R(Xp)];
    
    norm_RR(i+1) = norm(RR,2);
end


norm_ZZ = zeros(1,length(steps));
for i=1:n-1
   d_b_z = coef_b_z(i+1)-coef_b_z(i);
    d_c_z = coef_c_z(i+1)-coef_c_z(i);
    f_Z = @(x) d_b_z.*exp(1i.*x.*L)./(1-1i.*x.*a0) + d_c_z./sqrt(1-1i.*x)./(1-1i.*x.*a0);
    %     diff_Z = [f_Z(Xm), squeeze(Dif_Z(i+1,:)), f_Z(Xp)];
    
    diff_Z_p = squeeze(Dif_Z(i+1,i0:end));
    diff_Z_m = squeeze(Dif_Z(i+1,1:i0));
    
    d_Z_XP = [diff_Z_p(2:end) zeros(1,length(xp_big))];
    d_Z_XM = [zeros(1,length(xm_big)) diff_Z_m(1:end-1)];
    
    d_Z_p = interp1(XP,d_Z_XP,xx_big(xx_big>0),'spline');
    d_Z_p = [diff_Z_p(1) d_Z_p];
    d_Z_m = interp1(XM,d_Z_XM,xx_big(xx_big<0),'spline');
    d_Z_m = [d_Z_m diff_Z_p(end)];
    
    d_Z_p_psi = d_Z_p.*PSI_p(xx_big_p) + f_Z(xx_big_p).*(1-PSI_p(xx_big_p));
    d_Z_m_psi = d_Z_m.*PSI_m(xx_big_m) + f_Z(xx_big_m).*(1-PSI_m(xx_big_m));
    
    ZZ = [d_Z_m_psi d_Z_p_psi(2:end)];
%     RR = [f_R(Xm) RR f_R(Xp)];
    
    norm_ZZ(i+1) = norm(ZZ,2);
end


figure
ax1 = subplot(2,2,[1 2]);
semilogy(ax1,steps(2:N_END), discrepancy(2:N_END), 'b-v', 'LineWidth', 1)
grid on
xlabel('n','interpreter','LaTeX','FontSize',FS1)
title('$$\left\Vert F-F_n\right\Vert$$','interpreter','LaTeX','FontSize',FS1)
hold on
ax2 = subplot(2,2,3);
semilogy(ax2,steps(2:N_END), norm_RR(2:N_END),'b-v', 'LineWidth', 1)
grid on
xlabel('n','interpreter','LaTeX','FontSize',FS1)
title('$$\left\Vert\tilde{R}^-_{n+1}-\tilde{R}^-_{n}\right\Vert$$','interpreter','LaTeX','FontSize',FS1)
hold on
ax3 = subplot(2,2,4);
semilogy(ax3,steps(2:N_END), norm_ZZ(2:N_END), 'b-v', 'LineWidth', 1)
grid on
xlabel('n','interpreter','LaTeX','FontSize',FS1)
title('$$\left\Vert\tilde{Z}^+_{n+1}-\tilde{Z}^+_{n}\right\Vert$$','interpreter','LaTeX','FontSize',FS1)
grid on
hold off

function PP = load_P(tt,aa,bb)
PP = zeros(1,length(tt));
for j=1:length(tt)
    if abs(tt(j))<1e-14
        PP(j) = bb-aa;
    else
        PP(j) = 1i.*(exp(-bb.*1i.*tt(j))-exp(-aa.*1i.*tt(j)))./tt(j);
    end
end
end

function ll = log_p(z)
        ll = log(abs(z)) + 1i.*arg_p(z);
    end

    function tt = arg_p(z)
        TT = angle(z);
        tt = (TT>=-pi/2).*(TT<=pi).*TT + (TT>=-pi).*(TT<-pi/2).*(TT+2.*pi);
    end
    
        function ll = log_m(z)
        ll = log(abs(z)) + 1i.*arg_m(z);
    end

    function tt = arg_m(z)
        TT = angle(z);
        tt = (TT>=-pi).*(TT<=pi/2).*TT + (TT>pi/2).*(TT<=pi).*(TT-2.*pi);
    end
