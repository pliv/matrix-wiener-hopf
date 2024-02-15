
clear;
clc;

%% create plots for a given gamma
% note that before creating plots for this gamma, the results must be
% stored in results/gamma_{value_of_gamma} 
gamma = 1;
folderName = sprintf('results/gamma_%g', gamma); % This creates 'results/gamma_1' if gamma==1
fileName = 'results.mat';
filePath = fullfile(folderName, fileName);

if exist(filePath, 'file')
    data = load(filePath, 'params', 'struct');
    % Now 'data' is a structure with fields 'params' and 'struct'
    % You can access the variables as follows:
    params = data.params;
    structVar = data.struct; % Assuming 'struct' is the variable name; if it's a reserved keyword, use a different variable name.
else
    error('File does not exist: %s', filePath);
end


n_iter = 15;
ITERS = 0:n_iter;
J = 1;
discr_norm_1 = structVar(J).DISCR(1,:);
discr_norm_2 = structVar(J).DISCR(2,:);
diff_norm_U_m = structVar(J).DIFF(1,:);
diff_norm_U_p = structVar(J).DIFF(2,:);
diff_norm_V_m = structVar(J).DIFF(3,:);
diff_norm_V_p = structVar(J).DIFF(4,:);

DN_1_phi_001 = structVar(1).DISCR(1,:);
DN_1_phi_01 = structVar(2).DISCR(1,:);
DN_1_phi_1 = structVar(3).DISCR(1,:);


DN_2_phi_001 = structVar(1).DISCR(2,:);
DN_2_phi_01 = structVar(2).DISCR(2,:);
DN_2_phi_1 = structVar(3).DISCR(2,:);


DN_Um_phi_001 = structVar(1).DIFF(1,:);
DN_Um_phi_01 = structVar(2).DIFF(1,:);
DN_Um_phi_1 = structVar(3).DIFF(1,:);

DN_Up_phi_001 = structVar(1).DIFF(2,:);
DN_Up_phi_01 = structVar(2).DIFF(2,:);
DN_Up_phi_1 = structVar(3).DIFF(2,:);

DN_Vm_phi_001 = structVar(1).DIFF(3,:);
DN_Vm_phi_01 = structVar(2).DIFF(3,:);
DN_Vm_phi_1 = structVar(3).DIFF(3,:);

DN_Vp_phi_001 = structVar(1).DIFF(4,:);
DN_Vp_phi_01 = structVar(2).DIFF(4,:);
DN_Vp_phi_1 = structVar(3).DIFF(4,:);



FS1 = 16;
FS2 = 12;
color1 = 'r-*';
color2 = 'g-o';
color3 = 'b-v';
color4 = 'k-s';

x0 = 10;
y0 = 10;
width = 1000;
height = 400;

figure('Renderer', 'painters', 'Position', [x0 y0 width height])
ax1 = subplot(1,2,1);
semilogy(ax1,ITERS, DN_1_phi_001,color1,ITERS, DN_1_phi_01,color2,ITERS, DN_1_phi_1,color3)
xlabel('$n$','interpreter','latex','FontSize', FS1)
% title(['$\phi=$',num2str(phi(J))],'interpreter','latex')
ylabel('$$D_1$$','interpreter','LaTeX','FontSize',FS1, 'Rotation', 0)
%xticks(ITERS)
legend({'$\phi=0.01$','$\phi=0.1$','$\phi=1$','$\phi=10$'},'interpreter','latex', 'FontSize', FS2)
hold on
ax2 = subplot(1,2,2);
semilogy(ax2,ITERS, DN_2_phi_001,color1,ITERS, DN_2_phi_01,color2,ITERS, DN_2_phi_1,color3)
xlabel('$n$','interpreter','latex','FontSize', FS1)
legend({'$\phi=0.01$','$\phi=0.1$','$\phi=1$','$\phi=10$'},'interpreter','latex', 'FontSize', FS2)
ylabel('$$D_2$$','interpreter','LaTeX','FontSize',FS1, 'Rotation', 0)
%xticks(ITERS)


figure('Renderer', 'painters', 'Position', [x0 y0 width height])
ax1 = subplot(1,2,1);
semilogy(ax1,ITERS, DN_Um_phi_001,color1,ITERS, DN_Um_phi_01,color2,ITERS, DN_Um_phi_1,color3)
xlabel('$n$','interpreter','latex','FontSize', FS1)
ylabel('$$\Delta U^-$$','interpreter','LaTeX','FontSize',FS1, 'Rotation', 0)
%xticks(ITERS)
legend({'$\phi=0.01$','$\phi=0.1$','$\phi=1$','$\phi=10$'},'interpreter','latex', 'FontSize', FS2)
hold on
ax2 = subplot(1,2,2);
semilogy(ax2,ITERS, DN_Up_phi_001,color1,ITERS, DN_Up_phi_01,color2,ITERS, DN_Up_phi_1,color3)
xlabel('$n$','interpreter','latex','FontSize', FS1)
legend({'$\phi=0.01$','$\phi=0.1$','$\phi=1$','$\phi=10$'},'interpreter','latex', 'FontSize', FS2)
ylabel('$$\Delta U^+$$','interpreter','LaTeX','FontSize',FS1, 'Rotation', 0)
%xticks(ITERS)

figure('Renderer', 'painters', 'Position', [x0 y0 width height])
ax1 = subplot(1,2,1);
semilogy(ax1,ITERS, DN_Vm_phi_001,color1,ITERS, DN_Vm_phi_01,color2,ITERS, DN_Vm_phi_1,color3)
xlabel('$n$','interpreter','latex','FontSize', FS1)
ylabel('$$\Delta V^-$$','interpreter','LaTeX','FontSize',FS1, 'Rotation', 0)
%xticks(ITERS)
legend({'$\phi=0.01$','$\phi=0.1$','$\phi=1$','$\phi=10$'},'interpreter','latex', 'FontSize', FS2)
hold on
ax2 = subplot(1,2,2);
semilogy(ax2,ITERS, DN_Vp_phi_001,color1,ITERS, DN_Vp_phi_01,color2,ITERS, DN_Vp_phi_1,color3)
xlabel('$n$','interpreter','latex','FontSize', FS1)
legend({'$\phi=0.01$','$\phi=0.1$','$\phi=1$','$\phi=10$'},'interpreter','latex', 'FontSize', FS2)
ylabel('$$\Delta V^+$$','interpreter','LaTeX','FontSize',FS1, 'Rotation', 0)
%xticks(ITERS)

