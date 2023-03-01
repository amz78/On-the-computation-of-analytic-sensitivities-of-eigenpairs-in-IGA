clear;
close all;
addpath(genpath('')) %add path to the GeoPDEs installation
warning('off', 'geopdes:automaticGeneration')
warning('off', 'nrbderiv:SecondDerivative')

global DERC
global DERA
global NI
global NMAX
NMAX = 1; %highest order of the matrix derivatives to compute (in this repository 
% only the first derivative of the eigenpair can be computed, so we set NMAX = 1)
NI = 1; %running index of the current derivative, start with NI = 1

try 
    symbolicToolboxInstalled = any(any(contains(struct2cell(ver), 'Symbolic Math Toolbox')));
catch
    symbolicToolboxInstalled = false;
end
if symbolicToolboxInstalled
    DERC = der_Ct_At(NMAX, 'C');
    DERA = der_Ct_At(NMAX, 'A');
else
    warning('Symbolic Math Toolbox not found')
    load('der_Ct_At_funHandles.mat') %file with the first 10 function handles
    %to compute the derivatives of C(t) and A(t)
end

deg = 2;
sub = 2;
dispDetailed = false;

%% code for Sec. 4.1. Uncertainty quantification for a pillbox cavity of uncertain radius
% ==================================================================
fprintf('\nExample in Sec. 4.1. Uncertainty quantification for a pillbox cavity of uncertain radius: \n \n');

%Load the cavities for r = a, r = b, r = r_d and r = r_d + 1e-6 from files
[sys_ra.A, sys_ra.B, int_dofs, drchlt_dofs, geo_cavity_ra, sp_cavity_ra, msh_cavity_ra] = getMatrices_fromFile('pillbox_r_a.txt', deg, sub, dispDetailed);
[sys_rb.A, sys_rb.B, ~, ~, geo_cavity_rb, sp_cavity_rb, msh_cavity_rb] = getMatrices_fromFile('pillbox_r_b.txt', deg, sub, dispDetailed);
[sys_rd.A, sys_rd.B, ~, ~, geo_cavity_rd] = getMatrices_fromFile('pillbox_r_d.txt', deg, sub, dispDetailed);
[sys_rd_FD.A, sys_rd_FD.B] = getMatrices_fromFile('pillbox_r_d_FD.txt', deg, sub, dispDetailed);

%compare the derivatives of the system matrices with finite differences
dA_FD = (sys_rd_FD.A - sys_rd.A) / 1e-6;
dB_FD = (sys_rd_FD.B - sys_rd.B) / 1e-6;

[sys_rd.dA, sys_rd.dB] = getDerivativeMatrices(sp_cavity_ra, msh_cavity_ra, msh_cavity_rb, int_dofs, 0.5);

disp('Compare the shape derivatives of the system matrices to finite differences:')
dA_error = norm(sys_rd.dA{1}-dA_FD, 'fro') / norm(dA_FD, 'fro');
fprintf(' Relative error of the derivatives of the stiffness matrices measured in the frobenius norm is %1.2e \n', dA_error)
dB_error = norm(sys_rd.dB{1}-dB_FD, 'fro') / norm(dB_FD, 'fro');
fprintf(' Relative error of the derivatives of the mass matrices measured in the frobenius norm is %1.2e \n', dB_error)

%plot the cavity for r_d = 0.5m
figure(1)
for iptc = 1:numel(geo_cavity_rd)
    nrbplot(geo_cavity_rd(iptc).nurbs, [5, 5, 5])
    hold on
end
title('Pillbox cavity for r_d = 0.5m.')

%Compute the analytic and numerical results
G = 2.40482555769577; %first root of the Bessel function J_0(x)
c = 299792458;
f0 = @(x) ((G * c) ./ (2 * pi * x));
lambda0 = @(x) ((G^2) ./ (x.^2));

[eig_val, eig_vec] = solveMyEVP(sys_rd, 10, (2 * pi * f0(0.5))^2*1e-20);

mode_pb = 1; %consider the first eigenmode

%compute the derivatives of the eigenvector and the eigenvalue using (10)
dep = [sys_rd.A - eig_val(mode_pb) * sys_rd.B, -sys_rd.B * eig_vec(:, mode_pb); ...
    eig_vec(:, mode_pb)' * sys_rd.B, 0] \ [-sys_rd.dA{1} * eig_vec(:, mode_pb) + ...
    eig_val(mode_pb) * sys_rd.dB{1} * eig_vec(:, mode_pb); -eig_vec(:, mode_pb)' * sys_rd.dB{1} * eig_vec(:, mode_pb)];

%% Generate results as in Figure 3
a = 0.2;
b = 0.8;
radii = linspace(a, b, 50);

figure(2)
%plot the analytic result
plot(radii, lambda0(radii), 'k-')
hold on
taxis = linspace(-0.5, 0.5, 3);
%plot the Taylor series expansions
plot((taxis + 0.5)*(b - a)+a, (polyval([eig_val(mode_pb)], taxis))*1e20/c^2, 'color', [27, 158, 119]/256)
plot((taxis + 0.5)*(b - a)+a, (polyval([dep(end), eig_val(mode_pb)], taxis))*1e20/c^2, 'color', [217, 95, 2]/256)
ylim([-25, 160])
xlim([0.1, 0.9])
legend('Analytic', 'Order 0', 'Order 1', 'Location', 'northeast')
grid on
ylabel('\lambda')
xlabel('Radius r (m)')
title({'Eigenfrequencies of the pillbox cavity', 'estimated using the Taylor series expansion at rd = 0.5 m.'})

%% code for Sec. 4.2. Eigenvalue estimation for shape morphing
% ==================================================================
fprintf('\nExample in Sec. 4.2. Eigenvalue estimation for shape morphing: \n \n');

%Load the cavities from files
[sys_pb.A, sys_pb.B, int_dofs_pb, drchlet_dofs_pb, geometry_pb, space_pb, mesh_pb] = getMatrices_fromFile('shapeMorphing_pillbox.txt', deg, sub, dispDetailed);
[sys_tesla.A, sys_tesla.B, int_dofs_tesla, drchlet_dofs_tesla, geometry_tesla, space_tesla, mesh_tesla] = getMatrices_fromFile('shapeMorphing_tesla.txt', deg, sub, dispDetailed);
[sys_tesla_FD.A, sys_tesla_FD.B] = getMatrices_fromFile('shapeMorphing_tesla_FD.txt', deg, sub, dispDetailed);

[sys_tesla.dA, sys_tesla.dB] = getDerivativeMatrices(space_tesla, mesh_tesla, mesh_pb, int_dofs_tesla, 0);

%% plot the pillbox and the TESLA cavity
figure(3)
for iptc = 1:numel(geometry_pb)
    nrbplot(geometry_pb(iptc).nurbs, [5, 5, 5])
    hold on
end
title('Pillbox cavity for the shape morphing.')

figure(4)
for iptc = 1:numel(geometry_tesla)
    nrbplot(geometry_tesla(iptc).nurbs, [5, 5, 20])
    hold on
end
title('9-cell TESLA cavity for the shape morphing.')

%% compare the matrix derivatives to finite differences
tesla_dA = (sys_tesla_FD.A - sys_tesla.A) / 1e-6;
tesla_dB = (sys_tesla_FD.B - sys_tesla.B) / 1e-6;

disp('Compare the shape derivatives of the system matrices to finite differences:')
dA_error_tesla = norm(sys_tesla.dA{1}-tesla_dA, 'fro') / norm(tesla_dA, 'fro');
fprintf(' Relative error of the derivatives of the stiffness matrices measured in the frobenius norm is %1.2e \n', dA_error_tesla);

dB_error_tesla = norm(tesla_dB-sys_tesla.dB{1}, 'fro') / norm(tesla_dB, 'fro');
fprintf(' Relative error of the derivatives of the mass matrices measured in the frobenius norm is %1.2e \n', dB_error_tesla);

%%
[eig_val_tesla, eig_vec_tesla] = solveMyEVP(sys_tesla, 100, (2 * pi * 1.3e9)^2*1e-20);

mode_tesla = 9; %consider the accelerating mode of the 9-cell TESLA cavity

%compute the derivatives of the eigenvector and the eigenvalue using (10)
dep_pb = [sys_tesla.A - eig_val_tesla(mode_tesla) * sys_tesla.B, -sys_tesla.B * eig_vec_tesla(:, mode_tesla); ...
    eig_vec_tesla(:, mode_tesla)' * sys_tesla.B, 0] \ [-sys_tesla.dA{1} * eig_vec_tesla(:, mode_tesla) + ...
    eig_val_tesla(mode_tesla) * sys_tesla.dB{1} * eig_vec_tesla(:, mode_tesla); -eig_vec_tesla(:, mode_tesla)' * sys_tesla.dB{1} * eig_vec_tesla(:, mode_tesla)];

%% Generate results as in Figure 5
figure(5)
% load the results of the multistep method for comparison
Morphing = csvread('Morphing_Mode9_physic.csv',1,0);
tMorphing = Morphing(:, 1);
fMorphing = Morphing(:, 2);
plot(tMorphing, fMorphing, 'k-x');
hold on
taxis_tesla = linspace(0, 1, 20);
%plot the Taylor series expansions
plot(taxis_tesla, sqrt(abs((polyval([eig_val_tesla(mode_tesla)], taxis_tesla))*1e20))/2/pi, 'color', [27, 158, 119]/256)
plot(taxis_tesla, sqrt(abs((polyval([dep_pb(end), eig_val_tesla(mode_tesla)], taxis_tesla))*1e20))/2/pi, 'color', [217, 95, 2]/256)
ylim([1e9, 3.5e9])
xlim([-0.09, 1.09])
grid on
ylabel('f (Hz)')
xlabel('t')
title({'Tracking plot from TESLA cavity (t = 0) to pillbox cavity (t = 1).', ...
    'Taylor expansions at t = 0 compared to the tracking results.'})
legend('Tracking', 'Order 0', 'Order 1', 'Location', 'northwest')
