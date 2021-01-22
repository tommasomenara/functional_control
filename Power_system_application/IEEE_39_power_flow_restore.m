% this is the main script for the application of the weights adjustment 
% method on the IEEE 39 New England test case during a simulated fault

% NOTE: this code requires matpower7.1 (https://matpower.org/download/)
% and cvx (http://cvxr.com/cvx/)

% if you use this code, please cite the paper: 
% "T. Menara et al. (2021), Functional Control of Oscillator Networks"

% ------------------------------%
%       author: T. Menara       %
%              2021             %
% ------------------------------%

clear all
close all
clc

addpath('./functions/')

clear all
close all
clc

%% load parameters and solve nominal power flow

N = 39;

N_load = 29;

% generators inertia constants
M = [500 30.3 35.8 28.6 26 34.8 26.4 24.3 34.5 42];

% solve powerflow with MATPOWER
addpath('matpower7.1');
mpc = loadcase('case39');
results = runpf('case39');

% mechanical power injected at generators
P_m = results.gen(:,2);

% load power
P_l = results.bus(1:29,3);

% initial angles
theta_0 = results.bus(:,9)*pi/180;

% impedance, admittance and adjacency matrix (from [Dorfler & Bullo, 2010])
Z = zeros(39);
Y = zeros(39);
A = zeros(39);
for i = 1:46
    j = mpc.branch(i,1);
    k = mpc.branch(i,2);
    Z(j,k) = mpc.branch(i,3)+1i*mpc.branch(i,4);
    Z(k,j) = Z(j,k);
    Y(j,k) = -1/Z(j,k);
    Y(k,j) = Y(j,k);
    A(j,k) = imag(Y(j,k)) * 345 * (results.bus(j,8)) * (results.bus(k,8)); % a_ij = im(Y_ij)*|V_i|*|V_j| (345 converts p.u. to voltages)
    A(k,j) = A(j,k);
end

%% Kuramoto simulation only loads (N=29)

dt = .00005; % time step
T = 0.005; % final time
time = (0:dt:T)'; % time vector

figure, imagesc(A(1:29, 1:29)), colorbar, title('Loads adjacency matrix')

omega = P_l;

w = omega(1:29) - mean(omega(1:29)); % center at zero natural frequencies loads

Phases_evol = Kuramoto_fun(A(1:29, 1:29), 1, N_load, time, w, theta_0(1:29));

for t = 1:length(time)
    Phase_differences(t,:) = diff(Phases_evol(t,:));
end

figure, plot(time, Phase_differences), title('phase differences Kuramoto')

rho_matrix = zeros(N_load,N_load); % "correlation" matrix
c=0;
for i = 1:N_load-1
    for j = i+1:N_load
        rho_matrix(i,j) = mean(cos(Phases_evol(end,j)-Phases_evol(end,i)));
    end
end
rho_matrix = rho_matrix+rho_matrix'+eye(N_load);

figure, imagesc(rho_matrix), colorbar, title('functional pattern original')

%% Disconnect one branch

faulty = 25; % disconnect branch

% impedance, admittance and adjacency matrix (from [Dorfler & Bullo, 2010])
Z_fault = zeros(39);
Y_fault = zeros(39);
A_fault = zeros(39);
for i = [1:faulty-1 faulty+1:46] % do not load parameters from a branch 
    j = mpc.branch(i,1);
    k = mpc.branch(i,2);
    Z_fault(j,k) = mpc.branch(i,3)+1i*mpc.branch(i,4);
    Z_fault(k,j) = Z_fault(j,k);
    Y_fault(j,k) = -1/Z_fault(j,k);
    Y_fault(k,j) = Y_fault(j,k);
    A_fault(j,k) = imag(Y_fault(j,k)) * 345 * (results.bus(j,8)) * (results.bus(k,8)); % a_ij = im(Y_ij)*|V_i|*|V_j|
    A_fault(k,j) = A_fault(j,k);
end

Phases_evol_after_fault = Kuramoto_fun(A_fault(1:29, 1:29), 1, N_load, time, w, Phases_evol(end,:)');

for t = 1:length(time)
    Phase_differences_after_fault(t,:) = diff(Phases_evol_after_fault(t,:));
end

figure, plot(time, Phase_differences_after_fault), title('phase differences Kuramoto after fault')

rho_matrix_fault = zeros(N_load,N_load); % "correlation" matrix
c=0;
for i = 1:N_load-1
    for j = i+1:N_load
        rho_matrix_fault(i,j) = mean(cos(Phases_evol_after_fault(end,j)-Phases_evol_after_fault(end,i)));
    end
end
rho_matrix_fault = rho_matrix_fault+rho_matrix_fault'+eye(N_load);

figure, imagesc(rho_matrix_fault), colorbar, title('functional pattern faulty')

% difference between functional patterns before and after fault
figure, imagesc(rho_matrix-rho_matrix_fault), colorbar, title('functional pattern: original - faulty')

% difference between adjacency matrices before and after fault
figure, imagesc(A(1:29,1:29)-A_fault(1:29,1:29)), colorbar, title('adjacency matrix: original - faulty')

%% optimization

% define vector of original weights (except the one disconnected):
delta_A = [];
for i = 1:N_load-1
    for j = i+1:N_load
        if A_fault(i,j)>0
            delta_A(end+1,1) = A_fault(i,j);
        end
    end
end

B = (-adj2inc(sparse(triu(A_fault(1:29,1:29))>0)))'; % compute incidence matrix from adjacency matrix

A_UT = full(triu(A_fault(1:29,1:29)>0)); % upper triangular part of A_fault

x = []; % vector of phase differences for optimization
for i = 1:N_load-1
    for j = i+1:N_load
        x(end+1,1) = Phases_evol(end,j)-Phases_evol(end,i);
    end
end

c = 0; xx = [];
for i = 1:N_load-1
    for j = i+1:N_load
        c = c+1;
        if A_UT(i,j)==1
            xx(end+1,1) = x(c);
        end
    end
end

% diagonal matrix of sin(x_ij)
D = diag(sin(xx));

cvx_begin
variable delta_vec(34)
obj = norm(delta_vec-delta_A, 1);
minimize( 1*obj )
subject to
B*D*delta_vec == w;
delta_vec >= 0;
cvx_end

if size(cvx_status,2) == 6 && cvx_status(1) == 'S' % 'Solved' is 6 characters
    flag = 1;
        disp('Solution Found!');
        disp(datetime('now'));
        temp = zeros(N_load); c = 0;
        for i = 1:N_load-1
            for j = i+1:N_load
                if A_UT(i,j)==1
                    c = c+1;
                    temp(i,j) = delta_vec(c);
                end
            end
        end
        Solution = temp+temp';
else
    flag = 0;
end

%% integrate Kuramoto with new weights

if flag == 1
    Phases_evol_fixed = Kuramoto_fun(Solution, 1, N_load, time, w, Phases_evol_after_fault(end,:));
    
    for t = 1:length(time)
        Phase_differences_fixed(t,:) = diff(Phases_evol_fixed(t,:));
    end
    
    figure, plot(time, Phase_differences_fixed), title('phase differences Kuramoto after recovery')
    
    rho_matrix_fixed = zeros(N_load,N_load); % "correlation" matrix
    c=0;
    for i = 1:N_load-1
        for j = i+1:N_load
            rho_matrix_fixed(i,j) = (cos(Phases_evol_fixed(end,j)-Phases_evol_fixed(end,i)));
        end
    end
    rho_matrix_fixed = rho_matrix_fixed+rho_matrix_fixed'+eye(N_load);
    
    figure, imagesc(rho_matrix_fixed), colorbar, title('functional pattern after recovery')
    
    % difference between functional patterns before and after fault recovery
    figure, imagesc(rho_matrix-rho_matrix_fixed), colorbar, title('functional pattern: original - optimized')
    
    % difference between adjacency matrices before and after fault recovery
    figure, imagesc(A(1:29,1:29)-Solution), colorbar, title('adjacency matrix: original - optimized')
else
    error('infeasible optimization problem!')
end

%% fix colorbars

m = min([rho_matrix(:); rho_matrix_fault(:);rho_matrix_fixed(:)]);
figure(3), caxis([m 1]);
figure(5), caxis([m 1]);
figure(9), caxis([m 1]);