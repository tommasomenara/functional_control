% This is the main script for the example on multiple equilibria

% Please comment lines 32-35 and uncomment lines 39-117 to assign different equilibria

% if you use this code, please cite the paper: 
% "T. Menara et al. (2021), Functional Control of Oscillator Networks"

% NOTE: this code requires cvx (http://cvxr.com/cvx/) if you want to run
% your own optimization

% ------------------------------%
%       author: T. Menara       %
%              2021             %
% ------------------------------%

clear all
close all
clc

addpath('./functions/')

%% Kuramoto parameters

N = 7; % number of oscillators

A = 2*(ones(N) - eye(N)); % initialize adjacency matrix

rho_matrix = zeros(N,N); % initialize "correlation" matrix

% omega = rand(N,1); % random natural frequencies
% w = omega - mean(omega); % center at zero natural frequencies

%% load same parameters used in the paper example

load('w') % load natural frequencies
load('Solution') % load adjacency matrix
load('theta_0') % load initial phases
load('theta_0_perturbed') % load perturbed phases

%% optimization of network weights (uncomment to run the optimization from scratch)

% x_min = [pi/6 pi/6 pi/4 pi/4 pi/6 pi/4]; % first pattern
% x_desired1 = [0 x_min(1) x_min(2) x_min(3) x_min(4) x_min(5) x_min(6);...
%     0 0 x_min(2)-x_min(1) x_min(3)-x_min(1) x_min(4)-x_min(1) x_min(5)-x_min(1) x_min(6)-x_min(1);...
%     0 0 0 x_min(3)-x_min(2) x_min(4)-x_min(2) x_min(5)-x_min(2) x_min(6)-x_min(2);...
%     0 0 0 0 x_min(4)-x_min(3) x_min(5)-x_min(3) x_min(6)-x_min(3);...
%     0 0 0 0 0 x_min(5)-x_min(4) x_min(6)-x_min(4);...
%     0 0 0 0 0 0 x_min(6)-x_min(5)
%     0 0 0 0 0 0 0];
% 
% x_min2 = [pi/8 pi/3 pi/4 pi/4 pi/6 pi/4]; % second pattern
% x_desired2 = [0 x_min2(1) x_min2(2) x_min2(3) x_min2(4) x_min2(5) x_min2(6);...
%     0 0 x_min2(2)-x_min2(1) x_min2(3)-x_min2(1) x_min2(4)-x_min2(1) x_min2(5)-x_min2(1) x_min2(6)-x_min2(1);...
%     0 0 0 x_min2(3)-x_min2(2) x_min2(4)-x_min2(2) x_min2(5)-x_min2(2) x_min2(6)-x_min2(2);...
%     0 0 0 0 x_min2(4)-x_min2(3) x_min2(5)-x_min2(3) x_min2(6)-x_min2(3);...
%     0 0 0 0 0 x_min2(5)-x_min2(4) x_min2(6)-x_min2(4);...
%     0 0 0 0 0 0 x_min2(6)-x_min2(5)
%     0 0 0 0 0 0 0];
% 
% % define vector of original network weights:
% delta_A = [];
% for i = 1:N-1
%     for j = i+1:N
%         if A(i,j)>0
%             delta_A(end+1,1) = A(i,j);
%         end
%     end
% end
% 
% B = (-adj2inc(sparse(triu(A)>0)))'; % compute incidence matrix from adjacency matrix
% 
% A_UT = full(triu(A>0)); % upper triangular part of A
% 
% xx1 = [];
% for i = 1:N-1
%     for j = i+1:N
%         if A_UT(i,j)==1
%             xx1(end+1,1) = x_desired1(i,j);
%         end
%     end
% end
% xx2 = [];
% for i = 1:N-1
%     for j = i+1:N
%         if A_UT(i,j)==1
%             xx2(end+1,1) = x_desired2(i,j);
%         end
%     end
% end
% 
% % diagonal matrices of sin(x_ij)
% D1 = diag(sin(xx1));
% D2 = diag(sin(xx2));
% 
% cvx_begin
% variable delta_vec(size(delta_A))
% obj = norm(delta_vec-delta_A, 1);
% minimize( 1*obj )
% subject to
% [B*D1; B*D2]*delta_vec == [w;w];
% cvx_end
% 
% if size(cvx_status,2) == 6 && cvx_status(1) == 'S' % 'Solved' is 6 characters
%         disp('Solution Found!');
%         disp(datetime('now'));
%         temp = zeros(N); c = 0;
%         for i = 1:N-1
%             for j = i+1:N
%                 if A_UT(i,j)==1
%                     c = c+1;
%                     temp(i,j) = delta_vec(c);
%                 end
%             end
%         end
%         Solution = temp+temp';
% else
%     error('No solution found')
% end

%% integrate Kuramoto with ODE45

dt = .01; % time step
T = 150; % final time
time = (0:dt:T)'; % time vector

% theta_0 = [0 x_min]'+0.05*rand(N,1);

Phases_evol_temp1 = Kuramoto_fun(Solution, 1, N, time(1:5000), w, theta_0);

% theta_0_2 = [0 x_min]'+0.05*rand(N,1);

Phases_evol_temp2 = Kuramoto_fun(Solution, 1, N, time(5001:end), w, theta_0_2);

Phases_evol = [Phases_evol_temp1; Phases_evol_temp2];

for t = 1:length(time)
    for k = 2:7
        Phase_differences(t,k-1) = Phases_evol(t,k)-Phases_evol(t,1);
    end
end

figure, plot(time, Phase_differences), title('phase differences Kuramoto'), legend('x_{12}','x_{13}','x_{14}','x_{15}','x_{16}','x_{17}')