% This is the main script for the example used in Figure 1

% Please comment lines 40-43 and uncomment lines 51-109 and 140-197 to assign different equilibria

% NOTE: this code requires cvx (http://cvxr.com/cvx/) if you want to run
% your own optimization

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

%% Kuramoto parameters

N = 7; % number of oscillators

rho_matrix = zeros(N,N); % initialize "correlation" matrix

A = zeros(N); % initialize adjacency matrix
for i = 1:N-1
    for j = i+1:N
        if rand>0.1
            A(i,j) = 2;
        end
    end
end
A = A+A';

omega = rand(N,1); % random natural frequencies
w = omega - mean(omega); % center at zero natural frequencies

%% load same parameters as in Fig. 1

load('w')
load('b')
load('Solution')
load('Solution2')

%% optimization 1

x_min = [pi/8 pi/8 pi/6 pi/6 pi/3 2*pi/3]; % first pattern
% x_desired = [0 x_min(1) x_min(2) x_min(3) x_min(4) x_min(5) x_min(6);...
%     0 0 x_min(2)-x_min(1) x_min(3)-x_min(1) x_min(4)-x_min(1) x_min(5)-x_min(1) x_min(6)-x_min(1);...
%     0 0 0 x_min(3)-x_min(2) x_min(4)-x_min(2) x_min(5)-x_min(2) x_min(6)-x_min(2);...
%     0 0 0 0 x_min(4)-x_min(3) x_min(5)-x_min(3) x_min(6)-x_min(3);...
%     0 0 0 0 0 x_min(5)-x_min(4) x_min(6)-x_min(4);...
%     0 0 0 0 0 0 x_min(6)-x_min(5)
%     0 0 0 0 0 0 0];
% 
% % define vector of original weights (except the one disconnected):
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
% A_UT = full(triu(A>0)); % upper triangular part of A_fault
% 
% xx = [];
% for i = 1:N-1
%     for j = i+1:N
%         if A_UT(i,j)==1
%             xx(end+1,1) = x_desired(i,j);
%         end
%     end
% end
% 
% D = diag(sin(xx)); % diagonal matrix of sin(x_ij)
% 
% cvx_begin
% variables delta_vec(size(delta_A)) b(N)
% obj1 = norm(delta_vec, 2);
% obj2 = norm(b, 2);
% minimize( 1*obj1 + 1*obj2)
% subject to
% B*D*(delta_vec+delta_A) == omega+b;
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
%                     temp(i,j) = delta_A(c) + delta_vec(c);
%                 end
%             end
%         end
%         Solution = temp+temp';
%         w = omega+b;
% else
%     error('No solution found')
% end

%% integrate Kuramoto model

dt = .01; % time step
T = 20; % final time
time = (0:dt:T)'; % time vector

theta_0 = [0 x_min]' + rand(N,1);

Phases_evol = Kuramoto_fun(Solution, 1, N, time, w, theta_0);

for t = 1:length(time)
    Phase_differences(t,:) = diff(Phases_evol(t,:));
end

rho_matrix = zeros(N,N);
c=0;
for i = 1:N-1
    for j = i+1:N
        rho_matrix(i,j) = mean(cos(Phases_evol(end,j)-Phases_evol(end,i)));
    end
end
rho_matrix = rho_matrix+rho_matrix'+eye(N);

figure, imagesc(rho_matrix), colorbar, title('functional pattern original'), caxis([-1 1])
figure, plot(time, mod(Phases_evol,2*pi)), title('phases Kuramoto'), ylim([0 pi])

%% optimization 2

x_min = [2*pi/3 pi/3 pi/6 pi/6 pi/8 pi/8]; % second pattern
% x_desired = [0 x_min(1) x_min(2) x_min(3) x_min(4) x_min(5) x_min(6);...
%     0 0 x_min(2)-x_min(1) x_min(3)-x_min(1) x_min(4)-x_min(1) x_min(5)-x_min(1) x_min(6)-x_min(1);...
%     0 0 0 x_min(3)-x_min(2) x_min(4)-x_min(2) x_min(5)-x_min(2) x_min(6)-x_min(2);...
%     0 0 0 0 x_min(4)-x_min(3) x_min(5)-x_min(3) x_min(6)-x_min(3);...
%     0 0 0 0 0 x_min(5)-x_min(4) x_min(6)-x_min(4);...
%     0 0 0 0 0 0 x_min(6)-x_min(5)
%     0 0 0 0 0 0 0];
% 
% % define vector of original weights:
% delta_A = [];
% for i = 1:N-1
%     for j = i+1:N
%         if Solution(i,j)>0
%             delta_A(end+1,1) = Solution(i,j);
%         end
%     end
% end
% 
% B = (-adj2inc(sparse(triu(Solution)>0)))'; % compute incidence matrix from adjacency matrix
% 
% A_UT = full(triu(Solution>0)); % upper triangular part of A_fault
% 
% xx = [];
% for i = 1:N-1
%     for j = i+1:N
%         if A_UT(i,j)==1
%             xx(end+1,1) = x_desired(i,j);
%         end
%     end
% end
% 
% D = diag(sin(xx)); % diagonal matrix of sin(x_ij)
% 
% cvx_begin
% variables delta_vec(size(delta_A)) b(N)
% obj1 = norm(delta_vec, 1);
% obj2 = norm(b, 1);
% minimize( 1*obj1 + 1*obj2)
% subject to
% B*D*(delta_vec+delta_A) == w+b;
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
%                     temp(i,j) = delta_A(c) + delta_vec(c);
%                 end
%             end
%         end
%         Solution2 = temp+temp';
% else
%     error('No solution found')
% end


%% integrate Kuramoto model

theta_0 = Phases_evol(end,:)';

Phases_evol2 = Kuramoto_fun(Solution2, 1, N, time, w+b, theta_0);

for t = 1:length(time)
    Phase_differences2(t,:) = diff(Phases_evol2(t,:));
end

rho_matrix2 = zeros(N,N);
c=0;
for i = 1:N-1
    for j = i+1:N
        rho_matrix2(i,j) = mean(cos(Phases_evol2(end,j)-Phases_evol2(end,i)));
    end
end
rho_matrix2 = rho_matrix2+rho_matrix2'+eye(N);

figure, imagesc(rho_matrix2), colorbar, title('functional pattern original'), caxis([-1 1])
figure, plot(time, Phases_evol2), title('phases Kuramoto'), ylim([0 pi])
