% This is the main script for the example in the Supplementary Information
% section "A Heuristic Method To Promote Stability of Functional Patterns in Positive Networks"

% Please comment lines 43-46 and uncomment lines 51-132 and 165-250 to assign different equilibria

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

N = 7;

% A = 2*(ones(N) - eye(N));
% for i = 1:N-1
%     for j = i+1:N
%         if rand>0.7
%             A(i,j) = 0;
%             A(j,i) = 0;
%         end
%     end
% end
%  load('A') % adjacency matrix

% omega = rand(N,1); % natural frequencies
% w = omega - mean(omega); % center at zero natural frequencies

%% load parameters used in the paper example

load('w') % natural frequencies
load('Solution1') % initial matrix
load('Solution2') % tailored matrix for stability
load('theta_0') % initial phases

%% optimization 1

x_min = [pi/4 pi/6 pi/6 pi/8 pi/8 pi/3]; % phase differences
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
% xx = [];
% for i = 1:N-1
%     for j = i+1:N
%         if A_UT(i,j)==1
%             xx(end+1,1) = x_desired(i,j);
%         end
%     end
% end
% 
% % find weights associated to negative cos(x_ij)
% c = 0; index = zeros(N); neg_index =[]; pos_index=[]; zero_index=[];
% for i = 1:N-1
%     for j=i+1:N
%         if A_UT(i,j) == 1
%             c = c+1;
%             if cos(xx(c))<0
%                 index(i,j) = 1;
%                 neg_index(end+1) = c;
%             else
%                 pos_index(end+1) = c;
%             end
%         else
%             zero_index(end+1) = c;
%         end
%     end
% end
% 
% % diagonal matrices of sin(x_ij)
% D = diag(sin(xx));
% 
% cvx_begin
% variable delta_vec(size(delta_A))
% if isempty(neg_index)
%     obj1 = 0;
% else
%     obj1 = norm(delta_vec(neg_index)+delta_A(neg_index),1);
% end
% obj2 = norm(delta_vec(pos_index), 1);
% minimize( 10*obj1 + .1*obj2 )
% subject to
% B*D*(delta_A+delta_vec) == w;
% (delta_A+delta_vec) >= 0;
% cvx_end
% 
% if size(cvx_status,2) == 6 && cvx_status(1) == 'S' % 'Solved' is 6 characters
%     disp('Solution Found!');
%     disp(datetime('now'));
%     temp = zeros(N); c = 0;
%     for i = 1:N-1
%         for j = i+1:N
%             if A_UT(i,j)==1
%                 c = c+1;
%                 temp(i,j) = delta_A(c) + delta_vec(c);
%             end
%         end
%     end
%     Solution1 = temp+temp';
% else
%     error('no solution found')
% end

%% Kuramoto simulation

dt = .1; % time step
T = 30; % final time
time = (0:dt:T)'; % time vector

theta_0 = [0 x_min]'+0.1*rand(N,1);

Phases_evol = Kuramoto_fun(Solution1, 1, N, time, w, theta_0);

for t = 1:length(time)
    for k = 2:N
        Phase_differences(t,k-1) = Phases_evol(t,k)-Phases_evol(t,1);
    end
end

rho_matrix = zeros(N,N);
c=0;
for i = 1:N-1
    for j = i+1:N
        rho_matrix(i,j) = mean(cos(Phases_evol(end,j)-Phases_evol(end,i)));
    end
end
rho_matrix = rho_matrix+rho_matrix'+eye(N);

figure, imagesc(rho_matrix), colorbar, title('functional pattern 1'), caxis([-1 1])

figure, plot(time, Phase_differences), title('phase differences Kuramoto'), legend('x_{12}','x_{13}','x_{14}','x_{15}','x_{16}','x_{17}')

%% optimization 2

% x_min = [21*pi/32 pi/6 pi/6 pi/8 pi/8 pi/3]; % phase differences
% x_desired = [0 x_min(1) x_min(2) x_min(3) x_min(4) x_min(5) x_min(6);...
%     0 0 x_min(2)-x_min(1) x_min(3)-x_min(1) x_min(4)-x_min(1) x_min(5)-x_min(1) x_min(6)-x_min(1);...
%     0 0 0 x_min(3)-x_min(2) x_min(4)-x_min(2) x_min(5)-x_min(2) x_min(6)-x_min(2);...
%     0 0 0 0 x_min(4)-x_min(3) x_min(5)-x_min(3) x_min(6)-x_min(3);...
%     0 0 0 0 0 x_min(5)-x_min(4) x_min(6)-x_min(4);...
%     0 0 0 0 0 0 x_min(6)-x_min(5)
%     0 0 0 0 0 0 0];
% 
% Solution1 = Solution1.*full(double(Solution1>0.0001)); % smooth out spurious entries
% 
% % define vector of original weights:
% delta_A = [];
% for i = 1:N-1
%     for j = i+1:N
%         if Solution1(i,j)>0
%             delta_A(end+1,1) = Solution1(i,j);
%         end
%     end
% end
% 
% B = (-adj2inc(sparse(triu(Solution1)>0)))'; % compute incidence matrix from adjacency matrix
% 
% A_UT = full(triu(Solution1>0)); % upper triangular part of A
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
% % find weights associated to negative cos(x_ij)
% c = 0; index = zeros(N); neg_index =[]; pos_index=[]; zero_index=[];
% for i = 1:N-1
%     for j=i+1:N
%         if A_UT(i,j) == 1
%             c = c+1;
%             if cos(xx(c))<0
%                 index(i,j) = 1;
%                 neg_index(end+1) = c;
%             else
%                 pos_index(end+1) = c;
%             end
%         else
%             zero_index(end+1) = c;
%         end
%     end
% end
% 
% % diagonal matrices of sin(x_ij)
% D = diag(sin(xx));
% 
% cvx_begin
% variable delta_vec(size(delta_A))
% if isempty(neg_index)
%     obj1 = 0;
% else
%     obj1 = norm(delta_vec(neg_index)+delta_A(neg_index),1);
% end
% % obj2 = norm(delta_vec-delta_A, 1);
% obj2 = norm(delta_vec(pos_index), 1);
% minimize( 10*obj1 + .1*obj2 )
% subject to
% B*D*(delta_A+delta_vec) == w;
% (delta_A+delta_vec) >= 0;
% cvx_end
% 
% if size(cvx_status,2) == 6 && cvx_status(1) == 'S' % 'Solved' is 6 characters
%     disp('Solution Found!');
%     disp(datetime('now'));
%     temp = zeros(N); c = 0;
%     for i = 1:N-1
%         for j = i+1:N
%             if A_UT(i,j)==1
%                 c = c+1;
%                 temp(i,j) = delta_A(c) + delta_vec(c);
%             end
%         end
%     end
%     Solution2 = temp+temp';
% else
%     error('no solution found')
% end

%% integrate Kuramoto after network correction

T = 50; % final time
time = (0:dt:T)'; % time vector

Phases_evol_2 = Kuramoto_fun(Solution2, 1, N, time, w, theta_0);

for t = 1:length(time)
    for k = 2:N
        Phase_differences_2(t,k-1) = Phases_evol_2(t,k)-Phases_evol_2(t,1);
    end
end

rho_matrix2 = zeros(N,N);
c=0;
for i = 1:N-1
    for j = i+1:N
        rho_matrix2(i,j) = mean(cos(Phases_evol_2(end,j)-Phases_evol_2(end,i)));
    end
end
rho_matrix2 = rho_matrix2+rho_matrix2'+eye(N);

figure, imagesc(rho_matrix2), colorbar, title('functional pattern 2'), caxis([-1 1])

figure, plot(time, Phase_differences_2), title('phase differences Kuramoto'), legend('x_{12}','x_{13}','x_{14}','x_{15}','x_{16}','x_{17}')

%% This section is used to plot the Gerschgorin disks and the Jacobian's eigenvalues

% define vector of original weights:
% delta_A1 = [];
% for i = 1:N-1
%     for j = i+1:N
%         if Solution1(i,j)>0
%             delta_A1(end+1,1) = Solution1(i,j);
%         end
%     end
% end
% delta_A2 = [];
% for i = 1:N-1
%     for j = i+1:N
%         if Solution2(i,j)>0
%             delta_A2(end+1,1) = Solution2(i,j);
%         end
%     end
% end
% 
% n_points = 4;
% change1 = linspace(0.1706, 0, n_points);
% change2 = linspace(0.5796, 0.8948, n_points);
% change3 = linspace(1.3434, 0.4686, n_points);
% change4 = linspace(2, 2.9242, n_points);
% 
% figure;
% for k = 1:n_points
%     delta_A1(1) = change1(k);
%     delta_A1(2) = change2(k);
%     delta_A1(3) = change3(k);
%     delta_A1(5) = change4(k);
%     J = -B*diag(delta_A1.*cos(xx))*B';
%     ee = eig(J)
%     theta = linspace(0,2*pi);
%     for i = 1:N
%         center(i,k) = J(i,i);
%         radius(i,k) = sum(abs(J(i,[1:i-1,i+1:N])));
%         subplot(2,4,k)
%         hold on
%         r = 2;
%         xc = 4;
%         yc = 3;
%         x = radius(i,k)*cos(theta) + center(i,k);
%         y = radius(i,k)*sin(theta) + 0;
%         plot(x,y,'b')
%         plot([0,0],[-10,10],'--r')
%         grid on
%         plot(ee,zeros(N,1),'x','linewidth',2)
%         xlim([-14.5,0.5]), ylim([-7.5,7.5])
%         axis square
%         
%         subplot(2,4,n_points+k)
%         hold on
%         plot(ee,zeros(N,1),'x','linewidth',2)
%         grid on
%         plot([0,0],[-1,1],'--r')
%         xlim([-1,0.5]), ylim([-1,1])
%         title('zoomed in')
%     end
% end