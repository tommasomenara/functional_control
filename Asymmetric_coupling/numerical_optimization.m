% this is the main script for the weights adjustment method applied to
% directed networks (i.e., with asymmetric adjacency matrix)

% NOTE: this code requires cvx (http://cvxr.com/cvx/)

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

N = 7;

A = ones(N)-eye(N);
A(1,2) = 0; % disconnect one directed edge to generate asymmetry

omega = 0.1*[1 2 3 4 5 6 7]';
w = omega;

%% optimization

x_min = [pi/3 pi/4 pi/6 pi/8 pi/8 pi/6]; % phase differences
x_desired = [0 x_min(1) x_min(2) x_min(3) x_min(4) x_min(5) x_min(6);...
    0 0 x_min(2)-x_min(1) x_min(3)-x_min(1) x_min(4)-x_min(1) x_min(5)-x_min(1) x_min(6)-x_min(1);...
    0 0 0 x_min(3)-x_min(2) x_min(4)-x_min(2) x_min(5)-x_min(2) x_min(6)-x_min(2);...
    0 0 0 0 x_min(4)-x_min(3) x_min(5)-x_min(3) x_min(6)-x_min(3);...
    0 0 0 0 0 x_min(5)-x_min(4) x_min(6)-x_min(4);...
    0 0 0 0 0 0 x_min(6)-x_min(5)
    0 0 0 0 0 0 0];
x_desired = x_desired - x_desired';

% define vector of original weights (except the one disconnected):
delta_A = [];
for i = 1:N
    for j = 1:N
        if A(i,j)~=0 && i~=j
            delta_A(end+1,1) = A(i,j);
        end
    end
end

B = full(-adj2inc(sparse((A)>0)))'; % compute incidence matrix from adjacency matrix
B_sink = -double(B<0);
% A_UT = full(triu(A>0)); % upper triangular part of A_fault

xx = [];
for i = 1:N
    for j = 1:N
        if A(i,j)>0 % keep only differences between connected oscillators
            xx(end+1,1) = x_desired(i,j);
        end
    end
end

% diagonal matrix of sin(x_ij)
D = diag(sin(xx));

NN = (N^2-N)-1;

cvx_begin
variable delta_vec(NN)
obj = norm(delta_vec-delta_A, 1);
minimize( 1*obj )
subject to
w - B_sink*D*delta_vec == 1;
% delta_vec >= 0;
cvx_end

if size(cvx_status,2) == 6 && cvx_status(1) == 'S' % 'Solved' is 6 characters
    flag = 1;
        disp('Solution Found!');
        disp(datetime('now'));
        temp = zeros(N); c = 0;
        for i = 1:N
            for j = 1:N
                if A(i,j)>0
                    c = c+1;
                    temp(i,j) = delta_vec(c);
                end
            end
        end
        Solution = temp;
else
    flag = 0;
end

%% integrate Kuramoto with new weights

dt = .01; % time step
T = 3; % final time
time = (0:dt:T)'; % time vector

Theta_0 = [0 x_min]' +0.5*rand(N,1); % initial conditions

if flag == 1
    Phases_evol = Kuramoto_fun(Solution, 1, N, time, w, Theta_0);
    
    for j = [2 3 4 5 6 7]
        Phase_differences(:,j-1) = Phases_evol(:,j)-Phases_evol(:,1);
    end
    
    sum=0;
    for j = [2 3 4 5 6 7]
        sum = sum + A(1,j)*sin(Phases_evol(end,j)-Phases_evol(end,1));
    end
    frequency = w(1) + sum;
    
    figure, plot(time, Phases_evol), title('phases'), legend
    
    figure, plot(time, Phase_differences), title('phase differences'), legend
    
    rho_matrix = zeros(N,N); % "correlation" matrix
    c=0;
    for i = 1:N-1
        for j = i+1:N
            rho_matrix(i,j) = (cos(Phases_evol(end,j)-Phases_evol(end,i)));
        end
    end
    rho_matrix = rho_matrix+rho_matrix'+eye(N);
    
    figure, imagesc(rho_matrix), colorbar, title('functional pattern')
    
else
    error('infeasible optimization problem!')
end