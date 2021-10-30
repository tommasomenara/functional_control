% This script validates the usage of Kuramoto oscillators to map natural
% frequency changes to phase-locked fMRI time series

% NOTE: brain data are from the paper "A. Ponce-Alvarez et al., 2015 PLOS Comp Bio"

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

%% load data

subnetwork = [23 18 21 24 20 11 44 49 46 43 47 56]; % Cingulo Opercular regions of interest

N = length(subnetwork); % number of ROIs

subj = 18; % out of 24
run = 2; % out of 2
ts = dlmread(['Ponce-Alvarez_rs-fMRI_and_DTI/DataSet_rs-fMRI/subj', num2str(subj),'_block', num2str(run),'.txt']);

% reorder time series labels to match structural connectivity
load('Ponce-Alvarez_rs-fMRI_and_DTI/permutation');
timeseries = ts(:,permutation);

window = 249:269; % phase-locked time window

%% filter real BOLD signals

fs = 0.5; % sample rate (scanning TR = 2 seconds)

[temp_filtered,d] = bandpass(timeseries,[0.04  0.06],fs,'ImpulseResponse','iir','Steepness',0.95);

figure
subplot(2,1,1), plot(timeseries(:,:)),title('raw')
subplot(2,1,2), plot(temp_filtered(:,:), 'linewidth', 1.5),title('filtered')

filtered = temp_filtered(:, subnetwork);
filtered_window = temp_filtered(window, subnetwork);

%% Hilbert transform

z = hilbert(filtered);
for i = 1:N
    Hilbert_all(:,i) = atan2(real(z(:,i)),imag(z(:,i)));
end

BOLD_phases = Hilbert_all(window,:);

figure, plot(sin(BOLD_phases), 'Linewidth', 1.5), title('sin(BOLD phases)'), xlim([1 length(window)])

% remove bounds [-pi pi] on BOLD phases
BOLD_phases_no_mod = BOLD_phases;
for j = 1:size(BOLD_phases,2)
    for tt = 1:size(BOLD_phases,1)-1
        if abs(BOLD_phases(tt,j)-BOLD_phases(tt+1,j))>5
            BOLD_phases_no_mod(tt+1:end,j) = BOLD_phases_no_mod(tt+1:end,j)-2*pi;
        end
    end
end

figure, plot(BOLD_phases_no_mod, 'Linewidth', 1.5), title('BOLD phases'), xlim([1 length(window)])

%% Phase-locking values

PLV = zeros(N) + eye(N);
for i = 1:N-1
    for j = i+1:N
        temp = 0;
        for k = 1:length(window)
            temp = temp + exp(1i*(BOLD_phases(k,i)-BOLD_phases(k,j)));
        end
        PLV(i,j) = abs(1/(window(end)-window(1))*temp);
        PLV(j,i) = PLV(i,j);
    end
end

figure, imagesc(PLV), colorbar, caxis([0 1]), title('PLV matrix')

% total phase locking index
PL_index = 0;
for i = 1:N-1
    for j = i+1:N
        PL_index = PL_index + PLV(i,j);
    end
end
PL_index = PL_index/(N*(N-1)/2); % number between 0 (no PL) and 1 (full PL)
PL_tot(subj,run) = PL_index;

%%  SEED based correlation analysis: 
% https://fcp-indi.github.io/docs/user/sca.html

CORR_BOLD = corr(filtered_window);
CORR_BOLD_PHASES = corr(BOLD_phases);

for tt = 1:length(window)
    for i = 1:N
        for j = 1:N
            temp_cos(i,j,tt) = cos( BOLD_phases_no_mod(tt,j)-BOLD_phases_no_mod(tt,i) );
%             temp_diff(i,j,tt) = abs( BOLD_phases(tt,j)-BOLD_phases(tt,i) );
            temp_diff(i,j,tt) = ( BOLD_phases_no_mod(tt,j)-BOLD_phases_no_mod(tt,i) );
        end
    end
end
corr_cos = mean(temp_cos,3);
mean_phase_diff_matrix = mean(temp_diff,3);

figure
subplot(2,2,1), imagesc(CORR_BOLD), colorbar, caxis([-1 1])
title('Correlation BOLD signals')
subplot(2,2,2), imagesc(corr_cos), colorbar, caxis([-1 1])
title('$<\cos(\theta_j$ BOLD $ - $ $\theta_i$ BOLD$)>_{t}$', 'Interpreter','latex')
subplot(2,2,3), imagesc(mean_phase_diff_matrix), colorbar
title('$<\theta_j$ BOLD $ - $ $ \theta_i$ BOLD$>_{t}$', 'Interpreter','latex')
subplot(2,2,4), imagesc(cos(mean_phase_diff_matrix)), colorbar, caxis([-1 1])
title('$\cos(<\theta_j$ BOLD $ - $ $ \theta_i$ BOLD$>_{t})$', 'Interpreter','latex')

%% compute average phase differences and load parameters for optimization

C = zeros(N) + eye(N); % average phase differences as phasors
for i = 1:N-1
    for j = i+1:N
        C(i,j) = exp(1i*mean_phase_diff_matrix(i,j));
        C(j,i) = exp(-1i*mean_phase_diff_matrix(i,j));
    end
end

complex_phases = GPM_algo(C,N); % solve the nonlinear phase synchronization problem to estimate phases from C

for i = 1:N
    theta(i,1) = atan2(real(complex_phases(i)),imag(complex_phases(i)));
end

x = []; % vector of phase differences for optimization
for i = 1:N-1
    for j = i+1:N
        x(end+1,1) = theta(j)-theta(i);
    end
end

FC = zeros(N);
for i = 1:N
    for j = 1:N
        FC(i,j) = cos(theta(j)-theta(i));
    end
end

figure, imagesc(FC), colorbar, caxis([-1 1]), title('synthetic FC')

load('Ponce-Alvarez_rs-fMRI_and_DTI/connectivity'); % load structural connectivity data

SC = weights(subnetwork,subnetwork) - weights(subnetwork,subnetwork).*eye(N); % remove diagonal entries
SC = (SC+SC')/2; % make sure it is symmetric
figure, imagesc(SC), colorbar;

%% Kuramoto simulation

SC_binary = double(SC>0);

B_SC = (-adj2inc(sparse(triu(SC)>0)))'; % incidence matrix original SC matrix

A = double(SC>0);

A_UT = full(triu(A>0));

c = 0; xx = [];
for i = 1:N-1
    for j = i+1:N
        c = c+1;
        if A_UT(i,j)==1
            xx(end+1,1) = x(c);
        end
    end
end

% diagonal matrix of sin(x_ij)
D = diag(sin(xx));

% define vector of original weights:
delta_SC = [];
for i = 1:N-1
    for j = i+1:N
        if A_UT(i,j)==1
            delta_SC(end+1,1) = SC(i,j);
        end
    end
end

% spanning tree
[Tree, pred] = graphminspantree(sparse(SC_binary));
% incidence matrix of spanning tree
warning off, Bspan = (adj2inc(Tree,0))'; warning on

omega = B_SC*D*delta_SC;

dt = .05; % time step
T = 40; % final time
time = (0:dt:T)';

X0 = theta;
% X0 = X0 + .1*rand(N,1);

Phases_evol = Kuramoto_fun(SC, 1, N, time, omega+0.29, X0);

figure, plot(time, sin(Phases_evol)), title('Phases'), legend

for t = 1:length(time)
    Phase_differences(t,:) = Bspan'*Phases_evol(t,:)';
end

rho_matrix = zeros(N,N); % "correlation" matrix
c=0;
for i = 1:N-1
    for j = i+1:N
        rho_matrix(i,j) = mean(cos(Phases_evol(:,j)-Phases_evol(:,i)));
    end
end
rho_matrix = rho_matrix+rho_matrix'+eye(N);

figure
subplot(1,3,1), imagesc(FC), colorbar, title('$<\cos(\hat{\theta}_j-\hat{\theta}_i)>_{t}$', 'Interpreter','latex'), caxis([-1 1])
subplot(1,3,2), imagesc(rho_matrix), colorbar, title('Our synthetic FC'), caxis([-1 1])
subplot(1,3,3), imagesc(FC-rho_matrix), colorbar, title('Difference'), caxis([-1 1])
figure
subplot(1,3,1), imagesc(cos(mean_phase_diff_matrix)), colorbar, caxis([-1 1]), title('$\cos(<|$BOLD phase differences$|>_{t})$', 'Interpreter','latex')
subplot(1,3,2), imagesc(rho_matrix), colorbar, title('Our synthetic FC'), caxis([-1 1])
subplot(1,3,3), imagesc(cos(mean_phase_diff_matrix)-rho_matrix), colorbar, title('Difference')
figure
subplot(1,3,1), imagesc(CORR_BOLD), colorbar, caxis([-1 1])
title('Correlation BOLD signals')
subplot(1,3,2), imagesc(rho_matrix), colorbar, title('Our synthetic FC'), caxis([-1 1])
subplot(1,3,3), imagesc(CORR_BOLD-rho_matrix), colorbar, title('Difference')

%% assess basin of attraction manually

pert = .75*pi;
X0 = theta-(2*pert)*rand(N,1)+pert;

T = 5000; % final time
time = (0:dt:T)';

Phases_evol2 = Kuramoto_fun(SC, 1, N, time, omega+.29+100, X0);
for t = 1:length(time)
    Phase_differences2(t,:) = Bspan'*(Phases_evol2(t,:))';
end
for t = 1:length(Phase_differences2(:,1)) % modulus 2*pi
    for j = 1:11
        if Phase_differences2(t,j)>2*pi
            Phase_differences2(t:end,j) = Phase_differences2(t:end,j)-2*pi;
        elseif Phase_differences2(t,j)<-2*pi
            Phase_differences2(t:end,j) = Phase_differences2(t:end,j)+2*pi;
        end
    end
end

goal = Bspan'*theta;
figure, plot(time(1500/dt:end), Phase_differences2(1500/dt:end,:)), title('Phase differences')
hold on
for i = 1:11
    plot(time(end),goal(i),'rx','Linewidth', 1.25)
end

%% Compute distance from manifold

pert = .5*pi;
reduced_time = 1:500:length(time);

for loop = 1:100
    loop
    X0 = theta-(2*pert)*rand(N,1)+pert;
    
    Phases_evol2 = Kuramoto_fun(SC, 1, N, time, omega+.29+100, X0);
%     for i = 1:12
%         Phases_evol2(:,i) = mod(Phases_evol2(:,i),2*pi);
%     end
    Phases_evol_reduced = Phases_evol2(reduced_time,:); % reduce sampling time
    for t = 1:length(reduced_time)
        Phase_differences_reduced(t,:) = Bspan'*(Phases_evol_reduced(t,:))';
    end
    for t = 1:length(reduced_time)
        for j = 1:11
            phase_diff = Phase_differences_reduced(t,j);
            if phase_diff>=2*pi
                Phase_differences_reduced(t:end,j) = Phase_differences_reduced(t:end,j)-floor(phase_diff/(2*pi))*2*pi;
            elseif phase_diff<=-2*pi
                Phase_differences_reduced(t:end,j) = Phase_differences_reduced(t:end,j)+floor(abs(phase_diff)/(2*pi))*2*pi;
            end
        end
    end
    
    for tt = 1:length(reduced_time)
        Norm(loop,tt) = norm(Phase_differences_reduced(tt,:)-goal');
    end
    
end

M = mean(Norm);
M_min = min(Norm);
M_max = max(Norm);

figure
plot(reduced_time, Norm)
hold on
plot(reduced_time, M, 'k','Linewidth', 2)
plot(reduced_time, M_min, 'r','Linewidth', 2)
plot(reduced_time, M_max, 'b','Linewidth', 2)

lower_bound = M-M_min;
upper_bound = M_max-M;

dlmwrite('average_and_bounds.txt',[reduced_time', M', lower_bound', upper_bound'],'delimiter', ',', 'precision', 5);