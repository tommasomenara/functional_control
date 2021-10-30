% this is the main script for the simulation of the IEEE 39 test case on
% the third-order model (aka one-axis model) as in "K. Sharafutdinov et al (2018), Chaos"
%
% IEEE 39 test case parameters from:
% "Y. Susuki et al. (2011), J Nonlinear Sci"
% and
% "A. Moeini et al. (2015), Power Engineering Conference (UPEC)"
%
% NOTE1: values are adjusted to same Base MVA
%
% NOTE2: generators order is: 10 2 3 4 5 6 7 8 9 1
%
% if you use this code, please cite the paper: 
% "T. Menara et al. (2021), Functional Control of Oscillator Networks"
%
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

%% load parameters

N_gen = 10; % number of generators

% generators inertia constants (/pi /60Hz)
M = [0.2228 0.1607 0.1899 0.1517 0.1379 0.1846 0.1401 0.1289 0.1830 2.6526]/10;
% M = [4.728, 2.133, 1.003]';

% generators damping coefficients
D = 0.005*ones(N_gen,1);

% static reactances
X = [1.000 2.950 2.495 2.620 6.700 2.540 2.950 2.900 2.106 0.200]/10;

% transient reactances
X_prime = [0.310 0.697 0.531 0.436 1.320 0.500 0.490 0.570 0.570 0.060]/10;

% transient time constants
T = [10.20 6.560 5.700 5.690 5.400 7.300 5.660 6.700 4.790 7.000]/10;

% mechanical power injected at generators
% P_m = [ 1.2500 2.4943 3.2500 3.1600 2.5400 3.2500 2.8000 2.7000 4.1500 5.0000 ]'/5;
P_m = [250 520.81 650 632 508 650 560 540 830 1000]'/1000;

% internal voltage (or field flux)
V_f = abs([1.0568+0.0204i 0.9974+0.1770i 0.9842+0.1993i 0.9804+0.1809i 1.0855+0.3753i 1.0569+0.2160i 1.0423+0.2158i 0.9742+0.1817i 0.9352+0.3043i 1.0205-0.0431i])';

% impedance matrix G+B*i (pre fault) of (lossless) network-reduced power system model
Y =  [ 0.7239-15.0009i   0.2080+0.9484i   0.2536+1.2183i   0.2565+1.2209i   0.1033+0.4161i   0.2348+1.0950i   0.1950+0.9237i   0.0670+2.9064i   0.1961+1.5928i   0.6099+4.8881i
       0.2080+0.9484i    0.2519-9.1440i   0.2603+1.8170i   0.1657+0.6891i   0.0655+0.2346i   0.1513+0.6180i   0.1259+0.5214i   0.0916+0.5287i   0.1150+0.4400i   0.4159+2.7502i
       0.2536+1.2183i    0.2603+1.8170i   0.3870-10.9096i  0.2142+0.9763i   0.0857+0.3326i   0.1959+0.8756i   0.1629+0.7386i   0.1144+0.6848i   0.1471+0.5888i   0.4569+2.9961i
       0.2565+1.2209i    0.1657+0.6891i   0.2142+0.9763i   0.8131-12.0737i  0.2843+1.9774i   0.3178+1.7507i   0.2633+1.4766i   0.1608+0.7478i   0.2104+0.8320i   0.3469+1.6513i
       0.1033+0.4161i    0.0655+0.2346i   0.0857+0.3326i   0.2843+1.9774i   0.1964-5.5114i   0.1309+0.5973i   0.1088+0.5038i   0.0645+0.2548i   0.0826+0.2831i   0.1397+0.5628i
       0.2348+1.0950i    0.1513+0.6180i   0.1959+0.8756i   0.3178+1.7507i   0.1309+0.5973i   0.4550-11.1674i  0.3366+3.1985i   0.1471+0.6707i   0.1920+0.7461i   0.3175+1.4810i
       0.1950+0.9237i    0.1259+0.5214i   0.1629+0.7386i   0.2633+1.4766i   0.1088+0.5038i   0.3366+3.1985i   0.4039-9.6140i   0.1223+0.5657i   0.1599+0.6294i   0.2638+1.2493i
       0.0670+2.9064i    0.0916+0.5287i   0.1144+0.6848i   0.1608+0.7478i   0.0645+0.2548i   0.1471+0.6707i   0.1223+0.5657i   0.6650-10.0393i  0.3225+1.2618i   0.0991+2.5318i
       0.1961+1.5928i    0.1150+0.4400i   0.1471+0.5888i   0.2104+0.8320i   0.0826+0.2831i   0.1920+0.7461i   0.1599+0.6294i   0.3225+1.2618i   0.9403-7.5882i   0.2377+1.5792i
       0.6099+4.8881i    0.4159+2.7502i   0.4569+2.9961i   0.3469+1.6513i   0.1397+0.5628i   0.3175+1.4810i   0.2638+1.2493i   0.0991+2.5318i   0.2377+1.5792i   5.9222-18.6157i ];

%% Third-order simulation (N=10)

% initial values (from power flow computation)
% theta0 =[0.0193 0.1757 0.1998 0.1824 0.3329 0.2016 0.2042 0.1844 0.3145 0]';
theta0 =[0.0193 0.1757 0.1998 0.1824 0.3329 0.2016 0.2042 0.1844 0.3145]';
% V0 = abs([1.0568+0.0204i 0.9974+0.1770i 0.9842+0.1993i 0.9804+0.1809i 1.0855+0.3753i 1.0569+0.2160i 1.0423+0.2158i 0.9742+0.1817i 0.9352+0.3043i 1.0205-0.0431i])';
V0 = abs([1.0568+0.0204i 0.9974+0.1770i 0.9842+0.1993i 0.9804+0.1809i 1.0855+0.3753i 1.0569+0.2160i 1.0423+0.2158i 0.9742+0.1817i 0.9352+0.3043i])';
% omega0 = 60*2*pi*ones(N_gen,1);
omega0 = zeros(N_gen-1,1);

X0 = [theta0; omega0; V0];

dt = .01; % time step
T_final = 30; % final time
time = (0:dt:T_final)'; % time vector

options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

p.N_gen = N_gen;
p.M = M;
p.D = D;
p.P_m = P_m;
p.Y = Y;
p.X = X;
p.X_prime = X_prime;
p.V_f = V_f;
p.T = T;

[~, State_evol] = ode45(@(t,x) third_order(t,x,p), time, X0, options);

NN = N_gen-1; % number of generators integrated (gen 1 is constant because it is an infinite bus)

Phases_evol = [zeros(length(time),1), State_evol];
for t = 1:length(time)
    Phase_differences(t,:) = diff(Phases_evol(t,:));
end

figure(1), 
subplot(3,1,1); plot(time, Phase_differences), title('phase differences Kuramoto')
subplot(3,1,2); plot(time, State_evol(:,NN+1:2*NN)), title('omega')
subplot(3,1,3); plot(time, State_evol(:,2*NN+1:end)), title('voltage')

rho_matrix = zeros(N_gen,N_gen); % "correlation" matrix
c=0;
for i = 1:N_gen-1
    for j = i+1:N_gen
        rho_matrix(i,j) = mean(cos(Phases_evol(end,j)-Phases_evol(end,i)));
    end
end
rho_matrix = rho_matrix+rho_matrix'+eye(N_gen);

figure(2), 
subplot(1,3,1),
imagesc(rho_matrix), colorbar, title('functional pattern (pre-fault)')

%% fault: three-phase fault happens near bus 16
% 
% % impedance matrix (fault on)
% Y_fault = [ 0.5383-15.7638i   0.0901+0.5182i   0.0994+0.6084i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i  -0.0490+2.4392i   0.0472+1.0736i   0.3589+3.8563i
%         0.0901+0.5182i    0.1779-9.3864i   0.1628+1.4731i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0180+0.2653i   0.0219+0.1476i   0.2564+2.1683i
%         0.0994+0.6084i    0.1628+1.4731i   0.2591-11.3971i  0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0181+0.3113i   0.0241+0.1739i   0.2483+2.1712i
%         0.0000+0.0000i    0.0000+0.0000i   0.0000+0.0000i   0.4671-14.0254i  0.1411+1.3115i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i
%         0.0000+0.0000i    0.0000+0.0000i   0.0000+0.0000i   0.1411+1.3115i   0.1389-5.7383i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i
%         0.0000+0.0000i    0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.1633-12.7378i  0.0947+1.8739i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i
%         0.0000+0.0000i    0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0947+1.8739i   0.2035-10.7312i  0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i
%         -0.0490+2.4392i   0.0180+0.2653i   0.0181+0.3113i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.5925-10.3254i  0.2297+0.9440i   -0.0579+1.8999i
%         0.0472+1.0736i    0.0219+0.1476i   0.0241+0.1739i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.2297+0.9440i   0.8235-7.9409i   0.0363+0.8770i
%         0.3589+3.8563i    0.2564+2.1683i   0.2483+2.1712i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i  -0.0579+1.8999i   0.0363+0.8770i   5.5826-20.0113i ];
%    
% % update initial state and impedance matrix
% p.Y = Y_fault;
% X0 = State_evol(end,:);
% 
% dt = .005; % time step
% T_final = 1/3; % final time
% time = (0:dt:T_final)'; % time vector
% 
% [~, State_evol2] = ode45(@(t,x) third_order(t,x,p), time, X0, options);
% 
% Phases_evol2 = [zeros(length(time),1), State_evol2];
% for t = 1:length(time)
%     Phase_differences2(t,:) = diff(Phases_evol2(t,:));
% end
% 
% figure, 
% subplot(3,1,1); plot(time, Phase_differences2), title('phase differences Kuramoto')
% subplot(3,1,2); plot(time, State_evol2(:,NN+1:2*NN)), title('omega')
% subplot(3,1,3); plot(time, State_evol2(:,2*NN+1:end)), title('voltage')
% 
% rho_matrix2 = zeros(N_gen,N_gen); % "correlation" matrix
% c=0;
% for i = 1:N_gen-1
%     for j = i+1:N_gen
%         rho_matrix2(i,j) = mean(cos(Phases_evol2(end,j)-Phases_evol2(end,i)));
%     end
% end
% rho_matrix2 = rho_matrix2+rho_matrix2'+eye(N_gen);
% 
% figure, imagesc(rho_matrix2), colorbar, title('functional pattern (during fault)')

%% post-fault

% network impedance matrix (post fault)
Y_post = [ 0.8012-14.3511i  0.2163+0.9784i   0.2559+1.1997i   0.1629+0.5591i   0.0629+0.1900i   0.1483+0.5013i   0.1237+0.4230i   0.1385+3.3322i   0.3015+2.1485i   0.6576+5.3495i
        0.2163+0.9784i   0.2525-9.1427i   0.2603+1.8161i   0.1565+0.6587i   0.0619+0.2243i   0.1429+0.5908i   0.1189+0.4984i   0.0980+0.5482i   0.1239+0.4653i   0.4214+2.7715i
        0.2559+1.1997i   0.2603+1.8161i   0.3868-10.9091i  0.2124+0.9954i   0.0853+0.3392i   0.1943+0.8927i   0.1615+0.7530i   0.1153+0.6724i   0.1479+0.5726i   0.4586+2.9829i
        0.1629+0.5591i   0.1565+0.6587i   0.2124+0.9954i   0.9236-11.4000i  0.3306+2.2074i   0.4194+2.3551i   0.3474+1.9863i   0.0782+0.3146i   0.0903+0.2669i   0.2878+1.1812i
        0.0629+0.1900i   0.0619+0.2243i   0.0853+0.3392i   0.3306+2.2074i   0.2151-5.4330i   0.1734+0.8035i   0.1440+0.6778i   0.0308+0.1071i   0.0343+0.0905i   0.1135+0.4020i
        0.1483+0.5013i   0.1429+0.5908i   0.1943+0.8927i   0.4194+2.3551i   0.1734+0.8035i   0.5485-10.6253i  0.4139+3.6557i   0.0714+0.2821i   0.0820+0.2392i   0.2627+1.0592i
        0.1237+0.4230i   0.1189+0.4984i   0.1615+0.7530i   0.3474+1.9863i   0.1440+0.6778i   0.4139+3.6557i   0.4679-9.2284i   0.0594+0.2380i   0.0685+0.2019i   0.2187+0.8936i
        0.1385+3.3322i   0.0980+0.5482i   0.1153+0.6724i   0.0782+0.3146i   0.0308+0.1071i   0.0714+0.2821i   0.0594+0.2380i   0.7257-9.7609i   0.4096+1.6248i   0.1451+2.8344i
        0.3015+2.1485i   0.1239+0.4653i   0.1479+0.5726i   0.0903+0.2669i   0.0343+0.0905i   0.0820+0.2392i   0.0685+0.2019i   0.4096+1.6248i   1.0644-7.1152i   0.3063+1.9743i
        0.6576+5.3495i   0.4214+2.7715i   0.4586+2.9829i   0.2878+1.1812i   0.1135+0.4020i   0.2627+1.0592i   0.2187+0.8936i   0.1451+2.8344i   0.3063+1.9743i   5.9509-18.2881i ];

% Y_post = [ 0.5383-15.7638i   0.0901+0.5182i   0.0994+0.6084i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i  -0.0490+2.4392i   0.0472+1.0736i   0.3589+3.8563i
%         0.0901+0.5182i    0.1779-9.3864i   0.1628+1.4731i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0180+0.2653i   0.0219+0.1476i   0.2564+2.1683i
%         0.0994+0.6084i    0.1628+1.4731i   0.2591-11.3971i  0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0181+0.3113i   0.0241+0.1739i   0.2483+2.1712i
%         0.0000+0.0000i    0.0000+0.0000i   0.0000+0.0000i   0.4671-14.0254i  0.1411+1.3115i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i
%         0.0000+0.0000i    0.0000+0.0000i   0.0000+0.0000i   0.1411+1.3115i   0.1389-5.7383i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i
%         0.0000+0.0000i    0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.1633-12.7378i  0.0947+1.8739i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i
%         0.0000+0.0000i    0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0947+1.8739i   0.2035-10.7312i  0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i
%         -0.0490+2.4392i   0.0180+0.2653i   0.0181+0.3113i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.5925-10.3254i  0.2297+0.9440i   -0.0579+1.8999i
%         0.0472+1.0736i    0.0219+0.1476i   0.0241+0.1739i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.2297+0.9440i   0.8235-7.9409i   0.0363+0.8770i
%         0.3589+3.8563i    0.2564+2.1683i   0.2483+2.1712i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i   0.0000+0.0000i  -0.0579+1.8999i   0.0363+0.8770i   5.5826-20.0113i ];

p.Y = Y_post;
X0 = State_evol(end,:);

[~, State_evol3] = ode45(@(t,x) third_order(t,x,p), time, X0, options);

Phases_evol3 = [zeros(length(time),1), State_evol3];
for t = 1:length(time)
    Phase_differences3(t,:) = diff(Phases_evol3(t,:));
end

figure(4), 
subplot(3,1,1); plot(time, Phase_differences3), title('phase differences Kuramoto')
subplot(3,1,2); plot(time, State_evol3(:,NN+1:2*NN)), title('omega')
subplot(3,1,3); plot(time, State_evol3(:,2*NN+1:end)), title('voltage')

rho_matrix3 = zeros(N_gen,N_gen); % "correlation" matrix
c=0;
for i = 1:N_gen-1
    for j = i+1:N_gen
        rho_matrix3(i,j) = mean(cos(Phases_evol3(end,j)-Phases_evol3(end,i)));
    end
end
rho_matrix3 = rho_matrix3+rho_matrix3'+eye(N_gen);

figure(2), 
subplot(1,3,2),
imagesc(rho_matrix3), colorbar, title('functional pattern (fault)')


%% optimization of admittance matrix

% define vector of steady state voltages
VV = [State_evol3(end,2*NN+1:3*NN) abs(1.0205-0.0431i)];

% define vector of natural frequencies
bar_w = P_m - D.*[State_evol(end,NN+1:2*NN)'; 0];
% bar_w = P_m - D;
bar_w = bar_w-mean(bar_w);

% define vector of original weights:
delta_A = [];
for i = 1:N_gen-1
    for j = i+1:N_gen
        if Y(i,j)~=0
            delta_A(end+1,1) = imag(Y_post(i,j))*VV(i)*VV(j); % Im(Y_ij)*Vi*Vj
        end
    end
end

B = (-adj2inc(sparse(triu(Y-Y.*eye(N_gen))>0)))'; % compute incidence matrix from adjacency matrix

A_UT = full(triu(Y>0)); % upper triangular part of Y_post

x = []; % vector of phase differences for optimization
for i = 1:N_gen-1
    for j = i+1:N_gen
        x(end+1,1) = Phases_evol(end,j)-Phases_evol(end,i);
    end
end

c = 0; xx = [];
for i = 1:N_gen-1
    for j = i+1:N_gen
        c = c+1;
        if A_UT(i,j)==1
            xx(end+1,1) = x(c);
        end
    end
end

% diagonal matrix of sin(x_ij)
DD = diag(sin(xx));

% cvx_begin
% variable delta_vec(45)
% obj = norm(delta_vec-delta_A, 1);
% minimize( 1*obj )
% subject to
% B*DD*delta_vec == bar_w;
% delta_vec >= 0;
% cvx_end
% 
% if size(cvx_status,2) == 6 && cvx_status(1) == 'S' % 'Solved' is 6 characters
%     flag = 1;
%         disp('Solution Found!');
%         disp(datetime('now'));
%         temp = zeros(N_gen); c = 0;
%         for i = 1:N_gen-1
%             for j = i+1:N_gen
%                 if A_UT(i,j)==1
%                     c = c+1;
% %                     temp(i,j) = delta_vec(c)/VV(i)/VV(j); % new admittance matrix
%                     temp(i,j) = delta_vec(c); % new admittance matrix
%                 end
%             end
%         end
%         Solution = temp+temp';
% else
%     flag = 0;
% end
% 
% figure
% subplot(1,2,1), imagesc(imag(Y-Y.*eye(N_gen))), colorbar, title('Original admittance matrix'), caxis([0 max(delta_vec)])
% subplot(1,2,2), imagesc(Solution), colorbar, title('New admittance matrix'), caxis([0 max(delta_vec)])

%% natural frequency assignment

new_bar_w = B*DD*delta_A;
new_P_m = new_bar_w + D.*[State_evol(end,NN+1:2*NN)'; 0];

flag = 1;

%% integrate Kuramoto with new generators' mechanical inputs

new_P_m = new_P_m - min(new_P_m) + 0.25;
p.P_m = new_P_m/10; % rescale so that coupling strengths dominate nat. freq. heterogeneity
%     p.P_m = new_P_m/max(new_P_m);
%     p.V_f = State_evol(end,2*NN+1:end); % adjust internal voltages
X0 = State_evol3(end,:);

[~, State_evol4] = ode45(@(t,x) third_order(t,x,p), time, X0, options);

Phases_evol4 = [zeros(length(time),1), State_evol4];
for t = 1:length(time)
    Phase_differences4(t,:) = diff(Phases_evol4(t,:));
end

figure(6),
subplot(3,1,1); plot(time, Phase_differences4), title('phase differences Kuramoto')
subplot(3,1,2); plot(time, State_evol4(:,NN+1:2*NN)), title('omega')
subplot(3,1,3); plot(time, State_evol4(:,2*NN+1:end)), title('voltage')

rho_matrix4 = zeros(N_gen,N_gen); % "correlation" matrix
c=0;
for i = 1:N_gen-1
    for j = i+1:N_gen
        rho_matrix4(i,j) = mean(cos(Phases_evol4(end,j)-Phases_evol4(end,i)));
    end
end
rho_matrix4 = rho_matrix4+rho_matrix4'+eye(N_gen);

figure(2), subplot(1,3,3),
imagesc(rho_matrix4), colorbar, title('functional pattern (optimized)')

% difference between functional patterns before and after fault recovery
RR1 = abs(rho_matrix-rho_matrix3);
RR2 = abs(rho_matrix-rho_matrix4);
figure(8),
subplot(1,2,1), imagesc(RR1), colorbar, title('functional pattern: original - fault')
caxis([0, max([rho_matrix(:)-rho_matrix4(:);rho_matrix(:)-rho_matrix3(:)])])
subplot(1,2,2), imagesc(RR2), colorbar, title('functional pattern: original - optimized'),
caxis([0, max([rho_matrix(:)-rho_matrix4(:);rho_matrix(:)-rho_matrix3(:)])])
mm = min([rho_matrix(:); rho_matrix3(:); rho_matrix4(:)]);
figure(2), subplot(1,3,1), caxis([mm 1]);
figure(2), subplot(1,3,2), caxis([mm 1]);
figure(2), subplot(1,3,3), caxis([mm 1]);

disp('Frobenius norm original-fault:')
norm(rho_matrix-rho_matrix3,'fro')
disp('Frobenius norm original-recovered:')
norm(rho_matrix-rho_matrix4,'fro')