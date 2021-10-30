function dx = third_order(t,x,p)

N = p.N_gen; % number of generators
M = p.M; % inertia coefficients
D = p.D; % damping coefficients
P_m = p.P_m; % mechanical power
Y = p.Y; % matrix of susceptances
X = p.X; % static reactances
X_prime = p.X_prime; % transient reactances
V_f = p.V_f; % internal voltage or field flux
T = p.T; % transient time constants

NN = N-1; % number of generators integrated (gen 1 is constant because it is an infinite bus)

theta = x(1:NN);
theta(10) = 0;
omega= x(NN+1:2*NN);
omega(10) = 0;
V = x(2*NN+1:3*NN);
V(10) = abs(1.0205-0.0431i);

for i = 1:NN
    
    d_theta(i,1) = omega(i);
    
    coupling_omega = 0;
    coupling_V = 0;
    for j = 1:N
        coupling_omega = coupling_omega + imag(Y(i,j))*V(i)*V(j)*sin(theta(j)-theta(i));
        coupling_V = coupling_V + imag(Y(i,j))*V(j)*cos(theta(j)-theta(i));
    end
    
    d_omega(i,1)= 1/M(i) * (P_m(i) - D(i)*omega(i) + coupling_omega);
    
    d_V(i,1) = 1/T(i) * (V_f(i) - V(i) + (X(i) - X_prime(i)) * coupling_V);
    
end

dx = [d_theta; d_omega; d_V];
end