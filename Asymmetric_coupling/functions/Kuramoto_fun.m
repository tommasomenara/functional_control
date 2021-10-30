function Phases_evol = Kuramoto_fun(A, K, N, t, w, X0)

tspan = t;

p.A = A;
p.K = K;
p.N = N;
p.w = w;

options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

[~, Phases_evol] = ode45(@(t,x) Kur_45(t,x,p), tspan, X0, options);

end

%%

function phases = Kur_45(t, x, p)


A = p.A;
K = p.K;
w = p.w;
N = p.N;

% A: structural matrix 
% K: coupling strengths
% N: number of nodes
% w: natural frequencies

[Xj, Xi] = meshgrid(x(1:N));

coupling = K.*sum(A.*sin(Xj-Xi),2);

phases = w + coupling;

end