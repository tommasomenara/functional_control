function estimated_phases = GPM_algo(C,N)

for i = 1:N
    x0(i,1) = exp(1i * rand);
end

alpha = max([0, -min(eig(C))]);

I = eye(N);

C_tilde = C + alpha*I;

T_prec = zeros(N,1);
for i = 1:N
    if I(i,:)*C_tilde*x0 == 0
        T_prec(i,1) = x0(i);
    else
        T_prec(i,1) = (I(i,:)*C_tilde*x0)/abs(I(i,:)*C_tilde*x0);
    end
end

T_succ = zeros(N,1);
flag = 1;

while abs(flag)
        
        for i = 1:N
            if I(i,:)*C_tilde*T_prec == 0
                T_succ(i,1) = T_prec(i);
            else
                T_succ(i,1) = (I(i,:)*C_tilde*T_prec)/abs(I(i,:)*C_tilde*T_prec);
            end
        end
        
        if abs(T_prec-T_succ) < 1e-10
            flag = 0;
        end
        T_prec = T_succ; % update the estimate
        
end

estimated_phases = T_succ;

end