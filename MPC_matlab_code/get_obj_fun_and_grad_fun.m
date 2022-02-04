%% Check fmincon options for more details on f and g:

function [ f , g ] = get_obj_fun_and_grad_fun(du, x, v, eta_0, Ju, Px_Jv, lambda)

Np = length(du);

Qu  = zeros(3*Np, 3*Np);
eta = zeros(Np, 1);

uk = x(end-2);  
u = tril(ones(Np,Np)) * du + uk;

for i = 1 : Np
    alpha = u(i) * v(i);
    eta(i, 1) = (1/2)*(eta_0 + 1/eta_0) + (1/2)*(eta_0 - 1/eta_0)*tanh(10*alpha);
end

for i = 1 : (Np-1)
    Qu(3*(i-1)+2, 3*(i-1)+3) = eta(i);
    Qu(3*(i-1)+3, 3*(i-1)+2) = eta(i);
end
Qu(3*(Np-1)+2, 3*(Np-1)+3) = (1/2) * eta(Np);
Qu(3*(Np-1)+3, 3*(Np-1)+2) = (1/2) * eta(Np);

Hu = Ju' * Qu * Ju; 
Gu = Ju' * Qu * (Px_Jv);
Ku = (Px_Jv)' * Qu * (Px_Jv);

% Calculate objective function:
f = (1/2) * du' * Hu * du  ...
  + du' * Gu               ...
  + (1/2) * Ku             ...
  + lambda * (du' * du);

if nargout > 1 % gradient required
    
    d_eta = zeros(Np, 1);
    
    for i = 1 : Np
        alpha = u(i) * v(i);
        d_eta(i, 1) = 10*v(i)*(1/2)*(eta_0 - 1/eta_0)*(1 - (tanh(10*alpha))^2);
    end
    
    term_1 = zeros(Np, 1);
    term_2 = zeros(Np, 1);
    term_3 = zeros(Np, 1);
    term_4 = zeros(Np, 1);
    term_5 = zeros(Np, 1);
    term_6 = zeros(Np, 1);
    
    for i = 1 : Np
        d_Qu = zeros(3*Np, 3*Np);
        for j = i : (Np-1)
            d_Qu(3*(j-1)+2, 3*(j-1)+3) = d_eta(j);
            d_Qu(3*(j-1)+3, 3*(j-1)+2) = d_eta(j);
        end
        d_Qu(3*(Np-1)+2, 3*(Np-1)+3) = (1/2) * d_eta(Np);
        d_Qu(3*(Np-1)+3, 3*(Np-1)+2) = (1/2) * d_eta(Np);
        
        du_i = zeros(1, Np);
        du_i(1, i) = 1.0;
        
        term_1(i, 1) = du_i * Hu * du;
        term_2(i, 1) = (1/2) * du' * (Ju' * d_Qu * Ju) * du;
        term_3(i, 1) = du_i * Gu;
        term_4(i, 1) = du' * (Ju' * d_Qu * (Px_Jv));
        term_5(i, 1) = (1/2) * (Px_Jv)' * d_Qu * (Px_Jv);
        term_6(i, 1) = 2 * lambda * du(i, 1);
    end
    
    g = [];
    for i = 1 : Np
        d_Ju(i, 1) = term_1(i, 1) + term_2(i, 1) + term_3(i, 1)  ...
                   + term_4(i, 1) + term_5(i, 1) + term_6(i, 1);
        g = [ g  ;  d_Ju(i, 1)];
    end
    
end