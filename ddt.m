function dwdt = ddt(w, p, gamma, x, xc)
%Method 2 (talked about on meeting 061021)
global n

N = length(x);
int_num = length(xc);
dx = 1/N;
w1 = w(:,1); %rho
w2 = w(:,2); %rho*u
w3 = w(:,3); % E

u = w2 ./ w1;   
c = SoS_SG(gamma,w1, p, 0);

% Evaluate f(w)
F1 = w2;
F2 = w2.^2 ./ w1 + p;
F3 = (w3+ p) .* w2 ./ w1;
F = [F1 F2 F3];
    
nu = (abs(u) + c)/ 2;
nu_interface = max(nu(1:end-1), nu(2:end));
    
% Evaluate numerical flux  
F_interface = 0.5* (F(1:end-1,:) + F(2:end,:)) + nu_interface .* (w(1:end-1, :) - w(2:end, :));

% Calculate flux at material interface
k = 1;
F_interfacejm = zeros(int_num, 3);
F_interfacejp = zeros(int_num, 3);

for j = 1:N
    if any(xc==j)    
        gammal = gamma(j); p0l = 0; rhol = w1(j); vl = u(j); pl = p(j);
        gammar = gamma(j+1); p0r = 0; rhor = w1(j+1); vr = u(j+1); pr = p(j+1);

        if pl<0 || pr<0 || rhol<0 || rhor<0
            disp('Negative pressure or density')
            error_message;
        end
        %   EXACT RIEMANN SOLVER
        rsolve = rsol_func_mex(gammal, p0l, gammar, p0r, rhol, vl, pl, rhor, vr, pr);
        psr = rsolve(1); psl = rsolve(2); vsr = rsolve(3); vsl = rsolve(4); rhosr = rsolve(5); rhosl = rsolve(6);
        esl = e_SG(gammal,rhosl, psl, 0); esr = e_SG(gammar,rhosr, psr, 0);
        Esl = E_P2C(rhosl, vsl, esl); Esr = E_P2C(rhosr, vsr, esr);

        w1sl = rhosl; w2sl = rhosl*vsl; w3sl = Esl; wsl = [w1sl, w2sl, w3sl];
        w1sr = rhosr; w2sr = rhosr*vsr; w3sr = Esr; wsr = [w1sr, w2sr, w3sr];

        F1sl = w2sl; F2sl = w2sl^2/w1sl + psl; F3sl = (w3sl + psl) * w2sl / w1sl; Fsl = [F1sl, F2sl, F3sl];
        F1sr = w2sr; F2sr = w2sr^2/w1sr + psr; F3sr = (w3sr + psr) * w2sr / w1sr; Fsr = [F1sr, F2sr, F3sr];

        csl = SoS_SG(gammal,rhosl, psl, 0); csr = SoS_SG(gammar,rhosr, psr, 0);
        nusl = (max((abs(u(j)) + c(j)),(abs(vsl) + csl)))/2;
        nusr = (max((abs(u(j+1)) + c(j+1)),(abs(vsr) + csr)))/2;

        F_interfacejm(k, :) = 0.5*(F(j,:) + Fsl) + nusl*(w(j, :) - wsl);
        F_interfacejp(k, :) = 0.5*(Fsr + F(j+1,:)) + nusr*(wsr- w(j+1, :));
        k = k+1;
    end
end

% Apply BC
BC_left = [w(1,2), w(1,2)^2/w(1,1) + p(1),  (w(1,3) + p(1))* w(1,2)/w(1,1)];
BC_right = [w(end,2), w(end,2)^2/w(end,1) + p(end), (w(end,3) + p(end))* w(end,2)/w(end,1)];    
F_interface = [BC_left ; F_interface ; BC_right];

% Evaluate dwdt
dwdt = zeros(N,3);
k = 1;
for j = 1:N
    if any(xc==j)
        dwdt(j,:) = (F_interface(j, :) - F_interfacejm(k, :))/dx;
    elseif any(xc+1==j)
        dwdt(j,:) = (F_interfacejp(k, :) - F_interface(j+1, :))/dx;
        k = k+1;
    else
        dwdt(j,:) = (F_interface(j,:) - F_interface(j+1, :))/dx;
    end
end
n = n+1;
end