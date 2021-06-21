clear all; clc; close all;
global n 


% %Initial Conditions
% gammaIC = [1.4 1.4 1.2 1.2];
% rhoIC = [8 8 1 1];
% uIC = [0 0 0 0];
% pIC = [7 7 1 1];

% % Rallu IC 1
% gammaIC = [1.4 1.4 1.2 1.2];
% rhoIC = [1 1 0.125 0.125];
% uIC = [0 0 0 0];
% pIC = [1 1 0.1 0.1];

% % Rallu IC 2
gammaIC = [1.4 1.4 1.667 1.667];
rhoIC = [1 1 1 1];
uIC = [0 0 0 0];
pIC = [500 500 0.2 0.2];

% Discretize spatial domain
N = 1000;
x_in = 0; x_out = 1;
dx = 1/N;
x = (x_in + dx/2):dx:(x_out-dx/2);
% Temporal parameters
t = 0; dt = 1e-5; 
tf = 0.015; 
nsteps = round(tf/dt);
% nsteps = 100;
% Set interface locations
x_int = [0.25 0.5 0.75];
int_num = length(x_int);
mat_int = (1:int_num+1);

eIC = e_SG(gammaIC, rhoIC, pIC, 0);
EIC = E_P2C(rhoIC, uIC, eIC);

%Setup initial state matrix
x_change = zeros(1, int_num);
rho0 = rhoIC(end)*ones(N,1);
u0 = uIC(end)*ones(N,1);
E0 = EIC(end)*ones(N,1);
p0 = pIC(end)*ones(N,1);
for j = 1:int_num    
    x_change(j) = find(diff(sign(x - x_int(j))));
    xc = [0 x_change];
    rho0(((xc(j)+1):xc(j+1)) ) = rhoIC(j);
    u0(((xc(j)+1):xc(j+1)) ) = uIC(j);
    E0(((xc(j)+1):xc(j+1)) ) = EIC(j);
end
w = [rho0, rho0.*u0, E0];
gamma = gammaIC(end)*ones(N,1);


n = 1; 
W = zeros(3*N, nsteps+1); W(:, 1) = reshape(w, [3*N,1]);
p_mat = zeros(N, nsteps+1); p_mat(:, 1) = p0;
t_vec = zeros(1, nsteps+1);

X_int = zeros(nsteps, int_num); X_int(1,:) = x_int;
mat_id = (int_num+1) * ones(nsteps, N);


while n<=nsteps
    u = w(:,2)./w(:,1);
    e = e_C2P(w(:,1), w(:,2), w(:,3));
    p = p_SG(gamma, w(:,1), e, 0);
    %Interface Tracking Routine 
    x_brack = zeros(2, int_num);
    v_brack = zeros(2, int_num);
    v_int = zeros(1, int_num);
    x_int = X_int(n,:);
    mat_int = (1:int_num+1);
    mat_id(n, :) = mat_int(end)*ones(1,N);

    for j = 1:int_num
        x_int(j) = real(x_int(j));
        x_change(n+1, j) = find(diff(sign(x - x_int(j))));
        x_brack(1,j) = x(x_change(n+1, j)); x_brack(2,j) = x(x_change(n+1, j) +1);
        v_brack(1,j) = u(x_change(n+1, j)); v_brack(2,j) = u(x_change(n+1, j) +1);
        v_int(j) = ((x_int(j) - x_brack(1,j))*(v_brack(2,j) - v_brack(1,j))/(x_brack(2,j) - x_brack(1,j))) + v_brack(1,j);
        xc = [0, x_change(n+1,:)];
        mat_id(n, ((xc(j)+1):xc(j+1)) ) = mat_int(j);
        gamma(((xc(j)+1):xc(j+1)) ) = gammaIC(j);
    end    
    x_int = x_int + v_int*dt;
    X_int(n+1,:) = x_int;
    
%   Zeroth order nodal phase update
    if any(x_change(n+1, :)~=x_change(n, :))
%         disp(['Interface moved a node, state matrix updated at timestep ' num2str(n)])
        dxc = x_change(n+1, :) - x_change(n,:);
        for j = 1:int_num
            if dxc(j)==-1
                if j ==2
                s = [w(x_change(n, j), :); w(x_change(n,j)+1, :)];
                end
                w(x_change(n, j), :) = w(x_change(n,j)+1, :);
                gamma(x_change(n, j)) = gamma(x_change(n,j)+1);
%               disp('moved left')
            elseif dxc(j)==1
                if j ==2
                s = [w(x_change(n+1, j), :); w(x_change(n+1,j)-1, :)];
                end
                w(x_change(n+1, j),:) = w(x_change(n+1,j)-1, :);
                gamma(x_change(n+1, j)) = gamma(x_change(n+1,j)-1);
    %             disp('moved right')
            elseif abs(dxc(j))>1
                disp(['Interface moved more than one node, check at timestep ' num2str(n)])
            end
        end
    end
xc(1) = [];
     % ODE solver: Forward Euler
       dwdt = ddt(w, p, gamma, x, xc);
       w = w + dwdt*dt;
       
       % Save data and update loop
       W(:, n) = reshape(w, [3*N,1]);
       p_mat(:, n) = p;
       t = t+dt;
       t_vec(n) = t;
end

W1 = W(1:N,:);
W2 = W(N+1: 2*N, :);
W3 = W(2*N+1: 3*N, :);
v_mat = W2./W1;

%% Post processing (t = tf)
close all
lw = 2;

% % % Toggle comments to plot different variables % % %
% plot(x, W1(:,end), 'linewidth', lw); a = 1;
% plot(x, v_mat(:,end), 'linewidth', lw); a = 2;
plot(x, p_mat(:,end), 'linewidth', lw); a = 3;


for k = 1:int_num
    xline(X_int(end,k),'--', 'linewidth', lw-1)
end
txt = ['t = ' num2str(t_vec(end))];
text(0.1,0.1,txt, 'Units', 'normalized')
if a==1
    title('Density')
elseif a==2
    title('Velocity')
else
    title ('Pressure')
end
grid on

%% Post processing (animation)
close all
lw = 2;
ff = 20;% X times fast forward
for n = 1:ff:nsteps
% % % Toggle comments to plot different variables % % %
    plot(x, W1(:,n), 'linewidth', lw); a = 1;
%     plot(x, v_mat(:,n), 'linewidth', lw); a = 2;
%     plot(x, p_mat(:,n), 'linewidth', lw); a = 3;
    for k = 1:int_num
        xline(X_int(n,k),'--', 'linewidth', lw-1)
    end
    txt = ['t = ' num2str(t_vec(n))];
    text(0.1,0.1,txt, 'Units', 'normalized')
    if a==1
        title('Density')
    elseif a==2
        title('Velocity')
    elseif a==3
        title ('Pressure')
    end
    drawnow
end