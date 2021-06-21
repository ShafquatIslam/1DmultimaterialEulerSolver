function e = e_C2P(rho, rho_u, E)
% Performing a conservative to primitive operation
% Evaluate the internal energy per unit mass (e)
% Inputs: internal energy per unit density (rho), momentum (rho_u), total
% energy per unit volume (E)
u = rho_u./rho;
e = E./rho - 0.5*u.^2;
end