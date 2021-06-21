function e = e_SG(gamma,rho, p, p0)
% % % % EOS SPECIFIC FUNCTION: STIFFENED GAS % % % % %
% Output: internal energy per unit mass (e)
% Inputs: specific heat ratio (gamma), denisity (rho), pressure (p),
% pressure constant (p0->p0 = 0 for ideal gas EOS)
e = (p + gamma .* p0) ./ (rho .* (gamma-1));
end