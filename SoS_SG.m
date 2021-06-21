function c = SoS_SG(gamma,rho, p, p0)
% % % % EOS SPECIFIC FUNCTION: STIFFENED GAS % % % % %
% Output: Speed of sound (c)
% Inputs: specific heat ratio (gamma), denisity (rho), internal energy per unit mass (e),
% pressure constant (p0->p0 = 0 for ideal gas EOS)
c = sqrt(gamma .* (p+p0) ./ rho);
end