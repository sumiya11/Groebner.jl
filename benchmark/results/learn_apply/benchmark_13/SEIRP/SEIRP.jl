# SEIRP
#! format: off
using AbstractAlgebra, Groebner

ring, (S_8, I_8, E_7, S_7, I_7, E_6, S_6, I_6, E_5, S_5, I_5, E_4, S_4, I_4, E_3, S_3, I_3, E_2, S_2, I_2, E_1, S_1, I_1, E_0, S_0, I_0, z_aux, R_0, P_0, alpha_e_0, alpha_i_0, kappa_0, rho_0, beta_0, mu_0) = polynomial_ring(
    GF(1073741827), 
    ["S_8", "I_8", "E_7", "S_7", "I_7", "E_6", "S_6", "I_6", "E_5", "S_5", "I_5", "E_4", "S_4", "I_4", "E_3", "S_3", "I_3", "E_2", "S_2", "I_2", "E_1", "S_1", "I_1", "E_0", "S_0", "I_0", "z_aux", "R_0", "P_0", "alpha_e_0", "alpha_i_0", "kappa_0", "rho_0", "beta_0", "mu_0"], 
    ordering=:degrevlex
)
system = [
	1073741826*S_0 + 1073741826*I_0 + 564075071,
	E_0*S_0*alpha_e_0 + S_0*I_0*alpha_i_0 + S_1,
	1073741826*E_0*kappa_0 + I_0*beta_0 + I_0*mu_0 + I_1,
	1073741826*S_1 + 1073741826*I_1 + 343749228,
	S_1*E_0*alpha_e_0 + E_1*S_0*alpha_e_0 + I_1*S_0*alpha_i_0 + S_1*I_0*alpha_i_0 + S_2,
	1073741826*E_1*kappa_0 + I_1*beta_0 + I_1*mu_0 + I_2,
	1073741826*E_0*S_0*alpha_e_0 + 1073741826*S_0*I_0*alpha_i_0 + E_0*kappa_0 + E_0*rho_0 + E_1,
	1073741826*S_2 + 1073741826*I_2 + 748309925,
	1073741826*E_2*kappa_0 + I_2*beta_0 + I_2*mu_0 + I_3,
	2*E_1*S_1*alpha_e_0 + S_2*E_0*alpha_e_0 + E_2*S_0*alpha_e_0 + 2*S_1*I_1*alpha_i_0 + I_2*S_0*alpha_i_0 + S_2*I_0*alpha_i_0 + S_3,
	1073741826*S_1*E_0*alpha_e_0 + 1073741826*E_1*S_0*alpha_e_0 + 1073741826*I_1*S_0*alpha_i_0 + 1073741826*S_1*I_0*alpha_i_0 + E_1*kappa_0 + E_1*rho_0 + E_2,
	1073741826*S_3 + 1073741826*I_3 + 364851067,
	1073741826*E_3*kappa_0 + I_3*beta_0 + I_3*mu_0 + I_4,
	3*S_2*E_1*alpha_e_0 + 3*E_2*S_1*alpha_e_0 + S_3*E_0*alpha_e_0 + E_3*S_0*alpha_e_0 + 3*I_2*S_1*alpha_i_0 + 3*S_2*I_1*alpha_i_0 + I_3*S_0*alpha_i_0 + S_3*I_0*alpha_i_0 + S_4,
	1073741825*E_1*S_1*alpha_e_0 + 1073741826*S_2*E_0*alpha_e_0 + 1073741826*E_2*S_0*alpha_e_0 + 1073741825*S_1*I_1*alpha_i_0 + 1073741826*I_2*S_0*alpha_i_0 + 1073741826*S_2*I_0*alpha_i_0 + E_2*kappa_0 + E_2*rho_0 + E_3,
	1073741826*S_4 + 1073741826*I_4 + 916600629,
	1073741826*E_4*kappa_0 + I_4*beta_0 + I_4*mu_0 + I_5,
	6*E_2*S_2*alpha_e_0 + 4*S_3*E_1*alpha_e_0 + 4*E_3*S_1*alpha_e_0 + S_4*E_0*alpha_e_0 + E_4*S_0*alpha_e_0 + 6*S_2*I_2*alpha_i_0 + 4*I_3*S_1*alpha_i_0 + 4*S_3*I_1*alpha_i_0 + I_4*S_0*alpha_i_0 + S_4*I_0*alpha_i_0 + S_5,
	1073741824*S_2*E_1*alpha_e_0 + 1073741824*E_2*S_1*alpha_e_0 + 1073741826*S_3*E_0*alpha_e_0 + 1073741826*E_3*S_0*alpha_e_0 + 1073741824*I_2*S_1*alpha_i_0 + 1073741824*S_2*I_1*alpha_i_0 + 1073741826*I_3*S_0*alpha_i_0 + 1073741826*S_3*I_0*alpha_i_0 + E_3*kappa_0 + E_3*rho_0 + E_4,
	1073741826*S_5 + 1073741826*I_5 + 971489445,
	10*S_3*E_2*alpha_e_0 + 10*E_3*S_2*alpha_e_0 + 5*S_4*E_1*alpha_e_0 + 5*E_4*S_1*alpha_e_0 + S_5*E_0*alpha_e_0 + E_5*S_0*alpha_e_0 + 10*I_3*S_2*alpha_i_0 + 10*S_3*I_2*alpha_i_0 + 5*I_4*S_1*alpha_i_0 + 5*S_4*I_1*alpha_i_0 + I_5*S_0*alpha_i_0 + S_5*I_0*alpha_i_0 + S_6,
	1073741826*E_5*kappa_0 + I_5*beta_0 + I_5*mu_0 + I_6,
	1073741821*E_2*S_2*alpha_e_0 + 1073741823*S_3*E_1*alpha_e_0 + 1073741823*E_3*S_1*alpha_e_0 + 1073741826*S_4*E_0*alpha_e_0 + 1073741826*E_4*S_0*alpha_e_0 + 1073741821*S_2*I_2*alpha_i_0 + 1073741823*I_3*S_1*alpha_i_0 + 1073741823*S_3*I_1*alpha_i_0 + 1073741826*I_4*S_0*alpha_i_0 + 1073741826*S_4*I_0*alpha_i_0 + E_4*kappa_0 + E_4*rho_0 + E_5,
	1073741826*S_6 + 1073741826*I_6 + 677429293,
	20*E_3*S_3*alpha_e_0 + 15*S_4*E_2*alpha_e_0 + 15*E_4*S_2*alpha_e_0 + 6*S_5*E_1*alpha_e_0 + 6*E_5*S_1*alpha_e_0 + S_6*E_0*alpha_e_0 + E_6*S_0*alpha_e_0 + 20*S_3*I_3*alpha_i_0 + 15*I_4*S_2*alpha_i_0 + 15*S_4*I_2*alpha_i_0 + 6*I_5*S_1*alpha_i_0 + 6*S_5*I_1*alpha_i_0 + I_6*S_0*alpha_i_0 + S_6*I_0*alpha_i_0 + S_7,
	1073741826*E_6*kappa_0 + I_6*beta_0 + I_6*mu_0 + I_7,
	1073741817*S_3*E_2*alpha_e_0 + 1073741817*E_3*S_2*alpha_e_0 + 1073741822*S_4*E_1*alpha_e_0 + 1073741822*E_4*S_1*alpha_e_0 + 1073741826*S_5*E_0*alpha_e_0 + 1073741826*E_5*S_0*alpha_e_0 + 1073741817*I_3*S_2*alpha_i_0 + 1073741817*S_3*I_2*alpha_i_0 + 1073741822*I_4*S_1*alpha_i_0 + 1073741822*S_4*I_1*alpha_i_0 + 1073741826*I_5*S_0*alpha_i_0 + 1073741826*S_5*I_0*alpha_i_0 + E_5*kappa_0 + E_5*rho_0 + E_6,
	1073741826*S_7 + 1073741826*I_7 + 210406170,
	35*S_4*E_3*alpha_e_0 + 35*E_4*S_3*alpha_e_0 + 21*S_5*E_2*alpha_e_0 + 21*E_5*S_2*alpha_e_0 + 7*S_6*E_1*alpha_e_0 + 7*E_6*S_1*alpha_e_0 + S_7*E_0*alpha_e_0 + E_7*S_0*alpha_e_0 + 35*I_4*S_3*alpha_i_0 + 35*S_4*I_3*alpha_i_0 + 21*I_5*S_2*alpha_i_0 + 21*S_5*I_2*alpha_i_0 + 7*I_6*S_1*alpha_i_0 + 7*S_6*I_1*alpha_i_0 + I_7*S_0*alpha_i_0 + S_7*I_0*alpha_i_0 + S_8,
	1073741826*E_7*kappa_0 + I_7*beta_0 + I_7*mu_0 + I_8,
	1073741807*E_3*S_3*alpha_e_0 + 1073741812*S_4*E_2*alpha_e_0 + 1073741812*E_4*S_2*alpha_e_0 + 1073741821*S_5*E_1*alpha_e_0 + 1073741821*E_5*S_1*alpha_e_0 + 1073741826*S_6*E_0*alpha_e_0 + 1073741826*E_6*S_0*alpha_e_0 + 1073741807*S_3*I_3*alpha_i_0 + 1073741812*I_4*S_2*alpha_i_0 + 1073741812*S_4*I_2*alpha_i_0 + 1073741821*I_5*S_1*alpha_i_0 + 1073741821*S_5*I_1*alpha_i_0 + 1073741826*I_6*S_0*alpha_i_0 + 1073741826*S_6*I_0*alpha_i_0 + E_6*kappa_0 + E_6*rho_0 + E_7,
	1073741826*S_8 + 1073741826*I_8 + 86509150,
	z_aux + 1073741826
]

