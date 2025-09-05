#! format: off
#! source: https://github.com/alexeyovchinnikov/SIAN-Julia

import AbstractAlgebra

function SIR_R0_with_weights(; np=AbstractAlgebra, internal_ordering=:degrevlex, k=np.QQ)
    R, (In_4,S_3,In_3,S_2,In_2,aux_1,S_1,In_1,aux_0,S_0,In_0,z_aux,R_0,b_0,g_0) = np.polynomial_ring(k, [:In_4,:S_3,:In_3,:S_2,:In_2,:aux_1,:S_1,:In_1,:aux_0,:S_0,:In_0,:z_aux,:R_0,:b_0,:g_0], internal_ordering=internal_ordering)
    sys = [
    		-In_0 + 61627613,
		-S_0^2*In_0*b_0 + In_0*g_0 + In_1,
		-aux_0*g_0 - b_0 + k(713152382498549)//k(201304577)*g_0,
		aux_1,
		-In_1 + 742538890756902052683515,
		-In_1*S_0^2*b_0 - S_1^2*In_0*b_0 + In_1*g_0 + In_2,
		S_0^2*In_0*b_0 + S_1^2,
		-In_2 + 3506841728018291993077987240185403559981,
		-2*S_1^2*In_1*b_0 - In_2*S_0^2*b_0 - S_2^2*In_0*b_0 + In_2*g_0 + In_3,
		In_1*S_0^2*b_0 + S_1^2*In_0*b_0 + S_2^2,
		-In_3 - 114525628576740143026897344689006444933721316499484210709,
		-3*In_2*S_1^2*b_0 - 3*S_2^2*In_1*b_0 - In_3*S_0^2*b_0 - S_3^2*In_0*b_0 + In_3*g_0 + In_4,
		2*S_1^2*In_1*b_0 + In_2*S_0^2*b_0 + S_2^2*In_0*b_0 + S_3^2,
		-In_4 - 2398169492620307798398605936331379184030222780560926490178936617684529091,
		-aux_1*g_0,
		z_aux*g_0 - 1
    ]
end
