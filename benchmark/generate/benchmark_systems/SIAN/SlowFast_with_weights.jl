#! format: off
#! source: https://github.com/alexeyovchinnikov/SIAN-Julia

import AbstractAlgebra

function SlowFast_with_weights(; np=AbstractAlgebra, ordering=:degrevlex, k=np.QQ)
    R, (eC_3,xA_3,xB_3,xC_3,eA_3,eC_2,xA_2,xB_2,xC_2,eA_2,eC_1,xA_1,xB_1,xC_1,eA_1,eC_0,xA_0,xB_0,xC_0,eA_0,z_aux,k1_0,k2_0,eB_0) = np.PolynomialRing(k, [:eC_3,:xA_3,:xB_3,:xC_3,:eA_3,:eC_2,:xA_2,:xB_2,:xC_2,:eA_2,:eC_1,:xA_1,:xB_1,:xC_1,:eA_1,:eC_0,:xA_0,:xB_0,:xC_0,:eA_0,:z_aux,:k1_0,:k2_0,:eB_0], ordering=ordering)
    sys = [
    		-eA_0 + 238590141,
		eA_1,
		-xC_0 + 54769129,
		-xB_0*k2_0^2 + xC_1,
		-eC_0 + 582496235,
		eC_1,
		-eC_0*xC_0 - xA_0*eA_0 - xB_0*eB_0 + 388880808325572233,
		-xA_0*k1_0^2 + xB_0*k2_0^2 + xB_1,
		xA_0*k1_0^2 + xA_1,
		-xC_1 + 227880020762958234,
		-xB_1*k2_0^2 + xC_2,
		-xC_1*eC_0 - eA_1*xA_0 - eC_1*xC_0 - xA_1*eA_0 - xB_1*eB_0 + 27347063390464295018454846,
		eC_2,
		-xA_1*k1_0^2 + xB_1*k2_0^2 + xB_2,
		eA_2,
		xA_1*k1_0^2 + xA_2,
		-xC_2 - 71636259076367508137478078,
		-xB_2*k2_0^2 + xC_3,
		-2*eC_1*xC_1 - 2*xA_1*eA_1 - xC_2*eC_0 - eA_2*xA_0 - eC_2*xC_0 - xA_2*eA_0 - xB_2*eB_0 - 6930051459565503461629070659401882,
		xA_2*k1_0^2 + xA_3,
		-xA_2*k1_0^2 + xB_2*k2_0^2 + xB_3,
		eC_3,
		eA_3,
		-eA_1,
		-eA_2,
		-eA_3,
		-xC_3 + 24859703754859962821505780807346554,
		-eC_1,
		-eC_2,
		-eC_3,
		-3*xC_2*eC_1 - 3*eA_2*xA_1 - 3*eC_2*xC_1 - 3*xA_2*eA_1 - xC_3*eC_0 - eA_3*xA_0 - eC_3*xC_0 - xA_3*eA_0 - xB_3*eB_0 + 1971304241956446054427292096219317446982926,
		z_aux - 1
    ]
end
