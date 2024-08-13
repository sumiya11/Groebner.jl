#! format: off
#! source: https://github.com/alexeyovchinnikov/SIAN-Julia

import AbstractAlgebra

function SlowFast(; np=AbstractAlgebra, internal_ordering=:degrevlex, k=np.QQ)
    R, (eC_3,xA_3,xB_3,xC_3,eA_3,eC_2,xA_2,xB_2,xC_2,eA_2,eC_1,xA_1,xB_1,xC_1,eA_1,eC_0,xA_0,xB_0,xC_0,eA_0,z_aux,k1_0,k2_0,eB_0) = np.polynomial_ring(k, [:eC_3,:xA_3,:xB_3,:xC_3,:eA_3,:eC_2,:xA_2,:xB_2,:xC_2,:eA_2,:eC_1,:xA_1,:xB_1,:xC_1,:eA_1,:eC_0,:xA_0,:xB_0,:xC_0,:eA_0,:z_aux,:k1_0,:k2_0,:eB_0], internal_ordering=internal_ordering)
    sys = [
    		-eA_0 + 200461933,
		eA_1,
		-xC_0 + 594762172,
		-xB_0*k2_0 + xC_1,
		-eC_0 + 142222710,
		eC_1,
		-eC_0*xC_0 - xA_0*eA_0 - xB_0*eB_0 + 111787117091448305,
		-xA_0*k1_0 + xB_0*k2_0 + xB_1,
		xA_0*k1_0 + xA_1,
		-xC_1 + 12335904449409052,
		-xB_1*k2_0 + xC_2,
		-xC_1*eC_0 - eA_1*xA_0 - eC_1*xC_0 - xA_1*eA_0 - xB_1*eB_0 - 2178853312350972676271090,
		eC_2,
		-xA_1*k1_0 + xB_1*k2_0 + xB_2,
		eA_2,
		xA_1*k1_0 + xA_2,
		-xC_2 + 2584097121172474714496732,
		-xB_2*k2_0 + xC_3,
		-2*eC_1*xC_1 - 2*xA_1*eA_1 - xC_2*eC_0 - eA_2*xA_0 - eC_2*xC_0 - xA_2*eA_0 - xB_2*eB_0 + 1398478248631803968737625551837100,
		xA_2*k1_0 + xA_3,
		-xA_2*k1_0 + xB_2*k2_0 + xB_3,
		eC_3,
		eA_3,
		-eA_1,
		-eA_2,
		-eA_3,
		-xC_3 - 2740786820954960908989843079767008,
		-eC_1,
		-eC_2,
		-eC_3,
		-3*xC_2*eC_1 - 3*eA_2*xA_1 - 3*eC_2*xC_1 - 3*xA_2*eA_1 - xC_3*eC_0 - eA_3*xA_0 - eC_3*xC_0 - xA_3*eA_0 - xB_3*eB_0 - 878346468468028328602950295328224551413480,
		z_aux - 1
    ]
end
