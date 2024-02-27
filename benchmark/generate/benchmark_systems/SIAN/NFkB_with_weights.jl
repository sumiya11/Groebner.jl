#! format: off
#! source: https://github.com/alexeyovchinnikov/SIAN-Julia

import AbstractAlgebra

function NFkB_with_weights(; np=AbstractAlgebra, internal_ordering=:degrevlex, k=np.QQ)
    R, (x2_7,x3_6,x13_6,x10_6,x2_6,x5_6,x1_6,x8_6,x4_6,x3_5,x14_5,x13_5,x11_5,x10_5,x2_5,x12_5,x7_5,x5_5,x1_5,x8_5,x9_5,x6_5,x4_5,x3_4,x14_4,x13_4,x11_4,x10_4,x2_4,x12_4,x7_4,x5_4,x1_4,x8_4,x9_4,x6_4,x4_4,x3_3,x14_3,x13_3,x11_3,x10_3,x2_3,x12_3,x7_3,x5_3,x1_3,x8_3,x9_3,x6_3,x4_3,x3_2,x14_2,x13_2,x11_2,x10_2,x2_2,x12_2,x7_2,x5_2,x1_2,x8_2,x9_2,x6_2,x4_2,x3_1,x14_1,x13_1,x11_1,x10_1,x2_1,x12_1,x7_1,x5_1,x1_1,x8_1,x9_1,x6_1,x4_1,x3_0,x14_0,x13_0,x11_0,x10_0,x2_0,x12_0,x7_0,x5_0,x1_0,x8_0,x9_0,x6_0,x4_0,z_aux,k_prod_0,k_deg_0,k1_0,k3_0,t1_0,t2_0,k2_0,i1_0,c5_0,c_4a_0,i_1a_0,c_3a_0,e_2a_0) = np.PolynomialRing(k, [:x2_7,:x3_6,:x13_6,:x10_6,:x2_6,:x5_6,:x1_6,:x8_6,:x4_6,:x3_5,:x14_5,:x13_5,:x11_5,:x10_5,:x2_5,:x12_5,:x7_5,:x5_5,:x1_5,:x8_5,:x9_5,:x6_5,:x4_5,:x3_4,:x14_4,:x13_4,:x11_4,:x10_4,:x2_4,:x12_4,:x7_4,:x5_4,:x1_4,:x8_4,:x9_4,:x6_4,:x4_4,:x3_3,:x14_3,:x13_3,:x11_3,:x10_3,:x2_3,:x12_3,:x7_3,:x5_3,:x1_3,:x8_3,:x9_3,:x6_3,:x4_3,:x3_2,:x14_2,:x13_2,:x11_2,:x10_2,:x2_2,:x12_2,:x7_2,:x5_2,:x1_2,:x8_2,:x9_2,:x6_2,:x4_2,:x3_1,:x14_1,:x13_1,:x11_1,:x10_1,:x2_1,:x12_1,:x7_1,:x5_1,:x1_1,:x8_1,:x9_1,:x6_1,:x4_1,:x3_0,:x14_0,:x13_0,:x11_0,:x10_0,:x2_0,:x12_0,:x7_0,:x5_0,:x1_0,:x8_0,:x9_0,:x6_0,:x4_0,:z_aux,:k_prod_0,:k_deg_0,:k1_0,:k3_0,:t1_0,:t2_0,:k2_0,:i1_0,:c5_0,:c_4a_0,:i_1a_0,:c_3a_0,:e_2a_0], ordering=ordering)
    sys = [
    		-x3_0 - x2_0 - x1_0 + 2229734195589095673551583400052852,
		-630715481097022642713555616553466*x2_0*x8_0^2*k2_0 + x3_0*k_deg_0 - x2_0*k3_0 + x3_1,
		x1_0*k_deg_0 + 630715481097022642713555616553466*x1_0*k1_0 + x1_1 - k_prod_0,
		630715481097022642713555616553466*x2_0*x8_0^2*k2_0 - x4_0^2*t1_0 - x5_0^2*t2_0 + x13_0*x2_0 + k(1)//k(5)*x10_0*x2_0 + x2_0*k_deg_0 - 630715481097022642713555616553466*x1_0*k1_0 + x2_0*k3_0 + x2_1,
		-x12_0 + 987601179276526308473588913268612,
		x12_0*c_3a_0 + x12_1 - k(1)//k(2000000)*x7_0,
		-x2_0 + 536908731909698060886860269685544,
		-x9_0 + 1512790334540964598447372633631210,
		x9_1 - k(1)//k(2000000)*x7_0 + k(1)//k(2500)*x9_0,
		-x13_0 - x10_0 + 3060431816386813809911722164007348,
		-k(1)//k(2)*x10_0*x6_0^2 - x14_0^2*e_2a_0 + x13_0*x2_0 + x13_1 + k(1)//k(50000)*x13_0,
		k(1)//k(2)*x10_0*x6_0^2 - k(1)//k(2000)*x11_0^2 + k(1)//k(5)*x10_0*x2_0 - x12_0*c_4a_0 + x10_0*i_1a_0 + x10_1 + k(1)//k(10000)*x10_0,
		-x7_0 + 1541225865727688393785881037468687,
		k(1)//k(2)*x11_0^2*x7_0 - 5*x6_0^2*i1_0 + x7_1,
		-x3_1 - x2_1 - x1_1 + k(3857725759501414958956809419192070462155540509833775649770328393089)//k(5),
		x1_1*k_deg_0 + 630715481097022642713555616553466*x1_1*k1_0 + 1457708186765988809605264969354780*x1_0*k1_0 + x1_2,
		-630715481097022642713555616553466*x8_1^2*x2_0*k2_0 - 630715481097022642713555616553466*x2_1*x8_0^2*k2_0 - 1457708186765988809605264969354780*x2_0*x8_0^2*k2_0 + x3_1*k_deg_0 - x2_1*k3_0 + x3_2,
		630715481097022642713555616553466*x8_1^2*x2_0*k2_0 + 630715481097022642713555616553466*x2_1*x8_0^2*k2_0 + 1457708186765988809605264969354780*x2_0*x8_0^2*k2_0 - x4_1^2*t1_0 - x5_1^2*t2_0 + x2_1*x13_0 + k(1)//k(5)*x2_1*x10_0 + x13_1*x2_0 + k(1)//k(5)*x10_1*x2_0 + x2_1*k_deg_0 - 630715481097022642713555616553466*x1_1*k1_0 - 1457708186765988809605264969354780*x1_0*k1_0 + x2_1*k3_0 + x2_2,
		x4_0^2*t1_0 + x4_1^2 - k(1)//k(5)*x10_0*x2_0,
		x8_0^2*c5_0^3 + x8_1^2 - k(1)//k(2)*x9_0,
		x5_0^2*t2_0 + x5_1^2 - x13_0*x2_0,
		-x12_1 - k(2317620428069940318620746407112113600862443183503040438523361812042531313)//k(2000000),
		x12_1*c_3a_0 + x12_2 - k(1)//k(2000000)*x7_1,
		-x2_1 - k(666481218950821705084935524900992987667864399617002842083444560602994731303517846199636141667530994707143813446498374067122774090941)//k(5),
		-x13_1 - x10_1 - k(17397704328756113564134937847821824478043322995421960041360799281309059)//k(50000),
		-k(1)//k(2)*x6_1^2*x10_0 - k(1)//k(2)*x10_1*x6_0^2 - x14_1^2*e_2a_0 + x2_1*x13_0 + x13_1*x2_0 + x13_2 + k(1)//k(50000)*x13_1,
		k(1)//k(2)*x6_1^2*x10_0 + k(1)//k(2)*x10_1*x6_0^2 - k(1)//k(2000)*x11_1^2 + k(1)//k(5)*x2_1*x10_0 + k(1)//k(5)*x10_1*x2_0 - x12_1*c_4a_0 + x10_1*i_1a_0 + x10_2 + k(1)//k(10000)*x10_1,
		k(1)//k(2)*x10_0*x6_0^2 - x5_0^2*t2_0 + x6_0^2*i1_0 + x6_1^2 - k(1)//k(50000)*x13_0,
		k(1)//k(2)*x11_0^2*x7_0 + x11_1^2 + k(1)//k(400)*x11_0^2 - 5*x10_0*i_1a_0,
		-k(1)//k(2)*x11_0^2*x7_0 + 5*x14_0^2*e_2a_0 + x14_1^2,
		-x7_1 + k(13131052282828470497734952588949739257364463607936031516107869475833)//k(2),
		k(1)//k(2)*x7_1*x11_0^2 + k(1)//k(2)*x11_1^2*x7_0 - 5*x6_1^2*i1_0 + x7_2,
		-x3_2 - x2_2 - x1_2 + k(1681354817594246375258595073460506059893760142107898351557059598523886373382042584424956704852906174439778450181997237944292292627373003824702315154583439502626492907149)//k(6250),
		-1261430962194045285427111233106932*x2_1*x8_1^2*k2_0 - 630715481097022642713555616553466*x8_2^2*x2_0*k2_0 - 2915416373531977619210529938709560*x8_1^2*x2_0*k2_0 - 630715481097022642713555616553466*x2_2*x8_0^2*k2_0 - 2915416373531977619210529938709560*x2_1*x8_0^2*k2_0 - 45461424441375244328735581629804*x2_0*x8_0^2*k2_0 + x3_2*k_deg_0 - x2_2*k3_0 + x3_3,
		x1_2*k_deg_0 + 630715481097022642713555616553466*x1_2*k1_0 + 2915416373531977619210529938709560*x1_1*k1_0 + 45461424441375244328735581629804*x1_0*k1_0 + x1_3,
		1261430962194045285427111233106932*x2_1*x8_1^2*k2_0 + 630715481097022642713555616553466*x8_2^2*x2_0*k2_0 + 2915416373531977619210529938709560*x8_1^2*x2_0*k2_0 + 630715481097022642713555616553466*x2_2*x8_0^2*k2_0 + 2915416373531977619210529938709560*x2_1*x8_0^2*k2_0 + 45461424441375244328735581629804*x2_0*x8_0^2*k2_0 - x4_2^2*t1_0 - x5_2^2*t2_0 + 2*x13_1*x2_1 + k(2)//k(5)*x10_1*x2_1 + x2_2*x13_0 + k(1)//k(5)*x2_2*x10_0 + x13_2*x2_0 + k(1)//k(5)*x10_2*x2_0 + x2_2*k_deg_0 - 630715481097022642713555616553466*x1_2*k1_0 - 2915416373531977619210529938709560*x1_1*k1_0 - 45461424441375244328735581629804*x1_0*k1_0 + x2_2*k3_0 + x2_3,
		x4_1^2*t1_0 + x4_2^2 - k(1)//k(5)*x2_1*x10_0 - k(1)//k(5)*x10_1*x2_0,
		x5_1^2*t2_0 + x5_2^2 - x2_1*x13_0 - x13_1*x2_0,
		x8_1^2*c5_0^3 + x8_2^2 - k(1)//k(2)*x9_1,
		-x2_2 + k(206830876839269797023891266124144420029321262999849276543635972304257088540553837841129974811819303530979092801290784233355529751146625901026856897817861177386008256606783469427010613229918819172668025821838495191290672271432875382149)//k(6250),
		-x13_2 - x10_2 + k(336270963518849275051719014692101211978752028421579670311411919703742144313308810515178846500496204305450337034845033057945611295599681080841549168853808249924356170525995007)//k(1250000000),
		x10_1*x6_1^2 + k(1)//k(2)*x6_2^2*x10_0 + k(1)//k(2)*x10_2*x6_0^2 - k(1)//k(2000)*x11_2^2 + k(2)//k(5)*x10_1*x2_1 + k(1)//k(5)*x2_2*x10_0 + k(1)//k(5)*x10_2*x2_0 - x12_2*c_4a_0 + x10_2*i_1a_0 + x10_3 + k(1)//k(10000)*x10_2,
		-x10_1*x6_1^2 - k(1)//k(2)*x6_2^2*x10_0 - k(1)//k(2)*x10_2*x6_0^2 - x14_2^2*e_2a_0 + 2*x13_1*x2_1 + x2_2*x13_0 + x13_2*x2_0 + x13_3 + k(1)//k(50000)*x13_2,
		k(1)//k(2)*x6_1^2*x10_0 + k(1)//k(2)*x10_1*x6_0^2 - x5_1^2*t2_0 + x6_1^2*i1_0 + x6_2^2 - k(1)//k(50000)*x13_1,
		-k(1)//k(2)*x7_1*x11_0^2 - k(1)//k(2)*x11_1^2*x7_0 + 5*x14_1^2*e_2a_0 + x14_2^2,
		k(1)//k(2)*x7_1*x11_0^2 + k(1)//k(2)*x11_1^2*x7_0 + x11_2^2 + k(1)//k(400)*x11_1^2 - 5*x10_1*i_1a_0,
		-x7_2 - k(349518000993233712441686619013895343659429930634749016545363597772048523353382151414810570801548877002853)//k(20000),
		x11_1^2*x7_1 + k(1)//k(2)*x7_2*x11_0^2 + k(1)//k(2)*x11_2^2*x7_0 - 5*x6_2^2*i1_0 + x7_3,
		-x3_3 - x2_3 - x1_3 - k(10435585619303361711139234620368512483200110229785711734207725118702857427331929504598475062022396694222611985585684218988265089429847379761186555988084084149563158609572572753522438397276809399503685724106767070759217795472315547373247146189534304049646340371769556708881)//k(156250000),
		1892146443291067928140666849660398*x8_2^2*x2_1*k2_0 + 1892146443291067928140666849660398*x2_2*x8_1^2*k2_0 + 8746249120595932857631589816128680*x2_1*x8_1^2*k2_0 + 630715481097022642713555616553466*x8_3^2*x2_0*k2_0 + 4373124560297966428815794908064340*x8_2^2*x2_0*k2_0 + 136384273324125732986206744889412*x8_1^2*x2_0*k2_0 + 630715481097022642713555616553466*x2_3*x8_0^2*k2_0 + 4373124560297966428815794908064340*x2_2*x8_0^2*k2_0 + 136384273324125732986206744889412*x2_1*x8_0^2*k2_0 + 1758988856576644643919915561812650*x2_0*x8_0^2*k2_0 - x4_3^2*t1_0 - x5_3^2*t2_0 + 3*x2_2*x13_1 + k(3)//k(5)*x2_2*x10_1 + 3*x13_2*x2_1 + k(3)//k(5)*x10_2*x2_1 + x2_3*x13_0 + k(1)//k(5)*x2_3*x10_0 + x13_3*x2_0 + k(1)//k(5)*x10_3*x2_0 + x2_3*k_deg_0 - 630715481097022642713555616553466*x1_3*k1_0 - 4373124560297966428815794908064340*x1_2*k1_0 - 136384273324125732986206744889412*x1_1*k1_0 - 1758988856576644643919915561812650*x1_0*k1_0 + x2_3*k3_0 + x2_4,
		x1_3*k_deg_0 + 630715481097022642713555616553466*x1_3*k1_0 + 4373124560297966428815794908064340*x1_2*k1_0 + 136384273324125732986206744889412*x1_1*k1_0 + 1758988856576644643919915561812650*x1_0*k1_0 + x1_4,
		-1892146443291067928140666849660398*x8_2^2*x2_1*k2_0 - 1892146443291067928140666849660398*x2_2*x8_1^2*k2_0 - 8746249120595932857631589816128680*x2_1*x8_1^2*k2_0 - 630715481097022642713555616553466*x8_3^2*x2_0*k2_0 - 4373124560297966428815794908064340*x8_2^2*x2_0*k2_0 - 136384273324125732986206744889412*x8_1^2*x2_0*k2_0 - 630715481097022642713555616553466*x2_3*x8_0^2*k2_0 - 4373124560297966428815794908064340*x2_2*x8_0^2*k2_0 - 136384273324125732986206744889412*x2_1*x8_0^2*k2_0 - 1758988856576644643919915561812650*x2_0*x8_0^2*k2_0 + x3_3*k_deg_0 - x2_3*k3_0 + x3_4,
		x5_2^2*t2_0 + x5_3^2 - 2*x13_1*x2_1 - x2_2*x13_0 - x13_2*x2_0,
		x4_2^2*t1_0 + x4_3^2 - k(2)//k(5)*x10_1*x2_1 - k(1)//k(5)*x2_2*x10_0 - k(1)//k(5)*x10_2*x2_0,
		x8_2^2*c5_0^3 + x8_3^2 - k(1)//k(2)*x9_2,
		x9_2 - k(1)//k(2000000)*x7_1 + k(1)//k(2500)*x9_1,
		-x2_3 - k(1283727444906673887408376282422320619432104020668905676258724324407550829899329773845073090162139350145047628929223430522570856656854177121911389224862486170278757995113284975818509009253389556396985036844051697094085013595376830229107604198821457252061345600680005926437350000587258053902291065523691038772156419111640855300691819973881)//k(156250000),
		-x13_3 - x10_3 - k(8348468495442689368911387696294809986560088183828569387366180094931380930208989624024493475936468784096118308336818215150228063514137109304930245982837970002638689519045988978906409437744593849610338764419269515456350062809106516815621033131859215759291468227173119879323389789)//k(125000000000000),
		-k(3)//k(2)*x6_2^2*x10_1 - k(3)//k(2)*x10_2*x6_1^2 - k(1)//k(2)*x6_3^2*x10_0 - k(1)//k(2)*x10_3*x6_0^2 - x14_3^2*e_2a_0 + 3*x2_2*x13_1 + 3*x13_2*x2_1 + x2_3*x13_0 + x13_3*x2_0 + x13_4 + k(1)//k(50000)*x13_3,
		k(3)//k(2)*x6_2^2*x10_1 + k(3)//k(2)*x10_2*x6_1^2 + k(1)//k(2)*x6_3^2*x10_0 + k(1)//k(2)*x10_3*x6_0^2 - k(1)//k(2000)*x11_3^2 + k(3)//k(5)*x2_2*x10_1 + k(3)//k(5)*x10_2*x2_1 + k(1)//k(5)*x2_3*x10_0 + k(1)//k(5)*x10_3*x2_0 - x12_3*c_4a_0 + x10_3*i_1a_0 + x10_4 + k(1)//k(10000)*x10_3,
		x12_2*c_3a_0 + x12_3 - k(1)//k(2000000)*x7_2,
		-x11_1^2*x7_1 - k(1)//k(2)*x7_2*x11_0^2 - k(1)//k(2)*x11_2^2*x7_0 + 5*x14_2^2*e_2a_0 + x14_3^2,
		x10_1*x6_1^2 + k(1)//k(2)*x6_2^2*x10_0 + k(1)//k(2)*x10_2*x6_0^2 - x5_2^2*t2_0 + x6_2^2*i1_0 + x6_3^2 - k(1)//k(50000)*x13_2,
		x11_1^2*x7_1 + k(1)//k(2)*x7_2*x11_0^2 + k(1)//k(2)*x11_2^2*x7_0 + x11_3^2 + k(1)//k(400)*x11_2^2 - 5*x10_2*i_1a_0,
		-x7_3 + k(29022408155657476731624720212530331249928387815377397320529007642796296411152947625864459832951114995588737430261984092944522674216821234244153)//k(1000000000),
		k(3)//k(2)*x7_2*x11_1^2 + k(3)//k(2)*x11_2^2*x7_1 + k(1)//k(2)*x7_3*x11_0^2 + k(1)//k(2)*x11_3^2*x7_0 - 5*x6_3^2*i1_0 + x7_4,
		-x3_4 - x2_4 - x1_4 + k(259080227631266853148864360116307216875251128941799161540130582569659266597064841192837083821965588979737263156777645829068823969859570452250081426694709655928633907057663209817393612727315285021874599372074043394056530852155442254227244782901461130235144330800426333752867660503409942401903872974518625609152036251693735392482563340940191746407593880598504536129257548071577)//k(15625000000000),
		-3784292886582135856281333699320796*x2_2*x8_2^2*k2_0 - 2522861924388090570854222466213864*x8_3^2*x2_1*k2_0 - 17492498241191865715263179632257360*x8_2^2*x2_1*k2_0 - 2522861924388090570854222466213864*x2_3*x8_1^2*k2_0 - 17492498241191865715263179632257360*x2_2*x8_1^2*k2_0 - 545537093296502931944826979557648*x2_1*x8_1^2*k2_0 - 630715481097022642713555616553466*x8_4^2*x2_0*k2_0 - 5830832747063955238421059877419120*x8_3^2*x2_0*k2_0 - 272768546648251465972413489778824*x8_2^2*x2_0*k2_0 - 7035955426306578575679662247250600*x8_1^2*x2_0*k2_0 - 630715481097022642713555616553466*x2_4*x8_0^2*k2_0 - 5830832747063955238421059877419120*x2_3*x8_0^2*k2_0 - 272768546648251465972413489778824*x2_2*x8_0^2*k2_0 - 7035955426306578575679662247250600*x2_1*x8_0^2*k2_0 - 1592602613659573575745658548294093*x2_0*x8_0^2*k2_0 + x3_4*k_deg_0 - x2_4*k3_0 + x3_5,
		3784292886582135856281333699320796*x2_2*x8_2^2*k2_0 + 2522861924388090570854222466213864*x8_3^2*x2_1*k2_0 + 17492498241191865715263179632257360*x8_2^2*x2_1*k2_0 + 2522861924388090570854222466213864*x2_3*x8_1^2*k2_0 + 17492498241191865715263179632257360*x2_2*x8_1^2*k2_0 + 545537093296502931944826979557648*x2_1*x8_1^2*k2_0 + 630715481097022642713555616553466*x8_4^2*x2_0*k2_0 + 5830832747063955238421059877419120*x8_3^2*x2_0*k2_0 + 272768546648251465972413489778824*x8_2^2*x2_0*k2_0 + 7035955426306578575679662247250600*x8_1^2*x2_0*k2_0 + 630715481097022642713555616553466*x2_4*x8_0^2*k2_0 + 5830832747063955238421059877419120*x2_3*x8_0^2*k2_0 + 272768546648251465972413489778824*x2_2*x8_0^2*k2_0 + 7035955426306578575679662247250600*x2_1*x8_0^2*k2_0 + 1592602613659573575745658548294093*x2_0*x8_0^2*k2_0 - x4_4^2*t1_0 - x5_4^2*t2_0 + 6*x13_2*x2_2 + k(6)//k(5)*x10_2*x2_2 + 4*x2_3*x13_1 + k(4)//k(5)*x2_3*x10_1 + 4*x13_3*x2_1 + k(4)//k(5)*x10_3*x2_1 + x2_4*x13_0 + k(1)//k(5)*x2_4*x10_0 + x13_4*x2_0 + k(1)//k(5)*x10_4*x2_0 + x2_4*k_deg_0 - 630715481097022642713555616553466*x1_4*k1_0 - 5830832747063955238421059877419120*x1_3*k1_0 - 272768546648251465972413489778824*x1_2*k1_0 - 7035955426306578575679662247250600*x1_1*k1_0 - 1592602613659573575745658548294093*x1_0*k1_0 + x2_4*k3_0 + x2_5,
		x1_4*k_deg_0 + 630715481097022642713555616553466*x1_4*k1_0 + 5830832747063955238421059877419120*x1_3*k1_0 + 272768546648251465972413489778824*x1_2*k1_0 + 7035955426306578575679662247250600*x1_1*k1_0 + 1592602613659573575745658548294093*x1_0*k1_0 + x1_5,
		x5_3^2*t2_0 + x5_4^2 - 3*x2_2*x13_1 - 3*x13_2*x2_1 - x2_3*x13_0 - x13_3*x2_0,
		x4_3^2*t1_0 + x4_4^2 - k(3)//k(5)*x2_2*x10_1 - k(3)//k(5)*x10_2*x2_1 - k(1)//k(5)*x2_3*x10_0 - k(1)//k(5)*x10_3*x2_0,
		x8_3^2*c5_0^3 + x8_4^2 - k(1)//k(2)*x9_3,
		x9_3 - k(1)//k(2000000)*x7_2 + k(1)//k(2500)*x9_2,
		-x2_4 + k(31870602261910044185570276510195904805905789036208014598954956309624678863481937769761699231473405786989179325593221343849661428785358048726899449469164078341961383552011759481858885538149460607382647739586490713009027659007979823886823175505836005644024281895792664496579910457437284289483459902424773313957934263903473886756872877841847421084069946651683972527849564114079526074677416050250960182288443817894790882204738239597531860171577)//k(15625000000000),
		-x13_4 - x10_4 + k(25908022763126685314886436011630721687525112894179916154013058256870018320838129983852811756362654426780423444695378492849900165462857207377645905341364495513976765777897749829689051049994978022074551673117877414127824049151195769314283706529868793233239563785737843825421095236722400573610763347332391313147445636418936190308606968010129532564206794535370591304961011758795756541)//k(1562500000000000000),
		3*x10_2*x6_2^2 + 2*x6_3^2*x10_1 + 2*x10_3*x6_1^2 + k(1)//k(2)*x6_4^2*x10_0 + k(1)//k(2)*x10_4*x6_0^2 - k(1)//k(2000)*x11_4^2 + k(6)//k(5)*x10_2*x2_2 + k(4)//k(5)*x2_3*x10_1 + k(4)//k(5)*x10_3*x2_1 + k(1)//k(5)*x2_4*x10_0 + k(1)//k(5)*x10_4*x2_0 - x12_4*c_4a_0 + x10_4*i_1a_0 + x10_5 + k(1)//k(10000)*x10_4,
		-3*x10_2*x6_2^2 - 2*x6_3^2*x10_1 - 2*x10_3*x6_1^2 - k(1)//k(2)*x6_4^2*x10_0 - k(1)//k(2)*x10_4*x6_0^2 - x14_4^2*e_2a_0 + 6*x13_2*x2_2 + 4*x2_3*x13_1 + 4*x13_3*x2_1 + x2_4*x13_0 + x13_4*x2_0 + x13_5 + k(1)//k(50000)*x13_4,
		k(3)//k(2)*x6_2^2*x10_1 + k(3)//k(2)*x10_2*x6_1^2 + k(1)//k(2)*x6_3^2*x10_0 + k(1)//k(2)*x10_3*x6_0^2 - x5_3^2*t2_0 + x6_3^2*i1_0 + x6_4^2 - k(1)//k(50000)*x13_3,
		k(3)//k(2)*x7_2*x11_1^2 + k(3)//k(2)*x11_2^2*x7_1 + k(1)//k(2)*x7_3*x11_0^2 + k(1)//k(2)*x11_3^2*x7_0 + x11_4^2 + k(1)//k(400)*x11_3^2 - 5*x10_3*i_1a_0,
		x12_3*c_3a_0 + x12_4 - k(1)//k(2000000)*x7_3,
		-k(3)//k(2)*x7_2*x11_1^2 - k(3)//k(2)*x11_2^2*x7_1 - k(1)//k(2)*x7_3*x11_0^2 - k(1)//k(2)*x11_3^2*x7_0 + 5*x14_3^2*e_2a_0 + x14_4^2,
		-x7_4 - k(67974751019167478546295732078634114501575714235662183439355518940304144061155949729714481899127662570347021215865363302028764544243121769155060004231563107109571071365101997728193061889477304425896919826183684561631980387652204154583693099806653)//k(50000000000000),
		3*x11_2^2*x7_2 + 2*x7_3*x11_1^2 + 2*x11_3^2*x7_1 + k(1)//k(2)*x7_4*x11_0^2 + k(1)//k(2)*x11_4^2*x7_0 - 5*x6_4^2*i1_0 + x7_5,
		-x3_5 - x2_5 - x1_5 - k(804010512660020280244093676132037245100918987831826790436779180532011295308883740093256619198873648841569883586794508325854600084620529274135108529936131552951154624821246075494304676225992378805748382830107456361746296099561916478160307244644132483941192258878947158233565962379729392748224860617453322120339085782712804482080227090853941297234412904285832003923016872794634113706141655702095901338372651063549218376214226502518770067946302492418493032984338648056152964286413)//k(195312500000000000),
		6307154810970226427135556165534660*x8_3^2*x2_2*k2_0 + 6307154810970226427135556165534660*x2_3*x8_2^2*k2_0 + 43731245602979664288157949080643400*x2_2*x8_2^2*k2_0 + 3153577405485113213567778082767330*x8_4^2*x2_1*k2_0 + 29154163735319776192105299387095600*x8_3^2*x2_1*k2_0 + 1363842733241257329862067448894120*x8_2^2*x2_1*k2_0 + 3153577405485113213567778082767330*x2_4*x8_1^2*k2_0 + 29154163735319776192105299387095600*x2_3*x8_1^2*k2_0 + 1363842733241257329862067448894120*x2_2*x8_1^2*k2_0 + 35179777131532892878398311236253000*x2_1*x8_1^2*k2_0 + 630715481097022642713555616553466*x8_5^2*x2_0*k2_0 + 7288540933829944048026324846773900*x8_4^2*x2_0*k2_0 + 454614244413752443287355816298040*x8_3^2*x2_0*k2_0 + 17589888565766446439199155618126500*x8_2^2*x2_0*k2_0 + 7963013068297867878728292741470465*x8_1^2*x2_0*k2_0 + 630715481097022642713555616553466*x2_5*x8_0^2*k2_0 + 7288540933829944048026324846773900*x2_4*x8_0^2*k2_0 + 454614244413752443287355816298040*x2_3*x8_0^2*k2_0 + 17589888565766446439199155618126500*x2_2*x8_0^2*k2_0 + 7963013068297867878728292741470465*x2_1*x8_0^2*k2_0 + 1212419864798925874339613047207811*x2_0*x8_0^2*k2_0 - x4_5^2*t1_0 - x5_5^2*t2_0 + 10*x2_3*x13_2 + 2*x2_3*x10_2 + 10*x13_3*x2_2 + 2*x10_3*x2_2 + 5*x2_4*x13_1 + x2_4*x10_1 + 5*x13_4*x2_1 + x10_4*x2_1 + x2_5*x13_0 + k(1)//k(5)*x2_5*x10_0 + x13_5*x2_0 + k(1)//k(5)*x10_5*x2_0 + x2_5*k_deg_0 - 630715481097022642713555616553466*x1_5*k1_0 - 7288540933829944048026324846773900*x1_4*k1_0 - 454614244413752443287355816298040*x1_3*k1_0 - 17589888565766446439199155618126500*x1_2*k1_0 - 7963013068297867878728292741470465*x1_1*k1_0 - 1212419864798925874339613047207811*x1_0*k1_0 + x2_5*k3_0 + x2_6,
		x1_5*k_deg_0 + 630715481097022642713555616553466*x1_5*k1_0 + 7288540933829944048026324846773900*x1_4*k1_0 + 454614244413752443287355816298040*x1_3*k1_0 + 17589888565766446439199155618126500*x1_2*k1_0 + 7963013068297867878728292741470465*x1_1*k1_0 + 1212419864798925874339613047207811*x1_0*k1_0 + x1_6,
		-6307154810970226427135556165534660*x8_3^2*x2_2*k2_0 - 6307154810970226427135556165534660*x2_3*x8_2^2*k2_0 - 43731245602979664288157949080643400*x2_2*x8_2^2*k2_0 - 3153577405485113213567778082767330*x8_4^2*x2_1*k2_0 - 29154163735319776192105299387095600*x8_3^2*x2_1*k2_0 - 1363842733241257329862067448894120*x8_2^2*x2_1*k2_0 - 3153577405485113213567778082767330*x2_4*x8_1^2*k2_0 - 29154163735319776192105299387095600*x2_3*x8_1^2*k2_0 - 1363842733241257329862067448894120*x2_2*x8_1^2*k2_0 - 35179777131532892878398311236253000*x2_1*x8_1^2*k2_0 - 630715481097022642713555616553466*x8_5^2*x2_0*k2_0 - 7288540933829944048026324846773900*x8_4^2*x2_0*k2_0 - 454614244413752443287355816298040*x8_3^2*x2_0*k2_0 - 17589888565766446439199155618126500*x8_2^2*x2_0*k2_0 - 7963013068297867878728292741470465*x8_1^2*x2_0*k2_0 - 630715481097022642713555616553466*x2_5*x8_0^2*k2_0 - 7288540933829944048026324846773900*x2_4*x8_0^2*k2_0 - 454614244413752443287355816298040*x2_3*x8_0^2*k2_0 - 17589888565766446439199155618126500*x2_2*x8_0^2*k2_0 - 7963013068297867878728292741470465*x2_1*x8_0^2*k2_0 - 1212419864798925874339613047207811*x2_0*x8_0^2*k2_0 + x3_5*k_deg_0 - x2_5*k3_0 + x3_6,
		x4_4^2*t1_0 + x4_5^2 - k(6)//k(5)*x10_2*x2_2 - k(4)//k(5)*x2_3*x10_1 - k(4)//k(5)*x10_3*x2_1 - k(1)//k(5)*x2_4*x10_0 - k(1)//k(5)*x10_4*x2_0,
		x5_4^2*t2_0 + x5_5^2 - 6*x13_2*x2_2 - 4*x2_3*x13_1 - 4*x13_3*x2_1 - x2_4*x13_0 - x13_4*x2_0,
		x8_4^2*c5_0^3 + x8_5^2 - k(1)//k(2)*x9_4,
		x9_4 - k(1)//k(2000000)*x7_3 + k(1)//k(2500)*x9_3,
		-x2_5 - k(98904881694991424920200257492069182940086822763502057795774536588020508345202090063156180350606795102398333428649252184344471944683029832104481745811475689068001949914769971320455601158334918811088843058146272774055677671859862454732460720923478836326564821691788883022078653513840199403848051111601471656949070963117125450363805847084978232487187734968922517494279430861581971271682955817020149823572028679074267860213076802831601816180729805725593377487886447355553271506677286468664382749049630048543056435103052922210499965302580569536413)//k(195312500000000000),
		-x13_5 - x10_5 - k(1286416820256032448390549881811259592161470380530922864698846688846455914399360727868758976171156110401487671993854780918936812966575538119549140961325960945758241069205814755824365085931117835620822716306557357843853057548064892234693038742403534518237754037821884792394234015399207731963037216009678841485459379489970738935712500495949770807626625407785212918988761338843970705182386966865652902613217451110078991773226940509739997205932774475116060409043973828391192443015544845539)//k(312500000000000000000000),
		-5*x6_3^2*x10_2 - 5*x10_3*x6_2^2 - k(5)//k(2)*x6_4^2*x10_1 - k(5)//k(2)*x10_4*x6_1^2 - k(1)//k(2)*x6_5^2*x10_0 - k(1)//k(2)*x10_5*x6_0^2 - x14_5^2*e_2a_0 + 10*x2_3*x13_2 + 10*x13_3*x2_2 + 5*x2_4*x13_1 + 5*x13_4*x2_1 + x2_5*x13_0 + x13_5*x2_0 + x13_6 + k(1)//k(50000)*x13_5,
		5*x6_3^2*x10_2 + 5*x10_3*x6_2^2 + k(5)//k(2)*x6_4^2*x10_1 + k(5)//k(2)*x10_4*x6_1^2 + k(1)//k(2)*x6_5^2*x10_0 + k(1)//k(2)*x10_5*x6_0^2 - k(1)//k(2000)*x11_5^2 + 2*x2_3*x10_2 + 2*x10_3*x2_2 + x2_4*x10_1 + x10_4*x2_1 + k(1)//k(5)*x2_5*x10_0 + k(1)//k(5)*x10_5*x2_0 - x12_5*c_4a_0 + x10_5*i_1a_0 + x10_6 + k(1)//k(10000)*x10_5,
		3*x11_2^2*x7_2 + 2*x7_3*x11_1^2 + 2*x11_3^2*x7_1 + k(1)//k(2)*x7_4*x11_0^2 + k(1)//k(2)*x11_4^2*x7_0 + x11_5^2 + k(1)//k(400)*x11_4^2 - 5*x10_4*i_1a_0,
		-3*x11_2^2*x7_2 - 2*x7_3*x11_1^2 - 2*x11_3^2*x7_1 - k(1)//k(2)*x7_4*x11_0^2 - k(1)//k(2)*x11_4^2*x7_0 + 5*x14_4^2*e_2a_0 + x14_5^2,
		x12_4*c_3a_0 + x12_5 - k(1)//k(2000000)*x7_4,
		3*x10_2*x6_2^2 + 2*x6_3^2*x10_1 + 2*x10_3*x6_1^2 + k(1)//k(2)*x6_4^2*x10_0 + k(1)//k(2)*x10_4*x6_0^2 - x5_4^2*t2_0 + x6_4^2*i1_0 + x6_5^2 - k(1)//k(50000)*x13_4,
		-x2_6 + k(4910946102724319002785691532667847432895674410292330434007771871580131237876163120480328160985793792846598822599205542107985153250693041334748987858704267778556666461305204645103804170125653518143586958677875760782424743545599572651484001665635540725012397862435516605144513186886648815976384045338712633415980353775633802929360036360377688214237443837750157141454517900266308239909873746505729867621418841343610896828935041469076090775565201268489550918884982414514124436350870255622513543202164007422710423288357790614574521006252210920836277721014575290633126113356507242979515612778400567394404802834225361126599032168057735170275528726171027)//k(39062500000000000000000),
		12614309621940452854271112331069320*x2_3*x8_3^2*k2_0 + 9460732216455339640703334248301990*x8_4^2*x2_2*k2_0 + 87462491205959328576315898161286800*x8_3^2*x2_2*k2_0 + 9460732216455339640703334248301990*x2_4*x8_2^2*k2_0 + 87462491205959328576315898161286800*x2_3*x8_2^2*k2_0 + 4091528199723771989586202346682360*x2_2*x8_2^2*k2_0 + 3784292886582135856281333699320796*x8_5^2*x2_1*k2_0 + 43731245602979664288157949080643400*x8_4^2*x2_1*k2_0 + 2727685466482514659724134897788240*x8_3^2*x2_1*k2_0 + 105539331394598678635194933708759000*x8_2^2*x2_1*k2_0 + 3784292886582135856281333699320796*x2_5*x8_1^2*k2_0 + 43731245602979664288157949080643400*x2_4*x8_1^2*k2_0 + 2727685466482514659724134897788240*x2_3*x8_1^2*k2_0 + 105539331394598678635194933708759000*x2_2*x8_1^2*k2_0 + 47778078409787207272369756448822790*x2_1*x8_1^2*k2_0 + 630715481097022642713555616553466*x8_6^2*x2_0*k2_0 + 8746249120595932857631589816128680*x8_5^2*x2_0*k2_0 + 681921366620628664931033724447060*x8_4^2*x2_0*k2_0 + 35179777131532892878398311236253000*x8_3^2*x2_0*k2_0 + 23889039204893603636184878224411395*x8_2^2*x2_0*k2_0 + 7274519188793555246037678283246866*x8_1^2*x2_0*k2_0 + 630715481097022642713555616553466*x2_6*x8_0^2*k2_0 + 8746249120595932857631589816128680*x2_5*x8_0^2*k2_0 + 681921366620628664931033724447060*x2_4*x8_0^2*k2_0 + 35179777131532892878398311236253000*x2_3*x8_0^2*k2_0 + 23889039204893603636184878224411395*x2_2*x8_0^2*k2_0 + 7274519188793555246037678283246866*x2_1*x8_0^2*k2_0 + 1056756553613745686827527368684983*x2_0*x8_0^2*k2_0 - x4_6^2*t1_0 - x5_6^2*t2_0 + 20*x13_3*x2_3 + 4*x10_3*x2_3 + 15*x2_4*x13_2 + 3*x2_4*x10_2 + 15*x13_4*x2_2 + 3*x10_4*x2_2 + 6*x2_5*x13_1 + k(6)//k(5)*x2_5*x10_1 + 6*x13_5*x2_1 + k(6)//k(5)*x10_5*x2_1 + x2_6*x13_0 + k(1)//k(5)*x2_6*x10_0 + x13_6*x2_0 + k(1)//k(5)*x10_6*x2_0 + x2_6*k_deg_0 - 630715481097022642713555616553466*x1_6*k1_0 - 8746249120595932857631589816128680*x1_5*k1_0 - 681921366620628664931033724447060*x1_4*k1_0 - 35179777131532892878398311236253000*x1_3*k1_0 - 23889039204893603636184878224411395*x1_2*k1_0 - 7274519188793555246037678283246866*x1_1*k1_0 - 1056756553613745686827527368684983*x1_0*k1_0 + x2_6*k3_0 + x2_7,
		x4_5^2*t1_0 + x4_6^2 - 2*x2_3*x10_2 - 2*x10_3*x2_2 - x2_4*x10_1 - x10_4*x2_1 - k(1)//k(5)*x2_5*x10_0 - k(1)//k(5)*x10_5*x2_0,
		x5_5^2*t2_0 + x5_6^2 - 10*x2_3*x13_2 - 10*x13_3*x2_2 - 5*x2_4*x13_1 - 5*x13_4*x2_1 - x2_5*x13_0 - x13_5*x2_0,
		x8_5^2*c5_0^3 + x8_6^2 - k(1)//k(2)*x9_5,
		x9_5 - k(1)//k(2000000)*x7_4 + k(1)//k(2500)*x9_4,
		-x3_6 - x2_6 - x1_6 + k(39921712922862322535422086479229399452896317293310383818522157734413081865999362477995175268387940471028808932554374649955843653921701194312435745761970070855452005305759945479288342890525871983773386901637839916412519660071866903296104534063220823407633483230946052099114908654531729181845515418477926345286786043174603921524094392324973118856727356552851399455474304888657541464076023071559435674046293190315476026351917463538250313253485550731436328288298016103405608569917141327556504043404973210103942673569132654509538678159346118874623497878266649188783079207718752452171027)//k(39062500000000000000000),
		-x12_2 + k(5438799144146346099590196103702350017513853289130828432864294061647579030223326976471335829206434064913503)//k(4000000),
		-x12_3 - k(63816610718695209915078628982541606002098886162338754812654855138198161423055886142659136828531904741536771681389130990216970860296761975852853)//k(40000000000),
		-x12_4 + k(3743988052955298423195233024540679885001583764177489840985425157537393321864462777745870717759960599324900785768160896634655182004053139462349786876155727893628730527945663390994153)//k(2000000000000000),
		-x12_5 - k(67974751019167478546295732298286120191147161742083826426860731324114156885462722717598801924183089381514306649882358116328294968480361133176399246195627363237918256389730876341526165326390472781565185366679178373613460973419479709753526931556653)//k(100000000000000000000),
		-x2_7 - k(30480537475184718377924330893059838408457995978064571152560638777549457801593194188568669551925836653791905388687779677214453863700234357684056219910057861990515726072423007221173810346712251733303384761057197199567343245050807437900311220050699892605317264271477501598692722148476230279727315392893899242978878187777930408034520732700321392325489144991287313324993890016233830619139035776484476140446344593750810027416010382844788458273901502540753152913204955576272187241285613270473000134661007748880102915145795462040400493281369414926592725558198522095428286361909039396357072004600219020081006557313621709063745760358038218509917376710287310881939588649051450194183821379440552507524165737769957075441605500920492186187562126592531482816876451)//k(976562500000000000000000000),
		-x9_1 - k(1208691041767043990364112225867499313)//k(2000000),
		-x9_2 + k(16413815353535588122168690736187175280396621276964029759247062712290563)//k(5000000000),
		-x9_3 - k(218448750620771070276054136883684589819971337353789311585189630079904677656657087188184666269462173551364251)//k(25000000000000),
		-x9_4 + k(3627801019457184591453090026566291407114843479405258946170342502884275410674003802648214725459647894768210889409376360370803999354951348485976129)//k(250000000000000000),
		-x9_5 - k(1699368775479186963657393301965852862539392855891554585983887973507603601528898743242862047478191564273186734474462820916531473712343209857335874023410113462420646795664651585900841757829791512486014574727435671678304951174521101284397721439070841)//k(2500000000000000000000),
		-x13_6 - x10_6 + k(7984342584572464507084417295845879890579263458662076763704431546853059311233257534131911295169371692166072715334928778292290160640544658689849588125034269001635445224092270474555054143160080874678078730210327532740327385804670100139718156385163257334410806428157866960787055390052503386930355121486998467657762318903211429933051813256154966507579787187286314370108253313882803042253257820651903867650242847919266299422337187297086194636036614941947993066903952530646077558757970448467022699534128843978633376185419535791821951006010985855533707293590136601796399532991944034599569789957)//k(7812500000000000000000000000),
		-x7_5 + k(843791360144110130007822424528923579982999620972626228839578289643716915363999891278729082805237432628512746433347673577428644637453390338922517861645288271869682421385216071897670482991853252512596660615981613362854972531838331970676573446798355792367927403964294718967570451485930313180827629110964490804601418986081243505596666060910182603769153)//k(2500000000000000000),
		z_aux - 1
    ]
end