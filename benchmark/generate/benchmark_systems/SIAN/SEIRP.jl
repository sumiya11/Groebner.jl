#! format: off
#! source: https://github.com/alexeyovchinnikov/SIAN-Julia

import AbstractAlgebra

function SEIRP(; np=AbstractAlgebra, ordering=:degrevlex, k=np.QQ)
    R, (S_8,I_8,E_7,S_7,I_7,E_6,S_6,I_6,E_5,S_5,I_5,E_4,S_4,I_4,E_3,S_3,I_3,E_2,S_2,I_2,E_1,S_1,I_1,E_0,S_0,I_0,z_aux,R_0,P_0,alpha_e_0,alpha_i_0,kappa_0,rho_0,beta_0,mu_0) = np.PolynomialRing(k, [:S_8,:I_8,:E_7,:S_7,:I_7,:E_6,:S_6,:I_6,:E_5,:S_5,:I_5,:E_4,:S_4,:I_4,:E_3,:S_3,:I_3,:E_2,:S_2,:I_2,:E_1,:S_1,:I_1,:E_0,:S_0,:I_0,:z_aux,:R_0,:P_0,:alpha_e_0,:alpha_i_0,:kappa_0,:rho_0,:beta_0,:mu_0], ordering=ordering)
    sys = [
    		-S_0 - I_0 + 6444397146265595,
		E_0*S_0*alpha_e_0 + S_0*I_0*alpha_i_0 + S_1,
		-E_0*kappa_0 + I_0*beta_0 + I_0*mu_0 + I_1,
		-S_1 - I_1 - 44948294184682007389228261743674445967320528369,
		S_1*E_0*alpha_e_0 + E_1*S_0*alpha_e_0 + I_1*S_0*alpha_i_0 + S_1*I_0*alpha_i_0 + S_2,
		-E_1*kappa_0 + I_1*beta_0 + I_1*mu_0 + I_2,
		-E_0*S_0*alpha_e_0 - S_0*I_0*alpha_i_0 + E_0*kappa_0 + E_0*rho_0 + E_1,
		-S_2 - I_2 + 1432268801925020189298702697881794961864418626942147023934312288139114268265097,
		-E_2*kappa_0 + I_2*beta_0 + I_2*mu_0 + I_3,
		2*E_1*S_1*alpha_e_0 + S_2*E_0*alpha_e_0 + E_2*S_0*alpha_e_0 + 2*S_1*I_1*alpha_i_0 + I_2*S_0*alpha_i_0 + S_2*I_0*alpha_i_0 + S_3,
		-S_1*E_0*alpha_e_0 - E_1*S_0*alpha_e_0 - I_1*S_0*alpha_i_0 - S_1*I_0*alpha_i_0 + E_1*kappa_0 + E_1*rho_0 + E_2,
		-S_3 - I_3 - 32165859966425054626053827698820694059944055162234599498441224409106465764148851918149870804335974675679023531,
		-E_3*kappa_0 + I_3*beta_0 + I_3*mu_0 + I_4,
		3*S_2*E_1*alpha_e_0 + 3*E_2*S_1*alpha_e_0 + S_3*E_0*alpha_e_0 + E_3*S_0*alpha_e_0 + 3*I_2*S_1*alpha_i_0 + 3*S_2*I_1*alpha_i_0 + I_3*S_0*alpha_i_0 + S_3*I_0*alpha_i_0 + S_4,
		-2*E_1*S_1*alpha_e_0 - S_2*E_0*alpha_e_0 - E_2*S_0*alpha_e_0 - 2*S_1*I_1*alpha_i_0 - I_2*S_0*alpha_i_0 - S_2*I_0*alpha_i_0 + E_2*kappa_0 + E_2*rho_0 + E_3,
		-S_4 - I_4 - 262995399363993236909312192735093141418781666489672081867669841011453984098841208060062277843065863151990532340378233269780944662899383101115,
		-E_4*kappa_0 + I_4*beta_0 + I_4*mu_0 + I_5,
		6*E_2*S_2*alpha_e_0 + 4*S_3*E_1*alpha_e_0 + 4*E_3*S_1*alpha_e_0 + S_4*E_0*alpha_e_0 + E_4*S_0*alpha_e_0 + 6*S_2*I_2*alpha_i_0 + 4*I_3*S_1*alpha_i_0 + 4*S_3*I_1*alpha_i_0 + I_4*S_0*alpha_i_0 + S_4*I_0*alpha_i_0 + S_5,
		-3*S_2*E_1*alpha_e_0 - 3*E_2*S_1*alpha_e_0 - S_3*E_0*alpha_e_0 - E_3*S_0*alpha_e_0 - 3*I_2*S_1*alpha_i_0 - 3*S_2*I_1*alpha_i_0 - I_3*S_0*alpha_i_0 - S_3*I_0*alpha_i_0 + E_3*kappa_0 + E_3*rho_0 + E_4,
		-S_5 - I_5 + 87987184785914097528896216677640151030690819025902552819891673952857036211531147506913362262283487662906199838935688751635158374213377996339671565777095934178526804739708823,
		10*S_3*E_2*alpha_e_0 + 10*E_3*S_2*alpha_e_0 + 5*S_4*E_1*alpha_e_0 + 5*E_4*S_1*alpha_e_0 + S_5*E_0*alpha_e_0 + E_5*S_0*alpha_e_0 + 10*I_3*S_2*alpha_i_0 + 10*S_3*I_2*alpha_i_0 + 5*I_4*S_1*alpha_i_0 + 5*S_4*I_1*alpha_i_0 + I_5*S_0*alpha_i_0 + S_5*I_0*alpha_i_0 + S_6,
		-E_5*kappa_0 + I_5*beta_0 + I_5*mu_0 + I_6,
		-6*E_2*S_2*alpha_e_0 - 4*S_3*E_1*alpha_e_0 - 4*E_3*S_1*alpha_e_0 - S_4*E_0*alpha_e_0 - E_4*S_0*alpha_e_0 - 6*S_2*I_2*alpha_i_0 - 4*I_3*S_1*alpha_i_0 - 4*S_3*I_1*alpha_i_0 - I_4*S_0*alpha_i_0 - S_4*I_0*alpha_i_0 + E_4*kappa_0 + E_4*rho_0 + E_5,
		-S_6 - I_6 - 5481817266254453310584489995477985221677164093612788756079404243330385572632907710770069514651691442981797064285703646122661102252191757025912927472713845036503449301326336711240230552443968708415291875961,
		20*E_3*S_3*alpha_e_0 + 15*S_4*E_2*alpha_e_0 + 15*E_4*S_2*alpha_e_0 + 6*S_5*E_1*alpha_e_0 + 6*E_5*S_1*alpha_e_0 + S_6*E_0*alpha_e_0 + E_6*S_0*alpha_e_0 + 20*S_3*I_3*alpha_i_0 + 15*I_4*S_2*alpha_i_0 + 15*S_4*I_2*alpha_i_0 + 6*I_5*S_1*alpha_i_0 + 6*S_5*I_1*alpha_i_0 + I_6*S_0*alpha_i_0 + S_6*I_0*alpha_i_0 + S_7,
		-E_6*kappa_0 + I_6*beta_0 + I_6*mu_0 + I_7,
		-10*S_3*E_2*alpha_e_0 - 10*E_3*S_2*alpha_e_0 - 5*S_4*E_1*alpha_e_0 - 5*E_4*S_1*alpha_e_0 - S_5*E_0*alpha_e_0 - E_5*S_0*alpha_e_0 - 10*I_3*S_2*alpha_i_0 - 10*S_3*I_2*alpha_i_0 - 5*I_4*S_1*alpha_i_0 - 5*S_4*I_1*alpha_i_0 - I_5*S_0*alpha_i_0 - S_5*I_0*alpha_i_0 + E_5*kappa_0 + E_5*rho_0 + E_6,
		-S_7 - I_7 + 47751476404928385214547903573076773582782949046267202708483537670473257862225963459631627322608324568044126966677890554203481838108845014232730381281860505095341303245227005713217877724731242456838629163948275708833905675794517141039947,
		35*S_4*E_3*alpha_e_0 + 35*E_4*S_3*alpha_e_0 + 21*S_5*E_2*alpha_e_0 + 21*E_5*S_2*alpha_e_0 + 7*S_6*E_1*alpha_e_0 + 7*E_6*S_1*alpha_e_0 + S_7*E_0*alpha_e_0 + E_7*S_0*alpha_e_0 + 35*I_4*S_3*alpha_i_0 + 35*S_4*I_3*alpha_i_0 + 21*I_5*S_2*alpha_i_0 + 21*S_5*I_2*alpha_i_0 + 7*I_6*S_1*alpha_i_0 + 7*S_6*I_1*alpha_i_0 + I_7*S_0*alpha_i_0 + S_7*I_0*alpha_i_0 + S_8,
		-E_7*kappa_0 + I_7*beta_0 + I_7*mu_0 + I_8,
		-20*E_3*S_3*alpha_e_0 - 15*S_4*E_2*alpha_e_0 - 15*E_4*S_2*alpha_e_0 - 6*S_5*E_1*alpha_e_0 - 6*E_5*S_1*alpha_e_0 - S_6*E_0*alpha_e_0 - E_6*S_0*alpha_e_0 - 20*S_3*I_3*alpha_i_0 - 15*I_4*S_2*alpha_i_0 - 15*S_4*I_2*alpha_i_0 - 6*I_5*S_1*alpha_i_0 - 6*S_5*I_1*alpha_i_0 - I_6*S_0*alpha_i_0 - S_6*I_0*alpha_i_0 + E_6*kappa_0 + E_6*rho_0 + E_7,
		-S_8 - I_8 + 29603367900684556824607040982669409478753632988028947240836289703184329144446628101135585263451747383307217143792360133449465352860151164617407863296061727591660358468194878052538427901706188380468066911854802299216746457176086577506993575154973331800679584695696453681,
		z_aux - 1
    ]
end

