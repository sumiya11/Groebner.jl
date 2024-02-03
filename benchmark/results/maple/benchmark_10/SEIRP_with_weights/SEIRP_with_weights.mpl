# SEIRP_with_weights
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		-S_0 - I_0 + 13673625206296414,
		E_0^2*S_0*alpha_e_0 + S_0*I_0*alpha_i_0 + S_1,
		-E_0^2*kappa_0 + I_0*beta_0 + I_0*mu_0 + I_1,
		-S_1 - I_1 - 99357445294036280367546329428608401561777524420,
		S_1*E_0^2*alpha_e_0 + E_1^2*S_0*alpha_e_0 + I_1*S_0*alpha_i_0 + S_1*I_0*alpha_i_0 + S_2,
		-E_1^2*kappa_0 + I_1*beta_0 + I_1*mu_0 + I_2,
		E_0^2*rho_0^3 - E_0^2*S_0*alpha_e_0 - S_0*I_0*alpha_i_0 + E_0^2*kappa_0 + E_1^2,
		-S_2 - I_2 - 929567995557832077064050787673677227430956100579220917953774108769413773883684,
		-E_2^2*kappa_0 + I_2*beta_0 + I_2*mu_0 + I_3,
		2*E_1^2*S_1*alpha_e_0 + S_2*E_0^2*alpha_e_0 + E_2^2*S_0*alpha_e_0 + 2*S_1*I_1*alpha_i_0 + I_2*S_0*alpha_i_0 + S_2*I_0*alpha_i_0 + S_3,
		E_1^2*rho_0^3 - S_1*E_0^2*alpha_e_0 - E_1^2*S_0*alpha_e_0 - I_1*S_0*alpha_i_0 - S_1*I_0*alpha_i_0 + E_1^2*kappa_0 + E_2^2,
		-S_3 - I_3 + 64542696494011423001259456392269691443247919585653861435057914881184567659453893776525267004937536332666207676,
		-E_3^2*kappa_0 + I_3*beta_0 + I_3*mu_0 + I_4,
		3*S_2*E_1^2*alpha_e_0 + 3*E_2^2*S_1*alpha_e_0 + S_3*E_0^2*alpha_e_0 + E_3^2*S_0*alpha_e_0 + 3*I_2*S_1*alpha_i_0 + 3*S_2*I_1*alpha_i_0 + I_3*S_0*alpha_i_0 + S_3*I_0*alpha_i_0 + S_4,
		E_2^2*rho_0^3 - 2*E_1^2*S_1*alpha_e_0 - S_2*E_0^2*alpha_e_0 - E_2^2*S_0*alpha_e_0 - 2*S_1*I_1*alpha_i_0 - I_2*S_0*alpha_i_0 - S_2*I_0*alpha_i_0 + E_2^2*kappa_0 + E_3^2,
		-S_4 - I_4 + 2659491048644817986157377559256942536322862343689193368038638862045583793392179742254724870932887938273017967718564343675055693202002991201148,
		-E_4^2*kappa_0 + I_4*beta_0 + I_4*mu_0 + I_5,
		6*E_2^2*S_2*alpha_e_0 + 4*S_3*E_1^2*alpha_e_0 + 4*E_3^2*S_1*alpha_e_0 + S_4*E_0^2*alpha_e_0 + E_4^2*S_0*alpha_e_0 + 6*S_2*I_2*alpha_i_0 + 4*I_3*S_1*alpha_i_0 + 4*S_3*I_1*alpha_i_0 + I_4*S_0*alpha_i_0 + S_4*I_0*alpha_i_0 + S_5,
		E_3^2*rho_0^3 - 3*S_2*E_1^2*alpha_e_0 - 3*E_2^2*S_1*alpha_e_0 - S_3*E_0^2*alpha_e_0 - E_3^2*S_0*alpha_e_0 - 3*I_2*S_1*alpha_i_0 - 3*S_2*I_1*alpha_i_0 - I_3*S_0*alpha_i_0 - S_3*I_0*alpha_i_0 + E_3^2*kappa_0 + E_4^2,
		-S_5 - I_5 - 146192098734467153506935560270700961564715477333312532724999888151868936987328777554251152066424360388944972324072850746130345740913891458288635774519281530986705782103504260,
		10*S_3*E_2^2*alpha_e_0 + 10*E_3^2*S_2*alpha_e_0 + 5*S_4*E_1^2*alpha_e_0 + 5*E_4^2*S_1*alpha_e_0 + S_5*E_0^2*alpha_e_0 + E_5^2*S_0*alpha_e_0 + 10*I_3*S_2*alpha_i_0 + 10*S_3*I_2*alpha_i_0 + 5*I_4*S_1*alpha_i_0 + 5*S_4*I_1*alpha_i_0 + I_5*S_0*alpha_i_0 + S_5*I_0*alpha_i_0 + S_6,
		-E_5^2*kappa_0 + I_5*beta_0 + I_5*mu_0 + I_6,
		E_4^2*rho_0^3 - 6*E_2^2*S_2*alpha_e_0 - 4*S_3*E_1^2*alpha_e_0 - 4*E_3^2*S_1*alpha_e_0 - S_4*E_0^2*alpha_e_0 - E_4^2*S_0*alpha_e_0 - 6*S_2*I_2*alpha_i_0 - 4*I_3*S_1*alpha_i_0 - 4*S_3*I_1*alpha_i_0 - I_4*S_0*alpha_i_0 - S_4*I_0*alpha_i_0 + E_4^2*kappa_0 + E_5^2,
		-S_6 - I_6 - 15620880917431589260192029799780449047081925169104731969516096420102292860979390072095718035621848578662369930940834170256963219625513586410073911604180745733591146699955755967655922991079048958099538845444,
		20*E_3^2*S_3*alpha_e_0 + 15*S_4*E_2^2*alpha_e_0 + 15*E_4^2*S_2*alpha_e_0 + 6*S_5*E_1^2*alpha_e_0 + 6*E_5^2*S_1*alpha_e_0 + S_6*E_0^2*alpha_e_0 + E_6^2*S_0*alpha_e_0 + 20*S_3*I_3*alpha_i_0 + 15*I_4*S_2*alpha_i_0 + 15*S_4*I_2*alpha_i_0 + 6*I_5*S_1*alpha_i_0 + 6*S_5*I_1*alpha_i_0 + I_6*S_0*alpha_i_0 + S_6*I_0*alpha_i_0 + S_7,
		-E_6^2*kappa_0 + I_6*beta_0 + I_6*mu_0 + I_7,
		E_5^2*rho_0^3 - 10*S_3*E_2^2*alpha_e_0 - 10*E_3^2*S_2*alpha_e_0 - 5*S_4*E_1^2*alpha_e_0 - 5*E_4^2*S_1*alpha_e_0 - S_5*E_0^2*alpha_e_0 - E_5^2*S_0*alpha_e_0 - 10*I_3*S_2*alpha_i_0 - 10*S_3*I_2*alpha_i_0 - 5*I_4*S_1*alpha_i_0 - 5*S_4*I_1*alpha_i_0 - I_5*S_0*alpha_i_0 - S_5*I_0*alpha_i_0 + E_5^2*kappa_0 + E_6^2,
		-S_7 - I_7 + 534372767693344517649879922581492452696085044573547245228410867449666652631371434854544755388723625707318135550399286358811679850015084711554218055680414793992270776616238332173671855095780827268309814892675097886614345218038772276096636,
		35*S_4*E_3^2*alpha_e_0 + 35*E_4^2*S_3*alpha_e_0 + 21*S_5*E_2^2*alpha_e_0 + 21*E_5^2*S_2*alpha_e_0 + 7*S_6*E_1^2*alpha_e_0 + 7*E_6^2*S_1*alpha_e_0 + S_7*E_0^2*alpha_e_0 + E_7^2*S_0*alpha_e_0 + 35*I_4*S_3*alpha_i_0 + 35*S_4*I_3*alpha_i_0 + 21*I_5*S_2*alpha_i_0 + 21*S_5*I_2*alpha_i_0 + 7*I_6*S_1*alpha_i_0 + 7*S_6*I_1*alpha_i_0 + I_7*S_0*alpha_i_0 + S_7*I_0*alpha_i_0 + S_8,
		-E_7^2*kappa_0 + I_7*beta_0 + I_7*mu_0 + I_8,
		E_6^2*rho_0^3 - 20*E_3^2*S_3*alpha_e_0 - 15*S_4*E_2^2*alpha_e_0 - 15*E_4^2*S_2*alpha_e_0 - 6*S_5*E_1^2*alpha_e_0 - 6*E_5^2*S_1*alpha_e_0 - S_6*E_0^2*alpha_e_0 - E_6^2*S_0*alpha_e_0 - 20*S_3*I_3*alpha_i_0 - 15*I_4*S_2*alpha_i_0 - 15*S_4*I_2*alpha_i_0 - 6*I_5*S_1*alpha_i_0 - 6*S_5*I_1*alpha_i_0 - I_6*S_0*alpha_i_0 - S_6*I_0*alpha_i_0 + E_6^2*kappa_0 + E_7^2,
		-S_8 - I_8 + 151346014913786518775536769433542612972524399560281857293371870661348736536311856413904752800233356566113690210315518356411076707262604251010192052877962069158524837700706160088978323245083819477766421500967190497343699643434533190443047278890767385756340898665880232188,
		z_aux - 1
	]:
	print("Running SEIRP_with_weights");
	st := time[real]():
	G := Groebner[Basis](J, tdeg(S_8, I_8, E_7, S_7, I_7, E_6, S_6, I_6, E_5, S_5, I_5, E_4, S_4, I_4, E_3, S_3, I_3, E_2, S_2, I_2, E_1, S_1, I_1, E_0, S_0, I_0, z_aux, R_0, P_0, alpha_e_0, alpha_i_0, kappa_0, rho_0, beta_0, mu_0), method=fgb, characteristic=0):
	print("SEIRP_with_weights: ", time[real]() - st):
	runtime := min(runtime, time[real]() - st):
end do:

timings_fn := "/home/demin/Groebner.jl/benchmark/results/maple/benchmark_10/SEIRP_with_weights/timings":
FileTools[Text][WriteLine](timings_fn, "SEIRP_with_weights");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime))):

