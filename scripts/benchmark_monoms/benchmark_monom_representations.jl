using Pkg; Pkg.activate(temp=true);
Pkg.add(url="https://github.com/matthias314/Groebner.jl#m3/fixedvector", rev="cff7e7b191d92db162beba1dd2d290aaabccc482")
Pkg.add("AbstractAlgebra"); Pkg.add("BenchmarkTools"); Pkg.add("Primes");
Pkg.instantiate()

using Groebner, AbstractAlgebra, BenchmarkTools, Printf, Logging, InteractiveUtils, Primes
global_logger(SimpleLogger(stderr, Logging.Error))

function akt_system()
# akt_pathway
R, (Akt_0,Akt_1,Akt_2,Akt_3,Akt_4,Akt_5,EGF_EGFR_0,EGF_EGFR_1,EGF_EGFR_2,EGF_EGFR_3,EGF_EGFR_4,EGF_EGFR_5,S6_0,S6_1,S6_2,S6_3,S6_4,S6_5,S6_6,a2_0,a3_0,pAkt_0,pAkt_1,pAkt_2,pAkt_3,pAkt_4,pAkt_5,pAkt_6,pAkt_7,pAkt_S6_0,pAkt_S6_1,pAkt_S6_2,pAkt_S6_3,pAkt_S6_4,pAkt_S6_5,pAkt_S6_6,pAkt_S6_7,pEGFR_0,pEGFR_1,pEGFR_2,pEGFR_3,pEGFR_4,pEGFR_5,pEGFR_6,pEGFR_Akt_0,pEGFR_Akt_1,pEGFR_Akt_2,pEGFR_Akt_3,pEGFR_Akt_4,pEGFR_Akt_5,pEGFR_Akt_6,pS6_0,pS6_1,pS6_2,pS6_3,pS6_4,pS6_5,pS6_6,pS6_7,reaction_2_k1_0,reaction_2_k2_0,reaction_3_k1_0,reaction_4_k1_0,reaction_5_k1_0,reaction_5_k2_0,reaction_6_k1_0,reaction_7_k1_0,reaction_8_k1_0,reaction_9_k1_0) = polynomial_ring(QQ, ["Akt_0","Akt_1","Akt_2","Akt_3","Akt_4","Akt_5","EGF_EGFR_0","EGF_EGFR_1","EGF_EGFR_2","EGF_EGFR_3","EGF_EGFR_4","EGF_EGFR_5","S6_0","S6_1","S6_2","S6_3","S6_4","S6_5","S6_6","a2_0","a3_0","pAkt_0","pAkt_1","pAkt_2","pAkt_3","pAkt_4","pAkt_5","pAkt_6","pAkt_7","pAkt_S6_0","pAkt_S6_1","pAkt_S6_2","pAkt_S6_3","pAkt_S6_4","pAkt_S6_5","pAkt_S6_6","pAkt_S6_7","pEGFR_0","pEGFR_1","pEGFR_2","pEGFR_3","pEGFR_4","pEGFR_5","pEGFR_6","pEGFR_Akt_0","pEGFR_Akt_1","pEGFR_Akt_2","pEGFR_Akt_3","pEGFR_Akt_4","pEGFR_Akt_5","pEGFR_Akt_6","pS6_0","pS6_1","pS6_2","pS6_3","pS6_4","pS6_5","pS6_6","pS6_7","reaction_2_k1_0","reaction_2_k2_0","reaction_3_k1_0","reaction_4_k1_0","reaction_5_k1_0","reaction_5_k2_0","reaction_6_k1_0","reaction_7_k1_0","reaction_8_k1_0","reaction_9_k1_0"], internal_ordering=:degrevlex)

sys = [
-a2_0*pAkt_0 - a2_0*pAkt_S6_0 + 1870258352796715//4503599627370496,
S6_0*pAkt_0*reaction_5_k1_0 - pEGFR_Akt_0*reaction_3_k1_0 - pAkt_S6_0*reaction_5_k2_0 - pAkt_S6_0*reaction_6_k1_0 + pAkt_0*reaction_7_k1_0 + pAkt_1,
-S6_0*pAkt_0*reaction_5_k1_0 + pAkt_S6_0*reaction_5_k2_0 + pAkt_S6_0*reaction_6_k1_0 + pAkt_S6_1,
-a3_0*pS6_0 + 4258337759705301//18014398509481984,
-pAkt_S6_0*reaction_6_k1_0 + pS6_0*reaction_8_k1_0 + pS6_1,
-1785199063384347680030081822281357*pEGFR_0 - 1785199063384347680030081822281357*pEGFR_Akt_0 + 858895685354599//2251799813685248,
Akt_0*pEGFR_0*reaction_2_k1_0 - pEGFR_Akt_0*reaction_2_k2_0 - pEGFR_Akt_0*reaction_3_k1_0 + pEGFR_0*reaction_4_k1_0 - EGF_EGFR_0*reaction_9_k1_0 + pEGFR_1,
-Akt_0*pEGFR_0*reaction_2_k1_0 + pEGFR_Akt_0*reaction_2_k2_0 + pEGFR_Akt_0*reaction_3_k1_0 + pEGFR_Akt_1,
-a2_0*pAkt_1 - a2_0*pAkt_S6_1 - 7492296858826135//144115188075855872,
-S6_1*pAkt_0*reaction_5_k1_0 - S6_0*pAkt_1*reaction_5_k1_0 + pAkt_S6_1*reaction_5_k2_0 + pAkt_S6_1*reaction_6_k1_0 + pAkt_S6_2,
S6_1*pAkt_0*reaction_5_k1_0 + S6_0*pAkt_1*reaction_5_k1_0 - pEGFR_Akt_1*reaction_3_k1_0 - pAkt_S6_1*reaction_5_k2_0 - pAkt_S6_1*reaction_6_k1_0 + pAkt_1*reaction_7_k1_0 + pAkt_2,
S6_0*pAkt_0*reaction_5_k1_0 - pAkt_S6_0*reaction_5_k2_0 - pS6_0*reaction_8_k1_0 + S6_1,
-a3_0*pS6_1 - 8490510622091381//288230376151711744,
-pAkt_S6_1*reaction_6_k1_0 + pS6_1*reaction_8_k1_0 + pS6_2,
-1785199063384347680030081822281357*pEGFR_1 - 1785199063384347680030081822281357*pEGFR_Akt_1 - 8278789211907929//1152921504606846976,
-Akt_1*pEGFR_0*reaction_2_k1_0 - Akt_0*pEGFR_1*reaction_2_k1_0 + pEGFR_Akt_1*reaction_2_k2_0 + pEGFR_Akt_1*reaction_3_k1_0 + pEGFR_Akt_2,
Akt_1*pEGFR_0*reaction_2_k1_0 + Akt_0*pEGFR_1*reaction_2_k1_0 - pEGFR_Akt_1*reaction_2_k2_0 - pEGFR_Akt_1*reaction_3_k1_0 + pEGFR_1*reaction_4_k1_0 - EGF_EGFR_1*reaction_9_k1_0 + pEGFR_2,
EGF_EGFR_0*reaction_9_k1_0 + 553946985179424024508562363759772*EGF_EGFR_0 + EGF_EGFR_1,
Akt_0*pEGFR_0*reaction_2_k1_0 - pEGFR_Akt_0*reaction_2_k2_0 - pAkt_0*reaction_7_k1_0 + Akt_1,
-a2_0*pAkt_2 - a2_0*pAkt_S6_2 - 4623962009751873//288230376151711744,
-S6_2*pAkt_0*reaction_5_k1_0 - 2*S6_1*pAkt_1*reaction_5_k1_0 - S6_0*pAkt_2*reaction_5_k1_0 + pAkt_S6_2*reaction_5_k2_0 + pAkt_S6_2*reaction_6_k1_0 + pAkt_S6_3,
S6_2*pAkt_0*reaction_5_k1_0 + 2*S6_1*pAkt_1*reaction_5_k1_0 + S6_0*pAkt_2*reaction_5_k1_0 - pEGFR_Akt_2*reaction_3_k1_0 - pAkt_S6_2*reaction_5_k2_0 - pAkt_S6_2*reaction_6_k1_0 + pAkt_2*reaction_7_k1_0 + pAkt_3,
S6_1*pAkt_0*reaction_5_k1_0 + S6_0*pAkt_1*reaction_5_k1_0 - pAkt_S6_1*reaction_5_k2_0 - pS6_1*reaction_8_k1_0 + S6_2,
-a3_0*pS6_2 - 7833812891298975//576460752303423488,
-pAkt_S6_2*reaction_6_k1_0 + pS6_2*reaction_8_k1_0 + pS6_3,
-1785199063384347680030081822281357*pEGFR_2 - 1785199063384347680030081822281357*pEGFR_Akt_2 - 587731377280161//9007199254740992,
Akt_2*pEGFR_0*reaction_2_k1_0 + 2*Akt_1*pEGFR_1*reaction_2_k1_0 + Akt_0*pEGFR_2*reaction_2_k1_0 - pEGFR_Akt_2*reaction_2_k2_0 - pEGFR_Akt_2*reaction_3_k1_0 + pEGFR_2*reaction_4_k1_0 - EGF_EGFR_2*reaction_9_k1_0 + pEGFR_3,
-Akt_2*pEGFR_0*reaction_2_k1_0 - 2*Akt_1*pEGFR_1*reaction_2_k1_0 - Akt_0*pEGFR_2*reaction_2_k1_0 + pEGFR_Akt_2*reaction_2_k2_0 + pEGFR_Akt_2*reaction_3_k1_0 + pEGFR_Akt_3,
Akt_1*pEGFR_0*reaction_2_k1_0 + Akt_0*pEGFR_1*reaction_2_k1_0 - pEGFR_Akt_1*reaction_2_k2_0 - pAkt_1*reaction_7_k1_0 + Akt_2,
EGF_EGFR_1*reaction_9_k1_0 + 553946985179424024508562363759772*EGF_EGFR_1 + EGF_EGFR_2,
-a2_0*pAkt_3 - a2_0*pAkt_S6_3 + 6402045533440037//72057594037927936,
S6_3*pAkt_0*reaction_5_k1_0 + 3*S6_2*pAkt_1*reaction_5_k1_0 + 3*S6_1*pAkt_2*reaction_5_k1_0 + S6_0*pAkt_3*reaction_5_k1_0 - pEGFR_Akt_3*reaction_3_k1_0 - pAkt_S6_3*reaction_5_k2_0 - pAkt_S6_3*reaction_6_k1_0 + pAkt_3*reaction_7_k1_0 + pAkt_4,
-S6_3*pAkt_0*reaction_5_k1_0 - 3*S6_2*pAkt_1*reaction_5_k1_0 - 3*S6_1*pAkt_2*reaction_5_k1_0 - S6_0*pAkt_3*reaction_5_k1_0 + pAkt_S6_3*reaction_5_k2_0 + pAkt_S6_3*reaction_6_k1_0 + pAkt_S6_4,
S6_2*pAkt_0*reaction_5_k1_0 + 2*S6_1*pAkt_1*reaction_5_k1_0 + S6_0*pAkt_2*reaction_5_k1_0 - pAkt_S6_2*reaction_5_k2_0 - pS6_2*reaction_8_k1_0 + S6_3,
-a3_0*pS6_3 + 4338314774287177//72057594037927936,
-pAkt_S6_3*reaction_6_k1_0 + pS6_3*reaction_8_k1_0 + pS6_4,
-1785199063384347680030081822281357*pEGFR_3 - 1785199063384347680030081822281357*pEGFR_Akt_3 + 6516429318617311//72057594037927936,
Akt_3*pEGFR_0*reaction_2_k1_0 + 3*Akt_2*pEGFR_1*reaction_2_k1_0 + 3*Akt_1*pEGFR_2*reaction_2_k1_0 + Akt_0*pEGFR_3*reaction_2_k1_0 - pEGFR_Akt_3*reaction_2_k2_0 - pEGFR_Akt_3*reaction_3_k1_0 + pEGFR_3*reaction_4_k1_0 - EGF_EGFR_3*reaction_9_k1_0 + pEGFR_4,
-Akt_3*pEGFR_0*reaction_2_k1_0 - 3*Akt_2*pEGFR_1*reaction_2_k1_0 - 3*Akt_1*pEGFR_2*reaction_2_k1_0 - Akt_0*pEGFR_3*reaction_2_k1_0 + pEGFR_Akt_3*reaction_2_k2_0 + pEGFR_Akt_3*reaction_3_k1_0 + pEGFR_Akt_4,
EGF_EGFR_2*reaction_9_k1_0 + 553946985179424024508562363759772*EGF_EGFR_2 + EGF_EGFR_3,
Akt_2*pEGFR_0*reaction_2_k1_0 + 2*Akt_1*pEGFR_1*reaction_2_k1_0 + Akt_0*pEGFR_2*reaction_2_k1_0 - pEGFR_Akt_2*reaction_2_k2_0 - pAkt_2*reaction_7_k1_0 + Akt_3,
-a2_0*pAkt_4 - a2_0*pAkt_S6_4 - 7336860684178615//36028797018963968,
S6_4*pAkt_0*reaction_5_k1_0 + 4*S6_3*pAkt_1*reaction_5_k1_0 + 6*S6_2*pAkt_2*reaction_5_k1_0 + 4*S6_1*pAkt_3*reaction_5_k1_0 + S6_0*pAkt_4*reaction_5_k1_0 - pEGFR_Akt_4*reaction_3_k1_0 - pAkt_S6_4*reaction_5_k2_0 - pAkt_S6_4*reaction_6_k1_0 + pAkt_4*reaction_7_k1_0 + pAkt_5,
-S6_4*pAkt_0*reaction_5_k1_0 - 4*S6_3*pAkt_1*reaction_5_k1_0 - 6*S6_2*pAkt_2*reaction_5_k1_0 - 4*S6_1*pAkt_3*reaction_5_k1_0 - S6_0*pAkt_4*reaction_5_k1_0 + pAkt_S6_4*reaction_5_k2_0 + pAkt_S6_4*reaction_6_k1_0 + pAkt_S6_5,
S6_3*pAkt_0*reaction_5_k1_0 + 3*S6_2*pAkt_1*reaction_5_k1_0 + 3*S6_1*pAkt_2*reaction_5_k1_0 + S6_0*pAkt_3*reaction_5_k1_0 - pAkt_S6_3*reaction_5_k2_0 - pS6_3*reaction_8_k1_0 + S6_4,
-a3_0*pS6_4 - 641351474169837//4503599627370496,
-pAkt_S6_4*reaction_6_k1_0 + pS6_4*reaction_8_k1_0 + pS6_5,
-1785199063384347680030081822281357*pEGFR_4 - 1785199063384347680030081822281357*pEGFR_Akt_4 - 124452506710145//1125899906842624,
-Akt_4*pEGFR_0*reaction_2_k1_0 - 4*Akt_3*pEGFR_1*reaction_2_k1_0 - 6*Akt_2*pEGFR_2*reaction_2_k1_0 - 4*Akt_1*pEGFR_3*reaction_2_k1_0 - Akt_0*pEGFR_4*reaction_2_k1_0 + pEGFR_Akt_4*reaction_2_k2_0 + pEGFR_Akt_4*reaction_3_k1_0 + pEGFR_Akt_5,
Akt_4*pEGFR_0*reaction_2_k1_0 + 4*Akt_3*pEGFR_1*reaction_2_k1_0 + 6*Akt_2*pEGFR_2*reaction_2_k1_0 + 4*Akt_1*pEGFR_3*reaction_2_k1_0 + Akt_0*pEGFR_4*reaction_2_k1_0 - pEGFR_Akt_4*reaction_2_k2_0 - pEGFR_Akt_4*reaction_3_k1_0 + pEGFR_4*reaction_4_k1_0 - EGF_EGFR_4*reaction_9_k1_0 + pEGFR_5,
EGF_EGFR_3*reaction_9_k1_0 + 553946985179424024508562363759772*EGF_EGFR_3 + EGF_EGFR_4,
Akt_3*pEGFR_0*reaction_2_k1_0 + 3*Akt_2*pEGFR_1*reaction_2_k1_0 + 3*Akt_1*pEGFR_2*reaction_2_k1_0 + Akt_0*pEGFR_3*reaction_2_k1_0 - pEGFR_Akt_3*reaction_2_k2_0 - pAkt_3*reaction_7_k1_0 + Akt_4,
-a2_0*pAkt_5 - a2_0*pAkt_S6_5 + 6442496259164395//18014398509481984,
S6_5*pAkt_0*reaction_5_k1_0 + 5*S6_4*pAkt_1*reaction_5_k1_0 + 10*S6_3*pAkt_2*reaction_5_k1_0 + 10*S6_2*pAkt_3*reaction_5_k1_0 + 5*S6_1*pAkt_4*reaction_5_k1_0 + S6_0*pAkt_5*reaction_5_k1_0 - pEGFR_Akt_5*reaction_3_k1_0 - pAkt_S6_5*reaction_5_k2_0 - pAkt_S6_5*reaction_6_k1_0 + pAkt_5*reaction_7_k1_0 + pAkt_6,
-S6_5*pAkt_0*reaction_5_k1_0 - 5*S6_4*pAkt_1*reaction_5_k1_0 - 10*S6_3*pAkt_2*reaction_5_k1_0 - 10*S6_2*pAkt_3*reaction_5_k1_0 - 5*S6_1*pAkt_4*reaction_5_k1_0 - S6_0*pAkt_5*reaction_5_k1_0 + pAkt_S6_5*reaction_5_k2_0 + pAkt_S6_5*reaction_6_k1_0 + pAkt_S6_6,
S6_4*pAkt_0*reaction_5_k1_0 + 4*S6_3*pAkt_1*reaction_5_k1_0 + 6*S6_2*pAkt_2*reaction_5_k1_0 + 4*S6_1*pAkt_3*reaction_5_k1_0 + S6_0*pAkt_4*reaction_5_k1_0 - pAkt_S6_4*reaction_5_k2_0 - pS6_4*reaction_8_k1_0 + S6_5,
-a3_0*pS6_5 + 5106054353643969//18014398509481984,
-pAkt_S6_5*reaction_6_k1_0 + pS6_5*reaction_8_k1_0 + pS6_6,
-1785199063384347680030081822281357*pEGFR_5 - 1785199063384347680030081822281357*pEGFR_Akt_5 + 8291745550051749//72057594037927936,
Akt_5*pEGFR_0*reaction_2_k1_0 + 5*Akt_4*pEGFR_1*reaction_2_k1_0 + 10*Akt_3*pEGFR_2*reaction_2_k1_0 + 10*Akt_2*pEGFR_3*reaction_2_k1_0 + 5*Akt_1*pEGFR_4*reaction_2_k1_0 + Akt_0*pEGFR_5*reaction_2_k1_0 - pEGFR_Akt_5*reaction_2_k2_0 - pEGFR_Akt_5*reaction_3_k1_0 + pEGFR_5*reaction_4_k1_0 - EGF_EGFR_5*reaction_9_k1_0 + pEGFR_6,
-Akt_5*pEGFR_0*reaction_2_k1_0 - 5*Akt_4*pEGFR_1*reaction_2_k1_0 - 10*Akt_3*pEGFR_2*reaction_2_k1_0 - 10*Akt_2*pEGFR_3*reaction_2_k1_0 - 5*Akt_1*pEGFR_4*reaction_2_k1_0 - Akt_0*pEGFR_5*reaction_2_k1_0 + pEGFR_Akt_5*reaction_2_k2_0 + pEGFR_Akt_5*reaction_3_k1_0 + pEGFR_Akt_6,
EGF_EGFR_4*reaction_9_k1_0 + 553946985179424024508562363759772*EGF_EGFR_4 + EGF_EGFR_5,
Akt_4*pEGFR_0*reaction_2_k1_0 + 4*Akt_3*pEGFR_1*reaction_2_k1_0 + 6*Akt_2*pEGFR_2*reaction_2_k1_0 + 4*Akt_1*pEGFR_3*reaction_2_k1_0 + Akt_0*pEGFR_4*reaction_2_k1_0 - pEGFR_Akt_4*reaction_2_k2_0 - pAkt_4*reaction_7_k1_0 + Akt_5,
-a2_0*pAkt_6 - a2_0*pAkt_S6_6 - 7867619963948889//18014398509481984,
S6_6*pAkt_0*reaction_5_k1_0 + 6*S6_5*pAkt_1*reaction_5_k1_0 + 15*S6_4*pAkt_2*reaction_5_k1_0 + 20*S6_3*pAkt_3*reaction_5_k1_0 + 15*S6_2*pAkt_4*reaction_5_k1_0 + 6*S6_1*pAkt_5*reaction_5_k1_0 + S6_0*pAkt_6*reaction_5_k1_0 - pEGFR_Akt_6*reaction_3_k1_0 - pAkt_S6_6*reaction_5_k2_0 - pAkt_S6_6*reaction_6_k1_0 + pAkt_6*reaction_7_k1_0 + pAkt_7,
-S6_6*pAkt_0*reaction_5_k1_0 - 6*S6_5*pAkt_1*reaction_5_k1_0 - 15*S6_4*pAkt_2*reaction_5_k1_0 - 20*S6_3*pAkt_3*reaction_5_k1_0 - 15*S6_2*pAkt_4*reaction_5_k1_0 - 6*S6_1*pAkt_5*reaction_5_k1_0 - S6_0*pAkt_6*reaction_5_k1_0 + pAkt_S6_6*reaction_5_k2_0 + pAkt_S6_6*reaction_6_k1_0 + pAkt_S6_7,
S6_5*pAkt_0*reaction_5_k1_0 + 5*S6_4*pAkt_1*reaction_5_k1_0 + 10*S6_3*pAkt_2*reaction_5_k1_0 + 10*S6_2*pAkt_3*reaction_5_k1_0 + 5*S6_1*pAkt_4*reaction_5_k1_0 + S6_0*pAkt_5*reaction_5_k1_0 - pAkt_S6_5*reaction_5_k2_0 - pS6_5*reaction_8_k1_0 + S6_6,
-a3_0*pS6_6 - 7863926901127043//18014398509481984,
-pAkt_S6_6*reaction_6_k1_0 + pS6_6*reaction_8_k1_0 + pS6_7,
]

sys_zp = map(f -> map_coefficients(c -> GF(2^30+3)(numerator(c)) // GF(2^30+3)(denominator(c)), f), sys);
return sys_zp
end

function crauste_system()
# crauste-6
R, (E_0,E_1,E_2,E_3,E_4,E_5,M_0,M_1,M_2,M_3,M_4,M_5,N_0,N_1,N_2,N_3,N_4,P_0,P_1,P_2,P_3,P_4,P_5,P_6,S_0,S_1,S_2,S_3,S_4,S_5,delta_EL_0,delta_LM_0,delta_NE_0,mu_EE_0,mu_LE_0,mu_LL_0,mu_M_0,mu_N_0,mu_PE_0,mu_PL_0,mu_P_0,rho_E_0,rho_P_0) = polynomial_ring(QQ, ["E_0","E_1","E_2","E_3","E_4","E_5","M_0","M_1","M_2","M_3","M_4","M_5","N_0","N_1","N_2","N_3","N_4","P_0","P_1","P_2","P_3","P_4","P_5","P_6","S_0","S_1","S_2","S_3","S_4","S_5","delta_EL_0","delta_LM_0","delta_NE_0","mu_EE_0","mu_LE_0","mu_LL_0","mu_M_0","mu_N_0","mu_PE_0","mu_PL_0","mu_P_0","rho_E_0","rho_P_0"], internal_ordering=:degrevlex)

sys = [
-E_0 + 2843700732990995//4503599627370496,
E_0^2*mu_EE_0^2 - E_0*P_0*rho_E_0^2 - N_0*P_0*delta_NE_0 + E_0*delta_EL_0 + E_1,
-N_0 + 497546696075019//1125899906842624,
N_0*P_0*delta_NE_0^2 + N_0*mu_N_0^2 + N_1,
-P_0 + 5998952937064603//9007199254740992,
P_0*S_0*mu_PL_0^2 - P_0^2*rho_P_0^2 + E_0*P_0*mu_PE_0 + P_0*mu_P_0 + P_1,
-M_0 - S_0 + 2944392902821235//2251799813685248,
-S_0*delta_LM_0 + M_0*mu_M_0 + M_1,
E_0*S_0*mu_LE_0 + S_0^2*mu_LL_0 - S_0*delta_EL_0 + S_0*delta_LM_0 + S_1,
-E_1 - 6037301980703049//9007199254740992,
2*E_0*E_1*mu_EE_0^2 - E_1*P_0*rho_E_0^2 - E_0*P_1*rho_E_0^2 - N_1*P_0*delta_NE_0 - N_0*P_1*delta_NE_0 + E_1*delta_EL_0 + E_2,
-N_1 - 6631372290250479//9007199254740992,
N_1*P_0*delta_NE_0^2 + N_0*P_1*delta_NE_0^2 + N_1*mu_N_0^2 + N_2,
-P_1 - 4937411093259505//9007199254740992,
P_1*S_0*mu_PL_0^2 + P_0*S_1*mu_PL_0^2 - 2*P_0*P_1*rho_P_0^2 + E_1*P_0*mu_PE_0 + E_0*P_1*mu_PE_0 + P_1*mu_P_0 + P_2,
-M_1 - S_1 - 1175297854514141//1125899906842624,
E_1*S_0*mu_LE_0 + E_0*S_1*mu_LE_0 + 2*S_0*S_1*mu_LL_0 - S_1*delta_EL_0 + S_1*delta_LM_0 + S_2,
-S_1*delta_LM_0 + M_1*mu_M_0 + M_2,
-E_2 + 6541535684549225//9007199254740992,
2*E_1^2*mu_EE_0^2 + 2*E_0*E_2*mu_EE_0^2 - E_2*P_0*rho_E_0^2 - 2*E_1*P_1*rho_E_0^2 - E_0*P_2*rho_E_0^2 - N_2*P_0*delta_NE_0 - 2*N_1*P_1*delta_NE_0 - N_0*P_2*delta_NE_0 + E_2*delta_EL_0 + E_3,
-N_2 + 6614937622156439//4503599627370496,
N_2*P_0*delta_NE_0^2 + 2*N_1*P_1*delta_NE_0^2 + N_0*P_2*delta_NE_0^2 + N_2*mu_N_0^2 + N_3,
-P_2 + 4848135714772055//9007199254740992,
P_2*S_0*mu_PL_0^2 + 2*P_1*S_1*mu_PL_0^2 + P_0*S_2*mu_PL_0^2 - 2*P_1^2*rho_P_0^2 - 2*P_0*P_2*rho_P_0^2 + E_2*P_0*mu_PE_0 + 2*E_1*P_1*mu_PE_0 + E_0*P_2*mu_PE_0 + P_2*mu_P_0 + P_3,
-M_2 - S_2 + 5966318017369273//4503599627370496,
E_2*S_0*mu_LE_0 + 2*E_1*S_1*mu_LE_0 + E_0*S_2*mu_LE_0 + 2*S_1^2*mu_LL_0 + 2*S_0*S_2*mu_LL_0 - S_2*delta_EL_0 + S_2*delta_LM_0 + S_3,
-S_2*delta_LM_0 + M_2*mu_M_0 + M_3,
-E_3 - 4586571004957137//9007199254740992,
6*E_1*E_2*mu_EE_0^2 + 2*E_0*E_3*mu_EE_0^2 - E_3*P_0*rho_E_0^2 - 3*E_2*P_1*rho_E_0^2 - 3*E_1*P_2*rho_E_0^2 - E_0*P_3*rho_E_0^2 - N_3*P_0*delta_NE_0 - 3*N_2*P_1*delta_NE_0 - 3*N_1*P_2*delta_NE_0 - N_0*P_3*delta_NE_0 + E_3*delta_EL_0 + E_4,
-P_3 - 5619286159605997//4503599627370496,
P_3*S_0*mu_PL_0^2 + 3*P_2*S_1*mu_PL_0^2 + 3*P_1*S_2*mu_PL_0^2 + P_0*S_3*mu_PL_0^2 - 6*P_1*P_2*rho_P_0^2 - 2*P_0*P_3*rho_P_0^2 + E_3*P_0*mu_PE_0 + 3*E_2*P_1*mu_PE_0 + 3*E_1*P_2*mu_PE_0 + E_0*P_3*mu_PE_0 + P_3*mu_P_0 + P_4,
-M_3 - S_3 - 7007185265083143//2251799813685248,
E_3*S_0*mu_LE_0 + 3*E_2*S_1*mu_LE_0 + 3*E_1*S_2*mu_LE_0 + E_0*S_3*mu_LE_0 + 6*S_1*S_2*mu_LL_0 + 2*S_0*S_3*mu_LL_0 - S_3*delta_EL_0 + S_3*delta_LM_0 + S_4,
-S_3*delta_LM_0 + M_3*mu_M_0 + M_4,
-P_4 + 2764663660811609//562949953421312,
P_4*S_0*mu_PL_0^2 + 4*P_3*S_1*mu_PL_0^2 + 6*P_2*S_2*mu_PL_0^2 + 4*P_1*S_3*mu_PL_0^2 + P_0*S_4*mu_PL_0^2 - 6*P_2^2*rho_P_0^2 - 8*P_1*P_3*rho_P_0^2 - 2*P_0*P_4*rho_P_0^2 + E_4*P_0*mu_PE_0 + 4*E_3*P_1*mu_PE_0 + 6*E_2*P_2*mu_PE_0 + 4*E_1*P_3*mu_PE_0 + E_0*P_4*mu_PE_0 + P_4*mu_P_0 + P_5,
-M_4 - S_4 + 3601204812855993//281474976710656,
-S_4*delta_LM_0 + M_4*mu_M_0 + M_5,
E_4*S_0*mu_LE_0 + 4*E_3*S_1*mu_LE_0 + 6*E_2*S_2*mu_LE_0 + 4*E_1*S_3*mu_LE_0 + E_0*S_4*mu_LE_0 + 6*S_2^2*mu_LL_0 + 8*S_1*S_3*mu_LL_0 + 2*S_0*S_4*mu_LL_0 - S_4*delta_EL_0 + S_4*delta_LM_0 + S_5,
-P_5 - 3052340128689447//140737488355328,
P_5*S_0*mu_PL_0^2 + 5*P_4*S_1*mu_PL_0^2 + 10*P_3*S_2*mu_PL_0^2 + 10*P_2*S_3*mu_PL_0^2 + 5*P_1*S_4*mu_PL_0^2 + P_0*S_5*mu_PL_0^2 - 20*P_2*P_3*rho_P_0^2 - 10*P_1*P_4*rho_P_0^2 - 2*P_0*P_5*rho_P_0^2 + E_5*P_0*mu_PE_0 + 5*E_4*P_1*mu_PE_0 + 10*E_3*P_2*mu_PE_0 + 10*E_2*P_3*mu_PE_0 + 5*E_1*P_4*mu_PE_0 + E_0*P_5*mu_PE_0 + P_5*mu_P_0 + P_6,
6*E_2^2*mu_EE_0^2 + 8*E_1*E_3*mu_EE_0^2 + 2*E_0*E_4*mu_EE_0^2 - E_4*P_0*rho_E_0^2 - 4*E_3*P_1*rho_E_0^2 - 6*E_2*P_2*rho_E_0^2 - 4*E_1*P_3*rho_E_0^2 - E_0*P_4*rho_E_0^2 - N_4*P_0*delta_NE_0 - 4*N_3*P_1*delta_NE_0 - 6*N_2*P_2*delta_NE_0 - 4*N_1*P_3*delta_NE_0 - N_0*P_4*delta_NE_0 + E_4*delta_EL_0 + E_5,
N_3*P_0*delta_NE_0^2 + 3*N_2*P_1*delta_NE_0^2 + 3*N_1*P_2*delta_NE_0^2 + N_0*P_3*delta_NE_0^2 + N_3*mu_N_0^2 + N_4,
]

sys_zp = map(f -> map_coefficients(c -> GF(2^30+3)(numerator(c)) // GF(2^30+3)(denominator(c)), f), sys);
return sys_zp
end

function nfkb_reduced_system()
R, (i1_0,i1a_0,k_prod_0,t1_0,t2_0,u_0,u_1,u_2,u_3,x10_0,x10_1,x10_2,x10_3,x10_4,x11_0,x11_1,x11_2,x11_3,x11_4,x12_0,x12_1,x12_2,x12_3,x13_0,x13_1,x13_2,x13_3,x13_4,x14_0,x14_1,x14_2,x14_3,x1_0,x1_1,x1_2,x1_3,x1_4,x2_0,x2_1,x2_2,x2_3,x2_4,x3_0,x3_1,x3_2,x3_3,x3_4,x4_0,x4_1,x4_2,x4_3,x5_0,x5_1,x5_2,x5_3,x6_0,x6_1,x6_2,x6_3,x6_4,x7_0,x7_1,x7_2,x7_3,x7_4,x7_5,x8_0,x8_1,x8_2,x8_3,x9_0,x9_1,x9_2) = polynomial_ring(QQ, ["i1_0","i1a_0","k_prod_0","t1_0","t2_0","u_0","u_1","u_2","u_3","x10_0","x10_1","x10_2","x10_3","x10_4","x11_0","x11_1","x11_2","x11_3","x11_4","x12_0","x12_1","x12_2","x12_3","x13_0","x13_1","x13_2","x13_3","x13_4","x14_0","x14_1","x14_2","x14_3","x1_0","x1_1","x1_2","x1_3","x1_4","x2_0","x2_1","x2_2","x2_3","x2_4","x3_0","x3_1","x3_2","x3_3","x3_4","x4_0","x4_1","x4_2","x4_3","x5_0","x5_1","x5_2","x5_3","x6_0","x6_1","x6_2","x6_3","x6_4","x7_0","x7_1","x7_2","x7_3","x7_4","x7_5","x8_0","x8_1","x8_2","x8_3","x9_0","x9_1","x9_2"], internal_ordering=:degrevlex)

sys = [
-x7_0 + 1087917959349171//1125899906842624,
-1//10*i1_0*x6_0 + 1//10*x11_0*x7_0 + x7_1,
-x9_0 + 8617739573352685//9007199254740992,
1//10*x7_0 + 1//10*x9_0 + x9_1 - 1//10,
-x1_0 - x2_0 - x3_0 + 7627910363403261//2251799813685248,
1//10*u_0*x2_0*x8_0 - 1//10*u_0*x1_0 + 1//10*x10_0*x2_0 + 1//10*x13_0*x2_0 - t1_0*x4_0 - t2_0*x5_0 + 1//5*x2_0 + x2_1,
1//10*u_0*x1_0 - k_prod_0 + 1//10*x1_0 + x1_1,
-1//10*u_0*x2_0*x8_0 - 1//10*x2_0 + 1//10*x3_0 + x3_1,
-x10_0 - x13_0 + 2059383086697551//1125899906842624,
1//10*x13_0*x2_0 - 1//10*x10_0*x6_0 + 1//10*x13_0 + x13_1 - 1//10*x14_0,
i1a_0*x10_0 + 1//10*x10_0*x2_0 + 1//10*x10_0*x6_0 + 1//10*x10_0 + x10_1 - 1//10*x11_0 - 1//10*x12_0,
-x2_0 + 2801824555072859//2251799813685248,
-x12_0 + 4698329468064649//4503599627370496,
1//10*x12_0 + x12_1 - 1//10*x7_0 - 1//10,
-u_0 + 6530219459687117//4503599627370496,
u_1 - 1,
-x7_1 - 1265874459032373//18014398509481984,
-1//10*i1_0*x6_1 + 1//10*x11_1*x7_0 + 1//10*x11_0*x7_1 + x7_2,
-t2_0*x5_0 + i1_0*x6_0 + 1//10*x10_0*x6_0 - 1//10*x13_0 + x6_1,
-1//10*i1a_0*x10_0 + 1//10*x11_0*x7_0 + 27191544919973//2719154491997299*x11_0 + x11_1,
-x1_1 - x2_1 - x3_1 + 6653335799572365//9007199254740992,
1//10*u_1*x1_0 + 1//10*u_0*x1_1 + 1//10*x1_1 + x1_2,
1//10*u_1*x2_0*x8_0 + 1//10*u_0*x2_1*x8_0 + 1//10*u_0*x2_0*x8_1 - 1//10*u_1*x1_0 - 1//10*u_0*x1_1 + 1//10*x10_1*x2_0 + 1//10*x13_1*x2_0 + 1//10*x10_0*x2_1 + 1//10*x13_0*x2_1 - t1_0*x4_1 - t2_0*x5_1 + 1//5*x2_1 + x2_2,
-1//10*u_1*x2_0*x8_0 - 1//10*u_0*x2_1*x8_0 - 1//10*u_0*x2_0*x8_1 - 1//10*x2_1 + 1//10*x3_1 + x3_2,
1//10*x8_0 + x8_1 - 1//10*x9_0,
-1//10*x13_0*x2_0 + t2_0*x5_0 + x5_1,
-1//10*x10_0*x2_0 + t1_0*x4_0 + x4_1,
-x10_1 - x13_1 - 1604424493326049//4503599627370496,
i1a_0*x10_1 + 1//10*x10_1*x2_0 + 1//10*x10_0*x2_1 + 1//10*x10_1*x6_0 + 1//10*x10_0*x6_1 + 1//10*x10_1 + x10_2 - 1//10*x11_1 - 1//10*x12_1,
1//10*x13_1*x2_0 + 1//10*x13_0*x2_1 - 1//10*x10_1*x6_0 - 1//10*x10_0*x6_1 + 1//10*x13_1 + x13_2 - 1//10*x14_1,
-1//10*x11_0*x7_0 + 27191544919973//2719154491997299*x14_0 + x14_1,
-x2_1 + 7283230198371657//18014398509481984,
-x7_2 + 5615679340511373//288230376151711744,
-1//10*i1_0*x6_2 + 1//10*x11_2*x7_0 + 1//5*x11_1*x7_1 + 1//10*x11_0*x7_2 + x7_3,
-t2_0*x5_1 + 1//10*x10_1*x6_0 + i1_0*x6_1 + 1//10*x10_0*x6_1 - 1//10*x13_1 + x6_2,
-1//10*i1a_0*x10_1 + 1//10*x11_1*x7_0 + 1//10*x11_0*x7_1 + 27191544919973//2719154491997299*x11_1 + x11_2,
-x1_2 - x2_2 - x3_2 - 2154146894998627//4503599627370496,
-1//10*u_2*x2_0*x8_0 - 1//5*u_1*x2_1*x8_0 - 1//10*u_0*x2_2*x8_0 - 1//5*u_1*x2_0*x8_1 - 1//5*u_0*x2_1*x8_1 - 1//10*u_0*x2_0*x8_2 - 1//10*x2_2 + 1//10*x3_2 + x3_3,
1//10*u_2*x2_0*x8_0 + 1//5*u_1*x2_1*x8_0 + 1//10*u_0*x2_2*x8_0 + 1//5*u_1*x2_0*x8_1 + 1//5*u_0*x2_1*x8_1 + 1//10*u_0*x2_0*x8_2 - 1//10*u_2*x1_0 - 1//5*u_1*x1_1 - 1//10*u_0*x1_2 + 1//10*x10_2*x2_0 + 1//10*x13_2*x2_0 + 1//5*x10_1*x2_1 + 1//5*x13_1*x2_1 + 1//10*x10_0*x2_2 + 1//10*x13_0*x2_2 - t1_0*x4_2 - t2_0*x5_2 + 1//5*x2_2 + x2_3,
1//10*u_2*x1_0 + 1//5*u_1*x1_1 + 1//10*u_0*x1_2 + 1//10*x1_2 + x1_3,
u_2,
1//10*x8_1 + x8_2 - 1//10*x9_1,
-1//10*x13_1*x2_0 - 1//10*x13_0*x2_1 + t2_0*x5_1 + x5_2,
-1//10*x10_1*x2_0 - 1//10*x10_0*x2_1 + t1_0*x4_1 + x4_2,
-x10_2 - x13_2 + 2081185062807031//18014398509481984,
1//10*x13_2*x2_0 + 1//5*x13_1*x2_1 + 1//10*x13_0*x2_2 - 1//10*x10_2*x6_0 - 1//5*x10_1*x6_1 - 1//10*x10_0*x6_2 + 1//10*x13_2 + x13_3 - 1//10*x14_2,
i1a_0*x10_2 + 1//10*x10_2*x2_0 + 1//5*x10_1*x2_1 + 1//10*x10_0*x2_2 + 1//10*x10_2*x6_0 + 1//5*x10_1*x6_1 + 1//10*x10_0*x6_2 + 1//10*x10_2 + x10_3 - 1//10*x11_2 - 1//10*x12_2,
1//10*x12_1 + x12_2 - 1//10*x7_1,
-1//10*x11_1*x7_0 - 1//10*x11_0*x7_1 + 27191544919973//2719154491997299*x14_1 + x14_2,
-x2_2 - 76007005474747//140737488355328,
-x7_3 - 4931190233469493//576460752303423488,
-1//10*i1_0*x6_3 + 1//10*x11_3*x7_0 + 3//10*x11_2*x7_1 + 3//10*x11_1*x7_2 + 1//10*x11_0*x7_3 + x7_4,
-1//10*i1a_0*x10_2 + 1//10*x11_2*x7_0 + 1//5*x11_1*x7_1 + 1//10*x11_0*x7_2 + 27191544919973//2719154491997299*x11_2 + x11_3,
-t2_0*x5_2 + 1//10*x10_2*x6_0 + 1//5*x10_1*x6_1 + i1_0*x6_2 + 1//10*x10_0*x6_2 - 1//10*x13_2 + x6_3,
-x1_3 - x2_3 - x3_3 + 3500819284168905//9007199254740992,
-1//10*u_3*x2_0*x8_0 - 3//10*u_2*x2_1*x8_0 - 3//10*u_1*x2_2*x8_0 - 1//10*u_0*x2_3*x8_0 - 3//10*u_2*x2_0*x8_1 - 3//5*u_1*x2_1*x8_1 - 3//10*u_0*x2_2*x8_1 - 3//10*u_1*x2_0*x8_2 - 3//10*u_0*x2_1*x8_2 - 1//10*u_0*x2_0*x8_3 - 1//10*x2_3 + 1//10*x3_3 + x3_4,
1//10*u_3*x2_0*x8_0 + 3//10*u_2*x2_1*x8_0 + 3//10*u_1*x2_2*x8_0 + 1//10*u_0*x2_3*x8_0 + 3//10*u_2*x2_0*x8_1 + 3//5*u_1*x2_1*x8_1 + 3//10*u_0*x2_2*x8_1 + 3//10*u_1*x2_0*x8_2 + 3//10*u_0*x2_1*x8_2 + 1//10*u_0*x2_0*x8_3 - 1//10*u_3*x1_0 - 3//10*u_2*x1_1 - 3//10*u_1*x1_2 - 1//10*u_0*x1_3 + 1//10*x10_3*x2_0 + 1//10*x13_3*x2_0 + 3//10*x10_2*x2_1 + 3//10*x13_2*x2_1 + 3//10*x10_1*x2_2 + 3//10*x13_1*x2_2 + 1//10*x10_0*x2_3 + 1//10*x13_0*x2_3 - t1_0*x4_3 - t2_0*x5_3 + 1//5*x2_3 + x2_4,
1//10*u_3*x1_0 + 3//10*u_2*x1_1 + 3//10*u_1*x1_2 + 1//10*u_0*x1_3 + 1//10*x1_3 + x1_4,
u_3,
-1//10*x10_2*x2_0 - 1//5*x10_1*x2_1 - 1//10*x10_0*x2_2 + t1_0*x4_2 + x4_3,
-1//10*x13_2*x2_0 - 1//5*x13_1*x2_1 - 1//10*x13_0*x2_2 + t2_0*x5_2 + x5_3,
1//10*x8_2 + x8_3 - 1//10*x9_2,
1//10*x7_1 + 1//10*x9_1 + x9_2,
-x10_3 - x13_3 + 7413820825919623//144115188075855872,
1//10*x13_3*x2_0 + 3//10*x13_2*x2_1 + 3//10*x13_1*x2_2 + 1//10*x13_0*x2_3 - 1//10*x10_3*x6_0 - 3//10*x10_2*x6_1 - 3//10*x10_1*x6_2 - 1//10*x10_0*x6_3 + 1//10*x13_3 + x13_4 - 1//10*x14_3,
i1a_0*x10_3 + 1//10*x10_3*x2_0 + 3//10*x10_2*x2_1 + 3//10*x10_1*x2_2 + 1//10*x10_0*x2_3 + 1//10*x10_3*x6_0 + 3//10*x10_2*x6_1 + 3//10*x10_1*x6_2 + 1//10*x10_0*x6_3 + 1//10*x10_3 + x10_4 - 1//10*x11_3 - 1//10*x12_3,
1//10*x12_2 + x12_3 - 1//10*x7_2,
-1//10*x11_2*x7_0 - 1//5*x11_1*x7_1 - 1//10*x11_0*x7_2 + 27191544919973//2719154491997299*x14_2 + x14_3,
-x2_3 + 513728119802863//1125899906842624,
-x7_4 + 6517251002848763//1152921504606846976,
-1//10*i1_0*x6_4 + 1//10*x11_4*x7_0 + 2//5*x11_3*x7_1 + 3//5*x11_2*x7_2 + 2//5*x11_1*x7_3 + 1//10*x11_0*x7_4 + x7_5,
-t2_0*x5_3 + 1//10*x10_3*x6_0 + 3//10*x10_2*x6_1 + 3//10*x10_1*x6_2 + i1_0*x6_3 + 1//10*x10_0*x6_3 - 1//10*x13_3 + x6_4,
-1//10*i1a_0*x10_3 + 1//10*x11_3*x7_0 + 3//10*x11_2*x7_1 + 3//10*x11_1*x7_2 + 1//10*x11_0*x7_3 + 27191544919973//2719154491997299*x11_3 + x11_4,
]

sys_zp = map(f -> map_coefficients(c -> GF(2^30+3)(numerator(c)) // GF(2^30+3)(denominator(c)), f), sys);
return sys_zp
end

function nike4_system()
    R, (x2, x3, x4, a12, a13, a14, a22, a23, a24, Tag_1, Tag_16, Tag_17, Tag_18, Tag_20, Tag_21, Tag_22, Tag_2, Tag_4, Tag_6, Tag_7, Tag_14, Tag_15, Tag_19, Tag_8, Tag_9, Tag_10, Tag_11, Sat_1) = polynomial_ring(AbstractAlgebra.GF(2^30+3), [:x2, :x3, :x4, :a12, :a13, :a14, :a22, :a23, :a24, :Tag_1, :Tag_16, :Tag_17, :Tag_18, :Tag_20, :Tag_21, :Tag_22, :Tag_2, :Tag_4, :Tag_6, :Tag_7, :Tag_14, :Tag_15, :Tag_19, :Tag_8, :Tag_9, :Tag_10, :Tag_11, :Sat_1], internal_ordering=:degrevlex)

    sys = [x2*a12*Tag_17*Tag_22 - x2*a12*Tag_18*Tag_21 - x2*a13*Tag_16*Tag_22 + x2*a13*Tag_18*Tag_20 + x2*a14*Tag_16*Tag_21 - x2*a14*Tag_17*Tag_20 - x3*a12*a23*Tag_22 + x3*a12*a24*Tag_21 + x3*a13*a22*Tag_22 - x3*a13*a24*Tag_20 - x3*a14*a22*Tag_21 + x3*a14*a23*Tag_20 + x4*a12*a23*Tag_18 - x4*a12*a24*Tag_17 - x4*a13*a22*Tag_18 + x4*a13*a24*Tag_16 + x4*a14*a22*Tag_17 - x4*a14*a23*Tag_16 - Tag_1, x2*a12*Tag_17 + x2*a12*Tag_22 - x2*a13*Tag_16 - x2*a14*Tag_20 - x3*a12*a23 + x3*a13*a22 + x3*a13*Tag_22 - x3*a14*Tag_21 - x4*a12*a24 - x4*a13*Tag_18 + x4*a14*a22 + x4*a14*Tag_17 - Tag_2, x2*a12 + x3*a13 + x4*a14 - Tag_4, a22 - Tag_6 + Tag_17 + Tag_22, a12*Tag_14 + a13*Tag_15 + a14*Tag_19 - Tag_7, a22*Tag_17 + a22*Tag_22 - a23*Tag_16 - a24*Tag_20 - Tag_8 + Tag_17*Tag_22 - Tag_18*Tag_21, a22*Tag_17*Tag_22 - a22*Tag_18*Tag_21 - a23*Tag_16*Tag_22 + a23*Tag_18*Tag_20 + a24*Tag_16*Tag_21 - a24*Tag_17*Tag_20 - Tag_9, a12*a22*Tag_14 + a12*a23*Tag_15 + a12*a24*Tag_19 + a13*Tag_14*Tag_16 + a13*Tag_15*Tag_17 + a13*Tag_18*Tag_19 + a14*Tag_14*Tag_20 + a14*Tag_15*Tag_21 + a14*Tag_19*Tag_22 - Tag_10, -a12*a23*Tag_15*Tag_22 + a12*a23*Tag_18*Tag_19 + a12*a24*Tag_15*Tag_21 - a12*a24*Tag_17*Tag_19 + a12*Tag_14*Tag_17*Tag_22 - a12*Tag_14*Tag_18*Tag_21 + a13*a22*Tag_15*Tag_22 - a13*a22*Tag_18*Tag_19 - a13*a24*Tag_15*Tag_20 + a13*a24*Tag_16*Tag_19 - a13*Tag_14*Tag_16*Tag_22 + a13*Tag_14*Tag_18*Tag_20 - a14*a22*Tag_15*Tag_21 + a14*a22*Tag_17*Tag_19 + a14*a23*Tag_15*Tag_20 - a14*a23*Tag_16*Tag_19 + a14*Tag_14*Tag_16*Tag_21 - a14*Tag_14*Tag_17*Tag_20 - Tag_11, Sat_1 - 1]

    return sys
end

function nike4_ext_system()
    R, (x2, x3, x4, a12, a13, a14, a22, a23, a24, Tag_1, Tag_16, Tag_17, Tag_18, Tag_20, Tag_21, Tag_22, Tag_2, Tag_4, Tag_6, Tag_7, Tag_14, Tag_15, Tag_19, Tag_8, Tag_9, Tag_10, Tag_11, Sat_1, Dummy_1, Dummy_2, Dummy_3, Dummy_4, Dummy_5, Dummy_6, Dummy_7) = polynomial_ring(AbstractAlgebra.GF(2^30+3), [:x2, :x3, :x4, :a12, :a13, :a14, :a22, :a23, :a24, :Tag_1, :Tag_16, :Tag_17, :Tag_18, :Tag_20, :Tag_21, :Tag_22, :Tag_2, :Tag_4, :Tag_6, :Tag_7, :Tag_14, :Tag_15, :Tag_19, :Tag_8, :Tag_9, :Tag_10, :Tag_11, :Sat_1, :Dummy_1, :Dummy_2, :Dummy_3, :Dummy_4, :Dummy_5, :Dummy_6, :Dummy_7], internal_ordering=:degrevlex)

    sys = [x2*a12*Tag_17*Tag_22 - x2*a12*Tag_18*Tag_21 - x2*a13*Tag_16*Tag_22 + x2*a13*Tag_18*Tag_20 + x2*a14*Tag_16*Tag_21 - x2*a14*Tag_17*Tag_20 - x3*a12*a23*Tag_22 + x3*a12*a24*Tag_21 + x3*a13*a22*Tag_22 - x3*a13*a24*Tag_20 - x3*a14*a22*Tag_21 + x3*a14*a23*Tag_20 + x4*a12*a23*Tag_18 - x4*a12*a24*Tag_17 - x4*a13*a22*Tag_18 + x4*a13*a24*Tag_16 + x4*a14*a22*Tag_17 - x4*a14*a23*Tag_16 - Tag_1, x2*a12*Tag_17 + x2*a12*Tag_22 - x2*a13*Tag_16 - x2*a14*Tag_20 - x3*a12*a23 + x3*a13*a22 + x3*a13*Tag_22 - x3*a14*Tag_21 - x4*a12*a24 - x4*a13*Tag_18 + x4*a14*a22 + x4*a14*Tag_17 - Tag_2, x2*a12 + x3*a13 + x4*a14 - Tag_4, a22 - Tag_6 + Tag_17 + Tag_22, a12*Tag_14 + a13*Tag_15 + a14*Tag_19 - Tag_7, a22*Tag_17 + a22*Tag_22 - a23*Tag_16 - a24*Tag_20 - Tag_8 + Tag_17*Tag_22 - Tag_18*Tag_21, a22*Tag_17*Tag_22 - a22*Tag_18*Tag_21 - a23*Tag_16*Tag_22 + a23*Tag_18*Tag_20 + a24*Tag_16*Tag_21 - a24*Tag_17*Tag_20 - Tag_9, a12*a22*Tag_14 + a12*a23*Tag_15 + a12*a24*Tag_19 + a13*Tag_14*Tag_16 + a13*Tag_15*Tag_17 + a13*Tag_18*Tag_19 + a14*Tag_14*Tag_20 + a14*Tag_15*Tag_21 + a14*Tag_19*Tag_22 - Tag_10, -a12*a23*Tag_15*Tag_22 + a12*a23*Tag_18*Tag_19 + a12*a24*Tag_15*Tag_21 - a12*a24*Tag_17*Tag_19 + a12*Tag_14*Tag_17*Tag_22 - a12*Tag_14*Tag_18*Tag_21 + a13*a22*Tag_15*Tag_22 - a13*a22*Tag_18*Tag_19 - a13*a24*Tag_15*Tag_20 + a13*a24*Tag_16*Tag_19 - a13*Tag_14*Tag_16*Tag_22 + a13*Tag_14*Tag_18*Tag_20 - a14*a22*Tag_15*Tag_21 + a14*a22*Tag_17*Tag_19 + a14*a23*Tag_15*Tag_20 - a14*a23*Tag_16*Tag_19 + a14*Tag_14*Tag_16*Tag_21 - a14*Tag_14*Tag_17*Tag_20 - Tag_11, Sat_1 - 1]

    return sys
end

function multi(nbits; np=AbstractAlgebra)
    R, (x1, x2, x3, x4) =
        polynomial_ring(np.QQ, ["x1", "x2", "x3", "x4"], internal_ordering=:degrevlex)
    nbits_per_prime = 31
    nprimes = max(div(nbits, nbits_per_prime), 1)
    N = prod(map(BigInt, Primes.nextprimes(2^31 - 100, nprimes)))
    system = [
        x1 + x2 + x3 + x4,
        x1 * x2 + x1 * x3 + x1 * x4 + x2 * x3 + x2 * x4 + x3 * x4,
        x1 * x2 * x3 + x1 * x2 * x4 + x1 * x3 * x4 + x2 * x3 * x4,
        x1 * x2 * x3 * x4 + N
    ]
    system
end

function n_variable_set(n; internal_ordering=:degrevlex, k=GF(2^31 - 1))
    R, x = polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=internal_ordering)
    f = [sum(prod(x[i:(n - kk)], init=1) for i in 1:(kk + 1)) for kk in 0:(n - 1)]
    f
end

systems = Dict{String, Any}(
    "hexapod" => Groebner.Examples.hexapod(),
    "chandra9" => Groebner.Examples.chandran(9),
    "chandra10" => Groebner.Examples.chandran(10),
    "nvar-150" => n_variable_set(150),
    "nvar-200" => n_variable_set(200),
    "multi-10k" => multi(10_000),
    "multi-100k" => multi(100_000),
    "cyclic8" => Groebner.Examples.cyclicn(8, k=GF(2^30 + 3)),
    "cyclic9" => Groebner.Examples.cyclicn(9, k=GF(2^30 + 3)),
    "katsura10" => Groebner.Examples.katsuran(10, k=GF(2^30 + 3)),
    "katsura11" => Groebner.Examples.katsuran(11, k=GF(2^30 + 3)),
    "katsura12" => Groebner.Examples.katsuran(12, k=GF(2^30 + 3)),
    "noon9" => Groebner.Examples.noonn(9, k=GF(2^30 + 3)),
    "noon10" => Groebner.Examples.noonn(10, k=GF(2^30 + 3)),
    "eco13" => Groebner.Examples.econ(13, k=GF(2^30 + 3)),
    "eco14" => Groebner.Examples.econ(14, k=GF(2^30 + 3)),
    "jason210" => Groebner.Examples.jason210(k=GF(2^30 + 3)),
    "yang1" => Groebner.Examples.yang1(k=GF(2^30 + 3)),
    "bayes148" => Groebner.Examples.bayes148(k=GF(2^30 + 3)),
    "mayr42" => Groebner.Examples.mayr42(k=GF(2^30 + 3)),
    "nike4" => nike4_system(),
    "nike4_ext" => nike4_ext_system(),
    "crauste" => crauste_system(),
    "nfkb_reduced" => nfkb_reduced_system(),
    "akt" => akt_system(),
)

problems = [
    ("hexapod", (; kws...) -> @belapsed groebner($(systems["hexapod"]); $(kws)...)),
    ("nvar-200", (; kws...) -> @belapsed groebner($(systems["nvar-200"]); $(kws)...)),
    ("multi-10k", (; kws...) -> @belapsed groebner($(systems["multi-10k"]); $(kws)...)),
    ("multi-100k", (; kws...) -> @belapsed groebner($(systems["multi-100k"]); $(kws)...)),
    ("chandra9", (; kws...) -> @belapsed groebner($(systems["chandra9"]); $(kws)...)),
    ("chandra10", (; kws...) -> @elapsed groebner(systems["chandra10"]; kws...)),
    ("cyclic8", (; kws...) -> @belapsed groebner($(systems["cyclic8"]); $(kws)...)),
    ("cyclic9", (; kws...) -> @elapsed groebner(systems["cyclic9"]; kws...)),
    ("katsura11", (; kws...) -> @elapsed groebner(systems["katsura11"]; kws...)),
    ("katsura12", (; kws...) -> @elapsed groebner(systems["katsura12"]; kws...)),
    ("noon9", (; kws...) -> @elapsed groebner(systems["noon9"]; kws...)),
    ("noon10", (; kws...) -> @elapsed groebner(systems["noon10"]; kws...)),
    ("eco13", (; kws...) -> @elapsed groebner(systems["eco13"]; kws...)),
    ("eco14", (; kws...) -> @elapsed groebner(systems["eco14"]; kws...)),
    ("jason210", (; kws...) -> @belapsed groebner($(systems["jason210"]); $(kws)...)),
    ("yang1", (; kws...) -> @elapsed groebner(systems["yang1"]; kws...)),
    ("bayes148", (; kws...) -> @elapsed groebner(systems["bayes148"]; kws...)),
    ("mayr42", (; kws...) -> @elapsed groebner(systems["mayr42"]; kws...)),
    ("nike4", (; kws...) -> @elapsed groebner(systems["nike4"]; kws...)),
    ("nike4_ext", (; kws...) -> @elapsed groebner(systems["nike4_ext"]; kws...)),
    ("crauste", (; kws...) -> @belapsed groebner($(systems["crauste"]); $(kws)...)),
    ("nfkb_reduced", (; kws...) -> @belapsed groebner($(systems["nfkb_reduced"]); $(kws)...)),
    ("akt", (; kws...) -> @elapsed groebner(systems["akt"]; kws...)),
    ("katsura10+nf", (; kws...) -> begin gb = groebner(systems["katsura10"]); @belapsed (normalform($gb, $(systems["katsura10"]); $(kws)...); normalform($gb, $gb; $(kws)...)) end),
    ("katsura11+nf", (; kws...) -> begin gb = groebner(systems["katsura11"]); @belapsed (normalform($gb, $(systems["katsura11"]); $(kws)...); normalform($gb, $gb; $(kws)...)) end),
    ("katsura10+lr", (; kws...) -> begin t, _ = groebner_learn(systems["katsura10"]; kws...); @belapsed groebner_apply!($t, $(systems["katsura10"]); $(kws)...) end),
    ("katsura11+lr", (; kws...) -> begin t, _ = groebner_learn(systems["katsura11"]; kws...); @belapsed groebner_apply!($t, $(systems["katsura11"]); $(kws)...) end),
    ("akt+lr", (; kws...) -> begin t, _ = groebner_learn(systems["akt"]; (kws)...); @belapsed groebner_apply!($t, $(systems["akt"]); $(kws)...) end),
    ("akt+gt", (; kws...) -> begin gb = groebner(systems["akt"]); @belapsed isgroebner($gb; $(kws)...) end),
]

options = (
    (monoms=:dense,),
    (monoms=:packed,),
    (monoms=:fixed,),
    (monoms=:fixed2,),
    (monoms=:fixednodeg,),
    (monoms=:nibble,),
)

@printf("Timings in seconds\n")
@printf("%-16s %4s %10s %10s %10s %10s %12s %10s\n", "name", "N", ":dense", ":packed", ":fixed", ":fixed2", ":fixednodeg", ":nibble")
@printf("%-16s %4s %10s %10s %10s %10s %12s %10s\n", "-"^16, "-"^4, "-"^10, "-"^10, "-"^10, "-"^10, "-"^12, "-"^10)

for (name, problem) in problems
    sys = systems[split(name,"+")[1]]
    n = ngens(parent(sys[1]))
    @printf("%-16s %4d", name, n)

    for kwargs in options
        GC.gc()
        t = try
            problem(; kwargs...)
        catch e
            NaN
        end

        if isnan(t)
            @printf(" %10s", "ERR")
        else
            @printf(" %10.3f", t)
        end
    end

    @printf("\n")
end

versioninfo()
