# cyclic 8
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8,
		z1*z2 + z1*z8 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z7 + z7*z8,
		z1*z2*z3 + z1*z2*z8 + z1*z7*z8 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z7 + z6*z7*z8,
		z1*z2*z3*z4 + z1*z2*z3*z8 + z1*z2*z7*z8 + z1*z6*z7*z8 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z7 + z5*z6*z7*z8,
		z1*z2*z3*z4*z5 + z1*z2*z3*z4*z8 + z1*z2*z3*z7*z8 + z1*z2*z6*z7*z8 + z1*z5*z6*z7*z8 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z7 + z4*z5*z6*z7*z8,
		z1*z2*z3*z4*z5*z6 + z1*z2*z3*z4*z5*z8 + z1*z2*z3*z4*z7*z8 + z1*z2*z3*z6*z7*z8 + z1*z2*z5*z6*z7*z8 + z1*z4*z5*z6*z7*z8 + z2*z3*z4*z5*z6*z7 + z3*z4*z5*z6*z7*z8,
		z1*z2*z3*z4*z5*z6*z7 + z1*z2*z3*z4*z5*z6*z8 + z1*z2*z3*z4*z5*z7*z8 + z1*z2*z3*z4*z6*z7*z8 + z1*z2*z3*z5*z6*z7*z8 + z1*z2*z4*z5*z6*z7*z8 + z1*z3*z4*z5*z6*z7*z8 + z2*z3*z4*z5*z6*z7*z8,
		z1*z2*z3*z4*z5*z6*z7*z8 + 2147483646
	]:
	print("Running cyclic 8");
	st := time[real]():
	Groebner[Basis](J, tdeg(z1, z2, z3, z4, z5, z6, z7, z8), method=fgb, characteristic=2147483647):
	print("cyclic 8: ", time[real]() - st);
	runtime := min(runtime, time[real]() - st);
end do;

logs_fn := "/home/ademin/Groebner.jl/benchmark/arxiv_preprint/results/maple/benchmark_1/cyclic 8/logs":
FileTools[Text][WriteLine](logs_fn, "cyclic 8");
FileTools[Text][WriteLine](logs_fn, cat("total_time, ", String(runtime)));
