# cyclic 9
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9,
		z1*z2 + z1*z9 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z7 + z7*z8 + z8*z9,
		z1*z2*z3 + z1*z2*z9 + z1*z8*z9 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z7 + z6*z7*z8 + z7*z8*z9,
		z1*z2*z3*z4 + z1*z2*z3*z9 + z1*z2*z8*z9 + z1*z7*z8*z9 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z7 + z5*z6*z7*z8 + z6*z7*z8*z9,
		z1*z2*z3*z4*z5 + z1*z2*z3*z4*z9 + z1*z2*z3*z8*z9 + z1*z2*z7*z8*z9 + z1*z6*z7*z8*z9 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z7 + z4*z5*z6*z7*z8 + z5*z6*z7*z8*z9,
		z1*z2*z3*z4*z5*z6 + z1*z2*z3*z4*z5*z9 + z1*z2*z3*z4*z8*z9 + z1*z2*z3*z7*z8*z9 + z1*z2*z6*z7*z8*z9 + z1*z5*z6*z7*z8*z9 + z2*z3*z4*z5*z6*z7 + z3*z4*z5*z6*z7*z8 + z4*z5*z6*z7*z8*z9,
		z1*z2*z3*z4*z5*z6*z7 + z1*z2*z3*z4*z5*z6*z9 + z1*z2*z3*z4*z5*z8*z9 + z1*z2*z3*z4*z7*z8*z9 + z1*z2*z3*z6*z7*z8*z9 + z1*z2*z5*z6*z7*z8*z9 + z1*z4*z5*z6*z7*z8*z9 + z2*z3*z4*z5*z6*z7*z8 + z3*z4*z5*z6*z7*z8*z9,
		z1*z2*z3*z4*z5*z6*z7*z8 + z1*z2*z3*z4*z5*z6*z7*z9 + z1*z2*z3*z4*z5*z6*z8*z9 + z1*z2*z3*z4*z5*z7*z8*z9 + z1*z2*z3*z4*z6*z7*z8*z9 + z1*z2*z3*z5*z6*z7*z8*z9 + z1*z2*z4*z5*z6*z7*z8*z9 + z1*z3*z4*z5*z6*z7*z8*z9 + z2*z3*z4*z5*z6*z7*z8*z9,
		z1*z2*z3*z4*z5*z6*z7*z8*z9 + 2147483646
	]:
	print("Running cyclic 9");
	st := time[real]():
	Groebner[Basis](J, tdeg(z1, z2, z3, z4, z5, z6, z7, z8, z9), method=fgb, characteristic=2147483647):
	print("cyclic 9: ", time[real]() - st);
	runtime := min(runtime, time[real]() - st);
end do;

timings_fn := "/home/ademin/Groebner.jl/benchmark/arxiv_preprint/results/maple/benchmark_1/cyclic 9/timings":
FileTools[Text][WriteLine](timings_fn, "cyclic 9");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime)));
