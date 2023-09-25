# katsura 12
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		x0^2 + 2147483646*x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2 + 2*x8^2 + 2*x9^2 + 2*x10^2 + 2*x11^2 + 2*x12^2,
		2*x0*x1 + 2*x1*x2 + 2147483646*x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7 + 2*x7*x8 + 2*x8*x9 + 2*x9*x10 + 2*x10*x11 + 2*x11*x12,
		2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 + 2147483646*x2 + 2*x3*x5 + 2*x4*x6 + 2*x5*x7 + 2*x6*x8 + 2*x7*x9 + 2*x8*x10 + 2*x9*x11 + 2*x10*x12,
		2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 + 2*x3*x6 + 2147483646*x3 + 2*x4*x7 + 2*x5*x8 + 2*x6*x9 + 2*x7*x10 + 2*x8*x11 + 2*x9*x12,
		2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 + 2*x2*x6 + 2*x3*x7 + 2*x4*x8 + 2147483646*x4 + 2*x5*x9 + 2*x6*x10 + 2*x7*x11 + 2*x8*x12,
		2*x0*x5 + 2*x1*x4 + 2*x1*x6 + 2*x2*x3 + 2*x2*x7 + 2*x3*x8 + 2*x4*x9 + 2*x5*x10 + 2147483646*x5 + 2*x6*x11 + 2*x7*x12,
		2*x0*x6 + 2*x1*x5 + 2*x1*x7 + 2*x2*x4 + 2*x2*x8 + x3^2 + 2*x3*x9 + 2*x4*x10 + 2*x5*x11 + 2*x6*x12 + 2147483646*x6,
		2*x0*x7 + 2*x1*x6 + 2*x1*x8 + 2*x2*x5 + 2*x2*x9 + 2*x3*x4 + 2*x3*x10 + 2*x4*x11 + 2*x5*x12 + 2147483646*x7,
		2*x0*x8 + 2*x1*x7 + 2*x1*x9 + 2*x2*x6 + 2*x2*x10 + 2*x3*x5 + 2*x3*x11 + x4^2 + 2*x4*x12 + 2147483646*x8,
		2*x0*x9 + 2*x1*x8 + 2*x1*x10 + 2*x2*x7 + 2*x2*x11 + 2*x3*x6 + 2*x3*x12 + 2*x4*x5 + 2147483646*x9,
		2*x0*x10 + 2*x1*x9 + 2*x1*x11 + 2*x2*x8 + 2*x2*x12 + 2*x3*x7 + 2*x4*x6 + x5^2 + 2147483646*x10,
		2*x0*x11 + 2*x1*x10 + 2*x1*x12 + 2*x2*x9 + 2*x3*x8 + 2*x4*x7 + 2*x5*x6 + 2147483646*x11,
		x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 + 2*x8 + 2*x9 + 2*x10 + 2*x11 + 2*x12 + 2147483646
	]:
	print("Running katsura 12");
	st := time[real]():
	Groebner[Basis](J, tdeg(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12), method=fgb, characteristic=2147483647):
	print("katsura 12: ", time[real]() - st);
	runtime := min(runtime, time[real]() - st);
end do;

timings_fn := "/home/ademin/Groebner.jl/benchmark/arxiv_preprint/results/maple/benchmark_1/katsura 12/timings":
FileTools[Text][WriteLine](timings_fn, "katsura 12");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime)));
