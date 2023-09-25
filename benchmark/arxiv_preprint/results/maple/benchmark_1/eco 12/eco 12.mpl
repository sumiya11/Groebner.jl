# eco 12
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		x1*x2*x12 + x1*x12 + x2*x3*x12 + x3*x4*x12 + x4*x5*x12 + x5*x6*x12 + x6*x7*x12 + x7*x8*x12 + x8*x9*x12 + x9*x10*x12 + x10*x11*x12 + 2147483646,
		x1*x3*x12 + x2*x4*x12 + x2*x12 + x3*x5*x12 + x4*x6*x12 + x5*x7*x12 + x6*x8*x12 + x7*x9*x12 + x8*x10*x12 + x9*x11*x12 + 2147483645,
		x1*x4*x12 + x2*x5*x12 + x3*x6*x12 + x3*x12 + x4*x7*x12 + x5*x8*x12 + x6*x9*x12 + x7*x10*x12 + x8*x11*x12 + 2147483644,
		x1*x5*x12 + x2*x6*x12 + x3*x7*x12 + x4*x8*x12 + x4*x12 + x5*x9*x12 + x6*x10*x12 + x7*x11*x12 + 2147483643,
		x1*x6*x12 + x2*x7*x12 + x3*x8*x12 + x4*x9*x12 + x5*x10*x12 + x5*x12 + x6*x11*x12 + 2147483642,
		x1*x7*x12 + x2*x8*x12 + x3*x9*x12 + x4*x10*x12 + x5*x11*x12 + x6*x12 + 2147483641,
		x1*x8*x12 + x2*x9*x12 + x3*x10*x12 + x4*x11*x12 + x7*x12 + 2147483640,
		x1*x9*x12 + x2*x10*x12 + x3*x11*x12 + x8*x12 + 2147483639,
		x1*x10*x12 + x2*x11*x12 + x9*x12 + 2147483638,
		x1*x11*x12 + x10*x12 + 2147483637,
		x11*x12 + 2147483636,
		x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + 1
	]:
	print("Running eco 12");
	st := time[real]():
	Groebner[Basis](J, tdeg(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12), method=fgb, characteristic=2147483647):
	print("eco 12: ", time[real]() - st);
	runtime := min(runtime, time[real]() - st);
end do;

timings_fn := "/home/ademin/Groebner.jl/benchmark/arxiv_preprint/results/maple/benchmark_1/eco 12/timings":
FileTools[Text][WriteLine](timings_fn, "eco 12");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime)));
