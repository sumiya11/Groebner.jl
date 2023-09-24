# dummy
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		x^2 + y^2 + 2147483646,
		x + y
	]:
	print("Running dummy");
	st := time[real]():
	Groebner[Basis](J, tdeg(x, y), method=fgb, characteristic=2147483647):
	print("dummy: ", time[real]() - st);
	runtime := min(runtime, time[real]() - st);
end do;

timings_fn := "C:\data\projects\gbgb\Groebner.jl\benchmark\arxiv_preprint/results/maple/benchmark_1/dummy/timings":
FileTools[Text][WriteLine](timings_fn, "dummy");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime)));
