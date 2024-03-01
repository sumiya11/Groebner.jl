# dummy 2
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		x^2 + y^2 + 6,
		x + y
	]:
	print("Running dummy 2");
	st := time[real]():
	G := Groebner[Basis](J, tdeg(x, y), method=fgb, characteristic=7):
	print("dummy 2: ", time[real]() - st):
	runtime := min(runtime, time[real]() - st):
end do:

timings_fn := "C:\data\projects\gbgb\Groebner.jl\benchmark/results/maple/benchmark_0/dummy 2/timings":
FileTools[Text][WriteLine](timings_fn, "dummy 2");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime))):

