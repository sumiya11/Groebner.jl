# dummy 1
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		x^2 + y^2 + 6,
		x + y
	]:
	print("Running dummy 1");
	st := time[real]():
	G := Groebner[Basis](J, tdeg(x, y), method=fgb, characteristic=7):
	print("dummy 1: ", time[real]() - st):
	runtime := min(runtime, time[real]() - st):
end do:

timings_fn := "C:\data\projects\gbgb\Groebner.jl\benchmark/results/maple/benchmark_0/dummy 1/timings":
FileTools[Text][WriteLine](timings_fn, "dummy 1");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime))):

