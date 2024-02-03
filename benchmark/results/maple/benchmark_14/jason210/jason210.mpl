# jason210
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		x1^2*x3^4 + x2^2*x4^4 + x1*x2*x3^2*x5^2 + x1*x2*x4^2*x6^2 + x1*x2*x3*x4*x5*x7 + x1*x2*x3*x4*x6*x8,
		x2^6,
		x1^6
	]:
	print("Running jason210");
	st := time[real]():
	G := Groebner[Basis](J, tdeg(x1, x2, x3, x4, x5, x6, x7, x8), method=fgb, characteristic=1073741827):
	print("jason210: ", time[real]() - st):
	runtime := min(runtime, time[real]() - st):
end do:

timings_fn := "/home/demin/Groebner.jl/benchmark/results/maple/benchmark_14/jason210/timings":
FileTools[Text][WriteLine](timings_fn, "jason210");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime))):

output_fn := "/home/demin/Groebner.jl/benchmark/results/maple/benchmark_14/jason210//output":
FileTools[Text][WriteLine](output_fn, "x1, x2, x3, x4, x5, x6, x7, x8");
FileTools[Text][WriteLine](output_fn, "1073741827");
for poly in G do
    FileTools[Text][WriteLine](output_fn, cat(String(poly), ",")):
end do:


