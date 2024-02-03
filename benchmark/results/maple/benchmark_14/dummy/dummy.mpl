# dummy
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		x^2 + y^2 + 1073741826,
		x + y
	]:
	print("Running dummy");
	st := time[real]():
	G := Groebner[Basis](J, tdeg(x, y), method=fgb, characteristic=1073741827):
	print("dummy: ", time[real]() - st):
	runtime := min(runtime, time[real]() - st):
end do:

timings_fn := "/home/demin/Groebner.jl/benchmark/results/maple/benchmark_14/dummy/timings":
FileTools[Text][WriteLine](timings_fn, "dummy");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime))):

output_fn := "/home/demin/Groebner.jl/benchmark/results/maple/benchmark_14/dummy//output":
FileTools[Text][WriteLine](output_fn, "x, y");
FileTools[Text][WriteLine](output_fn, "1073741827");
for poly in G do
    FileTools[Text][WriteLine](output_fn, cat(String(poly), ",")):
end do:


