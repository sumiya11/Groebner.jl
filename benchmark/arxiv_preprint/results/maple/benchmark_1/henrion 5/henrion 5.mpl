# henrion 5
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		2*f1*f2*f3*f4*f5 + 2137660372,
		5*f1*f2*f3*f4 + 1288490193*f1*f2*f3*f5 + 858993463*f1*f2*f4*f5 + 858993462*f1*f3*f4*f5 + 1288490190*f2*f3*f4*f5 + 2143018522,
		4*f1*f2*f3 + 6*f1*f2*f4 + 429496733*f1*f2*f5 + 6*f1*f3*f4 + 1288490193*f1*f3*f5 + 1288490191*f1*f4*f5 + 4*f2*f3*f4 + 429496733*f2*f3*f5 + 1288490191*f2*f4*f5 + 429496731*f3*f4*f5 + 2147042161,
		3*f1*f2 + 4*f1*f3 + 4*f1*f4 + 1717986920*f1*f5 + 3*f2*f3 + 4*f2*f4 + 1717986920*f2*f5 + 3*f3*f4 + 1717986920*f3*f5 + 1717986919*f4*f5 + 2147468149,
		2*f1 + 2*f2 + 2*f3 + 2*f4 + 858993460*f5 + 2147483432,
		f1 + 2*f2 + 3*f3 + 4*f4 + 5*f5 + 6*t
	]:
	print("Running henrion 5");
	st := time[real]():
	Groebner[Basis](J, tdeg(f1, f2, f3, f4, f5, t), method=fgb, characteristic=2147483647):
	print("henrion 5: ", time[real]() - st);
	runtime := min(runtime, time[real]() - st);
end do;

timings_fn := "/home/ademin/Groebner.jl/benchmark/arxiv_preprint/results/maple/benchmark_1/henrion 5/timings":
FileTools[Text][WriteLine](timings_fn, "henrion 5");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime)));
