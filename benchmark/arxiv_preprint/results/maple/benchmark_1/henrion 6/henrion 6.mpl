# henrion 6
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		2*f1*f2*f3*f4*f5*f6 + 742755322,
		6*f1*f2*f3*f4*f5 + 357913947*f1*f2*f3*f4*f6 + 1431655770*f1*f2*f3*f5*f6 + 1073741828*f1*f2*f4*f5*f6 + 1431655768*f1*f3*f4*f5*f6 + 357913943*f2*f3*f4*f5*f6 + 1499147497,
		5*f1*f2*f3*f4 + 8*f1*f2*f3*f5 + 715827887*f1*f2*f3*f6 + 9*f1*f2*f4*f5 + 7*f1*f2*f4*f6 + 4*f1*f2*f5*f6 + 8*f1*f3*f4*f5 + 7*f1*f3*f4*f6 + 1431655770*f1*f3*f5*f6 + 3*f1*f4*f5*f6 + 5*f2*f3*f4*f5 + 715827887*f2*f3*f4*f6 + 4*f2*f3*f5*f6 + 3*f2*f4*f5*f6 + 715827884*f3*f4*f5*f6 + 2079886024,
		4*f1*f2*f3 + 6*f1*f2*f4 + 6*f1*f2*f5 + 1073741827*f1*f2*f6 + 6*f1*f3*f4 + 8*f1*f3*f5 + 715827887*f1*f3*f6 + 6*f1*f4*f5 + 715827887*f1*f4*f6 + 715827885*f1*f5*f6 + 4*f2*f3*f4 + 6*f2*f3*f5 + 1073741827*f2*f3*f6 + 6*f2*f4*f5 + 715827887*f2*f4*f6 + 715827885*f2*f5*f6 + 4*f3*f4*f5 + 1073741827*f3*f4*f6 + 715827885*f3*f5*f6 + 1073741825*f4*f5*f6 + 2144825947,
		3*f1*f2 + 4*f1*f3 + 4*f1*f4 + 4*f1*f5 + 1431655767*f1*f6 + 3*f2*f3 + 4*f2*f4 + 4*f2*f5 + 1431655767*f2*f6 + 3*f3*f4 + 4*f3*f5 + 1431655767*f3*f6 + 3*f4*f5 + 1431655767*f4*f6 + 1431655766*f5*f6 + 2147437404,
		2*f1 + 2*f2 + 2*f3 + 2*f4 + 2*f5 + 1789569707*f6 + 2147483289
	]:
	print("Running henrion 6");
	st := time[real]():
	Groebner[Basis](J, tdeg(f1, f2, f3, f4, f5, f6), method=fgb, characteristic=2147483647):
	print("henrion 6: ", time[real]() - st);
	runtime := min(runtime, time[real]() - st);
end do;

logs_fn := "C:\data\projects\gbgb\Groebner.jl\benchmark\arxiv_preprint/results/maple/benchmark_1/henrion 6/logs":
FileTools[Text][WriteLine](logs_fn, "henrion 6");
FileTools[Text][WriteLine](logs_fn, cat("total_time, ", String(runtime)));
