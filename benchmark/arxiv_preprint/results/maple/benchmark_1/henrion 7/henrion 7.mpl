# henrion 7
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		2*f1*f2*f3*f4*f5*f6*f7 + 955883441,
		7*f1*f2*f3*f4*f5*f6 + 306783385*f1*f2*f3*f4*f5*f7 + 1227133519*f1*f2*f3*f4*f6*f7 + 613566762*f1*f2*f3*f5*f6*f7 + 613566761*f1*f2*f4*f5*f6*f7 + 1227133516*f1*f3*f4*f5*f6*f7 + 306783380*f2*f3*f4*f5*f6*f7 + 1018741245,
		6*f1*f2*f3*f4*f5 + 10*f1*f2*f3*f4*f6 + 613566762*f1*f2*f3*f4*f7 + 12*f1*f2*f3*f5*f6 + 1840700278*f1*f2*f3*f5*f7 + 1840700274*f1*f2*f3*f6*f7 + 12*f1*f2*f4*f5*f6 + 1533916901*f1*f2*f4*f5*f7 + 613566764*f1*f2*f4*f6*f7 + 1533916895*f1*f2*f5*f6*f7 + 10*f1*f3*f4*f5*f6 + 1840700278*f1*f3*f4*f5*f7 + 613566764*f1*f3*f4*f6*f7 + 613566762*f1*f3*f5*f6*f7 + 1840700272*f1*f4*f5*f6*f7 + 6*f2*f3*f4*f5*f6 + 613566762*f2*f3*f4*f5*f7 + 1840700274*f2*f3*f4*f6*f7 + 1533916895*f2*f3*f5*f6*f7 + 1840700272*f2*f4*f5*f6*f7 + 613566758*f3*f4*f5*f6*f7 + 1202512894,
		5*f1*f2*f3*f4 + 8*f1*f2*f3*f5 + 8*f1*f2*f3*f6 + 920350139*f1*f2*f3*f7 + 9*f1*f2*f4*f5 + 12*f1*f2*f4*f6 + 306783385*f1*f2*f4*f7 + 9*f1*f2*f5*f6 + 306783385*f1*f2*f5*f7 + 306783382*f1*f2*f6*f7 + 8*f1*f3*f4*f5 + 12*f1*f3*f4*f6 + 306783385*f1*f3*f4*f7 + 12*f1*f3*f5*f6 + 1840700278*f1*f3*f5*f7 + 1840700274*f1*f3*f6*f7 + 8*f1*f4*f5*f6 + 306783385*f1*f4*f5*f7 + 1840700274*f1*f4*f6*f7 + 306783381*f1*f5*f6*f7 + 5*f2*f3*f4*f5 + 8*f2*f3*f4*f6 + 920350139*f2*f3*f4*f7 + 9*f2*f3*f5*f6 + 306783385*f2*f3*f5*f7 + 306783382*f2*f3*f6*f7 + 8*f2*f4*f5*f6 + 306783385*f2*f4*f5*f7 + 1840700274*f2*f4*f6*f7 + 306783381*f2*f5*f6*f7 + 5*f3*f4*f5*f6 + 920350139*f3*f4*f5*f7 + 306783382*f3*f4*f6*f7 + 306783381*f3*f5*f6*f7 + 920350136*f4*f5*f6*f7 + 1561634524,
		4*f1*f2*f3 + 6*f1*f2*f4 + 6*f1*f2*f5 + 6*f1*f2*f6 + 1227133516*f1*f2*f7 + 6*f1*f3*f4 + 8*f1*f3*f5 + 8*f1*f3*f6 + 920350139*f1*f3*f7 + 6*f1*f4*f5 + 8*f1*f4*f6 + 920350139*f1*f4*f7 + 6*f1*f5*f6 + 920350139*f1*f5*f7 + 920350137*f1*f6*f7 + 4*f2*f3*f4 + 6*f2*f3*f5 + 6*f2*f3*f6 + 1227133516*f2*f3*f7 + 6*f2*f4*f5 + 8*f2*f4*f6 + 920350139*f2*f4*f7 + 6*f2*f5*f6 + 920350139*f2*f5*f7 + 920350137*f2*f6*f7 + 4*f3*f4*f5 + 6*f3*f4*f6 + 1227133516*f3*f4*f7 + 6*f3*f5*f6 + 920350139*f3*f5*f7 + 920350137*f3*f6*f7 + 4*f4*f5*f6 + 1227133516*f4*f5*f7 + 920350137*f4*f6*f7 + 1227133514*f5*f6*f7 + 2135808562,
		3*f1*f2 + 4*f1*f3 + 4*f1*f4 + 4*f1*f5 + 4*f1*f6 + 1533916893*f1*f7 + 3*f2*f3 + 4*f2*f4 + 4*f2*f5 + 4*f2*f6 + 1533916893*f2*f7 + 3*f3*f4 + 4*f3*f5 + 4*f3*f6 + 1533916893*f3*f7 + 3*f4*f5 + 4*f4*f6 + 1533916893*f4*f7 + 3*f5*f6 + 1533916893*f5*f7 + 1533916892*f6*f7 + 2147367594,
		2*f1 + 2*f2 + 2*f3 + 2*f4 + 2*f5 + 2*f6 + 1840700270*f7 + 2147483094
	]:
	print("Running henrion 7");
	st := time[real]():
	Groebner[Basis](J, tdeg(f1, f2, f3, f4, f5, f6, f7), method=fgb, characteristic=2147483647):
	print("henrion 7: ", time[real]() - st);
	runtime := min(runtime, time[real]() - st);
end do;

logs_fn := "C:\data\projects\gbgb\Groebner.jl\benchmark\arxiv_preprint/results/maple/benchmark_1/henrion 7/logs":
FileTools[Text][WriteLine](logs_fn, "henrion 7");
FileTools[Text][WriteLine](logs_fn, cat("total_time, ", String(runtime)));
