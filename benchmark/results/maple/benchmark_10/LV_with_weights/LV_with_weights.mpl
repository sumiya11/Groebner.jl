# LV_with_weights
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		-x1_0 + 1722549543,
		x2_0^2*x1_0*b_0 - x1_0*a_0 + x1_1,
		-x1_1 - 1610979036514140676391781201,
		x1_1*x2_0^2*b_0 + x2_1^2*x1_0*b_0 - x1_1*a_0 + x1_2,
		-x2_0^2*x1_0*d_0^3 + x2_0^2*c_0^3 + x2_1^2,
		-x1_2 - 759191294294390104261578745130487712005195693,
		2*x2_1^2*x1_1*b_0 + x1_2*x2_0^2*b_0 + x2_2^2*x1_0*b_0 - x1_2*a_0 + x1_3,
		-x1_1*x2_0^2*d_0^3 - x2_1^2*x1_0*d_0^3 + x2_1^2*c_0^3 + x2_2^2,
		-x1_3 + 3880357818382056242632385002196133299683348488516066676543044291,
		3*x1_2*x2_1^2*b_0 + 3*x2_2^2*x1_1*b_0 + x1_3*x2_0^2*b_0 + x2_3^2*x1_0*b_0 - x1_3*a_0 + x1_4,
		-2*x2_1^2*x1_1*d_0^3 - x1_2*x2_0^2*d_0^3 - x2_2^2*x1_0*d_0^3 + x2_2^2*c_0^3 + x2_3^2,
		-x1_4 + 7820461320675201682731134250881663006061140640874866070160665368275210770857578643,
		6*x2_2^2*x1_2*b_0 + 4*x1_3*x2_1^2*b_0 + 4*x2_3^2*x1_1*b_0 + x1_4*x2_0^2*b_0 + x2_4^2*x1_0*b_0 - x1_4*a_0 + x1_5,
		-3*x1_2*x2_1^2*d_0^3 - 3*x2_2^2*x1_1*d_0^3 - x1_3*x2_0^2*d_0^3 - x2_3^2*x1_0*d_0^3 + x2_3^2*c_0^3 + x2_4^2,
		-x1_5 - 34324326740641627828145440935729833286402435741777464375584133668421521614636195398430649294133784801,
		z_aux - 1
	]:
	print("Running LV_with_weights");
	st := time[real]():
	G := Groebner[Basis](J, tdeg(x1_5, x2_4, x1_4, x2_3, x1_3, x2_2, x1_2, x2_1, x1_1, x2_0, x1_0, z_aux, a_0, b_0, c_0, d_0), method=fgb, characteristic=0):
	print("LV_with_weights: ", time[real]() - st):
	runtime := min(runtime, time[real]() - st):
end do:

timings_fn := "/home/demin/Groebner.jl/benchmark/results/maple/benchmark_10/LV_with_weights/timings":
FileTools[Text][WriteLine](timings_fn, "LV_with_weights");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime))):

