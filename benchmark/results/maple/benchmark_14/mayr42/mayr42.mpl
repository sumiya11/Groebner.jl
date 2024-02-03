# mayr42
with(Groebner):
with(PolynomialIdeals):
kernelopts(numcpus=1);

runtime := 2^1000:
for i from 1 by 1 to 1 do
	J := [
		x4*x49 + 1073741826*x10*x51,
		x3*x48 + 1073741826*x9*x51,
		x2*x47 + 1073741826*x8*x51,
		x1*x46 + 1073741826*x7*x51,
		x4*x44 + 1073741826*x9*x49,
		x3*x43 + 1073741826*x8*x48,
		x2*x42 + 1073741826*x7*x47,
		x1*x41 + 1073741826*x6*x46,
		x4*x39 + 1073741826*x9*x49,
		x3*x38 + 1073741826*x8*x48,
		x2*x37 + 1073741826*x7*x47,
		x1*x36 + 1073741826*x6*x46,
		x9*x34 + 1073741826*x9*x49,
		x4*x34 + 1073741826*x5*x51,
		x8*x33 + 1073741826*x8*x48,
		x3*x33 + 1073741826*x4*x51,
		x7*x32 + 1073741826*x7*x47,
		x2*x32 + 1073741826*x3*x51,
		x6*x31 + 1073741826*x6*x46,
		x1*x31 + 1073741826*x2*x51,
		x9*x14*x39 + 1073741826*x9*x29*x44,
		x8*x13*x38 + 1073741826*x8*x28*x43,
		x7*x12*x37 + 1073741826*x7*x27*x42,
		x6*x11*x36 + 1073741826*x6*x26*x41,
		x6*x26^2*x46 + 1073741826*x7*x51^3,
		x6*x11^2*x46 + 1073741826*x2*x51^3,
		x6*x21^2*x41 + 1073741826*x6*x46*x51^2,
		x6*x16^2*x36 + 1073741826*x6*x46*x51^2,
		x9*x24*x30*x39*x50 + 1073741826*x9*x29*x44*x50*x51,
		x8*x23*x29*x38*x49 + 1073741826*x8*x28*x43*x49*x51,
		x7*x22*x28*x37*x48 + 1073741826*x7*x27*x42*x48*x51,
		x6*x21*x27*x36*x47 + 1073741826*x6*x26*x41*x47*x51,
		x9*x24*x25*x39*x45 + 1073741826*x9*x29*x44*x45*x51,
		x8*x23*x24*x38*x44 + 1073741826*x8*x28*x43*x44*x51,
		x7*x22*x23*x37*x43 + 1073741826*x7*x27*x42*x43*x51,
		x6*x21*x22*x36*x42 + 1073741826*x6*x26*x41*x42*x51,
		x9*x20*x24*x39*x40 + 1073741826*x9*x29*x40*x44*x51,
		x8*x19*x23*x38*x39 + 1073741826*x8*x28*x39*x43*x51,
		x9*x15*x24*x35*x39 + 1073741826*x9*x29*x35*x44*x51,
		x7*x18*x22*x37*x38 + 1073741826*x7*x27*x38*x42*x51,
		x8*x14*x23*x34*x38 + 1073741826*x8*x28*x34*x43*x51,
		x6*x17*x21*x36*x37 + 1073741826*x6*x26*x37*x41*x51,
		x7*x13*x22*x33*x37 + 1073741826*x7*x27*x33*x42*x51,
		x6*x12*x21*x32*x36 + 1073741826*x6*x26*x32*x41*x51
	]:
	print("Running mayr42");
	st := time[real]():
	G := Groebner[Basis](J, tdeg(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51), method=fgb, characteristic=1073741827):
	print("mayr42: ", time[real]() - st):
	runtime := min(runtime, time[real]() - st):
end do:

timings_fn := "/home/demin/Groebner.jl/benchmark/results/maple/benchmark_14/mayr42/timings":
FileTools[Text][WriteLine](timings_fn, "mayr42");
FileTools[Text][WriteLine](timings_fn, cat("total_time, ", String(runtime))):

output_fn := "/home/demin/Groebner.jl/benchmark/results/maple/benchmark_14/mayr42//output":
FileTools[Text][WriteLine](output_fn, "x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51");
FileTools[Text][WriteLine](output_fn, "1073741827");
for poly in G do
    FileTools[Text][WriteLine](output_fn, cat(String(poly), ",")):
end do:


