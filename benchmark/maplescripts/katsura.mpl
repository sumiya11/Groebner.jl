with(Groebner):
with(PolynomialIdeals):

J := PolynomialIdeal({x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 + 2*x8 + 2*x9 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2 + 2*x8^2 + 2*x9^2, 2*x1*x2 + 2*x2*x3 - x2 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7 + 2*x7*x8 + 2*x8*x9, 2*x1*x3 + x2^2 + 2*x2*x4 + 2*x3*x5 - x3 + 2*x4*x6 + 2*x5*x7 + 2*x6*x8 + 2*x7*x9, 2*x1*x4 + 2*x2*x3 + 2*x2*x5 + 2*x3*x6 + 2*x4*x7 - x4 + 2*x5*x8 + 2*x6*x9, 2*x1*x5 + 2*x2*x4 + 2*x2*x6 + x3^2 + 2*x3*x7 + 2*x4*x8 + 2*x5*x9 - x5, 2*x1*x6 + 2*x2*x5 + 2*x2*x7 + 2*x3*x4 + 2*x3*x8 + 2*x4*x9 - x6, 2*x1*x7 + 2*x2*x6 + 2*x2*x8 + 2*x3*x5 + 2*x3*x9 + x4^2 - x7, 2*x1*x8 + 2*x2*x7 + 2*x2*x9 + 2*x3*x6 + 2*x4*x5 - x8}):
st := time[real]():
Basis(J, tdeg(x1,x2,x3,x4,x5,x6,x7,x8,x9), charactesistic=2^31-1):
time[real]() - st;

J := PolynomialIdeal({x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 + 2*x8 + 2*x9 + 2*x10 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2 + 2*x8^2 + 2*x9^2 + 2*x10^2, 2*x1*x2 + 2*x2*x3 - x2 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7 + 2*x7*x8 + 2*x8*x9 + 2*x9*x10, 2*x1*x3 + x2^2 + 2*x2*x4 + 2*x3*x5 - x3 + 2*x4*x6 + 2*x5*x7 + 2*x6*x8 + 2*x7*x9 + 2*x8*x10, 2*x1*x4 + 2*x2*x3 + 2*x2*x5 + 2*x3*x6 + 2*x4*x7 - x4 + 2*x5*x8 + 2*x6*x9 + 2*x7*x10, 2*x1*x5 + 2*x2*x4 + 2*x2*x6 + x3^2 + 2*x3*x7 + 2*x4*x8 + 2*x5*x9 - x5 + 2*x6*x10, 2*x1*x6 + 2*x2*x5 + 2*x2*x7 + 2*x3*x4 + 2*x3*x8 + 2*x4*x9 + 2*x5*x10 - x6, 2*x1*x7 + 2*x2*x6 + 2*x2*x8 + 2*x3*x5 + 2*x3*x9 + x4^2 + 2*x4*x10 - x7, 2*x1*x8 + 2*x2*x7 + 2*x2*x9 + 2*x3*x6 + 2*x3*x10 + 2*x4*x5 - x8, 2*x1*x9 + 2*x2*x8 + 2*x2*x10 + 2*x3*x7 + 2*x4*x6 + x5^2 - x9}):
st := time[real]():
Basis(J, tdeg(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10), characteristic=2^31-1):
time[real]() - st;        
