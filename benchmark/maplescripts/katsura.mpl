with(Groebner):
with(PolynomialIdeals):

J := PolynomialIdeal({x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2 + 2*x8^2,
2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7 + 2*x7*x8, 2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5 + 2*x4*x6 + 2*x5*x7 + 2*x6*x8,
2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 + 2*x3*x6 - x3 + 2*x4*x7 + 2*x5*x8,
2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 + 2*x2*x6 + 2*x3*x7 + 2*x4*x8 - x4,
2*x0*x5 + 2*x1*x4 + 2*x1*x6 + 2*x2*x3 + 2*x2*x7 + 2*x3*x8 - x5,
2*x0*x6 + 2*x1*x5 + 2*x1*x7 + 2*x2*x4 + 2*x2*x8 + x3^2 - x6,
2*x0*x7 + 2*x1*x6 + 2*x1*x8 + 2*x2*x5 + 2*x3*x4 - x7,
x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 + 2*x8 - 1}):
st := time[real]():
Basis(J, tdeg(x0,x1,x2,x3,x4,x5,x6)):
time[real]() - st;
