with(Groebner):
with(PolynomialIdeals):


J := PolynomialIdeal({z1 + z2 + z3 + z4 + z5 + z6 + z7, z1*z2 + z1*z7 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z7, z1*z2*z3 + z1*z2*z7 + z1*z6*z7 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z7, z1*z2*z3*z4 + z1*z2*z3*z7 + z1*z2*z6*z7 + z1*z5*z6*z7 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z7, z1*z2*z3*z4*z5 + z1*z2*z3*z4*z7 + z1*z2*z3*z6*z7 + z1*z2*z5*z6*z7 + z1*z4*z5*z6*z7 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z7, z1*z2*z3*z4*z5*z6 + z1*z2*z3*z4*z5*z7 + z1*z2*z3*z4*z6*z7 + z1*z2*z3*z5*z6*z7 + z1*z2*z4*z5*z6*z7 + z1*z3*z4*z5*z6*z7 + z2*z3*z4*z5*z6*z7, z1*z2*z3*z4*z5*z6*z7 - 1}, charactesistic=0):
print("Running cyclic 7");
st := time[real]():
Groebner[Basis](J, tdeg(z1, z2, z3, z4, z5, z6, z7), method=direct):
print("cyclic 7: ", time[real]() - st);

J := PolynomialIdeal({z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8, z1*z2 + z1*z8 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z7 + z7*z8, z1*z2*z3 + z1*z2*z8 + z1*z7*z8 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z7 + z6*z7*z8, z1*z2*z3*z4 + z1*z2*z3*z8 + z1*z2*z7*z8 + z1*z6*z7*z8 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z7 + z5*z6*z7*z8, z1*z2*z3*z4*z5 + z1*z2*z3*z4*z8 + z1*z2*z3*z7*z8 + z1*z2*z6*z7*z8 + z1*z5*z6*z7*z8 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z7 + z4*z5*z6*z7*z8, z1*z2*z3*z4*z5*z6 + z1*z2*z3*z4*z5*z8 + z1*z2*z3*z4*z7*z8 + z1*z2*z3*z6*z7*z8 + z1*z2*z5*z6*z7*z8 + z1*z4*z5*z6*z7*z8 + z2*z3*z4*z5*z6*z7 + z3*z4*z5*z6*z7*z8, z1*z2*z3*z4*z5*z6*z7 + z1*z2*z3*z4*z5*z6*z8 + z1*z2*z3*z4*z5*z7*z8 + z1*z2*z3*z4*z6*z7*z8 + z1*z2*z3*z5*z6*z7*z8 + z1*z2*z4*z5*z6*z7*z8 + z1*z3*z4*z5*z6*z7*z8 + z2*z3*z4*z5*z6*z7*z8, z1*z2*z3*z4*z5*z6*z7*z8 - 1}, charactesistic=0):
print("Running cyclic 8");
st := time[real]():
Groebner[Basis](J, tdeg(z1, z2, z3, z4, z5, z6, z7, z8), method=direct):
print("cyclic 8: ", time[real]() - st);

J := PolynomialIdeal({x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2 + 2*x8^2 + 2*x9^2, 2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7 + 2*x7*x8 + 2*x8*x9, 2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5 + 2*x4*x6 + 2*x5*x7 + 2*x6*x8 + 2*x7*x9, 2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 + 2*x3*x6 - x3 + 2*x4*x7 + 2*x5*x8 + 2*x6*x9, 2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 + 2*x2*x6 + 2*x3*x7 + 2*x4*x8 - x4 + 2*x5*x9, 2*x0*x5 + 2*x1*x4 + 2*x1*x6 + 2*x2*x3 + 2*x2*x7 + 2*x3*x8 + 2*x4*x9 - x5, 2*x0*x6 + 2*x1*x5 + 2*x1*x7 + 2*x2*x4 + 2*x2*x8 + x3^2 + 2*x3*x9 - x6, 2*x0*x7 + 2*x1*x6 + 2*x1*x8 + 2*x2*x5 + 2*x2*x9 + 2*x3*x4 - x7, 2*x0*x8 + 2*x1*x7 + 2*x1*x9 + 2*x2*x6 + 2*x3*x5 + x4^2 - x8, x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 + 2*x8 + 2*x9 - 1}, charactesistic=0):
print("Running katsura 9");
st := time[real]():
Groebner[Basis](J, tdeg(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9), method=direct):
print("katsura 9: ", time[real]() - st);

J := PolynomialIdeal({x0^2 - x0 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2 + 2*x8^2 + 2*x9^2 + 2*x10^2, 2*x0*x1 + 2*x1*x2 - x1 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7 + 2*x7*x8 + 2*x8*x9 + 2*x9*x10, 2*x0*x2 + x1^2 + 2*x1*x3 + 2*x2*x4 - x2 + 2*x3*x5 + 2*x4*x6 + 2*x5*x7 + 2*x6*x8 + 2*x7*x9 + 2*x8*x10, 2*x0*x3 + 2*x1*x2 + 2*x1*x4 + 2*x2*x5 + 2*x3*x6 - x3 + 2*x4*x7 + 2*x5*x8 + 2*x6*x9 + 2*x7*x10, 2*x0*x4 + 2*x1*x3 + 2*x1*x5 + x2^2 + 2*x2*x6 + 2*x3*x7 + 2*x4*x8 - x4 + 2*x5*x9 + 2*x6*x10, 2*x0*x5 + 2*x1*x4 + 2*x1*x6 + 2*x2*x3 + 2*x2*x7 + 2*x3*x8 + 2*x4*x9 + 2*x5*x10 - x5, 2*x0*x6 + 2*x1*x5 + 2*x1*x7 + 2*x2*x4 + 2*x2*x8 + x3^2 + 2*x3*x9 + 2*x4*x10 - x6, 2*x0*x7 + 2*x1*x6 + 2*x1*x8 + 2*x2*x5 + 2*x2*x9 + 2*x3*x4 + 2*x3*x10 - x7, 2*x0*x8 + 2*x1*x7 + 2*x1*x9 + 2*x2*x6 + 2*x2*x10 + 2*x3*x5 + x4^2 - x8, 2*x0*x9 + 2*x1*x8 + 2*x1*x10 + 2*x2*x7 + 2*x3*x6 + 2*x4*x5 - x9, x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 + 2*x8 + 2*x9 + 2*x10 - 1}, charactesistic=0):
print("Running katsura 10");
st := time[real]():
Groebner[Basis](J, tdeg(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10), method=direct):
print("katsura 10: ", time[real]() - st);

J := PolynomialIdeal({x1*x2*x10 + x1*x10 + x2*x3*x10 + x3*x4*x10 + x4*x5*x10 + x5*x6*x10 + x6*x7*x10 + x7*x8*x10 + x8*x9*x10 - 1, x1*x3*x10 + x2*x4*x10 + x2*x10 + x3*x5*x10 + x4*x6*x10 + x5*x7*x10 + x6*x8*x10 + x7*x9*x10 - 2, x1*x4*x10 + x2*x5*x10 + x3*x6*x10 + x3*x10 + x4*x7*x10 + x5*x8*x10 + x6*x9*x10 - 3, x1*x5*x10 + x2*x6*x10 + x3*x7*x10 + x4*x8*x10 + x4*x10 + x5*x9*x10 - 4, x1*x6*x10 + x2*x7*x10 + x3*x8*x10 + x4*x9*x10 + x5*x10 - 5, x1*x7*x10 + x2*x8*x10 + x3*x9*x10 + x6*x10 - 6, x1*x8*x10 + x2*x9*x10 + x7*x10 - 7, x1*x9*x10 + x8*x10 - 8, x9*x10 - 9, x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 1}, charactesistic=0):
print("Running eco 10");
st := time[real]():
Groebner[Basis](J, tdeg(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11), method=direct):
print("eco 10: ", time[real]() - st);

J := PolynomialIdeal({x1*x2*x11 + x1*x11 + x2*x3*x11 + x3*x4*x11 + x4*x5*x11 + x5*x6*x11 + x6*x7*x11 + x7*x8*x11 + x8*x9*x11 + x9*x10*x11 - 1, x1*x3*x11 + x2*x4*x11 + x2*x11 + x3*x5*x11 + x4*x6*x11 + x5*x7*x11 + x6*x8*x11 + x7*x9*x11 + x8*x10*x11 - 2, x1*x4*x11 + x2*x5*x11 + x3*x6*x11 + x3*x11 + x4*x7*x11 + x5*x8*x11 + x6*x9*x11 + x7*x10*x11 - 3, x1*x5*x11 + x2*x6*x11 + x3*x7*x11 + x4*x8*x11 + x4*x11 + x5*x9*x11 + x6*x10*x11 - 4, x1*x6*x11 + x2*x7*x11 + x3*x8*x11 + x4*x9*x11 + x5*x10*x11 + x5*x11 - 5, x1*x7*x11 + x2*x8*x11 + x3*x9*x11 + x4*x10*x11 + x6*x11 - 6, x1*x8*x11 + x2*x9*x11 + x3*x10*x11 + x7*x11 - 7, x1*x9*x11 + x2*x10*x11 + x8*x11 - 8, x1*x10*x11 + x9*x11 - 9, x10*x11 - 10, x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 1}, charactesistic=0):
print("Running eco 11");
st := time[real]():
Groebner[Basis](J, tdeg(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11), method=direct):
print("eco 11: ", time[real]() - st);

J := PolynomialIdeal({10*x1*x2^2 + 10*x1*x3^2 + 10*x1*x4^2 + 10*x1*x5^2 + 10*x1*x6^2 + 10*x1*x7^2 + 10*x1*x8^2 - 11*x1 + 10, 10*x1^2*x2 + 10*x2*x3^2 + 10*x2*x4^2 + 10*x2*x5^2 + 10*x2*x6^2 + 10*x2*x7^2 + 10*x2*x8^2 - 11*x2 + 10, 10*x1^2*x3 + 10*x2^2*x3 + 10*x3*x4^2 + 10*x3*x5^2 + 10*x3*x6^2 + 10*x3*x7^2 + 10*x3*x8^2 - 11*x3 + 10, 10*x1^2*x4 + 10*x2^2*x4 + 10*x3^2*x4 + 10*x4*x5^2 + 10*x4*x6^2 + 10*x4*x7^2 + 10*x4*x8^2 - 11*x4 + 10, 10*x1^2*x5 + 10*x2^2*x5 + 10*x3^2*x5 + 10*x4^2*x5 + 10*x5*x6^2 + 10*x5*x7^2 + 10*x5*x8^2 - 11*x5 + 10, 10*x1^2*x6 + 10*x2^2*x6 + 10*x3^2*x6 + 10*x4^2*x6 + 10*x5^2*x6 + 10*x6*x7^2 + 10*x6*x8^2 - 11*x6 + 10, 10*x1^2*x7 + 10*x2^2*x7 + 10*x3^2*x7 + 10*x4^2*x7 + 10*x5^2*x7 + 10*x6^2*x7 + 10*x7*x8^2 - 11*x7 + 10, 10*x1^2*x8 + 10*x2^2*x8 + 10*x3^2*x8 + 10*x4^2*x8 + 10*x5^2*x8 + 10*x6^2*x8 + 10*x7^2*x8 - 11*x8 + 10}, charactesistic=0):
print("Running noon 8");
st := time[real]():
Groebner[Basis](J, tdeg(x1, x2, x3, x4, x5, x6, x7, x8), method=direct):
print("noon 8: ", time[real]() - st);

J := PolynomialIdeal({10*x1*x2^2 + 10*x1*x3^2 + 10*x1*x4^2 + 10*x1*x5^2 + 10*x1*x6^2 + 10*x1*x7^2 + 10*x1*x8^2 + 10*x1*x9^2 - 11*x1 + 10, 10*x1^2*x2 + 10*x2*x3^2 + 10*x2*x4^2 + 10*x2*x5^2 + 10*x2*x6^2 + 10*x2*x7^2 + 10*x2*x8^2 + 10*x2*x9^2 - 11*x2 + 10, 10*x1^2*x3 + 10*x2^2*x3 + 10*x3*x4^2 + 10*x3*x5^2 + 10*x3*x6^2 + 10*x3*x7^2 + 10*x3*x8^2 + 10*x3*x9^2 - 11*x3 + 10, 10*x1^2*x4 + 10*x2^2*x4 + 10*x3^2*x4 + 10*x4*x5^2 + 10*x4*x6^2 + 10*x4*x7^2 + 10*x4*x8^2 + 10*x4*x9^2 - 11*x4 + 10, 10*x1^2*x5 + 10*x2^2*x5 + 10*x3^2*x5 + 10*x4^2*x5 + 10*x5*x6^2 + 10*x5*x7^2 + 10*x5*x8^2 + 10*x5*x9^2 - 11*x5 + 10, 10*x1^2*x6 + 10*x2^2*x6 + 10*x3^2*x6 + 10*x4^2*x6 + 10*x5^2*x6 + 10*x6*x7^2 + 10*x6*x8^2 + 10*x6*x9^2 - 11*x6 + 10, 10*x1^2*x7 + 10*x2^2*x7 + 10*x3^2*x7 + 10*x4^2*x7 + 10*x5^2*x7 + 10*x6^2*x7 + 10*x7*x8^2 + 10*x7*x9^2 - 11*x7 + 10, 10*x1^2*x8 + 10*x2^2*x8 + 10*x3^2*x8 + 10*x4^2*x8 + 10*x5^2*x8 + 10*x6^2*x8 + 10*x7^2*x8 + 10*x8*x9^2 - 11*x8 + 10, 10*x1^2*x9 + 10*x2^2*x9 + 10*x3^2*x9 + 10*x4^2*x9 + 10*x5^2*x9 + 10*x6^2*x9 + 10*x7^2*x9 + 10*x8^2*x9 - 11*x9 + 10}, charactesistic=0):
print("Running noon 9");
st := time[real]():
Groebner[Basis](J, tdeg(x1, x2, x3, x4, x5, x6, x7, x8, x9), method=direct):
print("noon 9: ", time[real]() - st);

J := PolynomialIdeal({2*f1*f2*f3*f4*f5*f6 - 1404728325, 6*f1*f2*f3*f4*f5 + 35//6*f1*f2*f3*f4*f6 + 16//3*f1*f2*f3*f5*f6 + 9//2*f1*f2*f4*f5*f6 + 10//3*f1*f3*f4*f5*f6 + 11//6*f2*f3*f4*f5*f6 - 648336150, 5*f1*f2*f3*f4 + 8*f1*f2*f3*f5 + 14//3*f1*f2*f3*f6 + 9*f1*f2*f4*f5 + 7*f1*f2*f4*f6 + 4*f1*f2*f5*f6 + 8*f1*f3*f4*f5 + 7*f1*f3*f4*f6 + 16//3*f1*f3*f5*f6 + 3*f1*f4*f5*f6 + 5*f2*f3*f4*f5 + 14//3*f2*f3*f4*f6 + 4*f2*f3*f5*f6 + 3*f2*f4*f5*f6 + 5//3*f3*f4*f5*f6 - 67597623, 4*f1*f2*f3 + 6*f1*f2*f4 + 6*f1*f2*f5 + 7//2*f1*f2*f6 + 6*f1*f3*f4 + 8*f1*f3*f5 + 14//3*f1*f3*f6 + 6*f1*f4*f5 + 14//3*f1*f4*f6 + 8//3*f1*f5*f6 + 4*f2*f3*f4 + 6*f2*f3*f5 + 7//2*f2*f3*f6 + 6*f2*f4*f5 + 14//3*f2*f4*f6 + 8//3*f2*f5*f6 + 4*f3*f4*f5 + 7//2*f3*f4*f6 + 8//3*f3*f5*f6 + 3//2*f4*f5*f6 - 2657700, 3*f1*f2 + 4*f1*f3 + 4*f1*f4 + 4*f1*f5 + 7//3*f1*f6 + 3*f2*f3 + 4*f2*f4 + 4*f2*f5 + 7//3*f2*f6 + 3*f3*f4 + 4*f3*f5 + 7//3*f3*f6 + 3*f4*f5 + 7//3*f4*f6 + 4//3*f5*f6 - 46243, 2*f1 + 2*f2 + 2*f3 + 2*f4 + 2*f5 + 7//6*f6 - 358}, charactesistic=0):
print("Running henrion 6");
st := time[real]():
Groebner[Basis](J, tdeg(f1, f2, f3, f4, f5, f6), method=direct):
print("henrion 6: ", time[real]() - st);

