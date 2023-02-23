using Oscar

R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])

@edit lex([x, w])
