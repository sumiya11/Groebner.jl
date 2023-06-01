using AbstractAlgebra
R1,
(
    a,
    b,
    c,
    d,
    e,
    f,
    g,
    h,
    i,
    j,
    k,
    l,
    m,
    n,
    o,
    p,
    q,
    r,
    s,
    t,
    u,
    v,
    w,
    x,
    y,
    z,
    A,
    B,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    J,
    K,
    L,
    M,
    N,
    O,
    P,
    Q,
    R,
    S,
    T,
    U,
    V
) = PolynomialRing(
    GF(101),
    [
        "a",
        "b",
        "c",
        "d",
        "e",
        "f",
        "g",
        "h",
        "i",
        "j",
        "k",
        "l",
        "m",
        "n",
        "o",
        "p",
        "q",
        "r",
        "s",
        "t",
        "u",
        "v",
        "w",
        "x",
        "y",
        "z",
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "J",
        "K",
        "L",
        "M",
        "N",
        "O",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "U",
        "V"
    ],
    ordering=:degrevlex
);
polys = [
    d * g * j * m - c * h * j * m - d * f * k * m + b * h * k * m + c * f * l * m -
    b * g * l * m - d * g * i * n +
    c * h * i * n +
    d * e * k * n - a * h * k * n - c * e * l * n +
    a * g * l * n +
    d * f * i * o - b * h * i * o - d * e * j * o +
    a * h * j * o +
    b * e * l * o - a * f * l * o - c * f * i * p +
    b * g * i * p +
    c * e * j * p - a * g * j * p - b * e * k * p + a * f * k * p,
    d * g * j * q - c * h * j * q - d * f * k * q + b * h * k * q + c * f * l * q -
    b * g * l * q - d * g * i * r +
    c * h * i * r +
    d * e * k * r - a * h * k * r - c * e * l * r +
    a * g * l * r +
    d * f * i * s - b * h * i * s - d * e * j * s +
    a * h * j * s +
    b * e * l * s - a * f * l * s - c * f * i * t +
    b * g * i * t +
    c * e * j * t - a * g * j * t - b * e * k * t + a * f * k * t,
    d * g * n * q - c * h * n * q - d * f * o * q + b * h * o * q + c * f * p * q -
    b * g * p * q - d * g * m * r +
    c * h * m * r +
    d * e * o * r - a * h * o * r - c * e * p * r +
    a * g * p * r +
    d * f * m * s - b * h * m * s - d * e * n * s +
    a * h * n * s +
    b * e * p * s - a * f * p * s - c * f * m * t +
    b * g * m * t +
    c * e * n * t - a * g * n * t - b * e * o * t + a * f * o * t,
    d * k * n * q - c * l * n * q - d * j * o * q + b * l * o * q + c * j * p * q -
    b * k * p * q - d * k * m * r +
    c * l * m * r +
    d * i * o * r - a * l * o * r - c * i * p * r +
    a * k * p * r +
    d * j * m * s - b * l * m * s - d * i * n * s +
    a * l * n * s +
    b * i * p * s - a * j * p * s - c * j * m * t +
    b * k * m * t +
    c * i * n * t - a * k * n * t - b * i * o * t + a * j * o * t,
    h * k * n * q - g * l * n * q - h * j * o * q + f * l * o * q + g * j * p * q -
    f * k * p * q - h * k * m * r +
    g * l * m * r +
    h * i * o * r - e * l * o * r - g * i * p * r +
    e * k * p * r +
    h * j * m * s - f * l * m * s - h * i * n * s +
    e * l * n * s +
    f * i * p * s - e * j * p * s - g * j * m * t +
    f * k * m * t +
    g * i * n * t - e * k * n * t - f * i * o * t + e * j * o * t,
    d * g * j * u - c * h * j * u - d * f * k * u + b * h * k * u + c * f * l * u -
    b * g * l * u - d * g * i * v +
    c * h * i * v +
    d * e * k * v - a * h * k * v - c * e * l * v +
    a * g * l * v +
    d * f * i * w - b * h * i * w - d * e * j * w +
    a * h * j * w +
    b * e * l * w - a * f * l * w - c * f * i * x +
    b * g * i * x +
    c * e * j * x - a * g * j * x - b * e * k * x + a * f * k * x,
    d * g * n * u - c * h * n * u - d * f * o * u + b * h * o * u + c * f * p * u -
    b * g * p * u - d * g * m * v +
    c * h * m * v +
    d * e * o * v - a * h * o * v - c * e * p * v +
    a * g * p * v +
    d * f * m * w - b * h * m * w - d * e * n * w +
    a * h * n * w +
    b * e * p * w - a * f * p * w - c * f * m * x +
    b * g * m * x +
    c * e * n * x - a * g * n * x - b * e * o * x + a * f * o * x,
    d * k * n * u - c * l * n * u - d * j * o * u + b * l * o * u + c * j * p * u -
    b * k * p * u - d * k * m * v +
    c * l * m * v +
    d * i * o * v - a * l * o * v - c * i * p * v +
    a * k * p * v +
    d * j * m * w - b * l * m * w - d * i * n * w +
    a * l * n * w +
    b * i * p * w - a * j * p * w - c * j * m * x +
    b * k * m * x +
    c * i * n * x - a * k * n * x - b * i * o * x + a * j * o * x,
    h * k * n * u - g * l * n * u - h * j * o * u + f * l * o * u + g * j * p * u -
    f * k * p * u - h * k * m * v +
    g * l * m * v +
    h * i * o * v - e * l * o * v - g * i * p * v +
    e * k * p * v +
    h * j * m * w - f * l * m * w - h * i * n * w +
    e * l * n * w +
    f * i * p * w - e * j * p * w - g * j * m * x +
    f * k * m * x +
    g * i * n * x - e * k * n * x - f * i * o * x + e * j * o * x,
    d * g * r * u - c * h * r * u - d * f * s * u + b * h * s * u + c * f * t * u -
    b * g * t * u - d * g * q * v +
    c * h * q * v +
    d * e * s * v - a * h * s * v - c * e * t * v +
    a * g * t * v +
    d * f * q * w - b * h * q * w - d * e * r * w +
    a * h * r * w +
    b * e * t * w - a * f * t * w - c * f * q * x +
    b * g * q * x +
    c * e * r * x - a * g * r * x - b * e * s * x + a * f * s * x,
    d * k * r * u - c * l * r * u - d * j * s * u + b * l * s * u + c * j * t * u -
    b * k * t * u - d * k * q * v +
    c * l * q * v +
    d * i * s * v - a * l * s * v - c * i * t * v +
    a * k * t * v +
    d * j * q * w - b * l * q * w - d * i * r * w +
    a * l * r * w +
    b * i * t * w - a * j * t * w - c * j * q * x +
    b * k * q * x +
    c * i * r * x - a * k * r * x - b * i * s * x + a * j * s * x,
    h * k * r * u - g * l * r * u - h * j * s * u + f * l * s * u + g * j * t * u -
    f * k * t * u - h * k * q * v +
    g * l * q * v +
    h * i * s * v - e * l * s * v - g * i * t * v +
    e * k * t * v +
    h * j * q * w - f * l * q * w - h * i * r * w +
    e * l * r * w +
    f * i * t * w - e * j * t * w - g * j * q * x +
    f * k * q * x +
    g * i * r * x - e * k * r * x - f * i * s * x + e * j * s * x,
    d * o * r * u - c * p * r * u - d * n * s * u + b * p * s * u + c * n * t * u -
    b * o * t * u - d * o * q * v +
    c * p * q * v +
    d * m * s * v - a * p * s * v - c * m * t * v +
    a * o * t * v +
    d * n * q * w - b * p * q * w - d * m * r * w +
    a * p * r * w +
    b * m * t * w - a * n * t * w - c * n * q * x +
    b * o * q * x +
    c * m * r * x - a * o * r * x - b * m * s * x + a * n * s * x,
    h * o * r * u - g * p * r * u - h * n * s * u + f * p * s * u + g * n * t * u -
    f * o * t * u - h * o * q * v +
    g * p * q * v +
    h * m * s * v - e * p * s * v - g * m * t * v +
    e * o * t * v +
    h * n * q * w - f * p * q * w - h * m * r * w +
    e * p * r * w +
    f * m * t * w - e * n * t * w - g * n * q * x +
    f * o * q * x +
    g * m * r * x - e * o * r * x - f * m * s * x + e * n * s * x,
    l * o * r * u - k * p * r * u - l * n * s * u + j * p * s * u + k * n * t * u -
    j * o * t * u - l * o * q * v +
    k * p * q * v +
    l * m * s * v - i * p * s * v - k * m * t * v +
    i * o * t * v +
    l * n * q * w - j * p * q * w - l * m * r * w +
    i * p * r * w +
    j * m * t * w - i * n * t * w - k * n * q * x +
    j * o * q * x +
    k * m * r * x - i * o * r * x - j * m * s * x + i * n * s * x,
    a * y + b * z + c * A + d * B,
    e * y + f * z + g * A + h * B,
    i * y + j * z + k * A + l * B,
    m * y + n * z + o * A + p * B,
    q * y + r * z + s * A + t * B,
    u * y + v * z + w * A + x * B,
    a * C + b * D + c * E + d * F,
    e * C + f * D + g * E + h * F,
    i * C + j * D + k * E + l * F,
    m * C + n * D + o * E + p * F,
    q * C + r * D + s * E + t * F,
    u * C + v * D + w * E + x * F,
    a * G + b * H + c * I + d * J,
    e * G + f * H + g * I + h * J,
    i * G + j * H + k * I + l * J,
    m * G + n * H + o * I + p * J,
    q * G + r * H + s * I + t * J,
    u * G + v * H + w * I + x * J,
    a * K + b * L + c * M + d * N,
    e * K + f * L + g * M + h * N,
    i * K + j * L + k * M + l * N,
    m * K + n * L + o * M + p * N,
    q * K + r * L + s * M + t * N,
    u * K + v * L + w * M + x * N,
    B * E * H * K - A * F * H * K - B * D * I * K + z * F * I * K + A * D * J * K -
    z * E * J * K - B * E * G * L +
    A * F * G * L +
    B * C * I * L - y * F * I * L - A * C * J * L +
    y * E * J * L +
    B * D * G * M - z * F * G * M - B * C * H * M +
    y * F * H * M +
    z * C * J * M - y * D * J * M - A * D * G * N +
    z * E * G * N +
    A * C * H * N - y * E * H * N - z * C * I * N + y * D * I * N,
    a * O + b * P + c * Q + d * R,
    e * O + f * P + g * Q + h * R,
    i * O + j * P + k * Q + l * R,
    m * O + n * P + o * Q + p * R,
    q * O + r * P + s * Q + t * R,
    u * O + v * P + w * Q + x * R,
    B * E * H * O - A * F * H * O - B * D * I * O + z * F * I * O + A * D * J * O -
    z * E * J * O - B * E * G * P +
    A * F * G * P +
    B * C * I * P - y * F * I * P - A * C * J * P +
    y * E * J * P +
    B * D * G * Q - z * F * G * Q - B * C * H * Q +
    y * F * H * Q +
    z * C * J * Q - y * D * J * Q - A * D * G * R +
    z * E * G * R +
    A * C * H * R - y * E * H * R - z * C * I * R + y * D * I * R,
    B * E * L * O - A * F * L * O - B * D * M * O + z * F * M * O + A * D * N * O -
    z * E * N * O - B * E * K * P +
    A * F * K * P +
    B * C * M * P - y * F * M * P - A * C * N * P +
    y * E * N * P +
    B * D * K * Q - z * F * K * Q - B * C * L * Q +
    y * F * L * Q +
    z * C * N * Q - y * D * N * Q - A * D * K * R +
    z * E * K * R +
    A * C * L * R - y * E * L * R - z * C * M * R + y * D * M * R,
    B * I * L * O - A * J * L * O - B * H * M * O + z * J * M * O + A * H * N * O -
    z * I * N * O - B * I * K * P +
    A * J * K * P +
    B * G * M * P - y * J * M * P - A * G * N * P +
    y * I * N * P +
    B * H * K * Q - z * J * K * Q - B * G * L * Q +
    y * J * L * Q +
    z * G * N * Q - y * H * N * Q - A * H * K * R +
    z * I * K * R +
    A * G * L * R - y * I * L * R - z * G * M * R + y * H * M * R,
    F * I * L * O - E * J * L * O - F * H * M * O + D * J * M * O + E * H * N * O -
    D * I * N * O - F * I * K * P +
    E * J * K * P +
    F * G * M * P - C * J * M * P - E * G * N * P +
    C * I * N * P +
    F * H * K * Q - D * J * K * Q - F * G * L * Q +
    C * J * L * Q +
    D * G * N * Q - C * H * N * Q - E * H * K * R +
    D * I * K * R +
    E * G * L * R - C * I * L * R - D * G * M * R + C * H * M * R,
    a * S + b * T + c * U + d * V,
    e * S + f * T + g * U + h * V,
    i * S + j * T + k * U + l * V,
    m * S + n * T + o * U + p * V,
    q * S + r * T + s * U + t * V,
    u * S + v * T + w * U + x * V,
    B * E * H * S - A * F * H * S - B * D * I * S + z * F * I * S + A * D * J * S -
    z * E * J * S - B * E * G * T +
    A * F * G * T +
    B * C * I * T - y * F * I * T - A * C * J * T +
    y * E * J * T +
    B * D * G * U - z * F * G * U - B * C * H * U +
    y * F * H * U +
    z * C * J * U - y * D * J * U - A * D * G * V +
    z * E * G * V +
    A * C * H * V - y * E * H * V - z * C * I * V + y * D * I * V,
    B * E * L * S - A * F * L * S - B * D * M * S + z * F * M * S + A * D * N * S -
    z * E * N * S - B * E * K * T +
    A * F * K * T +
    B * C * M * T - y * F * M * T - A * C * N * T +
    y * E * N * T +
    B * D * K * U - z * F * K * U - B * C * L * U +
    y * F * L * U +
    z * C * N * U - y * D * N * U - A * D * K * V +
    z * E * K * V +
    A * C * L * V - y * E * L * V - z * C * M * V + y * D * M * V,
    B * I * L * S - A * J * L * S - B * H * M * S + z * J * M * S + A * H * N * S -
    z * I * N * S - B * I * K * T +
    A * J * K * T +
    B * G * M * T - y * J * M * T - A * G * N * T +
    y * I * N * T +
    B * H * K * U - z * J * K * U - B * G * L * U +
    y * J * L * U +
    z * G * N * U - y * H * N * U - A * H * K * V +
    z * I * K * V +
    A * G * L * V - y * I * L * V - z * G * M * V + y * H * M * V,
    F * I * L * S - E * J * L * S - F * H * M * S + D * J * M * S + E * H * N * S -
    D * I * N * S - F * I * K * T +
    E * J * K * T +
    F * G * M * T - C * J * M * T - E * G * N * T +
    C * I * N * T +
    F * H * K * U - D * J * K * U - F * G * L * U +
    C * J * L * U +
    D * G * N * U - C * H * N * U - E * H * K * V +
    D * I * K * V +
    E * G * L * V - C * I * L * V - D * G * M * V + C * H * M * V,
    B * E * P * S - A * F * P * S - B * D * Q * S + z * F * Q * S + A * D * R * S -
    z * E * R * S - B * E * O * T +
    A * F * O * T +
    B * C * Q * T - y * F * Q * T - A * C * R * T +
    y * E * R * T +
    B * D * O * U - z * F * O * U - B * C * P * U +
    y * F * P * U +
    z * C * R * U - y * D * R * U - A * D * O * V +
    z * E * O * V +
    A * C * P * V - y * E * P * V - z * C * Q * V + y * D * Q * V,
    B * I * P * S - A * J * P * S - B * H * Q * S + z * J * Q * S + A * H * R * S -
    z * I * R * S - B * I * O * T +
    A * J * O * T +
    B * G * Q * T - y * J * Q * T - A * G * R * T +
    y * I * R * T +
    B * H * O * U - z * J * O * U - B * G * P * U +
    y * J * P * U +
    z * G * R * U - y * H * R * U - A * H * O * V +
    z * I * O * V +
    A * G * P * V - y * I * P * V - z * G * Q * V + y * H * Q * V,
    F * I * P * S - E * J * P * S - F * H * Q * S + D * J * Q * S + E * H * R * S -
    D * I * R * S - F * I * O * T +
    E * J * O * T +
    F * G * Q * T - C * J * Q * T - E * G * R * T +
    C * I * R * T +
    F * H * O * U - D * J * O * U - F * G * P * U +
    C * J * P * U +
    D * G * R * U - C * H * R * U - E * H * O * V +
    D * I * O * V +
    E * G * P * V - C * I * P * V - D * G * Q * V + C * H * Q * V,
    B * M * P * S - A * N * P * S - B * L * Q * S + z * N * Q * S + A * L * R * S -
    z * M * R * S - B * M * O * T +
    A * N * O * T +
    B * K * Q * T - y * N * Q * T - A * K * R * T +
    y * M * R * T +
    B * L * O * U - z * N * O * U - B * K * P * U +
    y * N * P * U +
    z * K * R * U - y * L * R * U - A * L * O * V +
    z * M * O * V +
    A * K * P * V - y * M * P * V - z * K * Q * V + y * L * Q * V,
    F * M * P * S - E * N * P * S - F * L * Q * S + D * N * Q * S + E * L * R * S -
    D * M * R * S - F * M * O * T +
    E * N * O * T +
    F * K * Q * T - C * N * Q * T - E * K * R * T +
    C * M * R * T +
    F * L * O * U - D * N * O * U - F * K * P * U +
    C * N * P * U +
    D * K * R * U - C * L * R * U - E * L * O * V +
    D * M * O * V +
    E * K * P * V - C * M * P * V - D * K * Q * V + C * L * Q * V,
    J * M * P * S - I * N * P * S - J * L * Q * S + H * N * Q * S + I * L * R * S -
    H * M * R * S - J * M * O * T +
    I * N * O * T +
    J * K * Q * T - G * N * Q * T - I * K * R * T +
    G * M * R * T +
    J * L * O * U - H * N * O * U - J * K * P * U +
    G * N * P * U +
    H * K * R * U - G * L * R * U - I * L * O * V +
    H * M * O * V +
    I * K * P * V - G * M * P * V - H * K * Q * V + G * L * Q * V
];
# @time G = groebner(polys)