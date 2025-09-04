using Nemo, StructuralIdentifiability

function siwr()
    ode = StructuralIdentifiability.@ODEmodel(
        S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
        I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
        W'(t) = xi * (I(t) - W(t)),
        R'(t) = gam * I(t) - (mu + a) * R(t),
        y(t) = k * I(t)
    )

    io_eqs = StructuralIdentifiability.find_ioequations(ode)
    id_funcs, bring = StructuralIdentifiability.extract_identifiable_functions_raw(io_eqs, ode, empty(ode.parameters), true)
    param_ring, _ = polynomial_ring(base_ring(bring), map(string, ode.parameters))
    id_funcs_no_states = map(polys -> map(poly -> StructuralIdentifiability.parent_ring_change(poly, param_ring), polys), id_funcs[:no_states])
    rff = StructuralIdentifiability.RationalFunctionField(id_funcs_no_states)

	# Output system
    StructuralIdentifiability.fractionfree_generators_raw(rff.mqs)[1]
end
