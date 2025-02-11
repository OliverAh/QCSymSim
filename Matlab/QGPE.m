clear all

global PX PY PZ
PX = struct([])

function X = pauli_x(id)
    name = "X"
    name = append(name, id)
    X = sym(name,[2 2])
    global PX
    PX(1).(id) = X
end
pauli_x("a")
pauli_x("b")
PX

out_xa_xb = kron(PX.a, PX.b)

phi = sym("phi", [2 1])

phiphi = kron(phi,phi)

ctranspose(out_xa_xb)*phiphi