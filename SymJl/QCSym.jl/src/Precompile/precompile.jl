import PrecompileTools
import SymJit

PrecompileTools.@compile_workload begin
    q1 = [0.5+0im, sqrt(3)/2+0im]
    q2 = [sqrt(5)+0im, 0.0+2.0im]
    statevec = QCSym.States.StateVector([q1, q2])
    gate = QCSym.Gates.CX_Gate
    expected_statevec_after = [q1[1]*q2[1],
                                q1[1]*q2[2],
                                q1[2]*q2[2],
                                q1[2]*q2[1]]
    qc = QCSym.Circuits.QuantumCircuit(name="TestCircuit")
    qreg = QCSym.Circuits.add_qreg(qc, "q_reg_1", 2)
    QCSym.Circuits.add_gate(qc, gate, qubits_t=[qreg[2]], qubits_c=[qreg[1]],step=1, is_treat_numeric_only=false)
    U = QCSym.Circuits.assemble_symbolic_unitary(qc, false, false)
    U_simplified = Symbolics.simplify.(U)
    U_statevec_numeric = U_simplified * statevec.vector_numeric
    statevec_after_numeric_subsd = QCSym.Circuits.substitute_numerics_from_gates(U_statevec_numeric, qc.gatecollection.collections[gate])
    #display(qc)
    vars = collect(Symbolics.get_variables(U_statevec_numeric))
    U_statevec_numeric_func = Symbolics.build_function(U_statevec_numeric, vars, expression=Val{false})
    U_statevec_numeric_func2 = Symbolics.build_function(U_statevec_numeric, vars, expression=Val{true})
    U_statevec_numeric_func_SymJit = SymJit.compile_func(vars, [U_statevec_numeric...], dtype=:complex)
end

