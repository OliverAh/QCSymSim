
module MottonenStatePreparation
import LinearAlgebra
import QCSym
import Symbolics

struct StatePreparator{T<:Union{ComplexF64, Symbolics.Num, Complex{Symbolics.Num}}}
    target_state::Array{T,1}
    is_complex::Bool
    state_size::Int
    num_qubits::Int
    alphas_y::Array{Union{Float64, Symbolics.Num},2}
    alphas_z::Array{Union{Float64, Symbolics.Num},2}
    theta_y_angles::Array{Union{Float64, Symbolics.Num},2}
    theta_z_angles::Array{Union{Float64, Symbolics.Num},2}
    control_qubits::Array{Int,2} #control qubits for Y rotations

    function StatePreparator(target_state::Array{T,1}) where {T}
        if T <: ComplexF64
            is_complex = !all(isapprox.(imag(target_state), 0.0; atol=floatmin(Float64)))
            TT = ComplexF64
        else
            if T <: Symbolics.Num
                is_complex = false
                TT = Symbolics.Num
            elseif T <: Complex{Symbolics.Num}
                is_complex = true
                TT = Symbolics.Num
            end
        end

        state_size = length(target_state)
        num_qubits = Int(log2(state_size))
        alphas_y =       zeros(TT, (2^(num_qubits-1), num_qubits))
        alphas_z =       zeros(TT, (2^(num_qubits-1), num_qubits))
        theta_y_angles = zeros(TT, (2^(num_qubits-1), num_qubits))
        theta_z_angles = zeros(TT, (2^(num_qubits-1), num_qubits))
        control_qubits = ones(Int, (2^(num_qubits-1), num_qubits-1))
        new{TT}(target_state, is_complex, state_size, num_qubits, alphas_y, alphas_z, theta_y_angles, theta_z_angles, control_qubits)
    end
end


function StatePreparator(target_state::Array{T}, compute_angles_and_controls::Bool) where {T<:Union{ComplexF64, Symbolics.Num, Complex{Symbolics.Num}}}
    sp = StatePreparator(target_state)
    if compute_angles_and_controls && sp.is_complex
        compute_y_angles(sp)
        compute_z_angles(sp)
        compute_theta_y_angles(sp)
        compute_theta_z_angles(sp)
        compute_control_qubits(sp)
    elseif compute_angles_and_controls && !sp.is_complex
        compute_y_angles(sp)
        compute_theta_y_angles(sp)
        compute_control_qubits(sp)
    end
    return sp
end


function compute_y_angles(sp::StatePreparator)
    """Compute the angles for the Y rotations."""
    for k in range(1, length=sp.num_qubits)
        for j in range(1, length=2^(sp.num_qubits-k))
            l_num = range(1, length=2^(k-1))
            l_den = range(1, length=2^(k))
            #println("k:", k, " j:", j, " l_num:", l_num, " l_den:", l_den)
            numerator =   sqrt(sum( abs.( (sp.target_state[(2*j-1)*(2^(k-1)) .+ l_num]) ).^2 ))
            denominator = sqrt(sum( abs.( (sp.target_state[(  j-1)*(2^(k  )) .+ l_den]) ).^2 ))
            #println("k:", k, " j:", j, " l_num:", l_num, " l_den:", l_den, " numerator:", numerator, " denominator:", denominator)
            
            sp.alphas_y[j, k] = 2*asin(numerator / denominator)
        end
    end
end

function compute_z_angles(sp::StatePreparator)
    """Compute the angles for the Z rotations."""
    for k in range(1, length=sp.num_qubits)
        for j in range(1, length=2^(sp.num_qubits-k))
            l1 = range(1, length=2^(k-1))
            l2 = l1
            summand1 = sum(angle.(sp.target_state[(2*j-1)*(2^(k-1)) .+ l1]))
            summand2 = sum(angle.(sp.target_state[(2*j-2)*(2^(k-1)) .+ l2]))
            #print('k:', k, 'j:', j, 'l_num:', l_num, 'l_den:', l_den, 'numerator:', numerator, 'denominator:', denominator)
            
            sp.alphas_z[j, k] = (summand1 - summand2) / 2^(k-1)
        end
    end
end

function _assemble_M_matrix_transposed(k::Int)
    """Assemble the M matrix for the Mottonen state preparation."""

    k2 = 2^k

    col = UInt.(range(0, length=k2))
    gray_code_col = xor.(col .>> 1, col)
    M = gray_code_col .& transpose(col) # vectors in Julia are already (n,1), so here we get (n x n)
    M = count_ones.(M) .% 2
    M = Int.(M) .* (-2) .+ 1
    return M 
end

function compute_theta_y_angles(sp::StatePreparator)
    """Compute the angles for the Y rotations."""
    k = Int(ceil(log2(size(sp.alphas_y, 1))))
    for j in range(0, length=size(sp.alphas_y, 2))
        first_n_elems_are_nonzero = 2^(k - j)
        _M = _assemble_M_matrix_transposed(k-j)
        sp.theta_y_angles[begin:first_n_elems_are_nonzero, j+1] = (_M * sp.alphas_y[begin:first_n_elems_are_nonzero, j+1]) / first_n_elems_are_nonzero
    end
end


function compute_theta_z_angles(sp::StatePreparator)
    """Compute the angles for the actual rotation z gates."""
    k = Int(ceil(log2(size(sp.alphas_z, 1))))
    for j in range(0, length=size(sp.alphas_z,2))
        first_n_elems_are_nonzero = 2^(k - j)
        _M = _assemble_M_matrix_transposed(k-j)
        sp.theta_z_angles[begin:first_n_elems_are_nonzero, j+1] = (_M * sp.alphas_z[begin:first_n_elems_are_nonzero, j+1]) / first_n_elems_are_nonzero
        
    end
end

function compute_control_qubits(sp::StatePreparator)
    """Compute the control qubits for the rotations."""
   num_control_qubits = sp.num_qubits - 1

    for k in (sp.num_qubits - 1):-1:0
        n = num_control_qubits - k
        code = [string((i >> 1) ⊻ i, base=2, pad=n) for i in 0:(2^n - 1)]

        if n==0 || code == ["0"]
            code = ["0", "1"]
        end

        num_selections = length(code)

        c_qubits = [
            Int(floor(log2(parse(Int, code[i], base=2) ⊻ parse(Int, code[mod1(i + 1, num_selections)], base=2))))
            for i in 1:num_selections
        ]

        c_qubits = [num_control_qubits - 1 - k - c for c in c_qubits]

        if k < sp.num_qubits - 1
            sp.control_qubits[1:length(c_qubits), k + 1] .= c_qubits .+ 1
        end
    end
end




function implement_state_prep!(qc::QCSym.Circuits.QuantumCircuit,
     sp::MottonenStatePreparation.StatePreparator,
     qreg::QCSym.BitsRegs.BitRegister{QCSym.BitsRegs.QBit}=qreg)
    
    num_qubits = sp.num_qubits
    num_control_qubits = num_qubits - 1
    
     step = 1
    for k in num_control_qubits:-1:0
        num_angles_nonzero = 2^(num_control_qubits-k)
        target_qubit = num_qubits - k
        for j in 0:(2^(num_control_qubits-k)-1)
            RY_gate = QCSym.Gates.RY_Gate
            QCSym.Circuits.add_gate(qc, RY_gate, qubits_t=[qreg[target_qubit]], step=step, param_values=Dict("θ" => sp.theta_y_angles[j+1, k+1]), is_treat_alt_only=true)
            step += 1

            if j < num_angles_nonzero-1
                CX_gate = QCSym.Gates.CX_Gate
                qubits_c = sp.control_qubits[j+1, k+1]
                QCSym.Circuits.add_gate(qc, CX_gate, qubits_t=[qreg[target_qubit]], qubits_c=[qreg[qubits_c]], step=step, is_treat_numeric_only=true)
                step += 1
            end
        end
        if k < num_qubits - 1
            CX_gate = QCSym.Gates.CX_Gate
            qubits_c = num_control_qubits - k
            QCSym.Circuits.add_gate(qc, CX_gate, qubits_t=[qreg[target_qubit]], qubits_c=[qreg[1]], step=step, is_treat_numeric_only=true)
            step += 1
        end
    end

    if sp.is_complex


        for k in num_control_qubits:-1:0
            num_angles_nonzero = 2^(num_control_qubits-k)
            target_qubit = num_qubits - k
            for j in 1:2^(num_control_qubits-k)
                RZ_gate = QCSym.Gates.RZ_Gate
                QCSym.Circuits.add_gate(qc, RZ_gate, qubits_t=[qreg[target_qubit]], step=step, param_values=Dict("θ" => sp.theta_z_angles[j, k+1]), is_treat_alt_only=true)
                step += 1

                if j < num_angles_nonzero
                    CX_gate = QCSym.Gates.CX_Gate
                    qubits_c = sp.control_qubits[j, k+1]
                    QCSym.Circuits.add_gate(qc, CX_gate, qubits_t=[qreg[target_qubit]], qubits_c=[qreg[qubits_c]], step=step, is_treat_numeric_only=true)
                    step += 1
                end
            end
            if k < num_qubits - 1
                CX_gate = QCSym.Gates.CX_Gate
                qubits_c = num_control_qubits - k
                QCSym.Circuits.add_gate(qc, CX_gate, qubits_t=[qreg[target_qubit]], qubits_c=[qreg[1]], step=step, is_treat_numeric_only=true)
                step += 1
            end
        end


        GP_gate = QCSym.Gates.GP_Gate
        QCSym.Circuits.add_gate(qc, GP_gate, qubits_t=[qreg[num_qubits]], step=step, param_values=Dict("γ" => sum(angle.(sp.target_state)/length(sp.target_state))), is_treat_alt_only=true)
        step += 1
    end
    return qc
end





end
