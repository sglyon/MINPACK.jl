module MINPACK

using Distances

using NLsolve: SolverState, SolverTrace
export fsolve

const cminpack = joinpath(dirname(dirname(@__FILE__)), "libcminpack.dylib")

# Just a testing function. Will delete soon...
function f!(x, fvec=similar(x))
    fvec[1] = (x[1]+3)*(x[2]^3-7)+18
    fvec[2] = sin(x[2]*exp(x[1])-1)
    fvec
end

type AlgoTrace
    f_calls::Int
    g_calls::Int
    show_trace::Bool
    x_old::Vector{Float64}
    trace::Vector{SolverState{Float64}}
    prev_iflag::Int  # allows us to track when we stop computing deriv

    function AlgoTrace(x_init::Vector{Float64}, verbose::Bool=false)
        x_old = verbose ? copy(x_init) : Array{Float64}(0)
        states = Array{SolverState{Float64}}(0)
        new(0, 0, verbose, x_old, states, 1)
    end
end

Base.unsafe_convert(::Type{Ptr{Void}}, o::AlgoTrace) = o

function Base.show(io::IO, trace::AlgoTrace)
    @printf io "Iter     f(x) inf-norm    Step 2-norm \n"
    @printf io "------   --------------   --------------\n"
    for state in trace.trace
        show(io, state)
    end
end

function Base.push!(trace::AlgoTrace, x::Vector{Float64}, fvec::Vector{Float64},
                    iflag::Cint)
    if iflag == 2  # computing derivative
        if trace.prev_iflag == 1
            # only increment if just starting to compute deriv
            trace.g_calls += 1
        end
    elseif iflag == 1  # computing function
        trace.f_calls += 1
        if trace.show_trace
            x_step = sqeuclidean(trace.x_old, x)
            f_norm = maximum(abs, fvec)
            ss = SolverState(trace.f_calls, f_norm, x_step, Dict())
            show(ss)
            push!(trace.trace, ss)
            copy!(trace.x_old, x)
        end
    end
    trace.prev_iflag = iflag
end

immutable SolverResults
    algo::String
    initial_x::Vector{Float64}
    x::Vector{Float64}
    f::Vector{Float64}
    return_code::Int
    converged::Bool
    msg::String
    trace::AlgoTrace
end

# NOTE: this method was adapted from NLsolve.jl
function Base.show(io::IO, s::SolverResults)
    @printf io "Results of Nonlinear Solver Algorithm\n"
    @printf io " * Algorithm: %s\n" s.algo
    @printf io " * Starting Point: %s\n" string(s.initial_x)
    @printf io " * Zero: %s\n" string(s.x)
    @printf io " * Inf-norm of residuals: %f\n" norm(s.f, Inf)
    # @printf io " * Iterations: %d\n" r.iterations
    @printf io " * Convergence: %s\n" s.converged
    @printf io " * Message: %s\n" s.msg
    @printf io " * Function Calls: %d\n" s.trace.f_calls
    @printf io " * Jacobian Calls (df/dx): %d" s.trace.g_calls
end

## Wrapping hybrd1 routine
const _hybrd1_func_ref = Ref{Function}()
function _hybrd1_func_wrapper(_p::Ptr{Void}, n::Cint, _x::Ptr{Cdouble},
                             _fvec::Ptr{Cdouble}, iflag::Cint)
    fvec = unsafe_wrap(Array, _fvec, n)
    x = unsafe_wrap(Array, _x, n)
    if iflag < 0
        print(fvec)
        return Cint(0)
    end
    _hybrd1_func_ref[](x, fvec)

    trace = unsafe_pointer_to_objref(_p)::AlgoTrace
    push!(trace, x, fvec, iflag)

    Cint(0)
end
const _hybrd1_cfunc = cfunction(_hybrd1_func_wrapper, Cint, (Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))

const _hybr_messages = Dict{Int,String}(
    0 => "improper input parameters",
    1 => string("algorithm estimates that the relative error between x and the ",
                "solution is at most tol"),
    2 => "maximum iterations has been exceeded",
    3 => "tol is too small, no further improvement in x is possible",
    4 => "iteration is not making good progress",
    -1 => "user terminated iterations with code "
)

function hybrd1(f!::Function, x0::Vector{Float64}; tol::Float64=1e-8,
                show_trace::Bool=false)
    x = copy(x0)
    fvec = similar(x)
    n = length(x)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    _hybrd1_func_ref[] = f!
    trace = AlgoTrace(x0, show_trace)

    if show_trace
        show(trace)
    end

    return_code = ccall(
        (:hybrd1, cminpack),
        Cint,
        (Ptr{Void}, Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cint),
        _hybrd1_cfunc, pointer_from_objref(trace), n, x, fvec, tol, wa, lwa
    )

    msg = _hybr_messages[max(-1, return_code)]
    if return_code < 0
        msg = msg * string(return_code)
    end

    coverged = return_code == 1

    SolverResults("Modified Powell", x0, x, fvec, return_code, coverged, msg, trace)
end

## Wrapping lmdif1 routine
const _lmdif1_func_ref = Ref{Function}()
function _lmdif1_func_wrapper(_p::Ptr{Void}, m::Cint, n::Cint, _x::Ptr{Cdouble},
                             _fvec::Ptr{Cdouble}, iflag::Cint)
    fvec = unsafe_wrap(Array, _fvec, m)
    x = unsafe_wrap(Array, _x, n)
    if iflag < 0
        print(fvec)
        return Cint(0)
    end
    _lmdif1_func_ref[](x, fvec)

    trace = unsafe_pointer_to_objref(_p)::AlgoTrace
    push!(trace, x, fvec, iflag)

    Cint(0)
end
const _lmdif1_cfunc = cfunction(_lmdif1_func_wrapper, Cint, (Ptr{Void}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))


const _lmdif1_messages = Dict{Int,String}(
    0 => "improper input parameters",
    1 => string("algorithm estimates that the relative error between x and the ",
                "solution is at most tol"),
    2 => string("algorithm estimates that the relative error between x and the ",
                "solution is at most tol"),
    3 => string("algorithm estimates that the relative error in the sum of ",
                "squares and the relative error between x and the solution is ",
                "at most tol"),
    4 =>  "fvec is orthogonal to the columns of the jacobian to machine precision.",
    5 => "maximum iterations has been exceeded",
    6 => "tol is too small, no further reduction of sum of squares is possible",
    7 => "tol is too small, no further improvement in x is possible",
    -1 => "user terminated iterations with code "
)

# NOTE: default doesn't always hold
function lmdif1(f!::Function, x0::Vector{Float64}, m::Int=length(x0); tol::Float64=1e-8,
                show_trace::Bool=false)
    x = copy(x0)
    n = length(x)
    if n > m
        msg = "Must have at least as many variables as equations"
        throw(ArgumentError(msg))
    end

    fvec = Array{Float64}(m)
    lwa = m*n+5*n+m
    iwa = Array{Int}(n)
    wa = Array{Float64}(lwa)
    _lmdif1_func_ref[] = f!
    trace = AlgoTrace(x0, show_trace)

    return_code = ccall(
        (:lmdif1,cminpack),
        Cint,
        (Ptr{Void}, Ptr{Void}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Cint),
        _lmdif1_cfunc, pointer_from_objref(trace), m, n, x, fvec, tol, iwa, wa, lwa
    )

    msg = _lmdif1_messages[max(-1, return_code)]
    if return_code < 0
        msg = msg * string(return_code)
    end
    converged = return_code in [1, 2, 3]

    SolverResults("Levenberg-Marquardt", x0, x, fvec, return_code, converged, msg, trace)
end

function fsolve(f!::Function, x0::Vector{Float64}, m::Int=length(x0); tol::Float64=1e-8,
                show_trace::Bool=false, method::Symbol=:hybr)
    if method == :hybr
        return hybrd1(f!, x0, tol=tol, show_trace=show_trace)
    elseif method == :lm
        return lmdif1(f!, x0, m, tol=tol, show_trace=show_trace)
    else
        error("unknown method $(method)")
    end
end

end  # module
