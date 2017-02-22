module MINPACK

using Distances

export fsolve

const cminpack = joinpath(dirname(dirname(@__FILE__)), "libcminpack.dylib")

# Just a testing function. Will delete soon...
function f!(x, fvec=similar(x))
    fvec[1] = (x[1]+3)*(x[2]^3-7)+18
    fvec[2] = sin(x[2]*exp(x[1])-1)
    fvec
end

immutable IterationState
    iteration::Int
    fnorm::Float64
    xnorm::Float64
    step_time::Float64
end

function Base.show(io::IO, t::IterationState)
    @printf io "%6d   %14e   %14e   %14e\n" t.iteration t.fnorm t.xnorm t.step_time
end

type AlgoTrace
    f_calls::Int
    g_calls::Int
    show_trace::Bool
    tracing::Bool
    x_old::Vector{Float64}
    trace::Vector{IterationState}
    prev_iflag::Int  # allows us to track when we stop computing deriv
    start_time::Float64
    last_feval_time::Float64
    tot_time::Float64
    io::IO

    function AlgoTrace(x_init::Vector{Float64}, verbose::Bool=false, tracing::Bool=false,
                       io::IO=STDOUT)
        if verbose
            tracing = true
        end
        x_old = tracing ? copy(x_init) : Array{Float64}(0)
        states = Array{IterationState}(0)
        new(0, 0, verbose, tracing, x_old, states, 1, time(), time(), NaN, io)
    end
end

Base.unsafe_convert(::Type{Ptr{Void}}, o::AlgoTrace) = o

function Base.show(io::IO, trace::AlgoTrace)
    @printf io "Iter     f(x) inf-norm    Step 2-norm      Step time\n"
    @printf io "------   --------------   --------------   --------------\n"
    for state in trace.trace
        show(io, state)
    end
end

function Base.push!(trace::AlgoTrace, x::Vector{Float64}, fvec::Vector{Float64},
                    iflag::Cint)
    if trace.tracing
        if iflag == 2  # computing derivative
            if trace.prev_iflag == 1
                # only increment if just starting to compute deriv
                trace.g_calls += 1
            end
        elseif iflag == 1  # computing function
            trace.f_calls += 1
            x_step = sqeuclidean(trace.x_old, x)
            f_norm = maximum(abs, fvec)
            elapsed = time() - trace.last_feval_time
            ss = IterationState(trace.f_calls, f_norm, x_step, elapsed)
            trace.show_trace && show(trace.io, ss)
            push!(trace.trace, ss)
            copy!(trace.x_old, x)
        end
        trace.prev_iflag = iflag
    end
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
    @printf io " * Convergence: %s\n" s.converged
    @printf io " * Message: %s\n" s.msg
    @printf io " * Total time: %f seconds\n" s.trace.tot_time
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

function hybrd1(f!::Function, x0::Vector{Float64}, tol::Float64,
                show_trace::Bool, tracing::Bool, io::IO)
    x = copy(x0)
    fvec = similar(x)
    n = length(x)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    _hybrd1_func_ref[] = f!
    trace = AlgoTrace(x0, show_trace, tracing, io)

    if show_trace
        show(io, trace)
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
    trace.tot_time = time() - trace.start_time

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
function lmdif1(f!::Function, x0::Vector{Float64}, m::Int, tol::Float64,
                show_trace::Bool, tracing::Bool, io::IO)
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
    trace = AlgoTrace(x0, show_trace, tracing, io)

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
    trace.tot_time = time() - trace.start_time

    SolverResults("Levenberg-Marquardt", x0, x, fvec, return_code, converged, msg, trace)
end

function fsolve(f!::Function, x0::Vector{Float64}, m::Int=length(x0); tol::Float64=1e-8,
                show_trace::Bool=false, tracing::Bool=false, method::Symbol=:hybr,
                io::IO=STDOUT)
    if method == :hybr
        return hybrd1(f!, x0, tol, show_trace, tracing, io)
    elseif method == :lm
        return lmdif1(f!, x0, m, tol, show_trace, tracing, io)
    else
        error("unknown method $(method)")
    end
end

end  # module
