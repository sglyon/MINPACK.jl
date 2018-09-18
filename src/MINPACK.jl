module MINPACK

using Distances
using Printf, LinearAlgebra

export fsolve

const _dl_ext = @static Sys.isapple() ? "dylib" : (Sys.iswindows() ? "dll" : "so")

const cminpack = joinpath(
    dirname(dirname(@__FILE__)), "deps", "libcminpack.$(_dl_ext)"
)

struct IterationState
    iteration::Int
    fnorm::Float64
    xnorm::Float64
    step_time::Float64
end

function Base.show(io::IO, t::IterationState)
    @printf io "%6d   %14e   %14e   %14f\n" t.iteration t.fnorm t.xnorm t.step_time
end

mutable struct AlgoTrace
    f_calls::Int
    g_calls::Int
    show_trace::Bool
    tracing::Bool
    maxit::Int
    x_old::Vector{Float64}
    trace::Vector{IterationState}
    prev_iflag::Int  # allows us to track when we stop computing deriv
    start_time::Float64
    last_feval_time::Float64
    tot_time::Float64

    function AlgoTrace(x_init::Vector{Float64}, verbose::Bool=false, tracing::Bool=false,
                       maxit::Int=typemax(Int))
        if verbose
            tracing = true
        end
        x_old = tracing ? copy(x_init) : Array{Float64}(undef, 0)
        states = Array{IterationState}(undef,0)
        new(0, 0, verbose, tracing, maxit, x_old, states, 1, time(), time(), NaN)
    end
end

Base.unsafe_convert(::Type{Ptr{Cvoid}}, o::AlgoTrace) = o

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
            now = time()
            elapsed = now - trace.last_feval_time
            trace.last_feval_time = now
            ss = IterationState(trace.f_calls, f_norm, x_step, elapsed)
            trace.show_trace && show(ss)
            push!(trace.trace, ss)
            copyto!(trace.x_old, x)
        end
        trace.prev_iflag = iflag
    end
end

struct SolverResults
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

function fsolve(f!::Function, x0::Vector{Float64}, m::Int=length(x0); tol::Float64=1e-8,
                show_trace::Bool=false, tracing::Bool=false, method::Symbol=:hybr,
                iterations::Int=typemax(Int), kwargs...)
    if method == :hybr
        return hybrd1(f!, x0, tol, show_trace, tracing, iterations)
    elseif method == :lm
        return lmdif1(f!, x0, m, tol, show_trace, tracing, iterations)
    elseif method == :lmdif
        return lmdif(f!, x0, m, tol, show_trace, tracing, iterations; kwargs...)
    elseif method == :hybrd
        return hybrd(f!, x0, tol, show_trace, tracing, iterations; kwargs...)
    else
        error("unknown method $(method)")
    end
end

function fsolve(f!::Function, g!::Function, x0::Vector{Float64}, m::Int=length(x0);
                tol::Float64=1e-8, show_trace::Bool=false, tracing::Bool=false,
                method::Symbol=:hybr, iterations::Int=typemax(Int),
                kwargs...)
    if method == :hybr
        return hybrj(f!, g!, x0, tol, show_trace, tracing, iterations; kwargs...)
    elseif method == :lm
        return lmder(f!, g!, x0, m, tol, show_trace, tracing, iterations; kwargs...)
    else
        error("unknown method $(method)")
    end
end

include("wrappers.jl")

end  # module
