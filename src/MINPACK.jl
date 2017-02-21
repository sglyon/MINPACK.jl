module MINPACK

const cminpack = joinpath(dirname(dirname(@__FILE__)), "libcminpack.dylib")

# Just a testing function. Will delete soon...
function f!(x, fvec=similar(x))
    fvec[1] = (x[1]+3)*(x[2]^3-7)+18
    fvec[2] = sin(x[2]*exp(x[1])-1)
    fvec
end


type AlgoData
    n_calls::Int
    show_trace::Bool

    AlgoData(verbose=false) = new(0, verbose)
end

Base.unsafe_convert(::Type{Ptr{Void}}, o::AlgoData) = o

immutable SolverResults
    algo::String
    initial_x::Vector{Float64}
    x::Vector{Float64}
    f::Vector{Float64}
    return_code::Int
    msg::String
    obj::AlgoData
end

# NOTE: this method was adapted from NLsolve.jl
function Base.show(io::IO, s::SolverResults)
    @printf io "Results of Nonlinear Solver Algorithm\n"
    @printf io " * Algorithm: %s\n" s.algo
    @printf io " * Starting Point: %s\n" string(s.initial_x)
    @printf io " * Zero: %s\n" string(s.x)
    @printf io " * Inf-norm of residuals: %f\n" norm(s.f, Inf)
    # @printf io " * Iterations: %d\n" r.iterations
    @printf io " * Convergence: %s\n" s.return_code==1
    @printf io " * Message: %s\n" s.msg
    # @printf io "   * |x - x'| < %.1e: %s\n" r.xtol r.x_converged
    # @printf io "   * |f(x)| < %.1e: %s\n" r.ftol r.f_converged
    @printf io " * Function Calls: %d\n" s.obj.n_calls
    # @printf io " * Jacobian Calls (df/dx): %d" r.g_calls
end

immutable ConvergenceError <: Exception
    msg::String
    return_code::Int
end

const _hybrd1_func_ref = Array(Function)
function _hybrd1_func_wrapper(_p::Ptr{Void}, n::Cint, _x::Ptr{Cdouble},
                             _fvec::Ptr{Cdouble}, iflag::Cint)
    fvec = unsafe_wrap(Array, _fvec, n)
    x = unsafe_wrap(Array, _x, n)
    if iflag < 0
        print(fvec)
        return Cint(0)
    end
    f!(x, fvec)

    obj = unsafe_pointer_to_objref(_p)::AlgoData
    obj.n_calls += 1
    if obj.show_trace
        println("it: $(obj.n_calls)\terr: $(norm(fvec, Inf))")
    end

    _hybrd1_func_ref[](x, fvec)
    Cint(0)
end
const _hybrd1_cfunc = cfunction(_hybrd1_func_wrapper, Cint, (Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))

const _hybrd1_messages = Dict{Int,String}(
    1 => "MINPACK says algorithm estimates that the relative error between x and the solution is at most tol",
    2 => "MINPACK says maximum iterations has been exceeded",
    3 => "MINPACK says tol is too small, no further improvement in x is possible",
    4 => "MINPACK says iteration is not making good progress",
    -1 => "MINPACK says user terminated iterations with code "
)

function hybrd1(f!::Function, x0::Vector{Float64}; tol::Float64=1e-8,
                show_trace::Bool=false)
    x = copy(x0)
    fvec = similar(x)
    n = length(x)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    _hybrd1_func_ref[] = f!
    obj = AlgoData(show_trace)

    return_code = ccall(
        (:hybrd1, cminpack),
        Cint,
        (Ptr{Void}, Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cint),
        _hybrd1_cfunc, pointer_from_objref(obj), n, x, fvec, tol, wa, lwa
    )

    msg = _hybrd1_messages[max(-1, return_code)]
    if return_code < 0
        msg = msg * string(return_code)
        throw(ConvergenceError(msg, return_code))
    end

    SolverResults("Modified Powell", x0, x, fvec, return_code, msg, obj)
end

const _lmdif1_func_ref = Array(Function)
function _lmdif1_func_wrapper(_p::Ptr{Void}, m::Cint, n::Cint, _x::Ptr{Cdouble},
                             _fvec::Ptr{Cdouble}, iflag::Cint)
    fvec = unsafe_wrap(Array, _fvec, m)
    x = unsafe_wrap(Array, _x, n)
    if iflag < 0
        print(fvec)
        return Cint(0)
    end
    f!(x, fvec)

    obj = unsafe_pointer_to_objref(_p)::AlgoData
    obj.n_calls += 1
    if obj.show_trace
        println("it: $(obj.n_calls)\terr: $(norm(fvec, Inf))")
    end

    _lmdif1_func_ref[](x, fvec)
    Cint(0)
end
const _lmdif1_cfunc = cfunction(_lmdif1_func_wrapper, Cint, (Ptr{Void}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))

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
    wa = ones(lwa)
    _lmdif1_func_ref[] = f!
    obj = AlgoData(show_trace)

    return_code = ccall(
        (:lmdif1,cminpack),
        Cint,
        (Ptr{Void}, Ptr{Void}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Cint),
        _lmdif1_cfunc, pointer_from_objref(obj), m, n, x, fvec, tol, iwa, wa, lwa
    )

    # handle return code
    if return_code == 0
        msg = "MINPACK says Improper input parameters"
        throw(ArgumentError(msg))
    elseif return_code == 1
        msg = string("MINPACK says algorithm estimates that the relative error",
                     " in the sum of squares is at most tol")
    elseif return_code == 2
        msg = string("MINPACK says algorithm estimates that the relative error",
                     " between x and the solution is at most tol")
    elseif return_code == 3
        string("MINPACK says algorithm estimates that the relative error in the",
               " sum of squares and the relative error between x and the",
               " solution is at most tol")
    elseif return_code == 4
        msg = string("MINPACK says fvec is orthogonal to the columns of ",
                     "the jacobian to machine precision.")
        warn(msg)
    elseif return_code == 5
        msg = "MINPACK says maximum iterations has been exceeded"
        throw(ConvergenceError(msg, return_code))
    elseif return_code == 6
        msg = string("MINPACK says tol is too small, no further reduction of ",
                     "sum of squares is possible")
        warn(msg)
    elseif return_code == 7
        msg = "MINPACK says tol is too small, no further improvement in x is possible"
        warn(msg)
    elseif return_code < 0
        msg = "MINPACK says user terminated iterations with code $(return_code)"
        throw(ConvergenceError(msg, return_code))
    end

    SolverResults("Levenberg-Marquardt", x0, x, fvec, return_code, msg, obj)
end

end  # module
