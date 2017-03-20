## Wrapping hybrd routine
const _hybrd_func_ref = Ref{Function}()
function _hybrd_func_wrapper(_p::Ptr{Void}, n::Cint, _x::Ptr{Cdouble},
                             _fvec::Ptr{Cdouble}, iflag::Cint)
    fvec = unsafe_wrap(Array, _fvec, n)
    x = unsafe_wrap(Array, _x, n)
    if iflag < 0
        print(fvec)
        return Cint(0)
    end
    _hybrd_func_ref[](x, fvec)

    trace = unsafe_pointer_to_objref(_p)::AlgoTrace
    push!(trace, x, fvec, iflag)

    trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
end
const _hybrd_cfunc = cfunction(_hybrd_func_wrapper, Cint, (Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))

const _hybrd_messages = Dict{Int, String}(
    0 => "improper input parameters",
    1 => string("algorithm estimates that the relative error between x and the ",
                "solution is at most tol"),
    2 => "maximum iterations has been exceeded",
    3 => "tol is too small, no further improvement in x is possible",
    4 => string("iteration is not making good progress, measured by ",
                "improvement from last 5 jacobian evaluations"),
    5 => string("iteration is not making good progress, measured by ",
                "improvement from last 10 iterations"),
    -1 => "exceeded user imposed number of iterations",
    -2 => "user terminated iterations with code "
)

function hybrd(f!::Function, x0::Vector{Float64}, tol::Float64,
               show_trace::Bool, tracing::Bool, maxit::Int, io::IO;
               _n::Int=length(x0), ml::Int=_n-1, mu::Int=_n-1,
               epsfcn::Float64=0.0, diag::Vector{Float64}=fill(1.0, _n),
               mode::Int=2, factor::Float64=100.0, nprint::Int=0,
               lr::Cint=ceil(Cint, _n*(_n+1)/2))
    n = length(x0)
    x = copy(x0)
    fvec = similar(x)
    fjac = Array{Float64}(n, n)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    _hybrd_func_ref[] = f!
    trace = AlgoTrace(x0, show_trace, tracing, maxit, io)

    r = Array{Float64}(lr)
    qtf = Array{Float64}(n)
    wa1 = Array{Float64}(n)
    wa2 = Array{Float64}(n)
    wa3 = Array{Float64}(n)
    wa4 = Array{Float64}(n)

    if show_trace
        show(io, trace)
    end

    return_code = ccall(
        (:hybrd, cminpack),
        Cint,
        (
            Ptr{Void}, Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble,
            Cint, Cint, Cint, Cdouble, Ptr{Cdouble}, Cint, Cdouble, Cint,
            Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble},
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}
        ),
        _hybrd_cfunc, pointer_from_objref(trace), n, x, fvec, tol, typemax(Cint),
        ml, mu, epsfcn, diag, mode, factor, nprint, [0], fjac, n, r, lr, qtf,
        wa1, wa2, wa3, wa4
    )

    msg = _hybrd_messages[max(-2, return_code)]
    if return_code < -1
        msg = msg * string(return_code)
    end

    coverged = return_code == 1
    trace.tot_time = time() - trace.start_time

    SolverResults("Modified Powell (Expert)", x0, x, fvec, return_code, coverged, msg, trace)
end

## Wrapping hybrj routine
const _hybrj_func_ref = Ref{Function}()
const _hybrj_jac_func_ref = Ref{Function}()
function _hybrj_func_wrapper(_p::Ptr{Void}, n::Cint, _x::Ptr{Cdouble},
                             _fvec::Ptr{Cdouble}, _fjac::Ptr{Cdouble},
                             ldfjac::Cint, iflag::Cint)
    fvec = unsafe_wrap(Array, _fvec, n)
    fjac_flat = unsafe_wrap(Array, _fjac, ldfjac*n)
    fjac = reshape(fjac_flat, Int(ldfjac), Int(n))
    x = unsafe_wrap(Array, _x, n)
    if iflag < 0
        print(fvec)
        return Cint(0)
    end

    if iflag == 1
        _hybrj_func_ref[](x, fvec)
    elseif iflag == 2
        _hybrj_jac_func_ref[](x, fjac)
    end

    trace = unsafe_pointer_to_objref(_p)::AlgoTrace
    push!(trace, x, fvec, iflag)

    trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
end
const _hybrj_cfunc = cfunction(_hybrj_func_wrapper, Cint, (Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint))

const _hybrj_messages = Dict{Int, String}(
    0 => "improper input parameters",
    1 => string("algorithm estimates that the relative error between x and the ",
                "solution is at most tol"),
    2 => "maximum iterations has been exceeded",
    3 => "tol is too small, no further improvement in x is possible",
    4 => string("iteration is not making good progress, measured by ",
                "improvement from last 5 jacobian evaluations"),
    5 => string("iteration is not making good progress, measured by ",
                "improvement from last 10 iterations"),
    -1 => "exceeded user imposed number of iterations",
    -2 => "user terminated iterations with code "
)

function hybrj(f!::Function, g!::Function, x0::Vector{Float64}, xtol::Float64,
               show_trace::Bool, tracing::Bool, maxit::Int, io::IO;
               _n::Int=length(x0), ml::Int=_n-1, mu::Int=_n-1,
               maxfev::Int=Int(typemax(Cint)),
               epsfcn::Float64=0.0, diag::Vector{Float64}=fill(1.0, _n),
               mode::Int=2, factor::Float64=100.0, nprint::Int=0,
               lr::Cint=ceil(Cint, _n*(_n+1)/2))
    n = length(x0)
    ldfjac = n
    x = copy(x0)
    fvec = similar(x)
    fjac = Array{Float64}(n, n)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    _hybrj_func_ref[] = f!
    _hybrj_jac_func_ref[] = g!
    trace = AlgoTrace(x0, show_trace, tracing, maxit, io)

    r = Array{Float64}(lr)
    qtf = Array{Float64}(n)
    wa1 = Array{Float64}(n)
    wa2 = Array{Float64}(n)
    wa3 = Array{Float64}(n)
    wa4 = Array{Float64}(n)

    if show_trace
        show(io, trace)
    end

    return_code = ccall(
        (:hybrj, cminpack),
        Cint,
        (
            Ptr{Void}, Ptr{Void}, Cint,
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble, Cint,
            Ptr{Cdouble}, Cint, Cdouble, Cint, Ptr{Cint}, Ptr{Cint},
            Ptr{Cdouble}, Cint,
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}
        ),
        _hybrj_cfunc, pointer_from_objref(trace), n,
        x, fvec, fjac, ldfjac, xtol, Cint(maxfev),
        diag, mode, factor, nprint, [0], [0],
        r, lr,
        qtf, wa1, wa2, wa3, wa4
    )

    msg = _hybrj_messages[max(-2, return_code)]
    if return_code < -1
        msg = msg * string(return_code)
    end

    coverged = return_code == 1
    trace.tot_time = time() - trace.start_time

    SolverResults("Modified Powell (User Jac, Expert)", x0, x, fvec, return_code, coverged, msg, trace)
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

    trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
end
const _hybrd1_cfunc = cfunction(_hybrd1_func_wrapper, Cint, (Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))

const _hybr_messages = Dict{Int, String}(
    0 => "improper input parameters",
    1 => string("algorithm estimates that the relative error between x and the ",
                "solution is at most tol"),
    2 => "maximum iterations has been exceeded",
    3 => "tol is too small, no further improvement in x is possible",
    4 => "iteration is not making good progress",
    -1 => "exceeded user imposed number of iterations",
    -2 => "user terminated iterations with code "
)

function hybrd1(f!::Function, x0::Vector{Float64}, tol::Float64,
                show_trace::Bool, tracing::Bool, maxit::Int, io::IO)
    x = copy(x0)
    fvec = similar(x)
    n = length(x)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    _hybrd1_func_ref[] = f!
    trace = AlgoTrace(x0, show_trace, tracing, maxit, io)

    if show_trace
        show(io, trace)
    end

    return_code = ccall(
        (:hybrd1, cminpack),
        Cint,
        (Ptr{Void}, Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cint),
        _hybrd1_cfunc, pointer_from_objref(trace), n, x, fvec, tol, wa, lwa
    )

    msg = _hybr_messages[max(-2, return_code)]
    if return_code < -1
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
    trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
end
const _lmdif1_cfunc = cfunction(_lmdif1_func_wrapper, Cint, (Ptr{Void}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))


const _lmdif1_messages = Dict{Int, String}(
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
    -1 => "exceeded user imposed number of iterations",
    -2 => "user terminated iterations with code "
)

# NOTE: default doesn't always hold
function lmdif1(f!::Function, x0::Vector{Float64}, m::Int, tol::Float64,
                show_trace::Bool, tracing::Bool, maxit::Int, io::IO)
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
    trace = AlgoTrace(x0, show_trace, tracing, maxit, io)

    if show_trace
        show(io, trace)
    end

    return_code = ccall(
        (:lmdif1, cminpack),
        Cint,
        (Ptr{Void}, Ptr{Void}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Cint),
        _lmdif1_cfunc, pointer_from_objref(trace), m, n, x, fvec, tol, iwa, wa, lwa
    )

    msg = _lmdif1_messages[max(-2, return_code)]
    if return_code < -1
        msg = msg * string(return_code)
    end
    converged = return_code in [1, 2, 3]
    trace.tot_time = time() - trace.start_time

    SolverResults("Levenberg-Marquardt", x0, x, fvec, return_code, converged, msg, trace)
end


## Wrapping fdjac1 routine
const _fdjac1_func_ref = Ref{Function}()
function _fdjac1_func_wrapper(_p::Ptr{Void}, n::Cint, _x::Ptr{Cdouble},
                              _fvec::Ptr{Cdouble}, iflag::Cint)
    fvec = unsafe_wrap(Array, _fvec, n)
    x = unsafe_wrap(Array, _x, n)
    _fdjac1_func_ref[](x, fvec)

    Cint(0)
end
const _fdjac1_cfunc = cfunction(_fdjac1_func_wrapper, Cint, (Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))

function fdjac1(f!::Function, x0::Vector{Float64};
               _n::Int=length(x0), ml::Int=_n-1, mu::Int=_n-1,
               epsfcn::Float64=0.0)
    n = length(x0)
    x = copy(x0)
    fvec = similar(x)
    fjac = Array{Float64}(n, n)
    ldjfac = n
    _fdjac1_func_ref[] = f!

    f!(x, fvec)

    wa1 = Array{Float64}(n)
    wa2 = Array{Float64}(n)

    return_code = ccall(
        (:fdjac1, cminpack),
        Cint,
        (
            Ptr{Void}, Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
            Cint, Cint, Cint, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}

        ),
        _fdjac1_cfunc, pointer_from_objref(nothing), n, x, fvec, fjac,
        ldjfac, ml, mu, epsfcn, wa1, wa2
    )

    return fjac
end
