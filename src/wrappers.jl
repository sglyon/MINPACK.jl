## Wrapping hybrd routine
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
               show_trace::Bool, tracing::Bool, maxit::Int;
               _n::Int=length(x0), ml::Int=_n-1, mu::Int=_n-1,
               epsfcn::Float64=0.0, diag::Vector{Float64}=fill(1.0, _n),
               mode::Int=2, factor::Float64=100.0, nprint::Int=0,
               lr::Cint=ceil(Cint, _n*(_n+1)/2))
    n = length(x0)
    x = copy(x0)
    fvec = similar(x)
    fjac = Array{Float64}(undef, n, n)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    trace = AlgoTrace(x0, show_trace, tracing, maxit)

    r = Array{Float64}(undef, lr)
    qtf = Array{Float64}(undef, n)
    wa1 = Array{Float64}(undef, n)
    wa2 = Array{Float64}(undef, n)
    wa3 = Array{Float64}(undef, n)
    wa4 = Array{Float64}(undef, n)

    if show_trace
        show(trace)
    end

    function _hybrd_func_wrapper(_p::Ptr{Cvoid}, n::Cint, _x::Ptr{Cdouble},
                                 _fvec::Ptr{Cdouble}, iflag::Cint)
        local fvec = unsafe_wrap(Array, _fvec, n)
        local x = unsafe_wrap(Array, _x, n)
        if iflag < 0
            print(fvec)
            return Cint(0)
        end
        f!(fvec, x)

        local trace = unsafe_pointer_to_objref(_p)::AlgoTrace
        push!(trace, x, fvec, iflag)

        trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
    end
    _hybrd_cfunc = @cfunction(
        $_hybrd_func_wrapper,
        Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint)
    )

    return_code = @ccall libcminpack.hybrd(
        _hybrd_cfunc::Ptr{Cvoid},
        pointer_from_objref(trace)::Ptr{Cvoid},
        n::Cint,
        x::Ptr{Cdouble},
        fvec::Ptr{Cdouble},
        tol::Cdouble,
        typemax(Cint)::Cint,
        ml::Cint,
        mu::Cint,
        epsfcn::Cdouble,
        diag::Ptr{Cdouble},
        mode::Cint,
        factor::Cdouble,
        nprint::Cint,
        [0]::Ptr{Cint},
        fjac::Ptr{Cdouble},
        n::Cint,
        r::Ptr{Cdouble},
        lr::Cint,
        qtf::Ptr{Cdouble},
        wa1::Ptr{Cdouble},
        wa2::Ptr{Cdouble},
        wa3::Ptr{Cdouble},
        wa4::Ptr{Cdouble},
    )::Cint

    msg = _hybrd_messages[max(-2, return_code)]
    if return_code < -1
        msg = msg * string(return_code)
    end

    coverged = return_code == 1
    trace.tot_time = time() - trace.start_time

    SolverResults("Modified Powell (Expert)", x0, x, fvec, return_code, coverged, msg, trace)
end

## Wrapping hybrj routine
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
               show_trace::Bool, tracing::Bool, maxit::Int;
               _n::Int=length(x0), ml::Int=_n-1, mu::Int=_n-1,
               maxfev::Int=Int(typemax(Cint)),
               epsfcn::Float64=0.0, diag::Vector{Float64}=fill(1.0, _n),
               mode::Int=2, factor::Float64=100.0, nprint::Int=0,
               lr::Cint=ceil(Cint, _n*(_n+1)/2))
    n = length(x0)
    ldfjac = n
    x = copy(x0)
    fvec = similar(x)
    fjac = Array{Float64}(undef, n, n)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    trace = AlgoTrace(x0, show_trace, tracing, maxit)

    r = Array{Float64}(undef, lr)
    qtf = Array{Float64}(undef, n)
    wa1 = Array{Float64}(undef, n)
    wa2 = Array{Float64}(undef, n)
    wa3 = Array{Float64}(undef, n)
    wa4 = Array{Float64}(undef, n)

    if show_trace
        show(trace)
    end

    function _hybrj_func_wrapper(_p::Ptr{Cvoid}, n::Cint, _x::Ptr{Cdouble},
                                 _fvec::Ptr{Cdouble}, _fjac::Ptr{Cdouble},
                                 ldfjac::Cint, iflag::Cint)
        local fvec = unsafe_wrap(Array, _fvec, n)
        fjac_flat = unsafe_wrap(Array, _fjac, ldfjac*n)
        fjac = reshape(fjac_flat, Int(ldfjac), Int(n))
        local x = unsafe_wrap(Array, _x, n)
        if iflag < 0
            print(fvec)
            return Cint(0)
        end

        if iflag == 1
            f!(fvec, x)
        elseif iflag == 2
            g!(fjac, x)
        end

        local trace = unsafe_pointer_to_objref(_p)::AlgoTrace
        push!(trace, x, fvec, iflag)

        trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
    end
    _hybrj_cfunc = @cfunction(
        $_hybrj_func_wrapper,
        Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint)
    )

    return_code = @ccall libcminpack.hybrj(
        _hybrj_cfunc::Ptr{Cvoid},
        pointer_from_objref(trace)::Ptr{Cvoid},
        n::Cint,
        x::Ptr{Cdouble},
        fvec::Ptr{Cdouble},
        fjac::Ptr{Cdouble},
        ldfjac::Cint,
        xtol::Cdouble,
        Cint(maxfev)::Cint,
        diag::Ptr{Cdouble},
        mode::Cint,
        factor::Cdouble,
        nprint::Cint,
        [0]::Ptr{Cint},
        [0]::Ptr{Cint},
        r::Ptr{Cdouble},
        lr::Cint,
        qtf::Ptr{Cdouble},
        wa1::Ptr{Cdouble},
        wa2::Ptr{Cdouble},
        wa3::Ptr{Cdouble},
        wa4::Ptr{Cdouble},
    )::Cint

    msg = _hybrj_messages[max(-2, return_code)]
    if return_code < -1
        msg = msg * string(return_code)
    end

    coverged = return_code == 1 || norm(fvec, Inf) <= xtol
    trace.tot_time = time() - trace.start_time

    SolverResults("Modified Powell (User Jac, Expert)", x0, x, fvec, return_code, coverged, msg, trace)
end

## Wrapping hybrd1 routine
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
                show_trace::Bool, tracing::Bool, maxit::Int)
    x = copy(x0)
    fvec = similar(x)
    n = length(x)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    trace = AlgoTrace(x0, show_trace, tracing, maxit)

    function _hybrd1_func_wrapper(_p::Ptr{Cvoid}, n::Cint, _x::Ptr{Cdouble},
                                 _fvec::Ptr{Cdouble}, iflag::Cint)
        local fvec = unsafe_wrap(Array, _fvec, n)
        local x = unsafe_wrap(Array, _x, n)
        if iflag < 0
            print(fvec)
            return Cint(0)
        end
        f!(fvec, x)

        local trace = unsafe_pointer_to_objref(_p)::AlgoTrace
        push!(trace, x, fvec, iflag)

        trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
    end

    _hybrd1_cfunc = @cfunction(
        $_hybrd1_func_wrapper,
        Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint)
    )

    if show_trace
        show(trace)
    end
    return_code = @ccall libcminpack.hybrd1(
        _hybrd1_cfunc::Ptr{Cvoid},
        pointer_from_objref(trace)::Ptr{Cvoid},
        n::Cint,
        x::Ptr{Cdouble},
        fvec::Ptr{Cdouble},
        tol::Cdouble,
        wa::Ptr{Cdouble},
        lwa::Cint,
    )::Cint

    msg = _hybr_messages[max(-2, return_code)]
    if return_code < -1
        msg = msg * string(return_code)
    end

    coverged = return_code == 1 || norm(fvec, Inf) <= tol
    trace.tot_time = time() - trace.start_time

    SolverResults("Modified Powell", x0, x, fvec, return_code, coverged, msg, trace)
end

## Wrapping lmdif1 routine
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
                show_trace::Bool, tracing::Bool, maxit::Int)
    x = copy(x0)
    n = length(x)
    if n > m
        msg = "Must have at least as many variables as equations"
        throw(ArgumentError(msg))
    end

    fvec = Array{Float64}(undef, m)
    lwa = m*n+5*n+m
    iwa = Array{Int}(undef, n)
    wa = Array{Float64}(undef, lwa)
    trace = AlgoTrace(x0, show_trace, tracing, maxit)

    if show_trace
        show(trace)
    end

    function _lmdif1_func_wrapper(_p::Ptr{Cvoid}, m::Cint, n::Cint, _x::Ptr{Cdouble},
                                 _fvec::Ptr{Cdouble}, iflag::Cint)
        local fvec = unsafe_wrap(Array, _fvec, m)
        local x = unsafe_wrap(Array, _x, n)
        if iflag < 0
            print(fvec)
            return Cint(0)
        end
        f!(fvec, x)

        local trace = unsafe_pointer_to_objref(_p)::AlgoTrace
        push!(trace, x, fvec, iflag)
        trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
    end
    _lmdif1_cfunc = @cfunction(
        $_lmdif1_func_wrapper,
        Cint,
        (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint)
    )

    return_code = @ccall libcminpack.lmdif1(
        _lmdif1_cfunc::Ptr{Cvoid},
        pointer_from_objref(trace)::Ptr{Cvoid},
        m::Cint,
        n::Cint,
        x::Ptr{Cdouble},
        fvec::Ptr{Cdouble},
        tol::Cdouble,
        iwa::Ptr{Cint},
        wa::Ptr{Cdouble},
        lwa::Cint,
    )::Cint

    msg = _lmdif1_messages[max(-2, return_code)]
    if return_code < -1
        msg = msg * string(return_code)
    end
    converged = return_code in [1, 2, 3] || norm(fvec, Inf) <= tol
    trace.tot_time = time() - trace.start_time

    SolverResults("Levenberg-Marquardt", x0, x, fvec, return_code, converged, msg, trace)
end

## Wrapping lmdif routine
const _lmdif_messages = Dict{Int, String}(
    0 => "improper input parameters",
    1 => string("Both actual and predicted relative errors in the sum of ",
                "squares are at most ftol"),
    2 => "relative error between two consecutive iterates is at most xtol",
    3 => string("Both actual and predicted relative errors in the sum of ",
                "squares are at most ftol",
                "\nAND relative error between two consecutive iterates is ",
                "at most xtol"),
    4 =>  string("the cosine of the angle between fvec and any column of the ",
                 "jacobian is at most gtol in absolute value"),
    5 => "number of calls to fcn has reached or exceeded maxfev",
    6 => "ftol is too small, no further reduction of sum of squares is possible",
    7 => "xtol is too small, no further improvement in x is possible",
    8 => string("gtol is too small. fvec is orthogonal to the columns of the ",
                "jacobian to machine precision"),
    -1 => "exceeded user imposed number of iterations",
    -2 => "user terminated iterations with code "
)

function lmdif(f!::Function, x0::Vector{Float64}, m::Int, tol::Float64,
                show_trace::Bool, tracing::Bool, maxit::Int;
                _n::Int = length(x0),
                gtol::Float64=0.0, ftol::Float64=tol, xtol::Float64=tol,
                epsfcn::Float64=0.0, mode::Int=1, nprint::Int=0,
                maxfev::Int=(_n+1)*200, factor::Float64=100.0,
                ldfjac::Int=m)
    x = copy(x0)
    n = length(x)
    if n > m
        msg = "Must have at least as many variables as equations"
        throw(ArgumentError(msg))
    end

    fvec = Array{Float64}(undef, m)
    lwa = m*n+5*n+m
    iwa = Array{Int}(undef, n)
    wa = Array{Float64}(undef, lwa)
    trace = AlgoTrace(x0, show_trace, tracing, maxit)

    diag = Array{Float64}(undef, n)
    nfev = [0]
    fjac = Array{Float64}(undef, m, n)
    ipvt = Array{Cint}(undef, n)
    qtf = Array{Float64}(undef, n)
    wa1 = Array{Float64}(undef, n)
    wa2 = Array{Float64}(undef, n)
    wa3 = Array{Float64}(undef, n)
    wa4 = Array{Float64}(undef, m)

    if show_trace
        show(trace)
    end

    function _lmdif_func_wrapper(_p::Ptr{Cvoid}, m::Cint, n::Cint, _x::Ptr{Cdouble},
                                 _fvec::Ptr{Cdouble}, iflag::Cint)
        local fvec = unsafe_wrap(Array, _fvec, m)
        local x = unsafe_wrap(Array, _x, n)
        if iflag < 0
            print(fvec)
            return Cint(0)
        end
        f!(fvec, x)

        local trace = unsafe_pointer_to_objref(_p)::AlgoTrace
        push!(trace, x, fvec, iflag)
        trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
    end
    _lmdif_cfunc = @cfunction(
        $_lmdif_func_wrapper,
        Cint,
        (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint)
    )

    return_code = @ccall libcminpack.lmdif(
        _lmdif_cfunc::Ptr{Cvoid},
        pointer_from_objref(trace)::Ptr{Cvoid},
        m::Cint,
        n::Cint,
        x::Ptr{Cdouble},
        fvec::Ptr{Cdouble},
        ftol::Cdouble,
        xtol::Cdouble,
        gtol::Cdouble,
        maxfev::Cint,
        epsfcn::Cdouble,
        diag::Ptr{Cdouble},
        mode::Cint,
        factor::Cdouble,
        nprint::Cint,
        nfev::Ptr{Cint},
        fjac::Ptr{Cdouble},
        ldfjac::Cint,
        ipvt::Ptr{Cint},
        qtf::Ptr{Cdouble},
        wa1::Ptr{Cdouble},
        wa2::Ptr{Cdouble},
        wa3::Ptr{Cdouble},
        wa4::Ptr{Cdouble},
    )::Cint

    msg = _lmdif_messages[max(-2, return_code)]
    if return_code < -1
        msg = msg * string(return_code)
    end
    converged = return_code in [1, 2, 3, 4] || norm(fvec, Inf) <= ftol
    trace.tot_time = time() - trace.start_time

    SolverResults("Levenberg-Marquardt (expert)", x0, x, fvec, return_code,
        converged, msg, trace)
end

## Wrapping lmder routine
const _lmder_messages = Dict{Int, String}(
    0 => "improper input parameters",
    1 => string("Both actual and predicted relative errors in the sum of ",
                "squares are at most ftol"),
    2 => "relative error between two consecutive iterates is at most xtol",
    3 => string("Both actual and predicted relative errors in the sum of ",
                "squares are at most ftol",
                "\nAND relative error between two consecutive iterates is ",
                "at most xtol"),
    4 =>  string("the cosine of the angle between fvec and any column of the ",
                 "jacobian is at most gtol in absolute value"),
    5 => "number of calls to fcn has reached or exceeded maxfev",
    6 => "ftol is too small, no further reduction of sum of squares is possible",
    7 => "xtol is too small, no further improvement in x is possible",
    8 => string("gtol is too small. fvec is orthogonal to the columns of the ",
                "jacobian to machine precision"),
    -1 => "exceeded user imposed number of iterations",
    -2 => "user terminated iterations with code "
)

function lmder(f!::Function, g!::Function, x0::Vector{Float64}, m::Int,
               tol::Float64, show_trace::Bool, tracing::Bool, maxit::Int;
               _n::Int=length(x0),
               gtol::Float64=0.0, ftol::Float64=tol, xtol::Float64=tol,
               epsfcn::Float64=0.0, mode::Int=1, nprint::Int=0,
               maxfev::Int=(_n+1)*200, factor::Float64=100.0,
               ldfjac::Int=m)
    x = copy(x0)
    n = length(x)
    if n > m
        msg = "Must have at least as many variables as equations"
        throw(ArgumentError(msg))
    end

    fvec = Array{Float64}(undef, m)
    lwa = m*n+5*n+m
    iwa = Array{Int}(undef, n)
    wa = Array{Float64}(undef, lwa)
    trace = AlgoTrace(x0, show_trace, tracing, maxit)

    diag = Array{Float64}(undef, n)
    nfev = [0]
    njev = [0]
    fjac = Array{Float64}(undef, m, n)
    ipvt = Array{Cint}(undef, n)
    qtf = Array{Float64}(undef, n)
    wa1 = Array{Float64}(undef, n)
    wa2 = Array{Float64}(undef, n)
    wa3 = Array{Float64}(undef, n)
    wa4 = Array{Float64}(undef, m)

    if show_trace
        show(trace)
    end

    function _lmder_func_wrapper(_p::Ptr{Cvoid}, m::Cint, n::Cint,
                                 _x::Ptr{Cdouble}, _fvec::Ptr{Cdouble},
                                 _fjac::Ptr{Cdouble}, ldfjac::Cint,
                                 iflag::Cint)
        local fvec = unsafe_wrap(Array, _fvec, m)
        fjac_flat = unsafe_wrap(Array, _fjac, ldfjac*n)
        fjac = reshape(fjac_flat, Int(ldfjac), Int(n))
        local x = unsafe_wrap(Array, _x, n)
        if iflag < 0
            print(fvec)
            return Cint(0)
        end

        if iflag == 1
            f!(fvec, x)
        elseif iflag == 2
            g!(fjac, x)
        end

        local trace = unsafe_pointer_to_objref(_p)::AlgoTrace
        push!(trace, x, fvec, iflag)

        trace.f_calls > trace.maxit ? Cint(-1) : Cint(0)
    end

    _lmder_cfunc = @cfunction(
        $_lmder_func_wrapper,
        Cint,
        (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint)
    )

    return_code = @ccall libcminpack.lmder(
        _lmder_cfunc::Ptr{Cvoid},
        pointer_from_objref(trace)::Ptr{Cvoid},
        m::Cint,
        n::Cint,
        x::Ptr{Cdouble},
        fvec::Ptr{Cdouble},
        fjac::Ptr{Cdouble},
        ldfjac::Cint,
        ftol::Cdouble,
        xtol::Cdouble,
        gtol::Cdouble,
        maxfev::Cint,
        diag::Ptr{Cdouble},
        mode::Cint,
        factor::Cdouble,
        nprint::Cint,
        nfev::Ptr{Cint},
        njev::Ptr{Cint},
        ipvt::Ptr{Cint},
        qtf::Ptr{Cdouble},
        wa1::Ptr{Cdouble},
        wa2::Ptr{Cdouble},
        wa3::Ptr{Cdouble},
        wa4::Ptr{Cdouble},
    )::Cint

    msg = _lmder_messages[max(-2, return_code)]
    if return_code < -1
        msg = msg * string(return_code)
    end
    converged = return_code in [1, 2, 3, 4, 8] || norm(fvec, Inf) <= ftol
    trace.tot_time = time() - trace.start_time

    SolverResults("Levenberg-Marquardt (expert)", x0, x, fvec, return_code,
        converged, msg, trace)
end

## Wrapping fdjac1 routine
function fdjac1(f!::Function, x0::Vector{Float64};
               _n::Int=length(x0), ml::Int=_n-1, mu::Int=_n-1,
               epsfcn::Float64=0.0)
    n = length(x0)
    x = copy(x0)
    fvec = similar(x)
    fjac = Array{Float64}(undef, n, n)
    ldjfac = n

    f!(x, fvec)

    wa1 = Array{Float64}(undef, n)
    wa2 = Array{Float64}(undef, n)

    function _fdjac1_func_wrapper(_p::Ptr{Cvoid}, n::Cint, _x::Ptr{Cdouble},
                                  _fvec::Ptr{Cdouble}, iflag::Cint)
        local fvec = unsafe_wrap(Array, _fvec, n)
        local x = unsafe_wrap(Array, _x, n)
        f!(fvec, x)

        Cint(0)
    end

    _fdjac1_cfunc = @cfunction(
        $_fdjac1_func_wrapper,
        Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint)
    )

    return_code = @ccall libcminpack.:fdjac1(
        _fdjac1_cfunc::Ptr{Cvoid},
        pointer_from_objref(nothing)::Ptr{Cvoid},
        n::Cint,
        x::Ptr{Cdouble},
        fvec::Ptr{Cdouble},
        fjac::Ptr{Cdouble},
        ldjfac::Cint,
        ml::Cint,
        mu::Cint,
        epsfcn::Cdouble,
        wa1::Ptr{Cdouble},
        wa2::Ptr{Cdouble},
    )::Cint

    return fjac
end
