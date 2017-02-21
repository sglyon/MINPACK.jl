module MINPACK

const cminpack = joinpath(dirname(dirname(@__FILE__)), "libcminpack.dylib")

# Just a testing function. Will delete soon...
function f!(x, fvec=similar(x))
    fvec[1] = (x[1]+3)*(x[2]^3-7)+18
    fvec[2] = sin(x[2]*exp(x[1])-1)
    fvec
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
    _hybrd1_func_ref[](x, fvec)
    Cint(0)
end
const _hybrd1_cfunc = cfunction(_hybrd1_func_wrapper, Cint, (Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))

function hybrd1(f!::Function, x0::Vector{Float64}; tol::Float64=1e-8)
    x = copy(x0)
    fvec = similar(x)
    n = length(x)
    lwa = ceil(Int, (n*(3*n+13))/2)
    wa = ones(lwa)
    _hybrd1_func_ref[] = f!

    return_code = ccall(
        (:hybrd1, cminpack),
        Cint,
        (Ptr{Void}, Ptr{Void}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cint),
        _hybrd1_cfunc, &nothing, n, x, fvec, tol, wa, lwa
    )

    # handle return code
    if return_code == 0
        msg = "MINPACK says Improper input parameters"
        throw(ArgumentError(msg))
    elseif return_code == 1
        # success
        nothing
    elseif return_code == 2
        msg = "MINPACK says maximum iterations has been exceeded"
        throw(ConvergenceError(msg, return_code))
    elseif return_code == 3
        msg = "MINPACK says tol is too small, no further improvement in x is possible"
        warn(msg)
    elseif return_code == 4
        msg = "MINPACK says iteration is not making good progress"
        warn(msg)
    elseif return_code < 0
        msg = "MINPACK says user terminated iterations with code $(return_code)"
        throw(ConvergenceError(msg, return_code))
    end

    x, fvec
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
    _lmdif1_func_ref[](x, fvec)
    Cint(0)
end
const _lmdif1_cfunc = cfunction(_lmdif1_func_wrapper, Cint, (Ptr{Void}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint))

# NOTE: default doesn't always hold
function lmdif1(f!::Function, x0::Vector{Float64}, m::Int=length(x0); tol::Float64=1e-8)
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

    return_code = ccall(
        (:lmdif1,cminpack),
        Cint,
        (Ptr{Void}, Ptr{Void}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Cint),
        _lmdif1_cfunc, &nothing, m, n, x, fvec, tol, iwa, wa, lwa
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

    x, fvec
end

end  # module
