# MINPACK

[![Build Status](https://travis-ci.org/sglyon/MINPACK.jl.svg?branch=master)](https://travis-ci.org/sglyon/MINPACK.jl)

[![Windows Build status](https://ci.appveyor.com/api/projects/status/hr1fjl9ldk62ql8v?svg=true)](https://ci.appveyor.com/project/spencerlyon2/minpack-jl)

[![Coverage Status](https://coveralls.io/repos/sglyon/MINPACK.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/sglyon/MINPACK.jl?branch=master)

[![codecov.io](http://codecov.io/github/sglyon/MINPACK.jl/coverage.svg?branch=master)](http://codecov.io/github/sglyon/MINPACK.jl?branch=master)

Julia interface to [cminpack](https://github.com/devernay/cminpack), a C/C++ rewrite of the MINPACK software (originally in fortran).

## Usage

Usage is quite simple, there are two main API methods:

```julia
fsolve(f!::Function, x0::Vector{Float64}, m::Int=length(x0); tol::Float64=1e-8,
       show_trace::Bool=false, tracing::Bool=false, method::Symbol=:hybr,
       iterations::Int=typemax(Int), io::IO=STDOUT, kwargs...)

fsolve(f!::Function, g!::Function, x0::Vector{Float64}, m::Int=length(x0);
       tol::Float64=1e-8, show_trace::Bool=false, tracing::Bool=false,
       method::Symbol=:hybr, iterations::Int=typemax(Int), io::IO=STDOUT,
       kwargs...)
```

The functions `f!` and `g!` should accept the current point (call it `x`) as the _second_ argument and fill the first argument with the function values and Jacobian matrix, repsectively. If no Jacobian is passed, one will be approximated using finite differences.

Example:

```julia
julia> using MINPACK

julia> function f!(fvec, x)
           fvec[1] = (x[1]+3)*(x[2]^3-7)+18
           fvec[2] = sin(x[2]*exp(x[1])-1)
           fvec
       end;

julia> function g!(fjac, x)
           fjac[1, 1] = x[2]^3 - 7
           fjac[1, 2] = 3 * (x[1] + 3) * x[2]*x[2]
           fjac[2, 1] = x[2] * exp(x[1]) * cos(x[2] * exp(x[1]) - 1)
           fjac[2, 2] = exp(x[1]) * cos(x[2] * exp(x[1]) - 1)
           fjac
       end
g! (generic function with 2 methods)

julia> res_jac = fsolve(f!, g!, ones(2))
Results of Nonlinear Solver Algorithm
 * Algorithm: Modified Powell (User Jac, Expert)
 * Starting Point: [1.0, 1.0]
 * Zero: [6.05177e-12, 1.0]
 * Inf-norm of residuals: 0.000000
 * Convergence: true
 * Message: algorithm estimates that the relative error between x and the solution is at most tol
 * Total time: 0.033416 seconds
 * Function Calls: 0
 * Jacobian Calls (df/dx): 0

julia> res_nojac = fsolve(f!, ones(2))
Results of Nonlinear Solver Algorithm
 * Algorithm: Modified Powell
 * Starting Point: [1.0, 1.0]
 * Zero: [6.05138e-12, 1.0]
 * Inf-norm of residuals: 0.000000
 * Convergence: true
 * Message: algorithm estimates that the relative error between x and the solution is at most tol
 * Total time: 0.000024 seconds
 * Function Calls: 0
 * Jacobian Calls (df/dx): 0
```

The additional available keyword arguments captured by `;kwargs...` vary by the method used.

The keyword argument `method` can take on different value depending on which method of `fsolve` you are calling.

Available methods for the version where only `f!` is pased are:

- `:hybr`: Modified version of Powell's algorithm. Uses MINPACK routine [`hybrd1`](https://github.com/devernay/cminpack/blob/d1f5f5a273862ca1bbcf58394e4ac060d9e22c76/hybrd1.c)
- `:lm`: Levenberg-Marquardt. Uses MINPACK routine [`lmdif1`](https://github.com/devernay/cminpack/blob/d1f5f5a273862ca1bbcf58394e4ac060d9e22c76/lmdif1.c)
- `:lmdif`: Advanced Levenberg-Marquardt (more options available with `;kwargs...`). See MINPACK routine [`lmdif`](https://github.com/devernay/cminpack/blob/d1f5f5a273862ca1bbcf58394e4ac060d9e22c76/lmdif.c) for more information
- `:hybrd`: Advacned modified version of Powell's algorithm (more options available with `;kwargs...`). See MINPACK routine [`hybrd`](https://github.com/devernay/cminpack/blob/d1f5f5a273862ca1bbcf58394e4ac060d9e22c76/hybrd.c) for more information

Available methods for the version where both `f!` and `g!` are passed are:

- `:hybr`: Advacned modified version of Powell's algorithm with user supplied Jacobian. Additional arguments are available via `;kwargs...`. See MINPACK routine [`hybrj`](https://github.com/devernay/cminpack/blob/d1f5f5a273862ca1bbcf58394e4ac060d9e22c76/hybrj.c) for more information
- `:lm`: Advanced Levenberg-Marquardt with user supplied Jacobian. Additional arguments are available via `;kwargs...`. See MINPACK routine [`lmder`](https://github.com/devernay/cminpack/blob/d1f5f5a273862ca1bbcf58394e4ac060d9e22c76/lmder.c) for more information
