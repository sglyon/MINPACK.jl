if Sys.isapple()
    url = "https://github.com/sglyon/MINPACK.jl/releases/download/v0.0.0/libcminpack.dylib"
    download(url, joinpath(dirname(@__FILE__), "libcminpack.dylib"))
elseif Sys.iswindows()
    if Sys.WORD_SIZE == 64
        dll_url = "https://github.com/sglyon/MINPACK.jl/releases/download/v0.0.0/libcminpack.dll"
    else
        dll_url = "https://github.com/sglyon/MINPACK.jl/releases/download/v0.0.0/libcminpack32.dll"
    end
    download(dll_url, joinpath(dirname(@__FILE__), "libcminpack.dll"))
elseif Sys.islinux()
    if Sys.WORD_SIZE == 64
        so_url = "https://github.com/sglyon/MINPACK.jl/releases/download/v0.0.0/libcminpack.so"
    else
        so_url = "https://github.com/sglyon/MINPACK.jl/releases/download/v0.0.0/libcminpack32.so"
    end
    download(so_url, joinpath(dirname(@__FILE__), "libcminpack.so"))
end
