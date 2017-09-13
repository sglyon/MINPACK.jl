url = "https://github.com/devernay/cminpack/archive/master.zip"

@static if is_unix()
    cd(joinpath(dirname(@__FILE__))) do
        download(url, "cminpack.zip")
        run(`rm -rf cminpack-master`)
        run(`unzip -a -q -u cminpack.zip`)
        rm("cminpack.zip")
        cd("cminpack-master") do
            run(`cmake -D BUILD_SHARED_LIBS=ON . -Wno-dev`)
            run(`make`)
        end
    end
elseif is_windows()
    if Sys.WORD_SIZE == 64
        dll_url = "https://github.com/sglyon/MINPACK.jl/releases/download/v0.0.0/libcminpack.dll"
    else
        dll_url = "https://github.com/sglyon/MINPACK.jl/releases/download/v0.0.0/libcminpack32.dll"
    end

    cd(joinpath(dirname(@__FILE__))) do
        mkdir("cminpack-master")
        cd("cminpack-master") do
            download(dll_url, "libcminpack.dll")
        end
    end
end
