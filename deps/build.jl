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
end
