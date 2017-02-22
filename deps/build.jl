url = "https://github.com/devernay/cminpack/archive/master.zip"

@static if is_apple()
    download(url, "cminpack.zip")
    run(`unzip -a -q cminpzck.zip`)
    run(`cd cminpack-master && cmake -D BUILD_SHARED_LIBS=ON . -Wno-dev && make`)
    rm("cminpack.zip")
end
