# Building dll for windows

I followed these steps

1. Use dockcross docker images to set up system for cross compiling to windows (I have `~/bin` on my PATH):
```bash
docker run --rm dockcross/windows-x86 > ./dockcross_windows32
chmod +x ./dockcross_windows32
mv ./dockcross_windows32 ~/bin
docker run --rm dockcross/windows-x64 > ./dockcross_windows64
chmod +x ./dockcross_windows64
mv ./dockcross_windows64 ~/bin
```
2. Download cminpack-master using the julia code in build.jl
```julia
url = "https://github.com/devernay/cminpack/archive/master.zip"
download(url, "cminpack.zip")
run(`rm -rf cminpack-master`)
run(`unzip -a -q -u cminpack.zip`)
rm("cminpack.zip")
```
3. Run `cmake` and `make` using one of the dockcross images (example with 32 windows32):
```bash
cd cminpack-master
dockcross_windows32 cmake -D BUILD_SHARED_LIBS=ON . -Wno-dev
dockcross_windows32 make
```
4. Move the .dll file to somewhere else and repeat step 3 with the other copiler:
```bash
rm -rf CMakeFiles
rm CMakeCache.txt
dockcross_windows64 cmake -D BUILD_SHARED_LIBS=ON . -Wno-dev
dockcross_windows64 make
```
