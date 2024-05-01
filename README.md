# how to use

Clone the reposity.
```bash
git clone https://github.com/Pupoin/dd4hep2fbx.git
```

add a line into parent's CMakelists.txt
```cmake
add_subdirectory(dd4hep2FBX)
```

In parent's folder
```bash
mkdir build
cd build 
cmake ..
make -j 20
# LD and bin
export LD_LIBRARY_PATH="/home/wln/STCFGeo/build/lib:$LD_LIBRARY_PATH" 
export PATH="$PATH:/home/wln/STCFGeo/build/bin"
```

# bugs
https://1drv.ms/p/s!At2Rqn8R-kGdhvdQU6wXxBG3SLWkXQ?e=jFxffn
