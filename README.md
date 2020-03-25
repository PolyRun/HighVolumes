# high-volumes

ASL group project

# dependencies

You must install: GLPK

# how to build:

```
mkdir build
cd build

cmake ..
make

#for debug mode:
cmake -DCMAKE_BUILD_TYPE=Debug ..

#make just your project with dependencies:
make peterem
```

never submit the build directory!

executables will end up in build/src

# use of C / C++

The volume library is compiled with C.
Currenly, peterem is compiled with C++.
