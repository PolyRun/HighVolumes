# high-volumes

ASL group project

[Sources we use for the Algorithm.](./papers/README.md)

# dependencies

You must install: GLPK and
[Armadillo](http://arma.sourceforge.net/download.html)


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

# Running Tests

After you built, in the build directory:
```
./test/test.py
```
Or for a single test:
```
./test/test.py preprocess/test_preprocess
```
