#!/bin/bash
(cd build/hemocell/ && cmake . && make -j4 && cp libhemocell_pre_all_deps.a libhemocell_local.a)
rm core/*.cpp
rm helper/*.cpp
rm mechanics/*.cpp
rm IO/*.cpp
rm config/*.cpp

patch -p1 < patch/make_precompiled.patch

git add -u
git add -f build/hemocell/libhemocell_local.a
