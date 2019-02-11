#!/bin/bash
for i in `grep -rhoP [^\\\\s]+\\\\.R 1204_patches/ | sort | uniq`; do
	git checkout -- "$i"
done
git pull
for i in `ls 1204_patches`; do
	patch -p0 < 1204_patches/$i
done
