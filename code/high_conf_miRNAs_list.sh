#!/usr/bin/env bash
mirbase_location="$@";
cat "$@/high_conf.dat" | grep 'AC  ' | awk '{split($0,a,"AC   "); split(a[2],b,";"); print b[1];}';
