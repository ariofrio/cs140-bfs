#!/bin/sh

tempfile=$(mktemp)
echo "./sequential 1 < sample.input" 
echo "./parallel 1 < sample.input"
./sequential 1 < sample.input &> $tempfile
./parallel 1 < sample.input 2>&1 | diff --side-by-side - $tempfile
ret=$?
rm $tempfile

exit $ret
