#!/bin/env bash

section() {
  echo
  tput setaf 8
  echo "$@"
  tput sgr0
}

indent() {
  c='s/^/  /'
  case $(uname) in
    Darwin) sed -l "$c";;
    *)      sed -u "$c";;
  esac
}

function compare()  {
  section "./parallel$1 1 < sample.input"
  if [ 1 ]; then
    ret=0
    tempfile=$(mktemp)
    ./sequential 1 < sample.input &> $tempfile
    ./parallel 1 < sample.input 2>&1 | diff --side-by-side - $tempfile
    ret=$?
    rm $tempfile
    return $ret
  fi > >(indent)
}

compare 0
compare 1

exit $ret
