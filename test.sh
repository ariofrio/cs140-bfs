#!/bin/env bash

indent() {
  c='s/^/  /'
  case $(uname) in
    Darwin) sed -l "$c";;
    *)      sed -u "$c";;
  esac
}

report_correct() {
  correct=$((correct + 1))
  echo -n '✓ '
  tput setaf 2;
  echo $test
  tput sgr0
  [ "$error" ] && echo "$error" | indent
  return 0
}

report_incorrect() {
  incorrect=$((incorrect + 1))
  tput setaf 1; tput bold
  echo -n '✘ '
  tput sgr0; tput setaf 1
  echo $test
  tput sgr0
  [ "$error" ] && echo "$error" | indent
  return 0
}

section() {
  echo
  tput setaf 8
  echo "$@"
  tput sgr0
}

incorrect=0
correct=0

function compare()  {
  section "$2"
  for file in tests/*; do
    [ "$file" != "${file%.cache}" ] && continue
    error=$(bash -c "./parallel$1 0 < $file" 2>&1 1>/dev/null)
    test="./parallel$1 0 < $file"
    [ "$error" ] && report_incorrect || report_correct

    if [ "$file" = tests/3 ]; then
      root=1011
      error=$(bash -c "./parallel$1 $root < $file" 2>&1 1>/dev/null)
      test="./parallel$1 $root < $file"
      [ "$error" ] && report_incorrect || report_correct
      root=711
      error=$(bash -c "./parallel$1 $root < $file" 2>&1 1>/dev/null)
      test="./parallel$1 $root < $file"
      [ "$error" ] && report_incorrect || report_correct
    fi

    line=$(( $RANDOM * ($(wc -l < $file) - 1) / 32767 + 2 ))
    root=$(sed -n "${line}p" < $file | awk "{print \$$(( $RANDOM * 2 / 32767 + 1 ))}")
    error=$(bash -c "./parallel$1 $root < $file" 2>&1 1>/dev/null)
    test="./parallel$1 $root < $file"
    [ "$error" ] && report_incorrect || report_correct
  done > >(indent)
}

compare 0 "Serial"
compare 1 "Naive linked list"
compare 2 "Smart linked list"

echo
if [ $incorrect -gt 0 ]; then
  tput setaf 1; tput bold
  echo -n '✘ '
  tput sgr0; tput setaf 1
  echo -n 'FAIL'
  tput sgr0
else
  echo -n '✓ '
  tput setaf 2
  echo -n 'OK'
  tput sgr0
fi

echo -n ' » '
if [ $correct -gt 0 ]; then
  tput bold
  echo -n $correct
  tput sgr0
  echo -n ' correct '
fi
if [ $incorrect -gt 0 ]; then
  tput bold
  echo -n $incorrect
  tput sgr0
  echo -n ' incorrect '
fi
echo

exit $(( $incorrect ))
