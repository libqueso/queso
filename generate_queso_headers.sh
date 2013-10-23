#!/bin/bash

# Find all teh headers
headers=`find src -name "*.h" -a -not -name queso.h -type f | grep -v 'ANN' | grep -v 'interface'`

# Find the thing we want to conduct wizardry on
quesoh=`find src -name "queso.h.in"`

# Nuke all previous includes
sed "/#include/d" $quesoh > tmpheader
mv tmpheader $quesoh

for header_with_path in $headers; do
  header=`basename $header_with_path`;

  # We need the backslash because this is going to be used with sed
  header="#include<queso\/$header>"

  # These go in backwards but whatever
  sed "/Insert magic header foo/a ${header}" $quesoh > tmpheader
  mv tmpheader $quesoh
done
