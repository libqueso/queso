#!/bin/bash
./test_uqEnvironmentOptionsPrint > tmp.inp
diff "$srcdir/test_uqEnvironmentOptions/test.inp" tmp.inp
retval=$?
rm tmp.inp
exit $retval
