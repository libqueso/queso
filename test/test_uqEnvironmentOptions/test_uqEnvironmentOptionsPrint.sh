#!/bin/bash
./test_uqEnvironmentOptionsPrint > tmp.inp
diff test.inp tmp.inp
retval=$?
rm tmp.inp
exit $retval
