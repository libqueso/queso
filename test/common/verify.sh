#!/bin/bash

function message_fail {
    echo " "
    echo "--------------------------------------------------------------------"
    echo '(rtest): ** ERROR:' $1
    echo "--------------------------------------------------------------------"
    echo " "
    exit 1
}

function message_passed {
    echo " "
    echo "--------------------------------------------------------------------"
    echo '(rtest): PASSED:' $1
    echo "--------------------------------------------------------------------"
    echo " "
}

function verify_file_exists {
    if [ ! -e "$1" ];then
	message_fail "$1 file does not exist"
    fi
}
