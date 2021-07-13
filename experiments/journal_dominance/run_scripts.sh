#!/bin/sh

set -e

if [ $# -ne 2 ] 
  then
    echo "usage: ./run_scripts.sh folder step"
    exit 1
fi

FOLDER=$1
STEP=$2

if [ ! -d $FOLDER ]
  then
    echo "first argument must be an existing folder"
    exit 1
fi

if [ "${FOLDER: -1}" != "/" ]
  then
    FOLDER=${FOLDER}"/"
fi

for file in "${FOLDER}"*
do
    if [ "${file: -3}" == ".py"  ]
      then
          python2 ${file} ${STEP}
    fi
done
