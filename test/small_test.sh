#!/bin/sh

# sanitise environment so all tools behave correctly
unset LC_ALL
unset LANG
LC_CTYPE="C"
LC_COLLATE="C"
LC_TIME="C"
LC_NUMERIC="C"
LC_MONETARY="C"
LC_MESSAGES="C"

MYTMP=/dev/NONEXISTANT

cleanup()
{
    echo "Removing ${MYTMP} ..."
    [ -d "${MYTMP}" ] && rm -r "${MYTMP}"
}

# catch interrupts and terms
trap 'cleanup' 0 1 2 3 15

# exit whenever a simple command returns non-zero
set -e

# error when reading from an undefined variable
set -u

# create MYTMP (THIS DOESNT WORK ON BSD AND MAC)
MYTMP=$(mktemp -d)

PROG="$1/lrcaller"
if [ ! -x "$PROG" ]; then
    echo "ERROR: lrcaller binary not found/executable." >&2
    echo "--------------------------------------" >&2
    exit 104
fi

DATADIR="$(realpath $(dirname $0))/small_data/"
if [ ! -d "$DATADIR" ]; then
    echo "ERROR: data directory not found/executable." >&2
    echo "--------------------------------------" >&2
    exit 105
fi

echo "Test start."

${PROG} -lsf 0.2 -w 100 -gtm ad -fa "${DATADIR}/chr1reg.fa" "${DATADIR}/reads.bam" "${DATADIR}/input.vcf" "${MYTMP}/output.vcf"

echo "Test done."

cd ${MYTMP}
ACTUAL=$(openssl md5 output.vcf)

cd ${DATADIR}
CONTROL=$(openssl md5 output.vcf)

if [ "$CONTROL" != "$ACTUAL" ]; then
    echo "Output files are not as expected."
    echo "GOT:"
    echo "-------"
    echo "$ACTUAL"
    echo "-------"
    echo "EXPECTED:"
    echo "-------"
    echo "$CONTROL"
    echo "DIFF:"
    diff -u "${MYTMP}/output.vcf" "${DATADIR}/output.vcf"
    exit 1
fi

