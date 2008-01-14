#!/bin/sh

topdir=..




files=windows/files

# delete everything
rm -rf ${files}

mkdir -p ${files}/doc
mkdir -p ${files}/examples

# copy all files
cp ${topdir}/bin/tibercad ${files}
cp ${topdir}/lib/*.dll ${files}
cp ${topdir}/share/tibercad.ico ${files}
cp ${topdir}/share/Copyright.txt ${files}
if test -e ${topdir}/manual/tiber_manual.pdf ; then
  cp ${topdir}/manual/tiber_manual.pdf ${files}/doc/manual.pdf
fi

# copy cygwin system libraries
test ! -x /usr/bin/cygcheck.exe && exit 1

for lib in `cygcheck windows/files/tibercad.exe | awk '/cyg/'` ; do
  cp `cygpath -u -a $lib` ${files}
done
