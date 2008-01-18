#!/bin/sh



make_windows_package () {

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

  return
}


prepare_linux_package () {

  mkdir -p ${files}/bin
  mkdir -p ${files}/lib/tibermodels
  mkdir -p ${files}/doc
  mkdir -p ${files}/examples
  mkdir -p ${files}/share

  cp ${topdir}/bin/tibercad ${files}/bin
  cp ${topdir}/lib/lib*.so* ${files}/lib
  find ${topdir}/lib/tibermodels -name "*.so" -exec cp {} ${files}/lib/tibermodels \;
  cp ${topdir}/share/tibercad.ico ${files}/share
  cp ${topdir}/share/Copyright.txt ${files}/share
  if test -e ${topdir}/manual/tiber_manual.pdf ; then
    cp ${topdir}/manual/tiber_manual.pdf ${files}/doc/manual.pdf
  fi

  return
}


make_tgz () {

  echo "Creating linux/${name}.tar.gz ..."
  tar zcf linux/${name}.tar.gz $name

  return
}

make_tbz () {

  echo "Creating linux/${name}.tar.bz2 ..."
  tar jcf linux/${name}.tar.bz2 $name

  return
}

make_deb () {

  files=debfiles/usr
  mkdir -m 0755 debfiles
  mkdir -p debfiles/usr/bin
  mkdir -p debfiles/usr/lib/tibermodels
  mkdir -p debfiles/usr/share/tibercad
  mkdir -p debfiles/usr/share/doc/tibercad

  cp ${topdir}/bin/tibercad ${files}/bin
  cp ${topdir}/lib/lib*.so* ${files}/lib
  find ${topdir}/lib/tibermodels -name "*.so" -exec cp {} ${files}/lib/tibermodels \;
  cp ${topdir}/share/tibercad.ico ${files}/share/tibercad
  cp ${topdir}/share/Copyright.txt ${files}/share/tibercad
  if test -e ${topdir}/manual/tiber_manual.pdf ; then
    cp ${topdir}/manual/tiber_manual.pdf ${files}/share/doc/tibercad/manual.pdf
  fi

  size=`du -s debfiles | awk '{print $1}'` 
  rm -rf debfiles/DEBIAN
  mkdir -m 0755 debfiles/DEBIAN
  sed "s/%computed_size%/$size/" linux/deb/control > debfiles/DEBIAN/control
  cd debfiles
  find usr/ -type f -exec md5sum {} >> DEBIAN/md5sums \;
  cd ..
  dpkg-deb --build debfiles linux
  rm -rf debfiles

  return
}





################### main ####################################

topdir=..
version=`awk '/TIBERVERSION/ {print $3}' \
  ${topdir}/include/base/common/tiber_config.h | sed 's/"//g'`

test $# -ge 1 || exit 1

files=$1/files

# delete everything
rm -rf ${files}

umask 022

case $1 in

  windows )
    make_windows_package ;;

  linux ) {
    name=tibercad-$version

    case $2 in

      deb )
        make_deb ;;

      tgz )
        files=$name
        prepare_linux_package
        make_tgz
        rm -rf $files ;;

      * )
        files=$name
        prepare_linux_package
        make_tbz
        rm -rf $files ;;
    esac

  } ;;

  * ) ;;

esac


