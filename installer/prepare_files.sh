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
  cp ${topdir}/lib/tibermodels/*.so ${files}/lib/tibermodels
  cp ${topdir}/share/tibercad.ico ${files}/share
  cp ${topdir}/share/Copyright.txt ${files}/share
  if test -e ${topdir}/manual/tiber_manual.pdf ; then
    cp ${topdir}/manual/tiber_manual.pdf ${files}/doc/manual.pdf
  fi

  return
}


make_tgz () {

  mv $files $name
  echo "Creating linux/${name}.tar.gz ..."
  tar zcf linux/${name}.tar.gz $name

  mv $name $files

  return
}

make_tbz () {

  mv $files $name
  echo "Creating linux/${name}.tar.bz2 ..."
  tar jcf linux/${name}.tar.bz2 $name

  mv $name $files

  return
}

make_deb () {

  mkdir -p debfiles/usr
  mkdir -p debfiles/usr/share/tibercad
  mkdir -p debfiles/usr/share/doc/tibercad
  mv $files/bin debfiles/usr/bin
  mv $files/lib debfiles/usr/lib/tibercad
  mv ${files}/share/* debfiles/usr/share/
  mv ${files}/share/Copyright.txt debfiles/usr/share/doc/
  mv ${files}/doc/* debfiles/usr/share/doc/tibercad/

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


case $1 in

  windows )
    make_windows_package ;;

  linux ) {
    name=tibercad-$version
    prepare_linux_package

    case $2 in

      deb )
        make_deb ;;

      tgz )
        make_tgz ;;

      * )
        make_tbz ;;
    esac

    rm -rf $files
  } ;;

  * ) ;;

esac


