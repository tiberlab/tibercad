#!/bin/sh



make_windows_package () {

  mkdir -m 0755 ${files}
  mkdir -p ${files}/doc
  mkdir -p ${files}/examples
  mkdir -p ${files}/materials

  # copy all files
  cp ${topdir}/bin/tibercad ${files}
  cp ${topdir}/lib/*.dll ${files}
  cp ${topdir}/share/tibercad.ico ${files}
  cp ${topdir}/share/Copyright.txt ${files}
  if test -e ${topdir}/manual/tiber_manual.pdf ; then
    cp ${topdir}/manual/tiber_manual.pdf ${files}/doc/manual.pdf
  fi
  cp -r ${topdir}/Tutorials/[^.]* ${files}/examples
  find ${files} -type d -name ".svn*" -exec rm -rf {} \;
  cp ${topdir}/materials/[^.]* ${files}/materials

  curdir=`pwd`
  cd ${files}
  find ./materials -type f -exec unix2dos {} \;
  find ./examples -type f -exec unix2dos {} \;
  unix2dos Copyright.txt

  #zip -r -l data.zip materials examples
  #zip -r  data.zip materials examples
  #cp /usr/bin/unzipsfx.exe data.exe
  #cat data.zip >> data.exe
  #zip -A data.exe
  #rm -f data.zip

  cd $curdir

  # copy cygwin system libraries
  test ! -x /usr/bin/cygcheck.exe && exit 1

  for lib in `cygcheck windows/files/tibercad.exe | awk '/cyg/'` ; do
    cp `cygpath -u -a $lib` ${files}
  done

  return
}


prepare_linux_package () {

  mkdir -m 0755 ${files}
  mkdir -p ${files}/bin
  mkdir -p ${files}/lib/tibermodels
  mkdir -p ${files}/doc
  mkdir -p ${files}/examples
  mkdir -p ${files}/share
  mkdir -p ${files}/license
  mkdir -p ${files}/materials

  cp ${topdir}/installer/linux/tar/install ${files}
  chmod a+x ${files}/install
  cp ${topdir}/installer/linux/tar/tibercad.sh ${files}/bin
  chmod a+x ${files}/bin/tibercad.sh

  # get some of the libraries
  libs=`ldd ${topdir}/lib/libtibercad.so | grep boost | awk '{print $3}'`
  libs="$libs `ldd ${topdir}/lib/libtibercad.so | grep lapack | awk '{print $3}'`"
  libs="$libs `ldd ${topdir}/lib/libtibercad.so | grep blas | awk '{print $3}'`"
  libs="$libs `ldd ${topdir}/lib/libtibercad.so | grep stdc | awk '{print $3}'`"
  libs="$libs `ldd ${topdir}/lib/libtibercad.so | grep gcc_s | awk '{print $3}'`"
  cp $libs ${files}/lib
  cp ${topdir}/bin/tibercad ${files}/bin/tibercad-$version
  chmod a+x ${files}/bin/tibercad-$version
  cp ${topdir}/lib/lib*.so* ${files}/lib
  find ${topdir}/lib/tibermodels -name "*.so" -exec cp {} ${files}/lib/tibermodels \;
  cp ${topdir}/share/tibercad.ico ${files}/share
  cp ${topdir}/share/Copyright.txt ${files}/share
  if test -e ${topdir}/manual/tiber_manual.pdf ; then
    cp ${topdir}/manual/tiber_manual.pdf ${files}/doc/manual.pdf
  fi
  cp -r ${topdir}/Tutorials/[^.]* ${files}/examples
  find ${files} -type d -name ".svn*" -exec rm -rf {} \;
  cp ${topdir}/materials/[^.]* ${files}/materials

  chmod -R a+r ${files}
  find ${files} -type d -exec chmod a+x {} \;

  su -c "chown -R 0:0 ${files} && tar $MODE $ARCHIVE $files && chmod a+w $ARCHIVE && rm -rf $files"

  return
}


make_deb () {

  bindir=debfiles/usr/bin
  files=debfiles/usr/share/tibercad-$version
  mkdir -m 0755 debfiles
  mkdir -p ${bindir}
  mkdir -p ${files}/lib/tibermodels
  mkdir -p ${files}/bin
  mkdir -p ${files}/doc
  mkdir -p ${files}/examples
  mkdir -p ${files}/license
  mkdir -p ${files}/materials

  cp ${topdir}/bin/tibercad ${files}/bin/tibercad-$version
  chmod a+x ${files}/bin/tibercad-$version
  cp ${topdir}/lib/lib*.so* ${files}/lib
  find ${topdir}/lib/tibermodels -name "*.so" -exec cp {} ${files}/lib/tibermodels \;
  cp ${topdir}/share/tibercad.ico ${files}
  cp ${topdir}/share/Copyright.txt ${files}
  if test -e ${topdir}/manual/tiber_manual.pdf ; then
    cp ${topdir}/manual/tiber_manual.pdf ${files}/doc/manual.pdf
  fi
  cp -r ${topdir}/Tutorials/[^.]* ${files}/examples
  find ${files} -type d -name ".svn*" -exec rm -rf {} \;
  cp ${topdir}/materials/[^.]* ${files}/materials

  chmod -R a+r debfiles
  find debfiles -type d -exec chmod a+x {} \;

  sed -e "s/<INSTALLDIR>/\/usr\/share\/tibercad-$version/g" ${topdir}/installer/linux/tar/tibercad.sh > ${files}/bin/tibercad.sh
  chmod a+x ${files}/bin/tibercad.sh
  cd ${bindir}
  ln -s ../share/tibercad-${version}/bin/tibercad.sh tibercad
  cd ../../..

  size=`du -s debfiles | awk '{print $1}'` 
  rm -rf debfiles/DEBIAN
  mkdir -m 0755 debfiles/DEBIAN
  sed "s/%computed_size%/$size/" linux/deb/control > debfiles/DEBIAN/control
  cd debfiles
  find usr/ -type f -exec md5sum {} >> DEBIAN/md5sums \;
  cd ..
  su -c "chown -R 0:0 debfiles && dpkg-deb --build debfiles linux && chmod a+w linux/*.deb && rm -rf debfiles"

  return
}





################### main ####################################

topdir=..
version=`cat $topdir/VERSION`

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
        make_deb
        ;;

      tgz )
        files=$name
        echo "Creating linux/${name}.tar.gz ..."
        ARCHIVE=linux/${name}.tar.gz
        MODE=zcf
        prepare_linux_package
        ;;

      * )
        files=$name
        echo "Creating linux/${name}.tar.bz2 ..."
        ARCHIVE=linux/${name}.tar.bz2
        MODE=jcf
        prepare_linux_package
        ;;
    esac

  } ;;

  * ) ;;

esac


