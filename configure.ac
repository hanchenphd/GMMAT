AC_INIT([GMMAT],[1.4.0],[Han.Chen.2@uth.tmc.edu])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([could not determine R_HOME])
fi

AC_LANG(C++)
AC_PROG_CXX

LIB_ZSTD=""
SOURCES_ZSTD=""
AC_CHECK_LIB([zstd], [ZSTD_decompress], [
  AC_CHECK_HEADERS([zstd.h], [
    LIB_ZSTD="-lzstd"
    AC_MSG_NOTICE([Will use system-wide zstd library])    
  ])                     
], [
  SOURCES_ZSTD=$(cat src/zstd.files)
  AC_MSG_NOTICE([Will build zstd library from the package])
])

LIB_LIBDEFLATE=""
SOURCES_LIBDEFLATE=""
AC_CHECK_LIB([deflate], [libdeflate_zlib_decompress], [
  AC_CHECK_HEADERS([libdeflate.h], [
    LIB_LIBDEFLATE="-ldeflate"
    AC_MSG_NOTICE([Will use system-wide libdeflate library])    
  ])                     
], [
  SOURCES_LIBDEFLATE=$(cat src/libdeflate.files)
  AC_MSG_NOTICE([Will build libdeflate library from the package])
])
             
AC_SUBST([PKG_CPPFLAGS], ["-Izstd/lib -Izstd/lib/common"])
AC_SUBST([PKG_LIBS], ["${LIB_ZSTD} ${LIB_LIBDEFLATE}"])
AC_SUBST([SOURCES_ZSTD])
AC_SUBST([SOURCES_LIBDEFLATE])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
