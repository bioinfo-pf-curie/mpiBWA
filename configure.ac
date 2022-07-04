# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([mpiBWA mpiBWAIdx], [1.5.3])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE([1.13 foreign -Wall])

# Checks for programs.
AC_PROG_CC([mpicc])
AC_PROG_RANLIB
AM_PROG_AR

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([z], [gzdopen])

# Checks for library functions.
AC_CHECK_FUNCS([dup2]) 
AC_CHECK_FUNCS([ftruncate]) 
AC_CHECK_FUNCS([gettimeofday]) 
AC_CHECK_FUNCS([memmove]) 
AC_CHECK_FUNCS([memset]) 
AC_CHECK_FUNCS([munmap]) 
AC_CHECK_FUNCS([select]) 
AC_CHECK_FUNCS([socket]) 
AC_CHECK_FUNCS([sqrt]) 
AC_CHECK_FUNCS([strdup]) 
AC_CHECK_FUNCS([strerror]) 
AC_CHECK_FUNCS([strstr]) 
AC_CHECK_FUNCS([strtol]) 
AC_CHECK_FUNCS([strtoull])

# Checks for headers
AC_CHECK_HEADERS([pthread.h], [ LIBS="-lpthread ${LIBS}"])
AC_CHECK_FUNCS([getcwd])
AC_CHECK_HEADERS([arpa/inet.h]) 
AC_CHECK_HEADERS([fcntl.h]) 
AC_CHECK_HEADERS([inttypes.h]) 
AC_CHECK_HEADERS([limits.h]) 
AC_CHECK_HEADERS([netdb.h]) 
AC_CHECK_HEADERS([stddef.h]) 
AC_CHECK_HEADERS([stdint.h]) 
AC_CHECK_HEADERS([stdlib.h]) 
AC_CHECK_HEADERS([string.h]) 
AC_CHECK_HEADERS([sys/socket.h]) 
AC_CHECK_HEADERS([sys/time.h]) 
AC_CHECK_HEADERS([unistd.h]) 
AC_CHECK_HEADER_STDBOOL

# Check for inline
AC_C_INLINE 

# Check for prototype
AC_FUNC_ALLOCA 
AC_FUNC_FORK 
AC_FUNC_MALLOC 
AC_FUNC_MMAP 
AC_FUNC_REALLOC 
AC_FUNC_STRTOD 

# Checks for types
AC_TYPE_INT16_T 
AC_TYPE_INT32_T 
AC_TYPE_INT64_T 
AC_TYPE_INT8_T 
AC_TYPE_OFF_T 
AC_TYPE_PID_T 
AC_TYPE_SIZE_T 
AC_TYPE_SSIZE_T 
AC_TYPE_UINT16_T 
AC_TYPE_UINT32_T 
AC_TYPE_UINT64_T 
AC_TYPE_UINT8_T 


# Check for large-file support.
AC_SYS_LARGEFILE

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
