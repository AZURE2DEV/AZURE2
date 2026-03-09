IF(NOT QT4_FOUND)
	FIND_PACKAGE( Qt4 REQUIRED QUIET )
ENDIF(NOT QT4_FOUND)

IF( QT4_FOUND )
	# Is Qwt6 installed? Look for header files
	FIND_PATH( Qwt6_INCLUDE_DIR qwt.h 
               PATHS ${QT_INCLUDE_DIR} /usr/local/qwt/include /usr/local/include /usr/include/qwt /usr/include /usr/local/qwt-6.0.1 /usr/local/qwt-6.0.0 $ENV{HOME}/local/qwt-6.0.0 $ENV{HOME}/local/qwt-6.0.1 /usr/local/qwt-6.1.0-svn /afs/crc.nd.edu/user/n/nsl/nuclear/x86_64_linux/src/qwt-6.1/src
               PATH_SUFFIXES qwt qwt6 qwt-qt4 qwt6-qt4 include qwt/include qwt6/include qwt-qt4/include qwt6-qt4/include lib/qwt.framework/Versions/6/Headers ENV PATH)
	
	# Find Qwt version
	IF( Qwt6_INCLUDE_DIR )
		FILE( READ ${Qwt6_INCLUDE_DIR}/qwt_global.h QWT_GLOBAL_H )
		STRING( REGEX MATCH "#define *QWT_VERSION *(0x06*)" QWT_IS_VERSION_6 ${QWT_GLOBAL_H})
		
		IF( QWT_IS_VERSION_6 )

		STRING(REGEX REPLACE ".*#define[\\t\\ ]+QWT_VERSION_STR[\\t\\ ]+\"([0-9]+\\.[0-9]+\\.[0-9]+)\".*" "\\1" Qwt_VERSION "${QWT_GLOBAL_H}")
		FIND_LIBRARY( Qwt6_Qt4_LIBRARY NAMES qwt6-qt4 qwt-qt4 qwt6 qwt PATHS $ENV{HOME}/local/qwt-6.0.0/lib $ENV{HOME}/local/qwt-6.0.1/lib /usr/local/qwt-6.0.0/lib /usr/local/qwt-6.0.1/lib /usr/local/qwt-6.1.0-svn/lib /usr/local/qwt/lib /usr/local/lib /usr/lib /afs/crc.nd.edu/user/n/nsl/nuclear/x86_64_linux/src/qwt-6.1/lib ${QT_LIBRARY_DIR})
						
		ENDIF( QWT_IS_VERSION_6 )		
	
	ENDIF( Qwt6_INCLUDE_DIR )

	include(FindPackageHandleStandardArgs)
	find_package_handle_standard_args(Qwt DEFAULT_MSG
       					  Qwt6_Qt4_LIBRARY Qwt6_INCLUDE_DIR)
        MARK_AS_ADVANCED( Qwt6_INCLUDE_DIR Qwt6_Qt4_LIBRARY )

ENDIF( QT4_FOUND )
