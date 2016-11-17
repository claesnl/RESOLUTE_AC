cmake_minimum_required(VERSION 2.8)

# name of the cpp project
project(RESOLUTE)

set(CMAKE_INSTALL_PREFIX /usr/local/bin)
message("CIP = ${CMAKE_INSTALL_PREFIX} (should be /usr/local/bin)")

# cpp flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# specify that we have a config file. This is required, since DCMTK package for Ubuntu 14.04 
# uses /usr/include/dcmtk/config/cfunix.h to setup dcmtk. The file is only loaded by
# DCMTK's /usr/include/dcmtk/config/osconfig.h when HAVE_CONFIG_H is set.
add_definitions(-DHAVE_CONFIG_H)


# define location of the source code. This is used so that we can refer in the example
# code to the DCMIMAGES folder relatively to the source folder.
add_definitions(-DSOURCE_CODE_LOCATION="${CMAKE_SOURCE_DIR}")

# our source file
set(SOURCE_FILES ../linux/main.cpp)

# search for DCMTK library and header files 
find_package(DCMTK REQUIRED)
find_package(Threads)

# specify DCMTK header include directories
include_directories(${DCMTK_INCLUDE_DIRS})

# set output executable of our program
add_executable(RESOLUTE ${SOURCE_FILES})

# link DCMTK library files 
target_link_libraries(RESOLUTE ${DCMTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
	include_directories( ${ZLIB_INCLUDE_DIRS} )
	target_link_libraries( RESOLUTE ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

install(TARGETS RESOLUTE DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ../linux/bin/insect_modified PERMISSIONS WORLD_EXECUTE WORLD_WRITE WORLD_READ DESTINATION /usr/local/share/RESOLUTE)

# <= 2 yrs
install(FILES ../linux/bin/nihpd_asym_17-21_sinus_and_nose2.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_17-21_base11.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_17-21_basesinus.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_17-21_nasal_septa_and_eithmoidal_sinuses.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_17-21_mask.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_17-21_t1w.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_17-21_t1w.mnc DESTINATION /usr/local/share/RESOLUTE)
# <= 4 yrs
install(FILES ../linux/bin/nihpd_asym_33-44_sinus_and_nose2.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_33-44_base11.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_33-44_basesinus.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_33-44_nasal_septa_and_eithmoidal_sinuses.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_33-44_mask.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_33-44_t1w.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_33-44_t1w.mnc DESTINATION /usr/local/share/RESOLUTE)
# < 8 yrs
install(FILES ../linux/bin/nihpd_asym_04.5-08.5_sinus_and_nose2.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_04.5-08.5_base11.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_04.5-08.5_basesinus.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_04.5-08.5_nasal_septa_and_eithmoidal_sinuses.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_04.5-08.5_mask.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_04.5-08.5_t1w.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_04.5-08.5_t1w.mnc DESTINATION /usr/local/share/RESOLUTE)
# <= 11 yrs
install(FILES ../linux/bin/nihpd_asym_07.0-11.0_sinus_and_nose2.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_07.0-11.0_base11.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_07.0-11.0_basesinus.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_07.0-11.0_nasal_septa_and_eithmoidal_sinuses.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_07.0-11.0_mask.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_07.0-11.0_t1w.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_07.0-11.0_t1w.mnc DESTINATION /usr/local/share/RESOLUTE)
# <= 18 yrs
install(FILES ../linux/bin/nihpd_asym_13.0-18.5_sinus_and_nose2.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_13.0-18.5_base11.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_13.0-18.5_basesinus.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_13.0-18.5_nasal_septa_and_eithmoidal_sinuses.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_13.0-18.5_mask.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_13.0-18.5_t1w.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nihpd_asym_13.0-18.5_t1w.mnc DESTINATION /usr/local/share/RESOLUTE)
# > 18 yrs
install(FILES ../linux/bin/sinus_and_nose2.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/base11.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/basesinus.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/nasal_septa_and_eithmoidal_sinuses.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/mni_icbm152_t1_tal_nlin_sym_09a.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/mni_icbm152_t1_tal_nlin_sym_09a_mask.nii DESTINATION /usr/local/share/RESOLUTE)
install(FILES ../linux/bin/mni_icbm152_t1_tal_nlin_sym_09a.mnc DESTINATION /usr/local/share/RESOLUTE)

# MNI-models and ctime.pl (only for CentOS)
install(FILES bin/mritotal.cfg DESTINATION /opt/minc/etc/mni_autoreg/)
install(FILES bin/ctime.pl DESTINATION /usr/share/perl5/)





