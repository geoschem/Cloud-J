# Cloud-J/src/Interfaces/CMakeLists.txt

# For now, always build standalone interface files
add_subdirectory(Standalone)

# Placeholder for other subdirectories, e.g. MAPL_ESMF
#if(CLOUDJ_EXTERNAL_CONFIG AND MAPL_ESMF)
#    add_subdirectory(MAPL_ESMF)
#endif()

# placeholder if org GEOS:
# option (BUILD_GEOS_INTERFACE FALSE)
#if (BUILD_GEOS_INTERFACE)
#    add_subdirectory(GEOS)
#endif()
