cmake_minimum_required(VERSION 3.8.0)

project("Preemptive response time analysis" VERSION 1.0.0 LANGUAGES CXX)

include_directories(include)
include_directories(lib/include)

option(PARALLEL_RUN "Enable parallel run" OFF)
option(USE_TBB_MALLOC "Use the Intel TBB scalable memory allocator" OFF)
option(USE_JE_MALLOC "Use the Facebook jemalloc scalable memory allocator" OFF)
option(COLLECT_SCHEDULE_GRAPHS "Enable the collection of schedule graphs (disables parallel)" OFF)
option(DEBUG "Enable debugging" OFF)

if (PARALLEL_RUN AND COLLECT_SCHEDULE_GRAPHS)
    message(FATAL_ERROR "Parallel run and schedule graph collection cannot be enabled at the same time")
endif ()

if (DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
    message(NOTICE "Debug build")
else ()
    set(CMAKE_BUILD_TYPE Release)
endif ()

if (PARALLEL_RUN OR USE_TBB_MALLOC)
    find_package(TBB REQUIRED)
endif()

if (COLLECT_SCHEDULE_GRAPHS)
    add_compile_definitions(CONFIG_COLLECT_SCHEDULE_GRAPH)
    if (PARALLEL_RUN)
        message(WARNING "Schedule graph collection is not supported in parallel mode")
    endif()
elseif (PARALLEL_RUN)
    set(TBB_LIB TBB::tbb)
    add_compile_definitions(CONFIG_PARALLEL)
endif ()

if (USE_JE_MALLOC)
    find_library(ALLOC_LIB NAMES jemalloc)
    message(NOTICE "Using Facebook jemalloc scalable memory allocator")
elseif (USE_TBB_MALLOC)
    set(ALLOC_LIB TBB::tbbmalloc_proxy)
    message(NOTICE "Using Intel TBB scalable memory allocator")
endif ()

set(CORE_LIBS ${TBB_LIB} ${ALLOC_LIB})

file(GLOB TEST_SOURCES "src/tests/*.cpp")
add_executable(runtests ${TEST_SOURCES} ${SOURCES})

set(ANALYSIS_SOURCES src/rt_analysis.cpp lib/src/OptionParser.cpp)
add_executable(rt_analysis ${ANALYSIS_SOURCES})

target_link_libraries(rt_analysis ${CORE_LIBS})
target_link_libraries(runtests ${CORE_LIBS})

target_compile_features(rt_analysis PUBLIC cxx_std_14)
target_compile_features(runtests PUBLIC cxx_std_14)



