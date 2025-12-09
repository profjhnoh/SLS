cmake_minimum_required(VERSION 3.10)

# project setting
SET ( PROJ_NAME Tbps_DLSLS )
PROJECT ( ${PROJ_NAME} CXX C )

# 실행 파일을 빌드 디렉토리에 생성 (증분 빌드를 위해)
# 빌드 후 루트로 복사하는 것은 POST_BUILD 명령으로 처리
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

# 실수 정밀도 옵션 (기본값: DOUBLE)
OPTION(USE_FLOAT "Use float precision instead of double" ON)

if(USE_FLOAT)
    ADD_DEFINITIONS(-DUSE_FLOAT)
    MESSAGE(STATUS "Using FLOAT precision")
else()
    MESSAGE(STATUS "Using DOUBLE precision (default)")
endif()

# C++ settings
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# 증분 빌드 최적화 설정
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # compile_commands.json 생성 (IDE 지원)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)  # 불필요한 재빌드 방지

# MinGW specific settings
if(MINGW)
    MESSAGE(STATUS "Building with MinGW")
    # Add MinGW specific flags
    # -Wa,-mbig-obj: Allow large object files (필수: 섹션이 많은 파일 컴파일)
    # -O2: 최적화 레벨 2 (릴리스 빌드용)
    # -march=native: CPU에 최적화된 명령어 사용 (AVX, SSE 등)
    # -mavx2: AVX2 명령어 활성화 (SIMD vectorization)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wa,-mbig-obj")

    # Debug 빌드: 최적화 없음, 디버그 정보 포함
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

    # Release 빌드: 최적화 활성화 + SIMD vectorization
    # -O3: 최고 수준 최적화 (loop unrolling, inline 등)
    # -march=native: 현재 CPU에 맞는 SIMD 명령어 자동 선택
    # -ffast-math: 부동소수점 연산 최적화 (IEEE 754 완전 준수 X)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -DNDEBUG -DEIGEN_NO_DEBUG")

    # 기본 빌드 타입이 설정되지 않았으면 Release로 설정
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    endif()
endif()

# Eigen3 설정
# 먼저 프로젝트 내부의 eigen3 디렉토리 확인
if(EXISTS "${CMAKE_SOURCE_DIR}/eigen3/Eigen/Core")
    set(EIGEN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/eigen3")
    MESSAGE(STATUS "Found local Eigen3 at: ${EIGEN3_INCLUDE_DIR}")
else()
    # 로컬에 없으면 시스템에서 찾기
    find_package(Eigen3 3.3 QUIET NO_MODULE)
    if(NOT Eigen3_FOUND)
        MESSAGE(STATUS "Eigen3 not found via find_package, searching manually...")
        if(MINGW)
            # MinGW/MSYS2 paths
            find_path(EIGEN3_INCLUDE_DIR
                NAMES Eigen/Core
                PATHS
                    $ENV{MINGW_PREFIX}/include/eigen3
                    /mingw64/include/eigen3
                    /mingw32/include/eigen3
                    C:/msys64/mingw64/include/eigen3
                    C:/msys64/mingw32/include/eigen3
                PATH_SUFFIXES eigen3
            )
        else()
            # macOS: Homebrew paths
            find_path(EIGEN3_INCLUDE_DIR
                NAMES Eigen/Core
                PATHS
                    /opt/homebrew/include/eigen3
                    /usr/local/include/eigen3
                PATH_SUFFIXES eigen3
            )
        endif()

        if(EIGEN3_INCLUDE_DIR)
            MESSAGE(STATUS "Found Eigen3 at: ${EIGEN3_INCLUDE_DIR}")
        else()
            MESSAGE(FATAL_ERROR "Eigen3 not found. Please download Eigen3 and extract to ${CMAKE_SOURCE_DIR}/eigen3")
        endif()
    endif()
endif()

# OpenMP 설정
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    MESSAGE(STATUS "Found OpenMP via find_package")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    MESSAGE(STATUS "OpenMP not found via find_package, trying manual configuration...")
    if(MINGW)
        # MinGW usually has OpenMP support built-in with -fopenmp
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        find_library(OpenMP_LIBRARY NAMES gomp pthread
            PATHS
                $ENV{MINGW_PREFIX}/lib
                /mingw64/lib
                /mingw32/lib
                C:/msys64/mingw64/lib
                C:/msys64/mingw32/lib
        )
        if(OpenMP_LIBRARY)
            MESSAGE(STATUS "Found OpenMP library at: ${OpenMP_LIBRARY}")
            set(OpenMP_omp_LIBRARY ${OpenMP_LIBRARY})
        else()
            MESSAGE(STATUS "Using MinGW built-in OpenMP support")
        endif()
    elseif(APPLE)
        # macOS: Homebrew libomp
        find_library(OpenMP_LIBRARY NAMES omp PATHS /opt/homebrew/lib NO_DEFAULT_PATH)
        find_path(OpenMP_INCLUDE_DIR NAMES omp.h PATHS /opt/homebrew/include NO_DEFAULT_PATH)

        if(OpenMP_LIBRARY AND OpenMP_INCLUDE_DIR)
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
            set(OpenMP_CXX_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY ${OpenMP_LIBRARY})
            MESSAGE(STATUS "Found OpenMP at: ${OpenMP_LIBRARY}")
        endif()
    endif()
endif()


# set source and headers
# set source and headers
SET (SRC_FILES
    src/main.cpp 
    src/allocateMemory.cpp 
    src/BS.cpp 
    src/channel_update.cpp 
    src/channel.cpp 
    src/generateLSP.cpp 
    src/generateSSP.cpp 
    src/Initiallization.cpp 
    src/Link.cpp 
    src/Logging_Point.cpp 
    src/measure.cpp 
    src/MS.cpp 
    src/scheduling.cpp 
    src/paramFastfading.cpp 
    src/PerdropStatistics.cpp 
    src/Print.cpp 
    src/receive_downlink.cpp 
    src/sector.cpp 
    src/setFastfadingParam.cpp 
    src/setSimulParam.cpp 
    src/standard.cpp
)

# Include directories
INCLUDE_DIRECTORIES(h)

# Platform-specific includes
if(MINGW)
    # MinGW/MSYS2 specific includes
    if(EIGEN3_INCLUDE_DIR)
        INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
        MESSAGE(STATUS "Adding Eigen3 include: ${EIGEN3_INCLUDE_DIR}")
    endif()
    # Check for MINGW_PREFIX environment variable
    if(DEFINED ENV{MINGW_PREFIX})
        INCLUDE_DIRECTORIES($ENV{MINGW_PREFIX}/include)
    endif()
else()
    # macOS specific paths
    INCLUDE_DIRECTORIES(/opt/homebrew/include)
    INCLUDE_DIRECTORIES(/opt/homebrew/include/eigen3/unsupported)
    INCLUDE_DIRECTORIES(/usr/local/include)
    if(EIGEN3_INCLUDE_DIR)
        INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
    endif()
endif()

# Eigen3 unsupported headers
if(EIGEN3_INCLUDE_DIR)
    get_filename_component(EIGEN3_ROOT ${EIGEN3_INCLUDE_DIR} DIRECTORY)
    if(EXISTS "${EIGEN3_INCLUDE_DIR}/unsupported")
        INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR}/unsupported)
        MESSAGE(STATUS "Adding Eigen3 unsupported: ${EIGEN3_INCLUDE_DIR}/unsupported")
    elseif(EXISTS "${EIGEN3_ROOT}/include/eigen3/unsupported")
        INCLUDE_DIRECTORIES(${EIGEN3_ROOT}/include/eigen3/unsupported)
        MESSAGE(STATUS "Adding Eigen3 unsupported: ${EIGEN3_ROOT}/include/eigen3/unsupported")
    elseif(EXISTS "${EIGEN3_INCLUDE_DIR}/../unsupported")
        INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR}/../unsupported)
        MESSAGE(STATUS "Adding Eigen3 unsupported: ${EIGEN3_INCLUDE_DIR}/../unsupported")
    endif()
endif()

# build
ADD_EXECUTABLE( ${PROJ_NAME} ${SRC_FILES} )

# 빌드 후 실행 파일을 루트 디렉토리로 복사 (선택적)
# 이렇게 하면 증분 빌드가 정상적으로 작동하면서도 루트에 실행 파일이 생성됨
add_custom_command(TARGET ${PROJ_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:${PROJ_NAME}>
        ${CMAKE_SOURCE_DIR}/$<TARGET_FILE_NAME:${PROJ_NAME}>
    COMMENT "Copying executable to project root directory"
)

# Link libraries
if(MINGW)
    # MinGW linking
    MESSAGE(STATUS "Configuring MinGW linking...")

    # MinGW OpenMP linking
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${PROJ_NAME} OpenMP::OpenMP_CXX)
        MESSAGE(STATUS "Linking OpenMP via find_package")
    elseif(OpenMP_omp_LIBRARY)
        target_link_libraries(${PROJ_NAME} ${OpenMP_omp_LIBRARY})
        MESSAGE(STATUS "Linking OpenMP: ${OpenMP_omp_LIBRARY}")
    else()
        # MinGW GCC usually has libgomp built-in with -fopenmp flag
        MESSAGE(STATUS "Using MinGW built-in OpenMP (libgomp)")
    endif()

else()
    # macOS/Linux linking
    if(OpenMP_omp_LIBRARY)
        target_link_libraries(${PROJ_NAME} ${OpenMP_omp_LIBRARY})
    elseif(OpenMP_CXX_FOUND)
        target_link_libraries(${PROJ_NAME} OpenMP::OpenMP_CXX)
    endif()
endif()

# Summary
MESSAGE(STATUS "")
MESSAGE(STATUS "========== Build Configuration Summary ==========")
MESSAGE(STATUS "Platform: ${CMAKE_SYSTEM_NAME}")
MESSAGE(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
if(USE_FLOAT)
    MESSAGE(STATUS "Precision: FLOAT")
else()
    MESSAGE(STATUS "Precision: DOUBLE")
endif()
if(EIGEN3_INCLUDE_DIR)
    MESSAGE(STATUS "Eigen3: Found")
else()
    MESSAGE(STATUS "Eigen3: Not found")
endif()
if(OpenMP_CXX_FOUND OR OpenMP_omp_LIBRARY OR MINGW)
    MESSAGE(STATUS "OpenMP: Enabled")
else()
    MESSAGE(STATUS "OpenMP: Disabled")
endif()
MESSAGE(STATUS "=================================================")
MESSAGE(STATUS "")
