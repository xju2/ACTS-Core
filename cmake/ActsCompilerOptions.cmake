# set Acts compiler flags
set(CUDA_PATH "${CUDA_ROOT}")
#set (ACTS_CXX_FLAGS "-fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=/usr/common/software/cuda/10.1.168  -ffp-contract=fast -fopenmp -Wall -Wextra -Wpedantic -Wshadow -Wunused-local-typedefs")
set (ACTS_CXX_FLAGS "-ffp-contract=fast -fopenmp -Wall -Wextra -Wpedantic -Wshadow -Wunused-local-typedefs -std=c++17")
#set (ACTS_CXX_FLAGS_DEBUG "--coverage")
set (ACTS_CXX_FLAGS_MINSIZEREL "")
set (ACTS_CXX_FLAGS_RELEASE "")
set (ACTS_CXX_FLAGS_RELWITHDEBINFO "")

#set_source_files_properties(SeedfinderTest.cpp PROPERTIES COMPILE_FLAGS -fopenmp-targets=nvptx64-nvidia-cuda)

# set Acts linker flags
#set (ACTS_EXE_LINKER_FLAGS_DEBUG "--coverage")
#set (ACTS_SHARED_LINKER_FLAGS_DEBUG "--coverage ")

# assign to global CXX flags
set (CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${ACTS_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} ${ACTS_CXX_FLAGS_DEBUG}")
set (CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${ACTS_CXX_FLAGS_MINSIZEREL}")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ACTS_CXX_FLAGS_RELEASE}")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${ACTS_CXX_FLAGS_RELWITHDEBINFO}")

# assign to global linker flags
set (CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${ACTS_EXE_LINKER_FLAGS_DEBUG}")
set (CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} ${ACTS_SHARED_LINKER_FLAGS_DEBUG}")

# silence warning about missing RPATH on Mac OSX
set (CMAKE_MACOSX_RPATH 1)

#ifeq ($(CXX),clang++)
#
#    CXXFLAGS = -O2 -ffast-math -ffp-contract=fast -fstrict-aliasing -Werror -Wall -Wno-unused-variable
#    CXXFLAGS += $(DEFINE)
#    CXXFLAGS += -std=c++11
#    CXXFLAGS += -lm
#    ifeq ($(OPENMP),y)
#        CXXFLAGS += -fopenmp
#    endif
#    ifeq ($(OPENMP_TARGET),y)
#        CXXFLAGS += -fopenmp
#        CXXFLAGS += -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=${CUDA_PATH} -ffp-contract=fast
#        CXXFLAGS += -D__NO_MATH_INLINES -U__SSE2_MATH__ -U__SSE_MATH__
#    endif
#endif
