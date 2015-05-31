#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/BlasrOutput.o \
	${OBJECTDIR}/DALIGN/DB.o \
	${OBJECTDIR}/DALIGN/QV.o \
	${OBJECTDIR}/DALIGN/align.o \
	${OBJECTDIR}/DalignWrapper.o \
	${OBJECTDIR}/Output.o \
	${OBJECTDIR}/SAMOutput.o \
	${OBJECTDIR}/SeedFinder.o \
	${OBJECTDIR}/SeedFinder_hashmap.o \
	${OBJECTDIR}/SeedFinder_hashmap_2bit.o \
	${OBJECTDIR}/SeedFinder_hashmap_2bit_1hashmap.o \
	${OBJECTDIR}/SeedFinder_test.o \
	${OBJECTDIR}/SeedFinder_test2.o \
	${OBJECTDIR}/Sequence.o \
	${OBJECTDIR}/Utility.o \
	${OBJECTDIR}/efficiency_test.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/newmain.o \
	${OBJECTDIR}/read_selector.o \
	${OBJECTDIR}/simple_test.o \
	${OBJECTDIR}/time_test.o \
	${OBJECTDIR}/trace_spacing_test.o


# C Compiler Flags
CFLAGS=-pg

# CC Compiler Flags
CCFLAGS=-pg
CXXFLAGS=-pg

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=`pkg-config --libs opencv`  

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}

${CND_DISTDIR}/${CND_CONF}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}
	g++ -o ${CND_DISTDIR}/${CND_CONF} ${OBJECTFILES} ${LDLIBSOPTIONS} -pg

${OBJECTDIR}/BlasrOutput.o: BlasrOutput.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BlasrOutput.o BlasrOutput.cpp

${OBJECTDIR}/DALIGN/DB.o: DALIGN/DB.c 
	${MKDIR} -p ${OBJECTDIR}/DALIGN
	${RM} "$@.d"
	$(COMPILE.c) -O2 `pkg-config --cflags opencv`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/DALIGN/DB.o DALIGN/DB.c

${OBJECTDIR}/DALIGN/QV.o: DALIGN/QV.c 
	${MKDIR} -p ${OBJECTDIR}/DALIGN
	${RM} "$@.d"
	$(COMPILE.c) -O2 `pkg-config --cflags opencv`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/DALIGN/QV.o DALIGN/QV.c

${OBJECTDIR}/DALIGN/align.o: DALIGN/align.c 
	${MKDIR} -p ${OBJECTDIR}/DALIGN
	${RM} "$@.d"
	$(COMPILE.c) -O2 `pkg-config --cflags opencv`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/DALIGN/align.o DALIGN/align.c

${OBJECTDIR}/DalignWrapper.o: DalignWrapper.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/DalignWrapper.o DalignWrapper.cpp

${OBJECTDIR}/Output.o: Output.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Output.o Output.cpp

${OBJECTDIR}/SAMOutput.o: SAMOutput.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SAMOutput.o SAMOutput.cpp

${OBJECTDIR}/SeedFinder.o: SeedFinder.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SeedFinder.o SeedFinder.cpp

${OBJECTDIR}/SeedFinder_hashmap.o: SeedFinder_hashmap.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SeedFinder_hashmap.o SeedFinder_hashmap.cpp

${OBJECTDIR}/SeedFinder_hashmap_2bit.o: SeedFinder_hashmap_2bit.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SeedFinder_hashmap_2bit.o SeedFinder_hashmap_2bit.cpp

${OBJECTDIR}/SeedFinder_hashmap_2bit_1hashmap.o: SeedFinder_hashmap_2bit_1hashmap.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SeedFinder_hashmap_2bit_1hashmap.o SeedFinder_hashmap_2bit_1hashmap.cpp

${OBJECTDIR}/SeedFinder_test.o: SeedFinder_test.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SeedFinder_test.o SeedFinder_test.cpp

${OBJECTDIR}/SeedFinder_test2.o: SeedFinder_test2.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SeedFinder_test2.o SeedFinder_test2.cpp

${OBJECTDIR}/Sequence.o: Sequence.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Sequence.o Sequence.cpp

${OBJECTDIR}/Utility.o: Utility.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Utility.o Utility.cpp

${OBJECTDIR}/efficiency_test.o: efficiency_test.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/efficiency_test.o efficiency_test.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/newmain.o: newmain.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/newmain.o newmain.cpp

${OBJECTDIR}/read_selector.o: read_selector.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/read_selector.o read_selector.cpp

${OBJECTDIR}/simple_test.o: simple_test.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/simple_test.o simple_test.cpp

${OBJECTDIR}/time_test.o: time_test.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/time_test.o time_test.cpp

${OBJECTDIR}/trace_spacing_test.o: trace_spacing_test.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags opencv` -std=c++11  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/trace_spacing_test.o trace_spacing_test.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
