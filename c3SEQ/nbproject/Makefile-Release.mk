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
	${OBJECTDIR}/3Seq.o \
	${OBJECTDIR}/CheckPRun.o \
	${OBJECTDIR}/FullRun.o \
	${OBJECTDIR}/GenPRun.o \
	${OBJECTDIR}/Genome.o \
	${OBJECTDIR}/GenomeSet.o \
	${OBJECTDIR}/HelpRun.o \
	${OBJECTDIR}/InfoRun.o \
	${OBJECTDIR}/Interface.o \
	${OBJECTDIR}/MatchRun.o \
	${OBJECTDIR}/Nucleotide.o \
	${OBJECTDIR}/PTable.o \
	${OBJECTDIR}/PTableFile.o \
	${OBJECTDIR}/Run.o \
	${OBJECTDIR}/Segment.o \
	${OBJECTDIR}/Sequence.o \
	${OBJECTDIR}/SequenceFile.o \
	${OBJECTDIR}/SingleBpRun.o \
	${OBJECTDIR}/String.o \
	${OBJECTDIR}/TextFile.o \
	${OBJECTDIR}/Triplet.o \
	${OBJECTDIR}/TripletRun.o \
	${OBJECTDIR}/Util.o \
	${OBJECTDIR}/VersionRun.o \
	${OBJECTDIR}/WriteRun.o \
	${OBJECTDIR}/stats.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk 3seq

3seq: ${OBJECTFILES}
	${LINK.cc} -o 3seq ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/3Seq.o: 3Seq.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/3Seq.o 3Seq.cpp

${OBJECTDIR}/CheckPRun.o: CheckPRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CheckPRun.o CheckPRun.cpp

${OBJECTDIR}/FullRun.o: FullRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FullRun.o FullRun.cpp

${OBJECTDIR}/GenPRun.o: GenPRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/GenPRun.o GenPRun.cpp

${OBJECTDIR}/Genome.o: Genome.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Genome.o Genome.cpp

${OBJECTDIR}/GenomeSet.o: GenomeSet.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/GenomeSet.o GenomeSet.cpp

${OBJECTDIR}/HelpRun.o: HelpRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HelpRun.o HelpRun.cpp

${OBJECTDIR}/InfoRun.o: InfoRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/InfoRun.o InfoRun.cpp

${OBJECTDIR}/Interface.o: Interface.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Interface.o Interface.cpp

${OBJECTDIR}/MatchRun.o: MatchRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MatchRun.o MatchRun.cpp

${OBJECTDIR}/Nucleotide.o: Nucleotide.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Nucleotide.o Nucleotide.cpp

${OBJECTDIR}/PTable.o: PTable.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/PTable.o PTable.cpp

${OBJECTDIR}/PTableFile.o: PTableFile.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/PTableFile.o PTableFile.cpp

${OBJECTDIR}/Run.o: Run.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Run.o Run.cpp

${OBJECTDIR}/Segment.o: Segment.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Segment.o Segment.cpp

${OBJECTDIR}/Sequence.o: Sequence.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Sequence.o Sequence.cpp

${OBJECTDIR}/SequenceFile.o: SequenceFile.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SequenceFile.o SequenceFile.cpp

${OBJECTDIR}/SingleBpRun.o: SingleBpRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SingleBpRun.o SingleBpRun.cpp

${OBJECTDIR}/String.o: String.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/String.o String.cpp

${OBJECTDIR}/TextFile.o: TextFile.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/TextFile.o TextFile.cpp

${OBJECTDIR}/Triplet.o: Triplet.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Triplet.o Triplet.cpp

${OBJECTDIR}/TripletRun.o: TripletRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/TripletRun.o TripletRun.cpp

${OBJECTDIR}/Util.o: Util.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Util.o Util.cpp

${OBJECTDIR}/VersionRun.o: VersionRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/VersionRun.o VersionRun.cpp

${OBJECTDIR}/WriteRun.o: WriteRun.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/WriteRun.o WriteRun.cpp

${OBJECTDIR}/stats.o: stats.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/stats.o stats.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} 3seq

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
