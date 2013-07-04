## DIRECTORY: /BRAIN/test
SET(EXTERN_INCLUDE
	D:/Development/boost/boost_1_49_0 # Boost.
)
INCLUDE_DIRECTORIES(. ${EXTERN_INCLUDE})
LINK_DIRECTORIES(D:/Development/boost/boost_1_49_0/stage/lib)

SET(PARAM
	ConfigLoader.h
	Param.h	Param.cpp
)
SET(VTFORMULAE
	Polynomials.h
	VietesFormulae.h	VietesFormulae.cpp
)
SET(NGFORMULAE
	${VTFORMULAE}
	Polynomials.h
	NewtonGirardFormulae.h	NewtonGirardFormulae.cpp
)
SET(PEAK
	Peak.h	Peak.cpp
	PeakList.h	PeakList.cpp
)
SET(PERIODICTABLE
	ConfigLoader.h	
	Element.h	Element.cpp
	PeriodicTable.h	PeriodicTable.cpp
)
SET(COMPOSITION
	${PERIODICTABLE} 
	Composition.h	Composition.cpp
)

SET(ISODIST
	${PEAK}
	${PARAM}
	${NGFORMULAE}
	${COMPOSITION}
	IsotopicDistribution.h IsotopicDistribution.cpp
	IsotopicConstants.h IsotopicConstants.cpp
)
SET(IO
	${ISODIST}
	CompoIO.h	CompoIO.cpp
)

#ADD_EXECUTABLE(param ${PARAM} ./test/ParamTest.cpp)
#ADD_EXECUTABLE(peak ${PEAK} ./test/PeakTest.cpp)
ADD_EXECUTABLE(vietesformulae ${VTFORMULAE} ./test/VietesFormulaeTest.cpp)
ADD_EXECUTABLE(isodist ${ISODIST} ./test/IsotopicDistributionTest.cpp)
ADD_EXECUTABLE(newtonformulae ${NGFORMULAE} ./test/NewtonGirardFormulaeTest.cpp)
#ADD_EXECUTABLE(periodictable ${PERIODICTABLE} ./test/PeriodicTableTest.cpp)
#ADD_EXECUTABLE(composition ${COMPOSITION} ./test/CompositionTest.cpp)
ADD_EXECUTABLE(useBRAIN ${IO} ./test/CompoIOTest.cpp)

INSTALL(TARGETS useBRAIN
		RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin
		CONFIGURATIONS Release
)

#ENABLE_TESTING()
#ADD_TEST(01_composition composition)