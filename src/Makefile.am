EXTRA_DIST = 

lib_LTLIBRARIES = libcopter.la
libcopter_la_CPPFLAGS = -DDATADIR='"$(datadir)/Copter"' $(BOOST_CPPFLAGS)
libcopter_la_CXXFLAGS = $(OPENMP_CFLAGS)
libcopter_la_LDFLAGS = $(OPENMP_CFLAGS)
libcopter_la_SOURCES = \
    Closure.cpp Common.cpp ConsistentSPT.cpp CorrelationFunction.cpp \
    Cosmology.cpp Datafile.cpp FlowingWithTime.cpp GrowthFunction.cpp \
    HigherSPT.cpp InterpolatedPS.cpp Kaiser.cpp \
    LagrangianResummation.cpp LargeN.cpp \
    LinearPS.cpp MonteCarlo.cpp NoWigglePS.cpp PowerSpectrum.cpp \
    Quadrature.cpp RPT.cpp SPT.cpp SimpleRG.cpp SpecialFunctions.cpp \
    Spline.cpp Timer.cpp \
    array.cpp cfg.c pstring.cpp rng.c

copterincdir = $(includedir)/Copter
copterinc_HEADERS = \
    Closure.h Common.h ConsistentSPT.h CorrelationFunction.h \
    Cosmology.h Datafile.h FlowingWithTime.h GrowthFunction.h \
    HigherSPT.h InterpolatedPS.h Kaiser.h \
    LagrangianResummation.h LargeN.h \
    LinearPS.h MonteCarlo.h NoWigglePS.h ODE.h ODE.inl PowerSpectrum.h \
    Quadrature.h Quadrature.inl RPT.h SPT.h SimpleRG.h SpecialFunctions.h \
    Spline.h Timer.h \
    array.h cfg.h pstring.h rng.h
