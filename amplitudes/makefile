ifndef COMPILER
COMPILER=gnu
endif

OBJSDIR = obj-$(COMPILER)
MODDIR = mod-$(COMPILER)

ifeq ("$(COMPILER)", "intel")
  FC = ifort
  SETFFOUT = -c -o
  SETLDOUT = -o
  LFLAGS = -lifcore -lirc -limf  -lpthread  -lifport -lifcoremt
  FFLAGS := $(FFLAGS) -132 -fPIC -diag-error-limit 1 -fpp -module $(MODDIR) $(OPT) $(DEBUG) $(FPE)
  ifdef DEBUG
	FFLAGS := $(FFLAGS) -check bounds
  endif
endif

ifeq ("$(COMPILER)", "gnu")
  FC = gfortran -I/$(PWD)/mod-gnu
  SETFFOUT = -c -o
  SETLDOUT = -o
  LFLAGS = -lgfortranbegin -lgfortran -lm
  FFLAGS := $(FFLAGS) -cpp  -ffree-line-length-none -fPIC -J $(MODDIR) $(OPT) $(DEBUG) $(FPE)
endif

VPATH = $(OBJSDIR) MCFMamps

##################################################################################################

# main dependencies

## general routines below

BasicDef = \
	mod_types.o \
	common_def.o \
	mod_consts_dp.o \
	set_parameters_in_amplitudes.o

AuxFunctions = \
	auxfunctions.o

# Opp = \
# 	mod_amplitudes.o \
# 	mod_superam.o \
# 	mod_MasterSubs.o \
# 	mod_resid4.o \
# 	mod_quadcut.o \
# 	mod_resid3.o \
# 	mod_tricut.o \
# 	mod_resid2.o \
# 	mod_doubcut.o \
# 	mod_checkgram.o

# Rational = \
# 	gggLLratbox.o \
# 	gggRRratbox.o \
# 	gggLLrattri3.o \
# 	gggRRrattri3.o \
# 	gggLLrattri4.o \
# 	gggRRrattri4.o \
# 	gggLLrattri5.o \
# 	gggRRrattri5.o \
# 	gggLLrattri6.o \
# 	gggRRrattri6.o \
# 	ratbox.o \
# 	rattri.o

# Amp = \
# 	ampLO_ggzz.o \
# 	ampNLO_ggzz.o \
# 	amplgggZZ_rat.o \
# 	amplgggZZ_sr.o \
# 	amplgggZZ.o \
# 	amplgggHZZ.o

Amp = \
	     func_for_h2vv_real.o \
       func_for_amp_1loop.o \
       func_for_h2vv.o \
       func_for_h2vv_2loop.o \
       mod_avec_mcfm.o \
       ampLO_gggh2zz.o \
       ampsqR_zz.o \
       ampLO_ggzz.o \
       ampLO_ggh2zz.o \
       ampNLO_ggh2zz.o \
       ampLO_ggzz_mass.o \
       DoubleTri.o \
       ampNLO_ggzz_mass.o	\
       heli_zz_heft.o \
       ampsqLO_zz.o \
       ampLO_ggww.o \
       ampLO_ggh2ww.o \
       ampsqLO_ww.o \
       Acalc.o \
       LRcalc.o \
       ZZmbc.o \
       Tri3masscoeff.o \
       ZZintegraleval.o \
       ZZD02x1x34LLmp.o \
       ZZD062x1x34LLmp.o \
       ZZC012x34LLmp.o \
       ZZC02x34LLmp.o \
       ZZC01x2LLmp.o \
       ZZC01x34LLmp.o \
       ZZbox2LL.o \
       ZZbox1LL.o \
       ZZtri1_2LL.o \
       ZZtri12_34LL.o \
       ZZtri1_34LL.o \
       ZZmassivebox.o \
       ZZmassivetri.o \
       ZZmassiveboxtri.o \
       ZZmassivebub.o \
       ggZZmassamp_new.o \
       WWd1six.o \
       WWd2six.o \
       WWd3six.o \
       WWd4six.o \
       WWbox1.o \
       WWbox3.o \
       WWbox5.o \
       WWtriangle6.o \
       WWtriangle7new.o \
       WWtriangle9new.o \
       WWtriangle11new.o \
       WWtriangle12.o \
       WWbub6.o \
       WWmbc.o \
       WWmassivebox6.o \
       WWmassivetri6.o \
       WWmassivebub.o \
       ff_DoubleTri.o \
       ffaxLO_ggzz_heft.o \
       ffaxNLO_ggzz_heft.o \
       ffvecLO_ggzz_heft.o \
       ffvecNLO_ggzz_heft.o \
       ampNLO_ggzz.o \
       ampsqNLOV_zz.o \
       ampNLO_ggww.o \
       ampNLO_ggh2ww.o \
       ampsqNLOV_ww.o \

#
OBJECTS_DP = $(BasicDef) $(AuxFunctions) $(Amp)
OBJECTS_QP = $(patsubst %.o, %_qp.o, $(Opp) $(Rational))
OBJECTS = $(OBJECTS_DP) $(OBJECTS_QP) cpolylog_main.o

#####################################################################

library: $(OBJECTS)
	$(MAKE) libamplitudes.so.1.0.0
	if [ ! -h $(OBJSDIR)/libamplitudes.so ] ; then ln -s libamplitudes.so.1.0.0  $(OBJSDIR)/libamplitudes.so; fi;

force_look :
	true

# target to generate object directory if it does not exist
$(OBJSDIR):
	mkdir -p $(OBJSDIR)

# target to generate modules directory if it does not exist
$(MODDIR):
	mkdir -p $(MODDIR)

%.o %.mod: %.f90 | $(OBJSDIR) $(MODDIR)
	$(FC) $(FFLAGS) $(SETFFOUT) $(OBJSDIR)/$@ $<

%.o: %.f | $(OBJSDIR) $(MODDIR)
	$(FC) $(FFLAGS) $(SETFFOUT) $(OBJSDIR)/$@ $<

%.o: %.F | $(OBJSDIR) $(MODDIR)
	$(FC) $(FFLAGS) $(SETFFOUT) $(OBJSDIR)/$@ $<


include amplitudes.dep

######################################################################


LLFLAGS := -Wl,-soname,libamplitudes.so.1.0.0 -Wl,--version-script,version_script

libamplitudes.so.1.0.0: $(OBJECTS)
	@echo " linking"
	@echo " library file is " $@
	$(FC) $(FFLAGS) -shared -o $(OBJSDIR)/$@ \
	$(patsubst %.o,$(OBJSDIR)/%.o,$(OBJECTS)) $(LLFLAGS)


clean:
	$(MAKE) remove
	rm -f *~
	rm -f *.so*

remove:
	rm -f $(MODDIR)/*
	rm -f $(OBJSDIR)/*

