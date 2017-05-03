# Note:  The environment variable MD_ARCH must be defined.  Use 
#        "setenv MD_ARCH `uname -s`" in your .bashrc or equivalent.
#
# Note:  To use the "make install" option, the environment vairable
#        MD_DIR should be defined to something like /usr/share/md/bin .
#        Use "setenv MD_DIR /usr/share/md/bin" in .bashrc or equivalent.

ANIMATINCS = implic genpar qpar anapar

ANIMATCOMS = $(ANIMATINCS) commons

ANIMATSRCS = animain.f qdynstuff.f Jstuff.f iostuff.f qdynmd.f ranstuff.f \
             geniostuff.f qpotstuff.f ewaldstuff.f mathstuff.f

ANIMATOBJS = ${ANIMATSRCS:.f=.o}

ANIMATLIBS = -lsunperf

FLUCQINCS = implic genpar

FLUCQCOMS = $(FLUCQINCS) commons qpar

FLUCQSRCS = qdyn.f qdynstuff.f qdynmd.f qpotstuff.f ewaldstuff.f Jstuff.f \
            ranstuff.f iostuff.f geniostuff.f mathstuff.f dgesv.f

FLUCQOBJS = ${FLUCQSRCS:.f=.o}

#FLUCQLIBS = -lsunperf
FLUCQLIBS = 

PROPSINCS = implic genpar

PROPSCOMS = $(PROPSINCS) commons qpar propar

PROPSSRCS = propsm.f qdynstuff.f qdynmd.f qpotstuff.f ewaldstuff.f Jstuff.f \
mathstuff.f ranstuff.f iostuff.f geniostuff.f dgesv.f

PROPSOBJS = ${PROPSSRCS:.f=.o}

#PROPSLIBS = -lsunperf
PROPSLIBS = 

# This makefile requires a VPATH-enabled make.  This includes both
# /usr/ccs/bin/make on Suns and GNU make.
VPATH = ..

# This should be a SVR4-compatible install:
INSTALL = /usr/sbin/install

.f.o:
	$(FC) $(FFLAGS) -c $<

flucq:	flucq_$(MD_ARCH)

dflucq:	dflucq_$(MD_ARCH)

props:	props_$(MD_ARCH)

animat: animat_$(MD_ARCH)

clean:
	-rm flucq_$(MD_ARCH) dflucq_$(MD_ARCH) props_$(MD_ARCH) animat_$(MD_ARCH)
	-cd $(MD_ARCH); rm *.o
	-cd $(MD_ARCH)_debug; rm *.o

flucq_objs:	$(FLUCQOBJS)
	$(FC) -o $(BIN) $(FLUCQOBJS) $(FLUCQLIBS)
	mv $(BIN) ..

props_objs:	$(PROPSOBJS)
	$(FC) -o $(BIN) $(PROPSOBJS) $(PROPSLIBS)
	mv $(BIN) ..

animat_objs:	$(ANIMATOBJS)
	$(FC) -o $(BIN) $(ANIMATOBJS) $(ANIMATLIBS)
	mv $(BIN) ..

# Be careful to preserve old executables in case anyone is running them.
# These datestamped files will need to be mopped up later by hand.
install:	flucq dflucq
	date=`date +%y%m%d`; \
	if [ ! -f $(MD_DIR)/flucq_$(MD_ARCH).$$date ] ; then \
	  mv $(MD_DIR)/flucq_$(MD_ARCH) $(MD_DIR)/flucq_$(MD_ARCH).$$date ; \
	else \
	  mv $(MD_DIR)/flucq_$(MD_ARCH) $(MD_DIR)/flucq_$(MD_ARCH).`date +%y%m%d%H%M%S` ; \
	fi
	$(INSTALL) -c $(MD_DIR) flucq_$(MD_ARCH)
	if [ ! -f $(MD_DIR)/dflucq_$(MD_ARCH).$$date ] ; then
	  mv $(MD_DIR)/dflucq_$(MD_ARCH) $(MD_DIR)/dflucq_$(MD_ARCH).$$date ; \
	else \
	  mv $(MD_DIR)/dflucq_$(MD_ARCH) $(MD_DIR)/dflucq_$(MD_ARCH).`date +%y%m%d%H%M%S` ; \
	fi
	$(INSTALL) -c $(MD_DIR) dflucq_$(MD_ARCH)

# Objects depend on include files
$(FLUCQOBJS):	$(FLUCQCOMS)
$(PROPSOBJS):	$(PROPSCOMS)
$(ANIMATOBJS):	$(ANIMATCOMS)

# What follows are the flags for each architecture.  The architectures are
# named (for convenience) according to the return of the `uname -s` command.

# Generic debug flags
FGFLAGS = -g

# THe R10K version
IRIX64_FLAGS = -vms_cc -64 -O2# -O3 is too aggressive
IRIX64_DEBUG_FLAGS = $(FGFLAGS)

# The non-R10K IRIX version
IRIX_FLAGS = -vms_cc -32 -O2 # -O3 is too aggressive (?) and -sopt breaks
IRIX_DEBUG_FLAGS = $(FGFLAGS)

# The Linux version
Linux_FLAGS = -Wall -O6
Linux_DEBUG_FLAGS = $(FGFLAGS)

# The Sun version
SunOS_FLAGS = -fast -fnonstd
SunOS_DEBUG_FLAGS = $(FGFLAGS)

# For the SGI people who have to put up with a defective make we
# do the following (ugly!)
#ARCH_FLAGS = $$($(MD_ARCH)_FLAGS)
ARCH_FLAGS=$(Linux_FLAGS)              # -eh
ARCH_DEBUG_FLAGS=$(Linux_DEBUG_FLAGS)  # -eh

# Here is where we actually build the executable
flucq_$(MD_ARCH): $(FLUCQSRCS) $(FLUCQCOMS)
	-mkdir $(MD_ARCH); cp makefile $(MD_ARCH)
	cd $(MD_ARCH); make flucq_objs FFLAGS="$(ARCH_FLAGS)" BIN=$@
	@echo "****Please look in the current directory for a new $@****"

ARCH_DEBUG_FLAGS = $$($(MD_ARCH)_DEBUG_FLAGS)
# How to build the debug version
dflucq_$(MD_ARCH): $(FLUCQSRCS) $(FLUCQCOMS)
	-mkdir $(MD_ARCH)_debug; cp makefile $(MD_ARCH)_debug
	cd $(MD_ARCH)_debug; make flucq_objs FFLAGS="$(ARCH_DEBUG_FLAGS)" BIN=$@
	@echo "****Please look in the current directory for a new $@****"

props_$(MD_ARCH):	$(PROPSSRCS) $(PROPSCOMS)
	-mkdir $(MD_ARCH); cp makefile $(MD_ARCH)
	cd $(MD_ARCH); make props_objs FFLAGS="$(ARCH_FLAGS)" BIN=$@
	@echo "****Please look in the current directory for a new $@****"

animat_$(MD_ARCH):	$(ANIMATSRCS) $(ANIMATCOMS)
	-mkdir $(MD_ARCH); cp makefile $(MD_ARCH)
	cd $(MD_ARCH); make animat_objs FFLAGS="$(ARCH_FLAGS)" BIN=$@
	@echo "****Please look in the current directory for a new $@****"
