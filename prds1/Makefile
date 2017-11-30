# Makefile created by mkmf $Id: mkmf,v 10.0 2003/10/27 23:31:04 arl Exp $ 

include ../install/f90_args


.DEFAULT:
	-touch $@
all: prds1.x 
compute_prds.o: compute_prds.f90 commvars.o
	$(FC) $(FFLAGS) -c	compute_prds.f90
gen_ui_uniform.o: gen_ui_uniform.f90 commvars.o
	$(FC) $(FFLAGS) -c	gen_ui_uniform.f90
gen_ui_lebedev.o: gen_ui_lebedev.f90 commvars.o
	$(FC) $(FFLAGS) -c	gen_ui_lebedev.f90
commvars.o: commvars.f90
	$(FC) $(FFLAGS) -c	commvars.f90
sphere_lebedev_rule.o: sphere_lebedev_rule.f90 commvars.o
	$(FC) $(FFLAGS) -c	sphere_lebedev_rule.f90
utilities.o: utilities.f90
	$(FC) $(FFLAGS) -c	utilities.f90
hyp_normal_vect.o: hyp_normal_vect.f90 commvars.o
	$(FC) $(FFLAGS) -c	hyp_normal_vect.f90
gen_ui_test.o: gen_ui_test.f90 commvars.o
	$(FC) $(FFLAGS) -c	gen_ui_test.f90
setup_rho.o: setup_rho.f90 commvars.o
	$(FC) $(FFLAGS) -c	setup_rho.f90
principle_radius.o: principle_radius.f90 commvars.o
	$(FC) $(FFLAGS) -c	principle_radius.f90
compute_prds_half.o: compute_prds_half.f90 commvars.o
	$(FC) $(FFLAGS) -c	compute_prds_half.f90
SRC = compute_prds.f90 gen_ui_uniform.f90 gen_ui_lebedev.f90 commvars.f90 sphere_lebedev_rule.f90 utilities.f90 hyp_normal_vect.f90 gen_ui_test.f90 setup_rho.f90 principle_radius.f90 compute_prds_half.f90
OBJ = compute_prds.o gen_ui_uniform.o gen_ui_lebedev.o commvars.o sphere_lebedev_rule.o utilities.o hyp_normal_vect.o gen_ui_test.o setup_rho.o principle_radius.o compute_prds_half.o
clean: neat
	-rm -f .cppdefs *.mod $(OBJ) prds1.x 
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
prds1.x: $(OBJ)
	$(LD) $(OBJ) -o prds1.x $(LDFLAGS)
