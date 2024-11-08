# FC = gfortran
FC = ifort

src_main = main.f90
src_mod1 = mod_parameter.f90  

obj_main = $(src_main:.f90=.o)
obj_mod1 = $(src_mod1:.f90=.o)

obj = $(notdir $(shell pwd))

$(obj):$(obj_main) $(obj_mod1)
	$(FC) -qopenmp -o $(obj) $(obj_main) $(obj_mod1)
	
$(obj_main):$(src_main) $(obj_mod1)
	$(FC) -qopenmp -c $(src_main)
	
$(obj_mod1):$(src_mod1)
	$(FC) -c $(src_mod1)
	
	
.PHONY: clean
clean:
	rm *.o
	rm *.mod
