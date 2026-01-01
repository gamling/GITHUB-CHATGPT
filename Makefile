
checkfiles= nca_U_bosons.f READIN.f utilities.f GINDERPOL.f selfenergies.f nullstelle.f\
	    nca_utilities.f suszi.f bubble.f selfenergies_2.f
obj= nca_U_bosons.o READIN.o utilities.o GINDERPOL.o selfenergies.o nullstelle.o\
	    nca_utilities.o suszi.o bubble.o selfenergies_2.o
	    FC = ifx -c -O3
	    linkF= ifx
	    librF= -lmkl_rt
.f.o:
	${FC} $<

nca:	$(obj)
	$(linkF) $@ $(obj) $(librF)
