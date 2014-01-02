DEST= $(R_LIBS_USER)/recombination.analyzer
EXTERNAL_CODE= $(DEST)/ext
CWD= $$(pwd)

all:	
	echo "current directory is "$(CWD)
	echo "destination directory for external code is "$(EXTERNAL_CODE)
	mkdir -p $(EXTERNAL_CODE)
	@echo "*********************************************** "
	@echo "********* BUILDING 3seq *********************** "
	@echo "*********************************************** "	
	make 3seq		
	@echo "*********************************************** "
	@echo "********* BUILDING standard-RAxML ************* "
	@echo "*********************************************** "	
	make ExaML-raxml	
	@echo "*********************************************** "
	@echo "********* BUILDING ExaML-parsimonator ********* "
	@echo "*********************************************** "	
	make EXAML-parsimonator
	@echo "*********************************************** "
	@echo "********* BUILDING ExaML ********************** "
	@echo "*********************************************** "	
	make examl	
	@echo "*********************************************** "
	@echo "********* MAKE DOCUMENTATION ****************** "
	@echo "*********************************************** "	
	make roxygen
	
3seq:
	echo $(CWD)	
	unzip ../inst/3seq_v1.10812_or.zip
	make -C 3seq_v1.10812_or PTFLAGS=-DPVT='\"'$(DEST)'/\"'
	mv $(CWD)/3seq_v1.10812_or/3seq $(EXTERNAL_CODE)/3seq
	rm -rf 3seq_v1.10812_or
	
ExaML-raxml:
	echo $(CWD)	
	git clone https://github.com/stamatak/standard-RAxML.git
	make -C standard-RAxML -f Makefile.SSE3.gcc
	mv $(CWD)/standard-RAxML/raxmlHPC-SSE3 $(EXTERNAL_CODE)/ExaML-raxml
	rm -rf standard-RAxML	
	
EXAML-parsimonator:	
	echo $(CWD)
	git clone https://github.com/stamatak/Parsimonator-1.0.2.git
	make -C Parsimonator-1.0.2 -f Makefile.SSE3.gcc
	mv $(CWD)/Parsimonator-1.0.2/parsimonator-SSE3 $(EXTERNAL_CODE)/ExaML-parsimonator 	
	rm -rf Parsimonator-1.0.2
	
examl:	
	echo $(CWD)	
	git clone https://github.com/stamatak/ExaML.git
	make -C ExaML/parser -f Makefile.SSE3.gcc
	mv $(CWD)/ExaML/parser/parser $(EXTERNAL_CODE)/ExaML-parser
	make -C ExaML/examl -f Makefile.SSE3.gcc
	mv $(CWD)/ExaML/examl/examl $(EXTERNAL_CODE)/ExaML-examl
	rm -rf ExaML
	mv $(EXTERNAL_CODE)/ExaML-examl $(EXTERNAL_CODE)/examl
	
make roxygen:
	echo $(CWD)
	../misc/3seq.startme.R -exe=ROXYGENIZE