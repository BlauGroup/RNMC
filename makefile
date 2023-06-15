GMC: 
	+$(MAKE) -C GMC

NPMC: 
	+$(MAKE) -C NPMC

LGMC: 
	+$(MAKE) -C LGMC

clean_LGMC:
	+$(MAKE) -C LGMC clean

clean-all:
	@for dir in LGMC NPMC GMC; do \
		$(MAKE) -C $$dir clean; \
	done

debug-all:
	@for dir in LGMC NPMC GMC; do \
		$(MAKE) -C $$dir debug; \
	done

make-all:
	@for dir in LGMC NPMC GMC; do \
		$(MAKE) -C $$dir; \
	done 

clean:
	find . -name "*.o" -type f -delete