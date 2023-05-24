# Clean all subdirectories
clean-all:
	@for dir in */; do \
		$(MAKE) -C $$dir clean; \
	done

GMC: 
	$(MAKE) -C GMC clean

(%_clean):
	+$(MAKE) -C LGMC clean


# Build all executables
all:
	+$(MAKE) -C LGMC