# The makefile calls for another makefile within Source_files directory
# which compiles and makes an executable TREKIS.x
# This file was written by N.Medvedev
# in 2023-2024
#-----------------------------------------------------

# To pass variables into the next make-file:
export

# Call makefile within the Source_files directory:
subsystem:
	cd Source_files && $(MAKE)

ifeq ($(mpi),y)
# Copy created executable into the parent directory:
	scp -r Source_files/TREKIS_MPI.x TREKIS_MPI.x
# Delete executable from the Source_files directory:
	rm -r Source_files/TREKIS_MPI.x
else
# Copy created executable into the parent directory:
	scp -r Source_files/TREKIS.x TREKIS.x
# Delete executable from the Source_files directory:
	rm -r Source_files/TREKIS.x
endif


# Clean all compiled files:
clean:
	rm -f Source_files/*.o
	rm -f Source_files/*.mod
	rm -f Source_files/*.obj
	rm -f Source_files/*.yaml
	rm -f Source_files/*.optrpt
	rm -f Source_files/*.x

	rm -f *.x

# Clean all results:
cleanresults:
	rm -f OUTPUT_*

# Clean all compiled files and the results:
veryclean:
	rm -f Source_files/*.o
	rm -f Source_files/*.mod
	rm -f Source_files/*.obj
	rm -f Source_files/*.yaml
	rm -f Source_files/*.optrpt
	rm -f Source_files/*.x

	rm -f *.x
	rm -f OUTPUT_*
