# Define variables ----------------------------------------------------------------------
# FORTRAN_COMPILER = gfortran
OBJECTS = mvrnorm.o main.o random.o
GFFLAGS = -fbounds-check -w

mvrnorm.app: $(OBJECTS) README.PDF
	gfortran $(GFFLAGS) $(OBJECTS) -o mvrnorm.app

main.o: main.f90
	gfortran -c main.f90

mvrnorm.o: mvrnorm.f90
	gfortran -c mvrnorm.f90

random.o: random.f90
	gfortran -c random.f90

README.PDF: README.md
	pandoc README.md -s -o README.PDF 


# Shorthand for target object files
#%.o: %.f90	
#	gfortran -c $< 

# Delete files resulting from compiling -------------------------------------------------
.PHONY: clean, run	
clean:
	rm ./mvrnorm.app

run: 
	@echo 
	@echo Running mvrnorm.app;
	./mvrnorm.app;

