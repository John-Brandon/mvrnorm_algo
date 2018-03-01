# Generate Multivariate Normal Random Variates

First iteration of the code was copied from (author John Burkardt):
http://people.sc.fsu.edu/~jburkardt/f_src/normal_dataset/normal_dataset.f90

As originally written, the code "creates a multivariate normal random dataset and writes it to a file". 

For a bit more description of the original code, see: 
http://people.sc.fsu.edu/~jburkardt/f_src/normal_dataset/normal_dataset.html

See also additional comments at the header of `program main`.

That code resides here in the `normal_dataset.f90` file. The `mvrnorm.f90` file takes
the random number procedures from the original code and wraps them in a module.
One of the purposes of this directory is to test the module of code, before
linking it into a larger application. 



