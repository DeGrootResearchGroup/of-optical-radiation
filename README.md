# README #

This is an implementation of the radiation library which excludes heat radiation and solves for multi-band, per-direction radiance fields - the idea is to have the same functionality as Ansys Fluent.

### How do I get set up? ###

1. This library can be compiled with wmake in OpenFoam. To keep things segregated from the main development, please change your OpenFoam bashrc file to include an appropriate directory to have your $FOAM_USER_SRC, then place the optical directory within this one. You will need to source your ~/.bashrc file once you modify it.

2. The library can be compiled with the command;

```
wmake libso

```

3. The most basic implementation can be seen in the "opticalFoam" directory. This includes the optical model in a SIMPLE implementation. 

4. Compile the solver by placing it in an appropriate directory (e.g. $FOAM_USER_APP -> <name>-<versionNumber>/applications) and calling
```
wmake

```
Any place can be defined for the destination of the binaries. It should be put in $FOAM_USER_APPBIN.

5. The solver, with all its functionality can be run

### Contribution guidelines ###