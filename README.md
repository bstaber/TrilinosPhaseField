## TrilinosPhaseField

### Phase field for brittle fracture implemented with the Trilinos Project.

### In progress.

* Three-dimensional
* Dirichlet conditions on displacement only
* No Dirichlet conditions on phase field
* Staggered algorithm
* Anisotropic model of Miehe
* Homogeneous, heterogeneous or random heterogeneous media

### Packages and libraries

* Built on top of Sandia's Trilinos Project [[Website](http://trilinos.org/) |
[Documentation](http://trilinos.org/about/documentation/) |
[Mailing List](https://trilinos.org/mailman/listinfo/trilinos-users) |
[Packages](http://trilinos.org/packages/) |
[GitHub](https://github.com/trilinos/Trilinos)] and uses, for instance, the packages:
	* [Epetra](https://trilinos.org/packages/epetra/) for linear algebra,
	* [Amesos](https://trilinos.org/packages/amesos/), [AztecOO](https://trilinos.org/packages/aztecoo/), [Stratimikos](https://trilinos.org/packages/stratimikos/) for direct and iterative linear solvers,
	* [Teuchos](https://trilinos.org/packages/teuchos/) for parsing XML lists of parameters,
	* [ML](https://trilinos.org/packages/ml/), [IfPack](https://trilinos.org/packages/ifpack/) for preconditioners.

* Additional required libraries: [GMSH](http://gmsh.info/), [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview), [Boost](https://www.boost.org/).
