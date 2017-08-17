
## Initial position generator

Q: Can the cells loll out of the domain?
A: Yes, cells are bound inside the domain by their centerpoints, so parts of the cell might reach over the domain. These cells are deleted during initialisation.

## Cell positioning

Q: What is the exact centerpoint of the domain, how big is my domain?
A: When you specify a domain of for example 10 cells, the exact
middle will be *BETWEEN* cell 5 and 6. imagine that dx is 0.5µm then the domain
is 9.5µm long, and the middle is 4.75µm. periodic boundaries of course add 1 dx
back to the length.
