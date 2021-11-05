# SimpleTreeMeshes

This code is intended to provide a quick and easy way to experiment
with algorithms on non-uniform "staggered grids". It does not provide
the performance or generality of the underlying [p4est](https::p4est.org) library (though relying on this
might allow for improvements in these directions in the future). Instead,
it aims to give the user quick access to data structures to quickly
define their algorithmic ideas, and plot them.

This is a research/experimental code, so one should follow the standard
practical assumption about such code: if it's not in the test suite, don't
believe it works!


```@autodocs
Modules = [SimpleTreeMeshes]
Order   = [:function, :type]
```
