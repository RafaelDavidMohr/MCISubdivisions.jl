# MCISubdivisions.jl

A julia package to compute tropical data associated to so-called
engineered complete intersections accompanying the paper TODO LINK.

# Installation

The package may be installed by adding it to your julia environment
directly from this git repository:

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/RafaelDavidMohr/MCISubdivisions.jl.git")
julia> using MCISubdivisions
```

# Usage

This package partially relies on the computer algebra system
[Oscar](https://www.oscar-system.org/).

An ECI $(V, A)$ is encoded by giving the underlying support set $A$ in
the form of a matrix of type `Matrix{Int}`, whose columns correspond
to the elements $A$, as well as the associated coefficient matrix of
type `Matrix{Int}` (encoding ECI's with rational coefficients) or
`Matrix{FqFieldElem}` (encoding ECI's with coefficients over a finite
field. Alternatively, an ECI can be encoded by directly giving a
polynomial system construced in Oscar in the form of a
`Vector{MPolyRingElem}`, in which case $A$ and $V$ are constructed
internally.

We refer to the file `src/MCISubdivisions.jl` and the docstrings of
the functions therein for details on the functionalities provided by
this package.
