<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_CHTML-full"></script>
# A Mini-App for teaching Accelerate

This is a simple application which solves a reaction-diffusion equation in
two-dimensions, using some common computational kernel to teach programming in
Accelerate.

The equation is [Fisher's equation](https://en.wikipedia.org/wiki/Fisher%27s_equation),
a simple non-linear reaction diffusion equation that is used to model simple
population dynamics.

\\[ \frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial^2 x} + R u (1 - u) \\]

Where the diffusion parameter \\(D\\) and growth rate \\(R\\) will be hard coded
for this example to give stable and interesting results.

The C implementation taken from: [bcummings/summer-school](https://github.com/bcumming/summer-school)


## Implementation notes

A 5-point Laplacian stencil is used for the spatial derivative. An implicit
temporal discritisation is employed, requiring the solution of a non-linear
system of equations at each time step.

## How to build

The Accelerate program can be built with the Haskell [`stack`](http://haskellstack.org) tool:

```sh
$ stack build
$ stack exec fisher-accelerate
```

## Visualising the results

The application outputs the final solution in the "brick of values" format,
which is stored in the two files __output.bin__ and __output.bov__. These can be
viewed using the popular visualization packages Paraview and Visit.

