# TimeInterpolatedGeoData

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cesaraustralia.github.io/TimeInterpolatedGeoData.jl/dev/)
[![Build Status](https://travis-ci.com/cesaraustralia/TimeInterpolatedGeoData.jl.svg?branch=master)](https://travis-ci.com/cesaraustralia/TimeInterpolatedGeoData.jl)
[![codecov.io](http://codecov.io/github/cesaraustralia/TimeInterpolatedGeoData.jl/coverage.svg?branch=master)](http://codecov.io/github/cesaraustralia/TimeInterpolatedGeoData.jl?branch=master)


Lazily interpolates spatial data in the time dimension. 

This can be useful when a model needs daily data but layers for a particular dataset are monthly.
It can also be useful to interpolate maximum and minimum values (such as
temperature) to hourly temperature values.

This package provides tools for doing both of these things.


To install:

```julia
] 
add https://github.com/cesaraustralia/SimpleSDMDataSources.jl
add https://github.com/cesaraustralia/TimeInterpolatedGeoData.jl
```
