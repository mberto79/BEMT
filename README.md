# BEMT.jl

## About this package

### Purpose

The `BEMT.jl` package provides a basic implementation of Blade Element Momentum Theory for the analysis of rotors and propellers. It has been designed to support the delivery of the Integrated Project (Level 2) as part of the Aerospace Engineering Programme at the University of Nottingham. Thus, the target users for this package are Aerospace students at Nottingham. The code has been released with an MIT license so anyone interested is welcome to use the package/code freely.

By design, the functionality provided is basic, but composable, to both encourage and allow users to extent/adopt the code to suit their needs.

The Julia programming language was chosen for implementing this package because it is a high-level programming language which is relatively straight forward to learn (especially for Matlab and Python users). Julia also offers excellent performance and it makes it almost trivial to share and deploy code.

### Reporting issues

University of Nottingham students should report any issues directly to Humberto Medina. Users outside the University are kindly asked to open an issue in Github.

Please include a minimal working example to illustrate the issue/bug when reporting or opening an issue.

### Feature requests and contributions

If you have a feature request, please open an issue explaining the feature, keeping in mind that this is not meant to be a production code. Ideally, feature request will include an educational argument.

Contributions from both University of Nottingham students or community users are, of course, welcome. Please submit contributions as a pull request.

Please note that updates to this package are likely to be sporadic (typically once per year around the summer months).

## Getting started

### Installing Julia

The first step is to make sure you install Julia locally. Simply visit the website for [Julia](https://julialang.org/) to download and install it. It is recommended that you choose the option to add Julia to the path during installation (by ticking the appropriate box in the installer).

### Using Julia

When Julia is installed, you will be able to find it in your list of installed programmes. When you launch Julia you will be presented with the Julia REPL (Read-Evaluate-Print Loop). The REPL is the simplest way to get started using Julia. However, as you start to develop more complex analysis or scripts, it is recommended to use a more complete IDE. Currently, **vscode** is the probably the leading IDE used in industry. You can learn more about setting up the vscode environment for Julia [on this link](https://code.visualstudio.com/docs/languages/julia). For information about how to use Julia and its syntax, you should refer to the [Julia documentation](https://docs.julialang.org/).

### Recommended packages

Once you have installed Julia, you can extent its functionality by adding packages from the Julia ecosystem. As a starting point, the following packages are highly recommended:

1. `Plots.jl` this package can be used to generate high quality plots and charts
2. `Revise.jl` this package enhances the development experience in Julia (reducing the need to reload the REPL when writing Julia code)

To add packages all you need to do is open a REPL window and press "]" on your keyboard (without the quotations). This will make the REPL enter *package mode*. To install a package simply type `add` followed by the package name (without the extension). For example, to install the package `Plot.jl`, you should type:

```julia
add Plots
```

The process above can be repeated to install (or add) any packages officially registered with the Julia package repository.

### Installing `BEMT.jl`

The `BEMT.jl` package is not meant to be a production package, thus, it has not been officially registered as a Julia package. To add `BEMT.jl` to your local Julia installation, you need to add it using the URL to this Github repository. To do this, enter *package mode* by pressing "]" in the Julia REPL and type:

```julia
add https://github.com/mberto79/BEMT.git
```


## Background

### Assumptions

### Blade element momentum theory

### Implementation limitations

## Basic usage example

