# Wildfire 2D

> **Archive notice:** This repository archives an old project from 2017. It is kept here for reference and is not actively maintained.
>
> 📺 **Explanation video:** [https://www.youtube.com/watch?v=KpPFPfFV0Sk](https://www.youtube.com/watch?v=KpPFPfFV0Sk)

A MATLAB simulation of wildfire spread across a 2D landscape. The model couples a **temperature** field with a **fuel** field and evolves them forward in time on a discrete mesh, producing an animation of the fire front as it spreads, gets pushed by wind, and burns through the available fuel.

## The model

The simulation solves a reaction–diffusion–advection system on the unit square:

```
∂u/∂t = k₀ ∇²u  −  v · ∇u  +  f(u, S)          (temperature)
∂S/∂t = −Cₛ · S · exp(−B / (u − Tₐ))            (fuel, only where u > Tₐ)
```

where

```
f(u, S) = A · ( S · exp(−B / (u − Tₐ))  −  C · (u − Tₐ) )
```

| Term | Meaning |
|------|---------|
| `k₀ ∇²u` | **Diffusion**: heat spreads to neighboring cells (thermal diffusivity `k₀`). |
| `v · ∇u` | **Advection**: wind with direction `θ` and magnitude `m` carries heat downwind. |
| `S · exp(−B / (u − Tₐ))` | **Combustion**: an Arrhenius-type burn rate. Heat is released only where there is fuel and the temperature is above the ignition temperature `Tₐ`. |
| `C · (u − Tₐ)` | **Convective cooling**: heat lost to the surrounding air. |
| `∂S/∂t` | **Fuel consumption**: fuel burns away at a rate set by `Cₛ` wherever the fire is active. |

## Numerical method

- **Spatial discretization:** an `n × n` grid (default `n = 100`) using a 5-point, second-order stencil for the Laplacian and centered differences for the wind (gradient) term. Both are assembled as sparse matrices.
- **Operator splitting (Godunov):** each time step is split into two parts:
  1. **Reaction step:** the nonlinear combustion/cooling term is applied explicitly, only in cells above the ignition temperature.
  2. **Diffusion + advection step:** the linear part is advanced with **backward Euler** (implicit), solving `(I − k(D − D₂)) uⁿ⁺¹ = u*`.
- **Prefactorization:** the implicit system matrix doesn't change between steps, so it is LU-factored once up front. Each step then only needs two cheap triangular solves.
- **Fuel update:** fuel is decreased semi-implicitly in burning cells and left unchanged elsewhere.

## Initial conditions

- **Spot fires:** `FireInitializer.m` creates the starting temperature field by placing circular hot spots (temperature `IC`) at given `[x, y, radius]` locations. `revamp2.m` starts six of them around the domain.
- **Random fuel:** the fuel field is randomized so that most cells have little fuel (0.1) and about 30% have a random, larger amount. This gives a patchy, uneven burn.
- **Fuel break (optional):** commented-out code in `revamp2.m` clears a band of fuel to show how a firebreak stops the spread.

## Running it

Requires MATLAB. From the project directory:

```matlab
revamp2
```

The script will:

1. Run for `T = 400` time steps and draw a filled contour plot of temperature at each step.
2. Write the animation to `revamped_fire.avi`.
3. Show figures of the **initial** and **final** fuel distributions, so you can see what burned.

Parameters such as grid size, wind direction and strength, diffusivity, burn rates, and spot-fire locations are all defined at the top of `revamp2.m` and can be changed.

## Files

| File | Description |
|------|-------------|
| `revamp2.m` | Main script: parameters, matrix assembly, time stepping, plotting, and video output. |
| `FireInitializer.m` | Helper that builds the initial temperature field from a list of circular spot fires. |
| `license.txt` | BSD license. |

## License

BSD license. See [`license.txt`](license.txt). © 2017 Cooper Maira.
