# variableExtinctionBox3D

Same geometry, BCs, and angular discretisation as
`absorbingScatteringBox3D`, but with `wideBandVariableExtinction`
instead of constant extinction. The total absorption and scattering
coefficients are computed at each cell from species concentration
fields:

    A_band = sum_i a_i,band * X_i      (absorbing species X1, X2)
    S_band = sum_i s_i,band * S_i      (scattering species S1, S2)

This is the form intended for use inside a coupled solver that
transports species concentrations alongside the radiation. It also
works standalone with `opticalRadiationFoam`: the model auto-loads any
species fields not already registered with the mesh.

With X1 = X2 = S1 = S2 = 0.5 uniform (the values shipped here) and
specific coefficients of 0.1 m^2/kg per band, the net coefficients
are A = 0.1 m^-1 and S = 0.1 m^-1 per band — same as
`absorbingScatteringBox3D`. Comparing the two cases is a useful
sanity check: the irradiance field `G` should be bit-for-bit
identical (Alltest enforces this).

## Files in `0.orig/`

- `I` — radiance [W/m^2/sr] (always required by opticalRadiationFoam).
- `X1`, `X2` — absorbing species concentrations [kg/m^3].
- `S1`, `S2` — scattering species concentrations [kg/m^3].

## Running

    ./Allrun        # mesh + solve
    paraFoam        # visualise

## Cleanup

    ./Allclean
