/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Application
    setFluenceRate

Description
    Writes a volScalarField G [W/m^2] = analytical radial fluence rate
    around an axial UV lamp, following the infinite-line-source model
    of Sozzi & Taghipour 2006 eq. (3):

        G(r) = P / (2 pi L_arc r) * exp(-sigma_w * (r - r_L))

    The lamp axis is assumed to lie along +x at (y, z) = (0, 0); r is
    the cylindrical radius sqrt(y^2 + z^2), clamped at r_L.

    The formula holds along the lamp's arc only. -xStart and -xEnd give
    the arc's axial extent; G is zero outside it, as in the radial model
    of Sozzi & Taghipour (2006), which has no lamp-end emission. Without
    them the formula is applied at every x, including any inlet pipe on
    the lamp axis, where r < r_L clamps G to its value at the sleeve.

    Defaults match the Sozzi 2006 25 GPM L-shape case (water, 70%
    transmissivity per cm, 35 W lamp, 80 cm arc).

    A timeSelector flag IS REQUIRED -- with no flags the utility
    defaults to writing the field into time 0/, which is rarely what
    you want (your converged flow lives at the latest time, not at
    t=0). Use -latestTime for the typical "drop G next to my converged
    flow snapshot" workflow, or -time <T> for a specific time
    directory.

Usage
    setFluenceRate -latestTime
    setFluenceRate -time 500 -P 35 -Larc 0.80 -rL 0.01 -sigmaW 35.67
    setFluenceRate -latestTime -xStart 0 -xEnd 0.80

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "mathematicalConstants.H"
#include "timeSelector.H"
#include "volFields.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption("P",      "scalar",
        "Total germicidal lamp power [W] (default 35)");
    argList::addOption("Larc",   "scalar",
        "Lamp arc length [m] (default 0.80)");
    argList::addOption("rL",     "scalar",
        "Lamp + sleeve outer radius [m] (default 0.01)");
    argList::addOption("sigmaW", "scalar",
        "Water absorption coefficient [1/m] (default 35.67, "
        "= -ln(0.7) per cm)");
    argList::addOption("xStart", "scalar",
        "Axial start of the lamp arc [m]; G = 0 for x below it "
        "(default: no limit)");
    argList::addOption("xEnd",   "scalar",
        "Axial end of the lamp arc [m]; G = 0 for x above it "
        "(default: no limit)");

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const scalar P      = args.optionLookupOrDefault<scalar>("P",      35.0);
    const scalar Larc   = args.optionLookupOrDefault<scalar>("Larc",   0.80);
    const scalar rL     = args.optionLookupOrDefault<scalar>("rL",     0.01);
    const scalar sigmaW = args.optionLookupOrDefault<scalar>("sigmaW", 35.67);
    const scalar xStart = args.optionLookupOrDefault<scalar>("xStart", -great);
    const scalar xEnd   = args.optionLookupOrDefault<scalar>("xEnd",    great);
    const scalar twoPi  = 2.0*constant::mathematical::pi;

    if (xEnd <= xStart)
    {
        FatalErrorInFunction
            << "-xEnd (" << xEnd << ") must exceed -xStart (" << xStart << ")"
            << exit(FatalError);
    }

    Info<< "setFluenceRate parameters:" << nl
        << "    P       = " << P      << " W" << nl
        << "    L_arc   = " << Larc   << " m" << nl
        << "    r_L     = " << rL     << " m" << nl
        << "    sigma_w = " << sigmaW << " 1/m" << nl;
    if (args.optionFound("xStart") || args.optionFound("xEnd"))
    {
        Info<< "    arc     = x in [" << xStart << ", " << xEnd << "] m"
            << " (G = 0 outside)" << nl;
    }
    Info<< endl;

    forAll(timeDirs, ti)
    {
        runTime.setTime(timeDirs[ti], ti);
        Info<< "Time = " << runTime.name() << endl;

        volScalarField G
        (
            IOobject
            (
                "G",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "G",
                dimensionSet(1, 0, -3, 0, 0, 0, 0),
                0.0
            )
        );

        auto Gat = [&](const scalar x, const scalar y, const scalar z) -> scalar
        {
            if (x < xStart || x > xEnd)
            {
                return 0;
            }
            const scalar r = std::max(std::sqrt(y*y + z*z), rL);
            return P/(twoPi*Larc*r)*std::exp(-sigmaW*(r - rL));
        };

        const volVectorField& C = mesh.C();
        forAll(G, i)
        {
            G[i] = Gat(C[i].x(), C[i].y(), C[i].z());
        }

        // The default 'calculated' BC has no evaluate() that fills from
        // the internal field, so face values stay at 0 unless we set
        // them explicitly. Vertex-averaged interpolators (e.g. the
        // radiationDose tracker) read these face values, so leaving
        // them at 0 biases interpolation low near every boundary.
        volScalarField::Boundary& Gb = G.boundaryFieldRef();
        forAll(Gb, patchi)
        {
            fvPatchScalarField& gp = Gb[patchi];
            const vectorField& cfp = mesh.boundary()[patchi].Cf();
            forAll(gp, facei)
            {
                gp[facei] = Gat(cfp[facei].x(), cfp[facei].y(), cfp[facei].z());
            }
        }

        Info<< "  G range: min = " << min(G).value()
            << " W/m^2, max = " << max(G).value() << " W/m^2" << nl
            << "  Writing G to " << runTime.name()/G.name() << endl;
        G.write();
    }

    return 0;
}


// ************************************************************************* //
