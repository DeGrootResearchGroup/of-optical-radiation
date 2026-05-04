/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
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

    Defaults match the Sozzi 2006 25 GPM L-shape case (water, 70%
    transmissivity per cm, 35 W lamp, 80 cm arc).

    Standard timeSelector flags (-time, -latestTime, ...) pick the
    time directory the field is written to.

Usage
    setFluenceRate -latestTime
    setFluenceRate -time 500 -P 35 -Larc 0.80 -rL 0.01 -sigmaW 35.67

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

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const scalar P      = args.optionLookupOrDefault<scalar>("P",      35.0);
    const scalar Larc   = args.optionLookupOrDefault<scalar>("Larc",   0.80);
    const scalar rL     = args.optionLookupOrDefault<scalar>("rL",     0.01);
    const scalar sigmaW = args.optionLookupOrDefault<scalar>("sigmaW", 35.67);
    const scalar twoPi  = 2.0*constant::mathematical::pi;

    Info<< "setFluenceRate parameters:" << nl
        << "    P       = " << P      << " W" << nl
        << "    L_arc   = " << Larc   << " m" << nl
        << "    r_L     = " << rL     << " m" << nl
        << "    sigma_w = " << sigmaW << " 1/m" << nl
        << endl;

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

        const volVectorField& C = mesh.C();
        forAll(G, i)
        {
            const scalar y = C[i].y();
            const scalar z = C[i].z();
            const scalar r = std::max(std::sqrt(y*y + z*z), rL);
            G[i] = P/(twoPi*Larc*r)*std::exp(-sigmaW*(r - rL));
        }
        G.correctBoundaryConditions();

        Info<< "  G range: min = " << min(G).value()
            << " W/m^2, max = " << max(G).value() << " W/m^2" << nl
            << "  Writing G to " << runTime.name()/G.name() << endl;
        G.write();
    }

    return 0;
}


// ************************************************************************* //
