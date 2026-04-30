    const wordList regionNames(rp["fluid"]);

    PtrList<fvMesh> calRegions(regionNames.size());

    forAll(regionNames, i)
    {
        Info<< "Create mesh for region " << regionNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        calRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    regionNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
