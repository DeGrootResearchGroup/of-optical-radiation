    // Initialise field pointer lists

    PtrList<volScalarField> I_cal(calRegions.size());
    PtrList<photoBio::photoBioModel> photoBioModels(calRegions.size());

    // Populate fluid field pointer lists
    forAll(calRegions, i)
    {
        Info<< "*** Reading  properties for region "
            << calRegions[i].name() << nl << endl;

        I_cal.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "I",
                    runTime.timeName(),
                    calRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                calRegions[i]
            )
        );

        photoBioModels.set
        (
            i,
            photoBio::photoBioModel::New(I_cal[i])
        );
    }
    
    
    
