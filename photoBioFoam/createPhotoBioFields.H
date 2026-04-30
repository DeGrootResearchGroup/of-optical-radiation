    // Initialise scalar fields

    volScalarField I_cal
    (
        IOobject
        (
            "I",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    // Initialise the photoBio model
    
    autoPtr<photoBio::photoBioModel> photoBioModel
    (
        photoBio::photoBioModel::New(I_cal)
    );
