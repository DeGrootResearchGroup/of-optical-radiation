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
    
    // Initialise the opticalRadiation model
    
    autoPtr<optical::radiationModel> radiationModel
    (
        optical::radiationModel::New(I_cal)
    );
