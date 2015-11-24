/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "turbulentPotential.H"
#include "addToRunTimeSelectionTable.H"
#include "backwardsCompatibilityWallFunctions.H"
#include "floatScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulentPotential, 0);
addToRunTimeSelectionTable(RASModel, turbulentPotential, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> turbulentPotential::Ts() const
{
    return max(k_/epsilon_, 6.0*sqrt(nu()/epsilon_));
}

tmp<volScalarField> turbulentPotential::TsEh() const
{
    return max(1.0/epsHat_, 6.0*sqrt(nu()/epsilon_));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentPotential::turbulentPotential
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),


    cEp1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEp1",
            coeffDict_,
            1.45
        )
    ),
    cEp2con_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEp2con",
            coeffDict_,
            1.83
        )
    ),
    cEp3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEp3",
            coeffDict_,
            0.15
        )
    ),
    cD1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cD1",
       	    coeffDict_,
            0.5
        )
    ),
    cD2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cD2",
       	    coeffDict_,
            1.88
        )
    ),
    cP1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP1",
            coeffDict_,
            2.0
        )
    ),
    cP2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP2",
            coeffDict_,
            0.6
        )
    ),
    cP3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP3",
            coeffDict_,
            0.12
        )
    ),
    cPphi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cPphi",
            coeffDict_,
            2.0
        )
    ),
    cMu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cMu",
            coeffDict_,
            0.21
        )
    ),

    cT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cT",
            coeffDict_,
            0.0033
        )
    ),

    cPr_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cPr",
            coeffDict_,
            1.0
        )
    ),

    cEhm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEhm",
            coeffDict_,
            10.0
        )
    ),

    cEhR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEhR",
            coeffDict_,
            1.0
        )
    ),

    sigmaKInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaKInit",
            coeffDict_,
            1.0
        )
    ),

    sigmaEpsInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEpsInit",
            coeffDict_,
            0.833
        )
    ),

    sigmaEpsVisc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEpsVisc",
            coeffDict_,
            1.0
        )
    ),

    sigmaPhiInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPhiInit",
            coeffDict_,
            0.33
        )
    ),

    sigmaPsiInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPsiInit",
            coeffDict_,
            1.0
        )
    ),

    nutScale_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "nutScale",
            coeffDict_,
            1.0
        )
    ),
    nutBlend_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "nutBlend",
            coeffDict_,
            1.0
        )
    ),
    psiNuFrac_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "psiNuFrac",
            coeffDict_,
            1.0
        )
    ),


   nutType_
   (
       coeffDict_.lookup("nutType")
   ),

   solveK_
   (
       coeffDict_.lookup("solveK")
   ),

   solveEps_
   (
       coeffDict_.lookup("solveEps")
   ),

   solvePsi_
   (
       coeffDict_.lookup("solvePsi")
   ),

   solvePhi_
   (
       coeffDict_.lookup("solvePhi")
   ),

   solveNut_
   (
       coeffDict_.lookup("solveNut")
   ),

   eqnSigmaK_
   (
       coeffDict_.lookup("eqnSigmaK")
   ),

   eqnSigmaEps_
   (
       coeffDict_.lookup("eqnSigmaEps")
   ),

   eqnSigmaPhi_
   (
       coeffDict_.lookup("eqnSigmaPhi")
   ),

   eqnSigmaPsi_
   (
       coeffDict_.lookup("eqnSigmaPsi")
   ),

   eqncEp2_
   (
       coeffDict_.lookup("eqncEp2")
   ),

   eqnEpsHat_
   (
       coeffDict_.lookup("eqnEpsHat")
   ),

   timeScaleEps_
   (
       coeffDict_.lookup("timeScaleEps")
   ),


    y_(mesh_),


    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    gradk_
    (
        IOobject
        (
            "gradk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(k_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    nutNorm_
    (
        IOobject
        (
            "nutNorm",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (nut_/max(nut_))
    ),
    tpphi_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    tpphisqrt_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt(tpphi_)
    ),
    vorticity_
    (
        IOobject
        (
            "vorticity",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(U_)
    ),
    vorticityTmp_
    (
        IOobject
        (
            "vorticityTmp",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::curl(U_)
    ),
    ivorticity_
    (
        IOobject
        (
            "ivorticity",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("iv", dimensionSet(0,0,1,0,0,0,0), vector(1.0,1.0,1.0))
    ),
    tppsi_
    (
        IOobject
        (
            "tppsi",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    uGrad_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(mag(U_))
    ),
	epsHat_
    (
        IOobject
        (
            "epsHat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsHat", dimensionSet(0,0,-1,0,0,0,0), 1.0)
    ),
    kol_
    (
        IOobject
        (
            "kol",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kol", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    ),
    kSafe_
    (
        IOobject
        (
            "kSafe",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (k_)
    ),
    kSqrt_
    (
        IOobject
        (
            "kSqrt",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    gradkSqrt_
    (
        IOobject
        (
            "gradkSqrt",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(kSqrt_)
    ),
    nutSafe_
    (
        IOobject
        (
            "nutSafe",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (nut_)
    ),
    epsilonSafe_
    (
        IOobject
        (
            "epsilonSafe",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (epsilon_)
    ),
    sigmaK_
    (
        IOobject
        (
            "sigmaK",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaKInit_
    ),
    sigmaEps_
    (
        IOobject
        (
            "sigmaEps",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaEpsInit_
    ),
    sigmaPhi_
    (
        IOobject
        (
            "sigmaPhi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaPhiInit_
    ),
    sigmaPsi_
    (
        IOobject
        (
            "sigmaPsi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaPsiInit_
    ),
    cEp2_
    (
        IOobject
        (
            "cEp2",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (cEp2con_ - 0.16*exp(-0.1*sqr(k_)/(nu()*epsilon_)))
    ),
    tpProd_
    (
        IOobject
        (
            "tpProd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cPr_*(tppsi_ & vorticity_)
    ),
    cP1eqn_
    (
        IOobject
        (
            "cP1eqn",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        2.0*(0.5+0.5*((tpProd_*k_)/epsilon_))
    ),
    dimRat_
    (
        IOobject
        (
            "dimRat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (psiReal() & psiReal())/(k_*phiReal())
    ),
    graddimRat_
    (
        IOobject
        (
            "graddimRat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(dimRat_)
    ),
    tpProdSqr_
    (
        IOobject
        (
            "tpProdSqr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(tppsi_ & vorticity_)
    ),
    tpProd3d_
    (
        IOobject
        (
            "tpProd3d",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (mag(psiReal() ^ vorticity_))
    )
{

    Info<< "Made it past constructors " << endl;

    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
		if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
		{
            nut_ = cMu_*k_*tpphi_*Ts();
            nut_.correctBoundaryConditions();
        }

        if(timeScaleEps_ == "epsHat")
		{
            nut_ = cMu_*k_*tpphi_*TsEh();
            nut_.correctBoundaryConditions();
        }

    }

    //Info<< "Made it past nut" << endl;

    kSafe_ = max(k_, dimensionedScalar("minK", k_.dimensions(), 1.0e-15));

    Info<< "solveK is: " <<solveK_ <<endl;
    Info<< "solveEps is: " <<solveEps_ <<endl;
    Info<< "solvePhi is: " <<solvePhi_ <<endl;
    Info<< "solvePsi is: " <<solvePsi_ <<endl;

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotential::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_))
        )
    );
}

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotential::devReff() const
{
    //Info<< "Using devReff" << endl;

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> turbulentPotential::divDevReff(volVectorField& U) const
{
    return
    (
       fvc::grad(phiReal())
     + fvc::curl(psiReal())
     + fvc::laplacian(nut_, U, "laplacian(nuEff,U)")
     - fvm::laplacian(nuEff(), U)
    );
}


tmp<fvVectorMatrix> turbulentPotential::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool turbulentPotential::read()
{
    if (RASModel::read())
    {
        Info<< "Updating Turbulence Coefficients" << endl;
        cEp1_.readIfPresent(coeffDict());
        cEp2con_.readIfPresent(coeffDict());
        cEp3_.readIfPresent(coeffDict());
        cP1_.readIfPresent(coeffDict());
        cP2_.readIfPresent(coeffDict());
        cP3_.readIfPresent(coeffDict());
        cMu_.readIfPresent(coeffDict());
		cPphi_.readIfPresent(coeffDict());
		cEhm_.readIfPresent(coeffDict());
		cEhR_.readIfPresent(coeffDict());
		cPr_.readIfPresent(coeffDict());
		cD1_.readIfPresent(coeffDict());
		cD2_.readIfPresent(coeffDict());
		cT_.readIfPresent(coeffDict());
		sigmaKInit_.readIfPresent(coeffDict());
        sigmaEpsInit_.readIfPresent(coeffDict());
        sigmaEpsVisc_.readIfPresent(coeffDict());
        sigmaPhiInit_.readIfPresent(coeffDict());
        sigmaPsiInit_.readIfPresent(coeffDict());
		nutBlend_.readIfPresent(coeffDict());
		nutScale_.readIfPresent(coeffDict());


        return true;
    }
    else
    {
        Info<< "NOT Updating Turbulence Coefficients" << endl;
        return false;
    }
}


void turbulentPotential::correct()
{

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
        bound(k_,dimensionedScalar("minK", k_.dimensions(), 1.0e-15));
        bound(epsilon_,dimensionedScalar("minEps", epsilon_.dimensions(), 1.0e-15));
    }

    // Vorticity
    vorticity_ = fvc::curl(U_);

    // Production
    tpProd_ = cPr_*(tppsi_ & vorticity_);
    tpProdSqr_ = sqr(tpProd_);
    tpProd3d_ = mag(psiReal() ^ vorticity_);

    // Epsilon-hat and Sigma Equations

        if(eqnEpsHat_ == "mod")
	{
            epsHat_ = (epsilon_)/(k_ + cEhm_*nu()*mag(gradkSqrt_));
            bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), 1.0e-15));
	}
	else if(eqnEpsHat_ == "dif")
	{
            epsHat_ = (epsilon_ - 2.0*nu()*sqr(mag(gradkSqrt_)))/k_;
            bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), 1.0e-15));
        }
	else if(eqnEpsHat_ == "rough")
	{
            epsHat_ = (epsilon_ - cEhR_*nu()*sqr(mag(gradkSqrt_)))/k_;
            bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), 1.0e-15));
	}
	else
	{
            Info<< "No EpsHat Model Chosen" <<endl;
	    epsHat_ = (epsilon_)/(k_ + cEhm_*nu()*mag(gradkSqrt_));
	    bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), 1.0e-15));
	}

    if(eqnSigmaK_ == "true")
    {
	    sigmaK_ = 0.67 + 0.33*(tpProd_/epsHat_);
	}

    if(eqnSigmaEps_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.5*(tpProd_/epsHat_);
    }

    if(eqnSigmaPhi_ == "true")
    {
	    sigmaPhi_ = 0.21 + 0.12*(tpProd_/epsHat_);
	}

    if(eqnSigmaPsi_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.4*(tpProd_/epsHat_);
    }

    epsilonSafe_ = epsilon_ + dimensionedScalar("minEps", epsilon_.dimensions(), ROOTVSMALL);

    if(eqncEp2_ == "true")
    {
        cEp2_ = cEp2con_ - 0.16*exp(-0.25*sqr(k_)/(nu()*epsilonSafe_));
    }
    else
    {
	    cEp2_ =  cEp2con_;
    }

    cP1eqn_ = 2.0*(0.5+0.5*((tpProd_*k_)/epsilonSafe_));


    //Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
       epsHat_*cEp1_*tpProd_*k_
     + fvm::Sp(-epsHat_*cEp2_,epsilon_)
     + epsHat_*cEp3_*tpProd3d_
    );

    if(solveEps_ == "true")
    {
    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_,dimensionedScalar("minEps", epsilon_.dimensions(), 1.0e-15));
    }


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (

        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        tpProd_*k_
      - fvm::Sp(epsilon_/k_,k_)
    );


    if(solveK_ == "true")
    {
    kEqn().relax();
    solve(kEqn);
    bound(k_,dimensionedScalar("minK", k_.dimensions(), 1.0e-15));
    }

    kSafe_ = max(k_, dimensionedScalar("minK", k_.dimensions(), 1.0e-15));
    kSqrt_ = sqrt(k_);
    kSqrt_.correctBoundaryConditions();

    gradk_ = fvc::grad(k_);
    gradkSqrt_ = fvc::grad(kSqrt_);


    // Phi/K equation
    tmp<fvScalarMatrix> tpphiEqn
    (
        fvm::ddt(tpphi_)
      + fvm::div(phi_, tpphi_)
      + fvm::SuSp(-fvc::div(phi_), tpphi_)
      - fvm::laplacian(DphiEff(), tpphi_)
      ==
        2.0*nutFrac()*epsHat_*(2.0*Alpha()-1.0)*tpphi_
      + fvm::Sp(-1.0*(1.0-cD2_)*tpProd_,tpphi_)
      + fvm::Sp(-1.0*Alpha()*(1.88*Alpha() - 1.0)*(epsilon_/kSafe_),tpphi_)
      + cT_*tpProd_*sqrt((nut_/nu()))
    );

    if(solvePhi_ == "true")
    {
    tpphiEqn().relax();
    solve(tpphiEqn);
    bound(tpphi_,dimensionedScalar("minTpphi", tpphi_.dimensions(), 1.0e-15));
    }


    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
		if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
		{
            nut_ = cMu_*k_*tpphi_*Ts();
            nut_.correctBoundaryConditions();
        }

        if(timeScaleEps_ == "epsHat")
		{
            nut_ = cMu_*k_*tpphi_*TsEh();
            nut_.correctBoundaryConditions();
        }

    }


    // Psi Equation
    tmp<fvVectorMatrix> tppsiEqn
    (
        fvm::ddt(tppsi_)
      + fvm::div(phi_, tppsi_)
      + fvm::Sp(-fvc::div(phi_), tppsi_)
      - fvm::laplacian(DpsiEff(), tppsi_)

      ==

        0.21*(2.0*Alpha()-1.0)*tpphi_*vorticity_
      + fvm::Sp(-1.0*(1.0 - cP2_)*tpProd_, tppsi_)
      + (1.0 - cP2_)*tpphi_*vorticity_
      + fvm::Sp(-cD1_*Alpha()*tpProd_,tppsi_)
      + fvm::Sp(-1.0*(cP1_*nutFrac()*(1.0-Alpha()))*epsHat_,tppsi_)
      + fvm::Sp(-0.12*Alpha()*(epsilon_/kSafe_),tppsi_)
      + cT_*sqrt((nut_/nu()))*vorticity_
    );

    if(solvePsi_ == "true")
    {
    tppsiEqn().relax();
    solve(tppsiEqn);
    }




}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
