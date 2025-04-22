/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "myWallHeatFlux.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(myWallHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, myWallHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::myWallHeatFlux::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "average");
    writeTabbed(os, "integral");
    os  << endl;
}


void Foam::functionObjects::myWallHeatFlux::calcHeatFlux
(
    const scalar alpha,
    const volScalarField T,
    volScalarField& myWallHeatFlux
)
{
    volScalarField::Boundary& myWallHeatFluxBf = myWallHeatFlux.boundaryFieldRef();

    const volScalarField::Boundary TBf = T.boundaryField();


    for (const label patchi : patchSet_)
    {
        myWallHeatFluxBf[patchi] = alpha*(TBf[patchi].snGrad());
    }


    // const auto* qrPtr = cfindObject<volScalarField>(qrName_);

    // if (qrPtr)
    // {
    //     const volScalarField::Boundary& radHeatFluxBf = qrPtr->boundaryField();

    //     for (const label patchi : patchSet_)
    //     {
    //         myWallHeatFluxBf[patchi] -= radHeatFluxBf[patchi];
    //     }
    // }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::myWallHeatFlux::myWallHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    patchSet_(),
    qrName_("qr")
{
    volScalarField* myWallHeatFluxPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName(typeName),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            // dimensionedScalar(dimMass/pow3(dimTime), Zero)
            // Dimensions: Heat flux/(rho*cp)
            dimensionedScalar(dimLength*dimTemperature/dimTime, Zero)
        )
    );

    mesh_.objectRegistry::store(myWallHeatFluxPtr);

    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::myWallHeatFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.getOrDefault<wordRes>("patches", wordRes())
        );

    dict.readIfPresent("qr", qrName_);

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::myWallHeatFlux::execute()
{
    auto& myWallHeatFlux = lookupObjectRef<volScalarField>(scopedName(typeName));

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh_.time().constant(),
            mesh_,
            // runTime.constant(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    const scalar alpha (transportProperties.getScalar("DT"));
    // const dimensionedScalar DT("DT", dimViscosity, transportProperties);

    const volScalarField& T = lookupObject<volScalarField>("T");

    calcHeatFlux(alpha, T, myWallHeatFlux);

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf = mesh_.magSf().boundaryField();

    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = myWallHeatFlux.boundaryField()[patchi];

        const scalar areaTot = gSum(magSf[patchi]);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);
        const scalar avgHfp = integralHfp/areaTot;

        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << avgHfp
                << token::TAB << integralHfp
                << endl;
        }

        Log << "    avg/integ(" << pp.name() << ") = "
            << avgHfp << ", "  << integralHfp << endl;

        this->setResult("avg(" + pp.name() + ")", avgHfp);
        this->setResult("int(" + pp.name() + ")", integralHfp);
    }


    return true;
}


bool Foam::functionObjects::myWallHeatFlux::write()
{
    const auto& myWallHeatFlux =
        lookupObject<volScalarField>(scopedName(typeName));

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << myWallHeatFlux.name() << endl;

    myWallHeatFlux.write();

    return true;
}


// ************************************************************************* //
