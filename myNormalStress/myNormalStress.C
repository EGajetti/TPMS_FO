/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "myNormalStress.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(myNormalStress, 0);
    addToRunTimeSelectionTable(functionObject, myNormalStress, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::myNormalStress::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Wall shear stress");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "average");
    writeTabbed(os, "integral");
    os  << endl;
}


void Foam::functionObjects::myNormalStress::calcShearStress
(
    const scalar nu,
    const volVectorField U,
    volScalarField& myShearStress
)
{
    volScalarField::Boundary& myShearStressBf = myShearStress.boundaryFieldRef();
    const volVectorField::Boundary UBf = U.boundaryField();

    for (const label patchi : patchSet_)
    {
        // const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        // const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        myShearStressBf[patchi] = nu*mag(UBf[patchi].snGrad());
        
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::myNormalStress::myNormalStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    patchSet_()
{
    read(dict);

    writeFileHeader(file());

    volScalarField* myNormalStressPtr
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
            dimensionedScalar(sqr(dimLength)/sqr(dimTime), Zero)
        )
    );

    mesh_.objectRegistry::store(myNormalStressPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::myNormalStress::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.getOrDefault<wordRes>("patches", wordRes())
        );

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
                    << "Requested wall shear stress on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::myNormalStress::execute()
{
    auto& myNormalStress =
        mesh_.lookupObjectRef<volScalarField>(scopedName(typeName));
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
    const scalar nu (transportProperties.getScalar("nu"));
    const volVectorField& U = lookupObject<volVectorField>("U");
    calcShearStress(nu, U, myNormalStress);
    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf = mesh_.magSf().boundaryField();
    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const polyPatch& patchTPMS = mesh_.boundaryMesh()[patchi];
        const scalarField& taufp = myNormalStress.boundaryField()[patchi];

        scalar sum_A_taui = 0.0;  
        scalar sum_A = 0.0;        
        forAll(patchTPMS, facei)
        {
            scalar Ai = mesh_.magSf().boundaryField()[patchi][facei];
            scalar taui = taufp[facei]; 
            sum_A_taui += Ai * taui;
            sum_A += Ai;
        }
        scalar tau = sum_A_taui / (sum_A);

        const scalar areaTot = gSum(magSf[patchi]);
        const scalar integralTaufp = gSum(magSf[patchi]*taufp);
        const scalar avgTaufp = gAverage(taufp);

        if (Pstream::master())
            {
                writeCurrentTime(file());

                file()
                    << token::TAB << pp.name()
                    << token::TAB << avgTaufp
                    << token::TAB << tau
                    << endl;
            }

            Log << "    avg/integ(" << pp.name() << ") = "
                << avgTaufp << ", "  << tau << endl;

            this->setResult("avg(" + pp.name() + ")", avgTaufp);
            this->setResult("int(" + pp.name() + ")", tau);
    }

    return true;
}


bool Foam::functionObjects::myNormalStress::write()
{
    const auto& myNormalStress =
        obr_.lookupObject<volScalarField>(scopedName(typeName));

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << myNormalStress.name() << endl;

    myNormalStress.write();

    return true;
}


// ************************************************************************* //
