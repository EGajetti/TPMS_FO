/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

Class
    Foam::functionObjects::tau_vector

Group
    grpForcesFunctionObjects

Description
    Computes the wall-shear stress at selected wall patches.

    \f[
        \vec \tau = \vec R \cdot \vec n
    \f]

    where
    \vartable
        \vec \tau   | Wall shear stress vector [m^2/s^2]
        \vec R | Shear-stress symmetric tensor (retrieved from turbulence model)
        \vec n | Patch normal vector (into the domain)
    \endvartable

    Operands:
    \table
      Operand        | Type              | Location
      input          | -                 | -
      output file    | dat | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<field\>
      output field   | volVectorField (only boundaryField) <!--
                 --> | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    tau_vector1
    {
        // Mandatory entries (unmodifiable)
        type        tau_vector;
        libs        (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        patches     (<patch1> ... <patchN>); // (wall1 "(wall2|wall3)");

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: tau_vector            | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      patches    | Names of operand patches   | wordList | no | all wall patches
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        <solver> -postProcess -func tau_vector
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::tau_vector

SourceFiles
    tau_vector.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_tau_vector_H
#define functionObjects_tau_vector_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class tau_vector Declaration
\*---------------------------------------------------------------------------*/

class tau_vector
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the shear-stress
        void calcShearStress
        (
            const scalar nu,
            const volVectorField U,
            volVectorField& myShearStress
        );


public:

    //- Runtime type information
    TypeName("tau_vector");


    // Constructors

        //- Construct from Time and dictionary
        tau_vector
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        tau_vector(const tau_vector&) = delete;

        //- No copy assignment
        void operator=(const tau_vector&) = delete;


    //- Destructor
    virtual ~tau_vector() = default;


    // Member Functions

        //- Read the tau_vector data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
