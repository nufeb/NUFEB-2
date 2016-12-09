/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "polyDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyDragModel> Foam::polyDragModel::New
(
    const dictionary& cloudDict,
    const IOdictionary& transDict,
    const scalarField& alpha,
    const scalarField& pd,
    const scalarField& aveD
)
{
    word polyDragModelType
    (
        cloudDict.lookup("polyDragModel")
    );


    Info<< "Selecting polyDragModel "
        << ": "
        << polyDragModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyDragModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyDragModel::New : " << endl
            << "    unknown polyDragModelType type "
            << polyDragModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid polyDragModel types are : " << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return cstrIter()(cloudDict, transDict, alpha, pd, aveD);
}


// ************************************************************************* //
