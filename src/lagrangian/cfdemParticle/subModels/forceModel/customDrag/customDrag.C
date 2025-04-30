/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Custom drag model for CFD-DEM coupling.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "customDrag.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(CustomDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    CustomDrag,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
CustomDrag::CustomDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict, sm),
    propsDict_(dict.subDict(typeName + "Props")),
    U_(sm.mesh().lookupObject<volVectorField>(propsDict_.lookup("velFieldName"))),
    scaleDia_(1.0),
    scaleDrag_(1.0)
{
    if (propsDict_.found("scale"))
        scaleDia_ = readScalar(propsDict_.lookup("scale"));
    if (propsDict_.found("scaleDrag"))
        scaleDrag_ = readScalar(propsDict_.lookup("scaleDrag"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

CustomDrag::~CustomDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CustomDrag::setForce() const
{
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    for (int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        vector drag(0,0,0);
        label cellI = particleCloud_.cellIDs()[index][0];

        if (cellI > -1) // Particle found in cell
        {
            vector Us = particleCloud_.velocity(index);
            vector Ur = U_[cellI] - Us;
            scalar ds = 2 * particleCloud_.radius(index);
            scalar ds_scaled = ds / scaleDia_;
            scalar nuf = nufField[cellI];
            scalar rho = rhoField[cellI];
            scalar magUr = mag(Ur);
            scalar Rep = 0;
            scalar Cd = 0;

            if (magUr > 0)
            {
                // Calculate particle Reynolds number
                Rep = ds_scaled * magUr / nuf;

                // Calculate drag coefficient using custom formula
                Cd = (24.0 / Rep) * (1.0 + 0.15 * pow(Rep, 0.681)) + (0.407 / (1.0 + (8710.0 / Rep)));

                // Calculate drag force
                drag = 0.125 * Cd * rho * M_PI * ds * ds * magUr * Ur * scaleDrag_;
            }

            // Write drag to particle data array
            forceSubM(0).partToArray(index, drag, vector::zero);
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
