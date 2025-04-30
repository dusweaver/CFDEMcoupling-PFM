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
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU ral Public License
    for more details.

    You should have received a copy of the GNU Geneublic License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "fuckVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fuckVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    fuckVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
fuckVoidFraction::fuckVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(propsDict_.found("verbose")),
      alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
      maxParticleSize_(readScalar(propsDict_.lookup("maxParticleSize"))),
      radiusFraction_(readScalar(propsDict_.lookup("radiusFraction"))),
      scaleVol_(weight()),
      scaleRadius_(cbrt(porosity())),
    interpolation_(false),
    cfdemUseOnly_(false)
{


        averagingRadius_ = radiusFraction_ * maxParticleSize_;

        // Generate uniform distribution offsets
        int m = 0;
        offsets[m++] = vector(0, 0, 0);
        for (scalar zeta = 0; zeta < 2.0 * M_PI; zeta += M_PI / 2.0)
        {
            for (scalar theta = 0; theta < M_PI; theta += M_PI / 2.0)
            {
                offsets[m++] = vector(
                    averagingRadius_ * Foam::sin(theta) * Foam::cos(zeta),
                    averagingRadius_ * Foam::sin(theta) * Foam::sin(zeta),
                    averagingRadius_ * Foam::cos(theta));
            }
        }

    Info << "UniformVoidFraction initialized with radius: " << averagingRadius_ << endl;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fuckVoidFraction::~fuckVoidFraction()
{}

void fuckVoidFraction::buildSuperCell(const label centralCellID, const scalar radius, labelHashSet& superCell) const {
    const vector centralPosition = particleCloud_.mesh().C()[centralCellID];
    const scalar searchRadius = radius;

    for (Foam::label cellI = 0; cellI < particleCloud_.mesh().nCells(); ++cellI) {
        const Foam::vector& cellPosition = particleCloud_.mesh().C()[cellI];

        // Calculate the distance from the central cell position to the current cell position
        if (mag(cellPosition - centralPosition) <= searchRadius) {
            superCell.insert(cellI);
        }
    }
}


void fuckVoidFraction::setvoidFraction(double** const& mask, double**& voidfractions, double**& particleWeights, double**& particleVolumes, double**& particleV)
{
    for (int index = 0; index < particleCloud_.numberOfParticles(); ++index) {
        if (!checkParticleType(index)) continue;

        // Initialize values for the current particle
        vector position = particleCloud_.position(index);
        label cellID = particleCloud_.cellIDs()[index][0];
        scalar radius = particleRadius(index);
        scalar volume = constant::mathematical::fourPiByThree * radius * radius * radius;

        // Set up super-cell based on the averaging radius
        labelHashSet superCell;
        buildSuperCell(cellID, averagingRadius_, superCell);

        // Calculate V_SC as the total volume of all cells in the super-cell
        scalar V_SC = 0.0;
        forEachIter(labelHashSet, superCell, it) {
            label scCellID = superCell.key(it);  // Extract the key as a label
            V_SC += particleCloud_.mesh().V()[scCellID];
        }

        // Ensure that V_SC is not zero to avoid division by zero errors
        if (V_SC > 0) {
            // Calculate the volume fraction contribution (alpha_pj) for this super-cell
            scalar alpha_pj = volume / V_SC;

            // Apply this volume fraction to all cells in the super-cell
            forEachIter(labelHashSet, superCell, it) {
                label scCellID = superCell.key(it);  // Extract the key as a label
                voidfractions[index][scCellID] = alpha_pj;
            }
        }
    }

    // Ensure boundary conditions are handled after setting void fractions
    voidfractionNext_.correctBoundaryConditions();

    // Bring voidfraction from Eulerian Field to particle array
    for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {
        for (int subcell = 0; subcell < cellsPerParticle()[index][0]; subcell++) {
            label cellID = particleCloud_.cellIDs()[index][subcell];
            if (cellID >= 0) {
                voidfractions[index][subcell] = voidfractionNext_[cellID];
            } else {
                voidfractions[index][subcell] = -1.0;
            }
        }
    }
}



inline double fuckVoidFraction::particleRadius(label index) const
{
    return particleCloud_.radius(index);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
