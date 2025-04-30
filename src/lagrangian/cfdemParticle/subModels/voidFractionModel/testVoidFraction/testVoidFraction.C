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
    along with CFDEMcoupling; if not, write to the Freftware Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "testVoidFraction.H"
#include "mathExtra.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(testVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    testVoidFraction,
    dictionary
);


// // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// // Construct from components
// testVoidFraction::testVoidFraction
// (
//     const dictionary& dict,
//     cfdemCloud& sm
// )
// :
//     voidFractionModel(dict,sm),
//     propsDict_(dict.subDict(typeName + "Props")),
//     alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
//     alphaLimited_(0)
// {
//     checkWeightNporosity(propsDict_);
//     if(porosity()!=1) FatalError << "porosity not used in testVoidFraction" << abort(FatalError);
// }


// // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// testVoidFraction::~testVoidFraction()
// {}


testVoidFraction::testVoidFraction(const dictionary& dict, cfdemCloud& sm)
    : voidFractionModel(dict, sm),
      propsDict_(dict.subDict(typeName + "Props")),
      alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
      maxParticleSize_(readScalar(propsDict_.lookup("maxParticleSize"))),
      radiusFraction_(readScalar(propsDict_.lookup("radiusFraction"))),
      scaleVol_(weight()),
      scaleRadius_(cbrt(porosity())),
    interpolation_(false),
    cfdemUseOnly_(false)
{
    averagingRadius_ = radiusFraction_ * maxParticleSize_;

    // Generate marker points
    int m = 0;
    offsets[m][0] = offsets[m][1] = offsets[m][2] = 0.0;
    m++;

    for (label ir = 0; ir <= 1; ir++)
    {
        scalar r = ir == 0 ? 0.5 : 1.0; // Uniform distribution scaling
        for (scalar zeta = 0; zeta < 2.0 * M_PI; zeta += M_PI / 2.0)
        {
            for (scalar theta = 0; theta < M_PI; theta += M_PI / 2.0)
            {
                offsets[m][0] = r * Foam::sin(theta) * Foam::cos(zeta);
                offsets[m][1] = r * Foam::sin(theta) * Foam::sin(zeta);
                offsets[m][2] = r * Foam::cos(theta);
                m++;
            }
        }
    }

    Info << "UniformVoidFraction initialized with radius: " << averagingRadius_ << endl;
}

testVoidFraction::~testVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void testVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV)
{
    // if (cfdemUseOnly_)
    //     reAllocArrays(particleCloud_.numberOfParticles());
    // else
    //     reAllocArrays();

}


// void testVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV)
// {
//     scalar radius(-1);
//     scalar volume(0);
//     scalar cellVol(0);
//     scalar scaleVol = weight();

//     for(int index=0; index< particleCloud_.numberOfParticles(); index++)
//     {
//         //if(mask[index][0])
//         //{
//             // reset
//             particleWeights[index][0]=0;
//             cellsPerParticle()[index][0]=1;

//             label cellI = particleCloud_.cellIDs()[index][0];

//             if (cellI >= 0)  // particel test is in domain
//             {
//                 if (multiWeights_) scaleVol = weight(index);
//                 cellVol = voidfractionNext_.mesh().V()[cellI];
//                 radius = particleCloud_.radius(index);
//                 volume = constant::mathematical::fourPiByThree*radius*radius*radius*scaleVol;

//                 // store volume for each particle
//                 particleVolumes[index][0] = volume;
//                 particleV[index][0] = volume;

//                 voidfractionNext_[cellI] -= volume/cellVol;

//                 if(voidfractionNext_[cellI] < alphaMin_ )
//                 {
//                     voidfractionNext_[cellI] = alphaMin_;
//                     alphaLimited_ = 1;
//                 }

//                 if(index==0 && alphaLimited_) Info<<"alpha limited to" <<alphaMin_<<endl;

//                 // store voidFraction for each particle
//                 voidfractions[index][0] = voidfractionNext_[cellI];

//                 // store cellweight for each particle  - this should not live here
//                 particleWeights[index][0] = 1;

//                 /*//OUTPUT
//                 if (index==0)
//                 {
//                     Info << "test cellI = " << cellI << endl;
//                     Info << "cellsPerParticle =" << cellsPerParticle()[index][0] << endl;

//                     for(int i=0;i<cellsPerParticle()[index][0];i++)
//                     {
//                        if(i==0)Info << "cellids, voidfractions, particleWeights, : \n";
//                        Info << particleCloud_.cellIDs()[index][i] << " ," << endl;
//                        Info << voidfractions[index][i] << " ," << endl;
//                        Info << particleWeights[index][i] << " ," << endl;
//                      }
//                 }*/
//             }
//         //}
//     }
//     voidfractionNext_.correctBoundaryConditions();

//     // bring voidfraction from Eulerian Field to particle array
//     for(int index=0; index< particleCloud_.numberOfParticles(); index++)
//     {
//         label cellID = particleCloud_.cellIDs()[index][0];

//         if(cellID >= 0)
//         {
//             voidfractions[index][0] = voidfractionNext_[cellID];
//         }
//         else
//         {
//             voidfractions[index][0] = -1.;
//         }
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
