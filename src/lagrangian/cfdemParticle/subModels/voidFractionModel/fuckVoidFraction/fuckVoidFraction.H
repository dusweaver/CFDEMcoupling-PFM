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
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
    fuckVoidFraction

SourceFiles
    fuckVoidFraction.C

\*---------------------------------------------------------------------------*/

#ifndef fuckVoidFraction_H
#define fuckVoidFraction_H

#include "mathExtra.H"
#include "voidFractionModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class fuckVoidFraction
:
    public voidFractionModel
{

private:
        static const int numberOfMarkerPoints = 29;

        dictionary propsDict_;

        const bool verbose_;

        const Switch procBoundaryCorrection_;

    const scalar alphaMin_;
    const scalar maxParticleSize_;
    const scalar radiusFraction_;
    scalar averagingRadius_;
        scalar scaleVol_;               // Volume scaling factor
        scalar scaleRadius_;            // Radius scaling factor
    

    bool cfdemUseOnly_;

    
        const bool interpolation_;

        vector offsets[numberOfMarkerPoints];

       inline scalar Vp(scalar radius, scalar scaleVol) const
        {
            return constant::mathematical::fourPiByThree * radius * radius * radius * scaleVol;
        }
        inline double particleRadius(label index) const;

        void buildSuperCell(const label centralCellID, const scalar radius, labelHashSet& superCell) const;

public:

    //- Runtime type information
    TypeName("fuck");


    // Constructors

        //- Construct from components
        fuckVoidFraction
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~fuckVoidFraction();


    // Member Functions
        void setvoidFraction(double** const& ,double**&, double**&, double**&, double**&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
