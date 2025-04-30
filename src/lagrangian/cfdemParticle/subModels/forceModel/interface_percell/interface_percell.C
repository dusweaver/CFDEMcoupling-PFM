/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    optimized interface_percell model with per-cell precompute
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "interface_percell.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorField.H"      // for vectorField type

namespace Foam
{

defineTypeNameAndDebug(interface_percell, 0);

addToRunTimeSelectionTable(
    forceModel,
    interface_percell,
    dictionary
);

// Constructor
interface_percell::interface_percell
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict, sm),
    propsDict_(dict.subDict(typeName + "Props")),
    VOFvoidfractionFieldName_(propsDict_.lookup("VOFvoidfractionFieldName")),
    alpha_(sm.mesh().lookupObject<volScalarField>(VOFvoidfractionFieldName_)),
    gradAlphaName_(propsDict_.lookup("gradAlphaName")),
    gradAlpha_(sm.mesh().lookupObject<volVectorField>(gradAlphaName_)),
    sigma_(readScalar(propsDict_.lookup("sigma"))),
    theta_(readScalar(propsDict_.lookup("theta"))),
    alphaThreshold_(readScalar(propsDict_.lookup("alphaThreshold"))),
    deltaAlphaIn_(readScalar(propsDict_.lookup("deltaAlphaIn"))),
    deltaAlphaOut_(readScalar(propsDict_.lookup("deltaAlphaOut"))),
    C_(1.0),
    interpolation_(false),
    alphaInterpolator_(interpolation<scalar>::New("cellPoint", alpha_)),
    gradAlphaInterpolator_(interpolation<vector>::New("cellPoint", gradAlpha_)),
    cellinterface_percellVec_(sm.mesh().nCells()),
    Fatt0_(
        mag(
            6.0 * sigma_
          * sin(M_PI - theta_/2.0)
          * sin(M_PI + theta_/2.0)
        )
      * M_PI
    )
{
    if (propsDict_.found("C")) C_ = readScalar(propsDict_.lookup("C"));
    if (propsDict_.found("interpolation")) interpolation_ = true;

    setForceSubModels(propsDict_);
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT, true);
    forceSubM(0).readSwitches();
    particleCloud_.checkCG(false);
}

// Destructor
interface_percell::~interface_percell() {}

// setForce: optimized version
void interface_percell::setForce() const
{
    // Precompute one interface_percell-vector per cell
    const fvMesh& mesh = alpha_.mesh();
    const scalarField& alphaFld = alpha_.internalField();
    const vectorField& gradFld  = gradAlpha_.internalField();
    label nCells = mesh.nCells();

    for (label c = 0; c < nCells; ++c)
    {
        scalar a = alphaFld[c];
        if (a <= alphaThreshold_ - deltaAlphaIn_
         || a >= alphaThreshold_ + deltaAlphaOut_)
        {
            cellinterface_percellVec_[c] = vector::zero;
        }
        else
        {
            vector nVec = gradFld[c] / max(mag(gradFld[c]), SMALL);
            scalar ramp = tanh(a - alphaThreshold_);
            cellinterface_percellVec_[c] = -nVec * (Fatt0_ * ramp * C_);
        }
    }

    // Apply cell vectors to each particle
    label nP = particleCloud_.numberOfParticles();
    for (label p = 0; p < nP; ++p)
    {
        label cellI = particleCloud_.cellIDs()[p][0];
        if (cellI < 0) continue;

        scalar dp = 2.0 * particleCloud_.radius(p);
        vector  F  = cellinterface_percellVec_[cellI] * dp;
        forceSubM(0).partToArray(p, F, vector::zero);
    }
}

} // End namespace Foam

// ************************************************************************* //
