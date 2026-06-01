IOdictionary fwhDict
(
	IOobject
	(
		"fwhDict",
		runTime.system(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

const dictionary& fwhSubDict       = fwhDict.subDict("FWH");
const dictionary& samplingSubDict  = fwhSubDict.subDict("sampling");
const dictionary& fourierTransform = fwhSubDict.subDict("fourierTransform");

// Acoustic constants
const scalar rho  = fwhSubDict.get<scalar>("rho");
const scalar c0   = fwhSubDict.get<scalar>("c0");
const scalar pRef = fwhSubDict.get<scalar>("pRef");

// Time controls
const scalar startTime      = fwhSubDict.get<scalar>("startTime");
const scalar endTime        = fwhSubDict.get<scalar>("endTime");
const label  timeStepStride = fwhSubDict.get<label>("timeStepStride");

// Observers
const List<vector> observers(fwhSubDict.lookup("observers"));

// Output controls
const Switch writeFT              = fourierTransform.get<Switch>("writeFT");
const Switch writeSignal          = fwhSubDict.get<Switch>("writeSignal");
const Switch writeSPL             = fwhSubDict.get<Switch>("writeSPL");
const Switch recomputeTimeHistory = fwhSubDict.get<Switch>("recomputeTimeHistory");
const Switch writeTimeHistory     = fwhSubDict.get<Switch>("writeTimeHistory");
const word   timeInterpolationScheme = fwhSubDict.get<word>("timeInterpolationScheme");
const word   fourierTransformType    = fourierTransform.get<word>("type");

// Sampling controls (used for validation and printing)
const word interpolationScheme = samplingSubDict.get<word>("faceInterpolationScheme");
const word surfaceFormat       = samplingSubDict.get<word>("surfaceFormat");
const wordList fields(samplingSubDict.lookup("fields"));

// * * * * * * * * * *  Volumetric field output switches  * * * * * * * * * * //
// Master switches; sub-dict entries are read lazily inside writeVolumeFields.H
// and writeOutputs.H when the corresponding switch is enabled.

const Switch writeVolumePPrime =
    fwhSubDict.getOrDefault<Switch>("writeVolumePPrime", false);

const Switch writeVolumeOASPL =
    fwhSubDict.getOrDefault<Switch>("writeVolumeOASPL", false);

// * * * * * * * * *  OASPL directivity (from observer arc) * * * * * * * * * //
const Switch writeOASPL =
    fwhSubDict.getOrDefault<Switch>("writeOASPL", false);

// OASPL geometry params read from OASPL{} sub-dict only when switch is on
// so that fwhDict without the sub-block still parses cleanly.
vector oasplCentre(vector::zero);
vector oasplRefDirection(vector(1, 0, 0));
word   oasplFileName("oaspl_directivity");

if (writeOASPL)
{
    const dictionary& oasplDict = fwhSubDict.subDict("OASPL");
    oasplCentre       = oasplDict.getOrDefault<vector>("centre",       vector::zero);
    oasplRefDirection = oasplDict.getOrDefault<vector>("refDirection", vector(1, 0, 0));
    oasplFileName     = oasplDict.getOrDefault<word>  ("fileName",     "oaspl_directivity");
}

autoPtr<sampledSurfaces> samplingPtr
(
	new sampledSurfaces
	(
		"fwhSampling",
		runTime,
		fwhSubDict.subDict("sampling")
	)
);
