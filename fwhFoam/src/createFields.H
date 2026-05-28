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

autoPtr<sampledSurfaces> samplingPtr
(
	new sampledSurfaces
	(
		"fwhSampling",
		runTime,
		fwhSubDict.subDict("sampling")
	)
);
