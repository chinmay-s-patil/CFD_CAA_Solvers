// createFWHFields.H
// Reads fwhDict and sets up FW-H sampling surfaces and history storage.
// Included once in pimpleFWHFoam.C before the time loop.

Info << "Reading fwhDict" << endl;

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

// Observers
const List<vector> observers(fwhSubDict.lookup("observers"));

// Output controls
const Switch writeFT              = fourierTransform.get<Switch>("writeFT");
const Switch writeSignal          = fwhSubDict.get<Switch>("writeSignal");
const Switch writeSPL             = fwhSubDict.get<Switch>("writeSPL");
const word   timeInterpolationScheme = fwhSubDict.get<word>("timeInterpolationScheme");
const word   fourierTransformType    = fourierTransform.get<word>("type");

// Sampling controls
const word interpolationScheme = samplingSubDict.get<word>("faceInterpolationScheme");
const word surfaceFormat       = samplingSubDict.get<word>("surfaceFormat");
const wordList fwhFields(samplingSubDict.lookup("fields"));

// Volumetric field output switches
const Switch writeVolumePPrime =
    fwhSubDict.getOrDefault<Switch>("writeVolumePPrime", false);

const Switch writeVolumeOASPL =
    fwhSubDict.getOrDefault<Switch>("writeVolumeOASPL", false);

// OASPL directivity
const Switch writeOASPL =
    fwhSubDict.getOrDefault<Switch>("writeOASPL", false);

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

autoPtr<sampledSurfaces> fwhSamplingPtr
(
	new sampledSurfaces
	(
		"fwhSampling",
		runTime,
		fwhSubDict.subDict("sampling")
	)
);

const label fwhNSurfaces = fwhSamplingPtr->size();

if (fwhNSurfaces == 0)
    FatalErrorInFunction
        << "No surfaces defined in fwhDict sampling sub-dict."
        << exit(FatalError);

const fileName fwhRootDir =
    runTime.rootPath() / runTime.globalCaseName()
    / "postProcessing" / "pimpleFWHFoam";

mkDir(fwhRootDir);

// Per-surface history storage (grows each time step)
List<DynamicList<scalar>>      fwhTimeValues(fwhNSurfaces);
List<DynamicList<scalarField>> fwhPHistory(fwhNSurfaces);
List<DynamicList<vectorField>> fwhUHistory(fwhNSurfaces);

// Per-surface geometry (grabbed on first sample)
List<pointField>  fwhFaceCentres(fwhNSurfaces);
List<vectorField> fwhFaceNormals(fwhNSurfaces);
List<scalarField> fwhFaceAreas(fwhNSurfaces);
List<bool>        fwhGeometryGrabbed(fwhNSurfaces, false);

Info << nl << "FW-H: " << fwhNSurfaces << " sampling surface(s) configured." << nl;
