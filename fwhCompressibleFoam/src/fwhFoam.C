/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
	\\  /    A nd           | www.openfoam.com
	 \\/     M anipulation  |
-------------------------------------------------------------------------------
	Copyright (C) 2026 AUTHOR,AFFILIATION
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

Application
	fwhCompressibleFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sampledSurfaces.H"
#include "interpolation.H"
#include "mathematicalConstants.H"

extern "C"
{
    #include <fftw3.h>
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Aeroacoustic solver implementing the Ffowcs Williams-Hawkings (FW-H)"
        " acoustic analogy to predict far-field sound from CFD flow data."
    );

	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"

	Info << "Reading fwhDict" << endl;

	#include "createFields.H"

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info << endl;
	Info << "Performing checks ..." << endl;

	#include "dictErrors.H"
	#include "printFWHSettings.H"
	#include "printSamplingSettings.H"
	#include "printVolumeFieldSettings.H"

	const label nSurfaces = samplingPtr->size();

	if (nSurfaces == 0)
		FatalErrorInFunction
			<< "No surfaces defined in fwhDict sampling sub-dict."
			<< exit(FatalError);

	const fileName rootDir =
		runTime.rootPath() / runTime.globalCaseName()
		/ "postProcessing" / "fwhCompressibleFoam";

	mkDir(rootDir);

	// Loop over every surface defined in the sampling block
	for (label surfI = 0; surfI < nSurfaces; ++surfI)
	{
		sampledSurface& surf = (*samplingPtr)[surfI];

		Info << nl
             << "======================================" << nl
             << "  Surface " << surfI + 1 << " / " << nSurfaces
             << ": " << surf.name() << nl
             << "======================================" << nl;

		#include "readAndSample.H"
		#include "fwhIntegrate.H"
		#include "writeVolumeFields.H"
		#include "writeOutputs.H"
	}

	Info << nl << nl;
	runTime.printExecutionTime(Info);
	Info << "End\n" << endl;

	return 0;
}


// ************************************************************************* //
