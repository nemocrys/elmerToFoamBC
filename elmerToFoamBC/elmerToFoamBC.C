/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 AUTHOR,AFFILIATION
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
    elmerToFoamBC

Description
    Maps boundary values from Elmer to an OpenFOAM mesh, including rotation to
    3D depending on options provided in elmerToFoamDict.

Usage
    \b elmerToFoamBC
    \table
        Property     | Description             | Required    | Default value
        coordinatePermut    | permutation of coordinates | no         | (0 1 2)
        coordStartColumn    | column of first coordinate in Elmer files | yes | -
        valueLabel          | column of value in elmer files | yes | -
        gradientLabel       | column of gradient in elmer files | yes | -
        gradientSize        | number of gradient components in elmer files | yes | -
        mapMethod           | interpolation method (interp1 or nearest) | no | interp1
        debug               | debug flag for additional output           | no | 0
        boundary            | dictionary with Elmer file names for OpenFOAM boundaries | yes | -
    \endtable

    Example of the boundary dictionary specification:
    \verbatim
    boundary 
    {
        <patchName>            "<fileName>";
    }
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    argList::addNote
    (
        "Reads in a elmer dat file and outputs a points file.\n"
        "Used in conjunction with elmerToFoamDict."
    );

    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    Info << "List of all available patches:" << endl
         << (mesh.boundary().size()) << endl
        << "(" << endl;
    // get mesh dimension
    label dim = 3;
    forAll (mesh.boundary(),patchI)
    {
        if
        (
            mesh.boundary()[patchI].type() == "wedge"
            ||
            mesh.boundary()[patchI].type() == "empty"
        )
        {
            dim = dim - 1;
            break;
        }
    }

    forAll (mesh.boundary(),patchI)
    {
        Info << mesh.boundary()[patchI].name() << endl;
    }
    
    Info << ")" << endl;

    Info << endl << "Reading elmerToFoamDict" << endl;

    IOdictionary elmerToFoamDict
    (
        IOobject
        (
            "elmerToFoamDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // INPUT
    // debug_ flag; currenty only two levels
    const label debug_ = elmerToFoamDict.getOrDefault<label>("debug", 1);
    // interpolation method
    const word method_ = elmerToFoamDict.getOrDefault<word>("mapMethod", "interp1");
    if (!(method_ == "nearest" ||method_ == "interp1"))
    {
         FatalErrorInFunction
            << "Interpolation method not supported. Available methods are:\n"
            << "nearest\ninterp1";
    }
    // permutation (x,y,z) -> (x,z,y)
    const vector coordinatePermut_ = elmerToFoamDict.getOrDefault<vector>("coordinatePermut", vector(0, 1, 2));
    // column of x
    const label coordStart_(elmerToFoamDict.get<label>("coordStartColumn"));
    // column of value to read
    const label valueLabel_(elmerToFoamDict.get<label>("valueLabel"));
    // column of gradient to read
    const label gradientLabel_(elmerToFoamDict.get<label>("gradientLabel"));
    // should the rotate 
    const label rot2D_ = elmerToFoamDict.getOrDefault<label>("rot2D", 0);
    // read boundary dictionary
    const dictionary& boundaryDict_(elmerToFoamDict.subDict("boundary"));
    word file(" ");

    if (debug_) {Info << elmerToFoamDict << endl;}

    forAll (boundaryDict_.keys(), boundaryI)
    {
        const label& currPatchID = 
                mesh.boundary().findPatchID(boundaryDict_.keys()[boundaryI]);

        file = boundaryDict_.get<word>(boundaryDict_.keys()[boundaryI]);

        if(currPatchID==-1)
        {
            FatalErrorInFunction
                << "Patch "
                << boundaryDict_.keys()[boundaryI] 
                << " not found." << endl 
                << abort(FatalError);
        }

        Info << endl << "Reading file " 
             << file
             << " for boundary " 
             << mesh.boundary()[currPatchID].name()
             << endl;


        fileName posName(file);
        IFstream inFile(posName);

        // find max iteration
        scalar maxIter(Zero);

        while (inFile.good())
        {
            string line;
            inFile.getLine(line);
            scalar dd;
            if (line.size()>1)
            {
                // Set up up a stream from this line
                std::stringstream lineStr(line);

                DynamicList<scalar> data;
                string item;
                while (lineStr >> dd)
                {
                    data.append(dd);
                }
                if (data[0] > maxIter)
                {
                    maxIter = data[0];
                }
            }
        }

        Info << "Only reading Elmer values of maxIter = " << maxIter << endl;
        // Rewind the stream so that it may be read again.
        inFile.rewind();
        
        DynamicList<vector> coordinates;
        DynamicList<vector> gradient;
        DynamicList<scalar> value;

        while (inFile.good())
        {
            string line;
            inFile.getLine(line);
            scalar dd;
            if (line.size()>1)
            {
                // Set up up a stream from this line
                std::stringstream lineStr(line);

                DynamicList<scalar> data;
                string item;
                while (lineStr >> dd)
                {
                    data.append(dd);
                }

                // write only last elmer iteration
                if (data[0] == maxIter)
                {
                    // read elmer coordinates
                    // permutation (x,y,z) -> (x,z,y)
                    vector coord(Zero);
                    forAll (coord, xi)
                    {
                        coord[coordinatePermut_[xi]] = data[coordStart_+xi];
                    }            
                    
                    // read elmer gradient
                    vector field(Zero);

                    field[coordinatePermut_[0]] = data[gradientLabel_+0];
                    field[coordinatePermut_[1]] = data[gradientLabel_+1];
                    
                    coordinates.append(coord);
                    gradient.append(field);
                    // read elmer value
                    value.append(data[valueLabel_]);
                }
            }
        }

        Info<< "Reading field T" << endl;
        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        scalarField gradT(T.boundaryField()[currPatchID].size(),Zero); 
        vectorField gradTComp(T.boundaryField()[currPatchID].size(),Zero);  

        
        vectorField Nf = mesh.boundary()[currPatchID].nf(); 
        scalar maxdistall = -1000;
        
        // map elmer values to OpenFOAM faces
        Info << "Mapping file " 
             << file
             << " to boundary " 
             << T.boundaryField()[currPatchID].patch().name()
             << " with type " 
             << T.boundaryField()[currPatchID].type()
             << endl;
                
        Info<< "Interpolation method: "
            << method_ << endl;

        // TODO: Combine interpolation methods and simplify for readability
        // Interpolation method with nearest neighbor
        if (method_ == "nearest")
        {
            forAll(T.boundaryField()[currPatchID], facesi) 
            {
                vector p = mesh.Cf().boundaryField()[currPatchID][facesi];
                vector por = p;

                // project to xz plane
                scalar rr = Foam::sqrt(
                             p.x()*p.x() + p.z()*p.z() 
                             );
                p.x() = rr;
                p.y() = p.y();
                p.z() = 0;
                // Info << "; projected mesh point = " << p;

                // find nearest neighbour
                scalar mindist = 1000;
                label minpoint = 0;
                for (label i=0; i<coordinates.size(); i++) 
                {
                    scalar dist = mag(p-coordinates[i]);
                    if (dist < mindist) 
                    { 
                        mindist = dist;  
                        minpoint = i; 
                    }
                }

                if (mindist>maxdistall) 
                { 
                    maxdistall = mindist; 
                }

                if
                (
                    T.boundaryField()[currPatchID].type() 
                    == "fixedValue"
                    ||
                    T.boundaryField()[currPatchID].type() 
                    == "timeVaryingMappedFixedValue"
                )
                {
                    // write nearest value to current patch
                    T.boundaryFieldRef()[currPatchID][facesi] = value[minpoint];
                }
                if
                (
                    T.boundaryField()[currPatchID].type() 
                    == "fixedGradient"
                )
                {
                    // write nearest value to current patch
                    T.boundaryFieldRef()[currPatchID][facesi] = value[minpoint];

                    vector vecRot(Zero);
                    if (p.x()>1e-10) 
                    {
                        vecRot.x() = gradient[minpoint].x()*por.x()/(p.x()+SMALL);
                        vecRot.z() = gradient[minpoint].x()*por.z()/(p.x()+SMALL); 
                    }
                    else 
                    {
                        vecRot.x() = 0;
                        vecRot.z() = 0;
                    }
                    
                    vecRot.y() = gradient[minpoint].y();
                    
                    vector N = Nf[facesi];
                    gradT[facesi] = ( N & vecRot ); // grad in normal direction
                    if (debug_) {gradTComp[facesi] = vecRot;}
                }
            } // forAll
        Info<< "Max distance between mapping = " << maxdistall << endl;
        } // if nearest

        // Interpolation method linear between two nearest points
        if (method_ == "interp1")
        {
            forAll(T.boundaryField()[currPatchID], facesi) 
            {
                vector p = mesh.Cf().boundaryField()[currPatchID][facesi];
                vector por = p;
                scalar rr = Foam::sqrt( p.x()*p.x() + p.z()*p.z() );
                p.x() = rr;
                p.y() = p.y();
                p.z() = 0;
                // find 2 nearest points to p
                label minpoint1=0, minpoint2=0;
                scalar mindist = 1000;
                for (label i=0; i<coordinates.size(); i++) 
                {
                    scalar dist = mag(p-coordinates[i]);
                    if (dist < mindist) { mindist = dist;  minpoint1 = i; }
                }
                if (mindist>maxdistall) { maxdistall = mindist; }

                mindist = 1000;
                for (label i=0; i<coordinates.size(); i++) 
                    {
                    scalar dist = mag(p-coordinates[i]);
                    if (dist < mindist && i!=minpoint1) { mindist = dist;  minpoint2 = i; }
                    }
                if (mindist>maxdistall) { maxdistall = mindist; }

                // weighted average of 2 nearest points
                scalar dist1 = mag(p-coordinates[minpoint1]);
                scalar dist2 = mag(p-coordinates[minpoint2]);
                // write nearest value to current patch
                scalar w = dist1/(dist1+dist2);
                scalar interpValue = (1-w)*value[minpoint1] + w*value[minpoint2];

                if
                (
                    T.boundaryField()[currPatchID].type() 
                    == "fixedValue"
                    ||
                    T.boundaryField()[currPatchID].type() 
                    == "timeVaryingMappedFixedValue"
                )
                {
                    // write nearest value to current patch
                    T.boundaryFieldRef()[currPatchID][facesi] = interpValue;
                }
                if
                (
                    T.boundaryField()[currPatchID].type() 
                    == "fixedGradient"
                )
                {
                    // write nearest value to current patch
                    T.boundaryFieldRef()[currPatchID][facesi] = interpValue;
                    vector interpValue = (1-w)*gradient[minpoint1] + w*gradient[minpoint2];
                    vector vecRot(Zero);
                    if (p.x()>1e-10) 
                    {
                        vecRot.x() = interpValue.x()*por.x()/(p.x()+SMALL);
                        vecRot.z() = interpValue.x()*por.z()/(p.x()+SMALL); 
                    }
                    else 
                    {
                        vecRot.x() = 0;
                        vecRot.z() = 0;
                    }
                    
                    vecRot.y() = interpValue.y();
                    
                    vector N = Nf[facesi];
                    gradT[facesi] = ( N & vecRot ); // grad in normal direction
                    if(debug_){gradTComp[facesi] = vecRot;} // grad components

                }
            } // forAll
        } //if

        // write gradient value for boundaries of type fixedGradient
        if
        (
            T.boundaryField()[currPatchID].type() 
            == "fixedGradient"
        )
        {
            Info<< "Writing gradient" << endl;
            fixedGradientFvPatchScalarField& gradTPatch =
                refCast<fixedGradientFvPatchScalarField>
                (
                    const_cast<fvPatchScalarField&>
                    (
                        T.boundaryField()[currPatchID]
                    )
                );
            gradTPatch.gradient() = gradT;
        }
        
        // Write T field
        Info<< "Writing value" << endl;
        T.write();

        if (debug_)
        {
            // write out data for testing
            const word testDir = "debugElmerToFoamBC";

            if (!exists(testDir))
            {
                mkDir(
                    testDir
                    );
            }

            OFstream osTesting
            (
                testDir/boundaryDict_.keys()[boundaryI]
            );
            // write header
            osTesting << "x y z T gradT gradTx gradTy gradTz" << nl;
            // write values
            forAll (T.boundaryField()[currPatchID], facesi)
            {
                vector p = mesh.Cf().boundaryField()[currPatchID][facesi];
                osTesting  << float(p.x()) << ' '
                << float(p.y()) << ' '
                << float(p.z()) << ' '
                << T.boundaryField()[currPatchID][facesi] << ' '
                << gradT[facesi] << ' ' 
                << gradTComp[facesi].x() << ' ' 
                << gradTComp[facesi].y() << ' ' 
                << gradTComp[facesi].z()
                << nl;
            }

            // Write out points to constant/boundaryData/*
            const word outDir = runTime.constant()/
                                "boundaryData"/
                                boundaryDict_.keys()[boundaryI];

            if (!exists(outDir/"0"))
            {
                mkDir(
                    outDir/"0"
                    );
            }

            scalarField valueOnMesh(T.boundaryFieldRef()[currPatchID]);
            scalarField gradientOnMesh(gradT);
            vectorField coordinatesOnMesh(mesh.Cf().boundaryField()[currPatchID]);

            OFstream os
            (
                outDir/"points"
            );
            Info<< "Writing points to " << os.name() << nl;
            os  << "// Points" << nl;
            os << coordinates << nl;

            OFstream osVector
            (
                outDir/"0/T"
            );
            Info<< "Writing vector field to " << osVector.name() << nl;
            osVector  << "// Data on points"  << nl;
            osVector << gradient << nl;
            
            OFstream osScalar
            (
                outDir/"0/gradT"
            );
            Info<< "Writing scalar field to " << osScalar.name() << nl;
            osScalar  << "// Data on points"  << nl;
            osScalar << value << nl;
        }
    }

    Info << endl << "End\n" << endl;
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //