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
    Info << "dim = " << dim << endl;
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
    // size of gradient in Elmer 
    const label gradientSize_(elmerToFoamDict.get<label>("gradientSize"));
    // should the rotate 
    const bool rot2D_ = elmerToFoamDict.getOrDefault<label>("rot2D", 0);
    // read boundary dictionary
    const dictionary& boundaryDict_(elmerToFoamDict.subDict("boundary"));
    word file(" ");

    const vector axis_ = elmerToFoamDict.getOrDefault<vector>("axis",vector(0,0,0));
    const vector rotDir_ = elmerToFoamDict.getOrDefault<vector>("rotDir",vector(0,0,0));
    
    // vector axis_(0,0,1);
    // vector rotDir_(0,1,0);
    

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

                    for (int i = 0; i<gradientSize_; i=i+1)
                    {
                        field[coordinatePermut_[i]] = data[gradientLabel_+i];
                    }
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
        label rrLabel;
        label rotLabel;
        label axisLabel;
        if (rot2D_)
        {
            rrLabel = (vector(1,1,1)-axis_-rotDir_) & vector(0,1,2);
            rotLabel = rotDir_ & vector(0,1,2);
            axisLabel = axis_ & vector(0,1,2);    
            Info << rrLabel << ' '<< rotLabel << ' '<< axisLabel << ' ' << endl;
        }
        forAll(T.boundaryField()[currPatchID], facesi) 
        {
            vector p = mesh.Cf().boundaryField()[currPatchID][facesi];
            vector por = p;
            // calculte r
            scalar rr = 0;
            // calculate radius rr from axis_
            if(rot2D_)
            {
                forAll (p, xi)
                {
                    rr = rr + (1-axis_[xi])*p[xi]*p[xi];
                }
                rr = Foam::sqrt(rr);
                // flatten p
                p[rrLabel] = rr;
                p[axisLabel] = p[axisLabel];
                p[rotLabel] = 0;
            }
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
            // Info << "por = " << por << "; p = " << p << "; minpoint = " << coordinates[minpoint] << endl; 
            if (mindist>maxdistall) 
            { 
                maxdistall = mindist; 
            }     

            scalar interpValue;
            vector interpGradient;
            if (method_ == "nearest")
            {
                interpValue = value[minpoint];
                interpGradient = gradient[minpoint];
            }
            else if (method_ == "interp1")
            {
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
                interpValue = (1-w)*value[minpoint1] + w*value[minpoint2];
                interpGradient = (1-w)*gradient[minpoint1] + w*gradient[minpoint2];
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

                vector vecRot(Zero);
                if(rot2D_)
                {
                    if (p[rrLabel]>1e-10) 
                    {
                        vecRot[rrLabel] = interpGradient[rrLabel]*por[rrLabel]/(p[rrLabel]+SMALL);
                        vecRot[rotLabel] = interpGradient[rrLabel]*por[rotLabel]/(p[rrLabel]+SMALL); 
                    }
                    else 
                    {
                        vecRot[rrLabel] = 0;
                        vecRot[rotLabel] = 0;
                    }
                    
                    vecRot[axisLabel] = interpGradient[axisLabel];
                }
                else
                {
                    vecRot = interpGradient;
                }
                vector N = Nf[facesi];
                gradT[facesi] = ( N & vecRot ); // grad in normal direction
                if(debug_){gradTComp[facesi] = vecRot;} // grad components
            } // patch if-statement
        } //forAll

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
        }
    }

    Info << endl << "End\n" << endl;
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
