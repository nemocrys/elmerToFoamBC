/*

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------

*/

#include "argList.H"
#include "IFstream.H"
#include "fvCFD.H"

using namespace Foam;

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
    Info << elmerToFoamDict << endl;
    
    // INPUT
    // interpolation method
    const word method = elmerToFoamDict.getOrDefault<word>("method", "interp1");
    // permutation (x,y,z) -> (x,z,y)
    const vector coordinatePermut(elmerToFoamDict.lookup("coordinatePermut"));
    // column of x
    const label coordStart(elmerToFoamDict.get<label>("coordStartColumn"));
    // column of value to read
    const label valueLabel(elmerToFoamDict.get<label>("valueLabel"));
    // column of gradient to read
    const label gradientLabel(elmerToFoamDict.get<label>("gradientLabel"));
    // size of gradient in elmer (with permutation)
    const label gradientSize(elmerToFoamDict.get<label>("gradientSize"));
    // name of value
    const word sName(elmerToFoamDict.lookup("valueName"));
    // name of gradient
    const word vName(elmerToFoamDict.lookup("gradientName"));

    const dictionary& boundaryDict(elmerToFoamDict.subDict("boundary"));
    word file(" ");


    forAll (boundaryDict.keys(), boundaryI)
    {
        Info << "here"<<endl;
        const label& currPatchID = 
                mesh.boundary().findPatchID(boundaryDict.keys()[boundaryI]);

        file = boundaryDict.get<word>(boundaryDict.keys()[boundaryI]);

        if(currPatchID==-1)
        {
            FatalErrorInFunction
                << "Patch "
                << boundaryDict.keys()[boundaryI] 
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
                        coord[coordinatePermut[xi]] = data[coordStart+xi];
                    }            
                    
                    // read elmer gradient
                    vector field(Zero);
                    for (int i=0;i<gradientSize;i++)
                    {
                        field[coordinatePermut[0]] = data[gradientLabel+0];
                    }
                    field[coordinatePermut[1]] = data[gradientLabel+1];
                    
                    coordinates.append(coord);
                    gradient.append(field);
                    // read elmer value
                    value.append(data[valueLabel]);
                }
            }
        }

        // ---
        Info<< "\nMapping data\n" << endl;

        Info<< "Reading field T\n" << endl;
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
        vectorField gradTFull(T.boundaryField()[currPatchID].size(),Zero);  
        vectorField Nf = mesh.boundary()[currPatchID].nf(); 
        
        scalar maxdistall = -1000;
        
        // map elmer values to OpenFOAM faces
        Info << endl << "Mapping file " 
             << file
             << " to boundary " 
             << T.boundaryField()[currPatchID].patch().name()
             << " with type " 
             << T.boundaryField()[currPatchID].type()
             << endl;
        
        forAll(T.boundaryField()[currPatchID], facesi) 
        {
            vector p = mesh.Cf().boundaryField()[currPatchID][facesi];
            // Info << "current mesh point = " << p;
            vector por = p;

            // project to xz plane
            scalar rr = Foam::sqrt( p.x()*p.x() + p.y()*p.y() );
            p.x() = rr;
            p.y() = 0;
            p.z() = p.z();
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
                    vecRot.y() = gradient[minpoint].x()*por.y()/(p.x()+SMALL); 
                }
                else 
                {
                    vecRot.x() = 0;
                    vecRot.y() = 0;
                }
                
                vecRot.z() = gradient[minpoint].z();
                
                vector N = Nf[facesi];
                gradT[facesi] = ( N & vecRot ); // grad in normal direction
                gradTFull[facesi] = vecRot;
                // Info << "gradient[minpoint] = "  << gradient[minpoint] << endl;
                // Info << "vecRot = " <<  vecRot << endl;
                // Info << "N = " << N << endl;
                // Info <<  "gradT[facesi] = "<<  gradT[facesi] << endl;
            }
            
            // Info << "; nearest elmer neighbor = " << coordinates[minpoint] << endl;

        } // forAll
        
        Info<< "Max distance between mapping = " << maxdistall << endl;

        // Other mapping method
        forAll(T.boundaryField()[currPatchID], facesi) 
        {
            vector p = mesh.Cf().boundaryField()[currPatchID][facesi];
            vector por = p;
            // find 3 nearest points to p
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

            // mindist = 1000;
            // for (label i=0; i<coordinates.size(); i++) 
            //     {
            //     scalar dist = mag(p-coordinates[i]);
            //     if (dist < mindist && i!=minpoint1 && i!=minpoint2) { mindist = dist;  minpoint3 = i; }
            //     }
            // if (mindist>maxdistall) { maxdistall = mindist; }

            // Info << "minpoint1 = " << minpoint1
            //      << "; minpoint2 = " << minpoint2
            //      << "; minpoint3 = " << minpoint3 
            //      << endl;
            // Info << "p = " << p << endl;
            // Info << "coordinates1 = " << coordinates[minpoint1]
            //      << "; coordinates2 = " << coordinates[minpoint2]
            //      << "; coordinates3 = " << coordinates[minpoint3]
            //      << endl;

            // weighted average of 2 nearest points
            scalar dist1 = mag(p-coordinates[minpoint1]);
            scalar dist2 = mag(p-coordinates[minpoint2]);
            // write nearest value to current patch
            // Info << dist1/(dist1+dist2) << endl;
            scalar w = dist1/(dist1+dist2);
            scalar interpValue = (1-w)*value[minpoint1] + w*value[minpoint2];
            // Info << "dist1 = " << dist1 << "; dist2 = " << dist2 << "; w =" << w << endl;
            // Info << "value[minpoint1]" << value[minpoint1] 
            //     << "value[minpoint2]" << value[minpoint2] 
            //     << "interpValue = " << interpValue << endl;
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
                    vecRot.y() = interpValue.x()*por.y()/(p.x()+SMALL); 
                }
                else 
                {
                    vecRot.x() = 0;
                    vecRot.y() = 0;
                }
                
                vecRot.z() = interpValue.z();
                
                vector N = Nf[facesi];
                gradT[facesi] = ( N & vecRot ); // grad in normal direction
                gradTFull[facesi] = vecRot;
                // Info << "interpValue = "  << interpValue << endl;
                // Info << "vecRot = " <<  vecRot << endl;
                // Info << "N = " << N << endl;
                // Info <<  "gradT[facesi] = "<<  gradT[facesi] << endl;
            }
        }

        // write gradient value for boundaries of type fixedGradient
        if
        (
            T.boundaryField()[currPatchID].type() 
            == "fixedGradient"
        )
        {
            fixedGradientFvPatchScalarField& gradTPatch =
                refCast<fixedGradientFvPatchScalarField>
                (
                    const_cast<fvPatchScalarField&>
                    (
                        T.boundaryField()[currPatchID]
                    )
                );
            gradTPatch.gradient() = gradT;
            // Info << gradTPatch.gradient()<< endl;
        }
        
        // Write T field
        T.write();

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
            testDir/boundaryDict.keys()[boundaryI]
        );

        forAll (T.boundaryField()[currPatchID], facesi)
        {
            vector p = mesh.Cf().boundaryField()[currPatchID][facesi];
            osTesting  << float(p.x()) << ' '
            << float(p.y()) << ' '
            << float(p.z()) << ' '
            << T.boundaryField()[currPatchID][facesi] << ' '
            << gradT[facesi] << ' ' 
            << gradTFull[facesi].x() << ' ' 
            << gradTFull[facesi].y() << ' ' 
            << gradTFull[facesi].z()
            << nl;
        }

        // Write out points to constant/boundaryData/*
        const word outDir = runTime.constant()/
                            "boundaryData"/
                            boundaryDict.keys()[boundaryI];

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
            outDir/"0"/vName
        );
        Info<< "Writing vector field to " << osVector.name() << nl;
        osVector  << "// Data on points"  << nl;
        osVector << gradient << nl;
        
        OFstream osScalar
        (
            outDir/"0"/sName
        );
        Info<< "Writing scalar field to " << osScalar.name() << nl;
        osScalar  << "// Data on points"  << nl;
        osScalar << value << nl;

    }

    Info << endl << "End\n" << endl;
    return 0;
}
