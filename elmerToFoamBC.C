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
        const label& currPatchID = 
                mesh.boundary().findPatchID(boundaryDict.keys()[boundaryI]);

        file = boundaryDict.get<word>(boundaryDict.keys()[boundaryI]);
        Info << endl << "Reading file " 
             << file
             << " for boundary " 
             << mesh.boundary()[currPatchID].name()
             << endl;


        fileName posName(file);
        IFstream inFile(posName);

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

                // permutation (x,y,z) -> (x,z,y)
                vector coord(Zero);
                forAll (coord, xi)
                {
                    coord[coordinatePermut[xi]] = data[coordStart+xi];
                }            
                

                // read gradient
                vector field(Zero);
                for (int i=0;i<gradientSize;i++)
                {
                    field[coordinatePermut[0]] = data[gradientLabel+0];
                }
                field[coordinatePermut[1]] = data[gradientLabel+1];
                coordinates.append(coord);
                gradient.append(field);
                value.append(data[valueLabel]);
                }
        }

        Info<< "  coordinates.size() = " << coordinates.size() << endl;
        Info<< "  gradient.size() = " << gradient.size() << endl;
        Info<< "  value.size() = " << value.size() << endl;
        Info<< "  coordinates[points.size()-1] = " << coordinates[coordinates.size()-1] << endl;
        Info<< "  gradient[gradient.size()-1] = " << gradient[gradient.size()-1] << endl;
        Info<< "  value[value.size()-1] = " << value[value.size()-1] << endl;

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
        vectorField Nf = mesh.boundary()[currPatchID].nf(); 
        
        scalar maxdistall = -1000;
        
        // map elmer values to OpenFOAM faces
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
                // Info << T.boundaryField()[currPatchID].snGrad() << endl;
                // Info << gradT << endl;

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
                
                vecRot.z() = gradient[minpoint].y();
                
                vector N = Nf[facesi];
                gradT[facesi] = ( N & vecRot ); // grad in normal direction
                // Info << gradT[facesi] << endl;
            }
            
            // Info << "; nearest elmer neighbor = " << coordinates[minpoint] << endl;

        } // forAll
        
        Info<< ";  maxdistall = " << maxdistall << endl;
        
        // THIS IS NOT WRITTEN...
        T.boundaryField()[currPatchID].snGrad() = gradT;
        
        // Info << T.boundaryField()[currPatchID].gradient() << endl;
        
        T.write();

        // write out data for testing
        OFstream osTesting
        (
            boundaryDict.keys()[boundaryI]
        );
        forAll (T.boundaryField()[currPatchID], facesi)
        {
            vector p = mesh.Cf().boundaryField()[currPatchID][facesi];
            osTesting  << float(p.x()) << ' '
            << float(p.y()) << ' '
            << float(p.z()) << ' '
            << T.boundaryField()[currPatchID][facesi] 
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
        os << coordinatesOnMesh << nl;

        OFstream osVector
        (
            outDir/"0"/vName
        );
        Info<< "Writing vector field to " << osVector.name() << nl;
        osVector  << "// Data on points"  << nl;
        osVector << gradientOnMesh << nl;
        
        OFstream osScalar
        (
            outDir/"0"/sName
        );
        Info<< "Writing scalar field to " << osScalar.name() << nl;
        osScalar  << "// Data on points"  << nl;
        osScalar << valueOnMesh << nl;
    
    }

    Info << endl << "End\n" << endl;
    return 0;
}
