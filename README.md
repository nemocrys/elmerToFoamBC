# elmerToFoamBC
Maps boundary values from Elmer to an OpenFOAM mesh, including rotation to
3D depending on options provided in elmerToFoamDict. Developed for OpenFOAM-v2112.

## Compilation 
    wmake elmerToFoamBC

## Application
    elmerToFoamBC

## Description
    Maps boundary values from Elmer to an OpenFOAM mesh, including rotation to
    3D depending on options provided in elmerToFoamDict.

## Usage
    elmerToFoamBC

| Property         | Description                | Required | Default value |
|------------------|----------------------------|----------|---------------|
| coordinatePermut | permutation of coordinates | no       | (0 1 2)       |
| coordStartColumn    | column of first coordinate in Elmer files | yes | - |
| valueLabel          | column of value in elmer files | yes | - |
| gradientLabel       | column of gradient in elmer files | yes | - |
| gradientSize        | number of gradient components in elmer files | yes | - |
| mapMethod           | interpolation method (interp1 or nearest) | no | interp1 |
| debug               | debug flag for additional output           | no | 0 |
| boundary            | dictionary with Elmer file names for OpenFOAM boundaries | yes | - |

Example of the boundary dictionary specification:

    boundary 
    {
        <patchName>            "<fileName>";
    }

## License

nemoblock is published under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.html).

## Acknowledgements

[This project](https://www.researchgate.net/project/NEMOCRYS-Next-Generation-Multiphysical-Models-for-Crystal-Growth-Processes) has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 851768).

<img src="https://raw.githubusercontent.com/nemocrys/nemoblock/master/EU-ERC.png">

## Contribution

Any help to improve this package is very welcome!
