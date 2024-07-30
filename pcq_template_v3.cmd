Number of uncertainty variables 3

## Input parameter ranges
PARAMETER x1-min 0.5
PARAMETER x1-max 4.5
PARAMETER x2-min 0.2
PARAMETER x2-max 0.5
PARAMETER x3-min 1000
PARAMETER x3-max 10000

## Parameters for PCQ
## "NUMBER OF QUADRATURE POINTS" & "NUMBER OF SECONDARY SAMPLING POINTs" indicate those for each parameter
MAXIMUM DEGREE OF EXPANSION 10
NUMBER OF QUADRATURE POINTS 10
NUMBER OF SECONDARY SAMPLING POINTs 100

## Threshold value
cri 0.5

## Path for output asc file (numerical result) & save folder
./RUN_RES/
./PCQ_RES/
