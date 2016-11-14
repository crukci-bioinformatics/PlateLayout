#!/bin/bash

# locate test script directory

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
TestDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
TestDIR=`readlink -f $TestDIR`


# export binary file location
BinaryDIR=${TestDIR/test/bin}
export PATH=$BinaryDIR:$PATH

# tests
echo "Run: metadata_template_TEST_11x4.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_11x4.csv -o MyPlateLayout1 -b Sex,ExtractionInformation -r 3

echo "Run: metadata_template_TEST_12-12-14.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_12-12-14.csv -o MyPlateLayout -b Sex,ExtractionInformation,PassageNumber

echo "Run: metadata_template_TEST_12x3.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_12x3.csv -o MyPlateLayout2 -b Sex,ExtractionInformation -r 2

echo "Run: metadata_template_TEST_12x4.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_12x4.csv

echo "Run: metadata_template_TEST_12x8.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_12x8.csv -o MyPlateLayout3 -G

echo "Run: metadata_template_TEST_21x8.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_21x8.csv -b Sex,PassageNumber

echo "Run: metadata_template_TEST_23x8.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_23x8.csv -o MyPlateLayout5 -b Sex,ExtractionInformation,PassageNumber

echo "Run: metadata_template_TEST_24x6.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_24x6.csv -o MyPlateLayout6 -b Sex,ExtractionInformation,PassageNumber

echo "Run: metadata_template_TEST_24x7.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_24x7.csv -b Sex,ExtractionInformation,PassageNumber -r 2

echo "Run: metadata_template_TEST_24x8.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_24x8.csv -o MyPlateLayout8

echo "Run: metadata_template_TEST_5x3.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_5x3.csv

echo "Run: metadata_template_TEST_100.csv"
PlateLayoutRandomistation -d $TestDIR/metadata_template_TEST_100.csv -o MyPlateLayout10 -b ExtractionInformation,PassageNumber
