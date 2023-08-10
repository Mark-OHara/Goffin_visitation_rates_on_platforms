# Goffin_visitation_rates_on_platforms
Datasets and R-code for the publication: "A novel platform to attract wild Tanimbar corellas (Cacatua goffiniana) for behavioural research"
Mark O’Hara, Alice M.I. Auersperg, Dewi M. Prawiradilaga, Ludwig Huber and Berenika Mioduszewska

This README.txt file was generated on 2023-08-10 by Mark O'Hara

GENERAL INFORMATION

Title of Dataset: Datasets and R-code for the publication: "A novel platform to attract wild Tanimbar corellas (Cacatua goffiniana) for behavioural research"

Author Information A. Principal Investigator Contact Information Name: Mark O'Hara Institution: Messerli Research Institute, University of Veterinary Medicine Vienna, Austria Address: Veterinärplatz 1, 1210 Vienna, Austria Email: mark.ohara@vetmeduni.ac.at

Date of data collection: 2/2019 to 3/2019

Geographic location of data collection: Tanimbar Islands, Indonesia (S7.81321°,E131.37395°)

Information about funding sources that supported the collection of the data: This research was funded in part by the Austrian Science Fund (FWF) grants awarded to M.O. (FWF project: J4169) and Alice M.I. Auersperg (FWF START project: Y01309).

SHARING/ACCESS INFORMATION

Licenses/restrictions placed on the data: For the purpose of open access, the authors have applied a CC BY public copyright license

Links to publications that cite or use the data: ...

Recommended citation for this dataset: ...

DATA & FILE OVERVIEW

File List:

Data.txt Full data set used for analysis

Weather.txt Full data set used for environmental data during visitations

Functions.R Additional used functions in "A novel platform to attract wild Tanimbar corellas (Cacatua goffiniana) for behavioural research" provided and written by Roger Mundry

Visitations.RData image file containing stored objects from the analysis

Visitation rate of Wild Goffin’s cockatoos.R Annotated R-script to carry out the full analysis and plotting of figures

METHODOLOGICAL INFORMATION

Description of methods used for collection/generation of data: Experimental design and protocols used for data collection can be found in ...

Methods for processing the data: Submitted data were generated from the raw or collected data as described in ...

DATA-SPECIFIC INFORMATION FOR: Data.tx

Number of variables: 16
Number of cases/rows: 7967
Variable List:
 $ Location   : Factor with 8 levels
 $ Bout       : Factor with 2 levels
 $ Filename   : Factor with 7967 levels
 $ Date       : Factor with 62 levels
 $ Time       : Factor wiht 753 levels
 $ Manualcount: integer
 $ datetime   : Factor with 7932 levels
 $ timemin    : integer
 $ Playback   : Factor wihth 2 levels
 $ daysbait   : integer
 $ Temp       : numerical
 $ Humidity   : integer
 $ Lux        : numerical
 $ UVI        : integer
 $ Rain       : numerical
 $ Wind       : numumerical


DATA-SPECIFIC INFORMATION FOR: Weather.txt

Number of variables: 8
Number of cases/rows: 19660
Variable List:
 $ Temp    : numerical
 $ Hum     : integer
 $ Lux     : numerical
 $ UVI     : integer
 $ Rain    : numerical
 $ Wind    : numerical
 $ timemin : integer
