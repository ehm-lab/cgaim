# 1.0.0

## Resubmission
This is a resubmission. In this version I have:

* Added a \donttest statement for R CMD check to ignore the most computationally expensive example.
* Added the reference to the methodological paper in the DESCRIPTION file.
* Replaced the T and F that previously escaped me by TRUE and FALSE.
* Added \value fields to the documentation of exported functions missing one.

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE on winbuilder only:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Pierre Masselot <pierre.masselot@lshtm.ac.uk>'

New submission

Possibly misspelled words in DESCRIPTION:
  Campagna (10:182)
  Chebana (10:173)
  Gosselin (10:209)
  Groupwise (2:20)
  Lavigne (10:192)
  Masselot (10:163)
  Ouarda (10:201)
  groupwise (10:31, 10:238)

I believe these notes are just flagging the fact it is a new submission.
Regarding misspelled words, this is how these words are spelled and used in associated publications and are therefore correct (groupwise) and the names of authors in the citation.

## Downstream dependencies
There are currently no downstream dependencies for this package.
