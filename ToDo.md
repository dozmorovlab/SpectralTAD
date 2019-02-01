Note 80 characters per line limit, read history https://en.wikipedia.org/wiki/Characters_per_line

- Add an example how-to run to the help section of each function.

- Fix things noted by John

# John's review

- warning messages when installing package from github:

  Warning: C:/Users/stansfieldjc/AppData/Local/Temp/RtmpspSFFX/R.INSTALL213c52a63ae5/SpectralTAD/man/SpectralTAD.Rd:30: unexpected section header '\description'
  Warning: C:/Users/stansfieldjc/AppData/Local/Temp/RtmpspSFFX/R.INSTALL213c52a63ae5/SpectralTAD/man/SpectralTAD.Rd:33: unexpected section header '\details'
  Warning: C:/Users/stansfieldjc/AppData/Local/Temp/RtmpspSFFX/R.INSTALL213c52a63ae5/SpectralTAD/man/SpectralTAD.Rd:36: unexpected END_OF_INPUT '
  '
  Warning: C:/Users/stansfieldjc/AppData/Local/Temp/RtmpspSFFX/R.INSTALL213c52a63ae5/SpectralTAD/man/SpectralTAD_Par.Rd:34: unexpected section header '\description'
  Warning: C:/Users/stansfieldjc/AppData/Local/Temp/RtmpspSFFX/R.INSTALL213c52a63ae5/SpectralTAD/man/SpectralTAD_Par.Rd:37: unexpected section header '\details'
  Warning: C:/Users/stansfieldjc/AppData/Local/Temp/RtmpspSFFX/R.INSTALL213c52a63ae5/SpectralTAD/man/SpectralTAD_Par.Rd:40: unexpected END_OF_INPUT '
  Warning messages:
  1: In untar2(tarfile, files, list, exdir) :
    skipping pax global extended headers
  2: In untar2(tarfile, files, list, exdir) :
    skipping pax global extended headers
    
- math text isn't displaying correctly in README.md
- link to vignette in README should probably be changed to seen by `browseVignettes("SpectralTAD")` since it just links to the raw version of the vignette on github.
- in README maybe have some introduction - explain what TADs are and what levels are
+ why are there 2 license files? should only be 1
- need a line in README for installing dependancies


- vignette references aren't working

- bioconductor is going to make you add some tests; see hadley wickams test section in the package building guide

- functions need an #' @return statement in the help for the object that is returned

- need examples in the ?SpectralTAD help

- lines too long in code

- gap_threshold error in SpectralTAD_Par
