## Test environments

* local Windows 11, R 4.2.3
* win-builder.r-project.org, Windows Server 2022 x64, R 4.3.0 beta
* local ubuntu 22.04.2 LTS , R 4.2.3

## R CMD check results:

* There were no ERRORs or WARNINGs or NOTEs 

## Downstream dependencies

There are currently no downstream dependencies.

## Additional comment:

Used higher tolerance in expect_equal() in the tests that failed on M1mac. 
