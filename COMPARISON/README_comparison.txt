This is our comparison program, 
which will sort out and merge the solved results in "/NETWORK/result" into "/COMPARISON/result", 
and calculate the average error of each pattern according to the error formula given by the competition, 
the maximum error, the minimum error, the standard error , 2sigma region and the error of each item, 
and the results are stored in the "/COMPARISON/result/error" directory. 
Of course, you can also view the calculated table file and the actual table file separately.

There is a clean method in the makefile, which can clear unnecessary files, 
and if you want to recalculate other files, you can enter the "make" command or ./check.sh.

"/COMPARISON/bin" directory is used to store executable files (script calls, no active calls).
"/COMPARISON/Case" directory is used to store the data provided by the Matchmaker.
"/COMPARISON/result" is used to store error files, real tabor files and calculated table files.
"/COMPARISON/src" directory is used to store the source files.
"/COMPARISON/tmp" is used to store temporary documents.