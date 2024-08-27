These are our solution capacitance programs.

Usage:
        [Calculate]
        cd NETWORK/bin
        ./solver -i ../Cases/pattern/file -o path/to/save/your/file

        [Calculate example]
        cd NETWORK/bin
        ./solver -i ../Cases/input_diag/BEM_INPUT_11_386073.txt -o a.txt

        [Comparison]
        cd COMPARISON
        ./check.sh
        cd result

NETWORK is used to solve capacitance, the specific method of use is to pass in an input file, and an output file.
The calculation results will be automatically stored in the output file, using the c++17 standard.

COMPARISON is used to compare the results of our program calculation with the field solver solution results provided by the competition.
Using the C++11 standard.
The specific operation is to use the results of our program and the formula given by the competition party to solve the error file of each table file,
including the average error, the maximum error, the minimum error, the standard error, 2sigma region and the error of each item,
and the results are stored in the COMPARISON/result/eror directory

OTHER_METHODS is used to store other methods that have been tried. Each of these methods has its own advantages and disadvantages. 
Some methods show significant potential, but due to time constraints, we have only implemented a part of them. 
Nonetheless, they remain highly relevant and can be used as reference for future solutions.
