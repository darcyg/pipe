program: fusepipe  
version: 0.0.0  
updated: 15:28:35 Jun  4 2019  
Usage: fusepipe [OPTIONS]  

|  Options                    | Explanations
|-----------------------------|---------------------------------
|  -h,--help                  | Print this help message and exit
|  -s,--slist FILE REQUIRED   | sample list file
|  -r,--ref FILE REQUIRED     | reference file
|  -b,--bed FILE REQUIRED     | bed region file
|  -t,--gset FILE             | gene list
|  -v,--vread TEXT            | fastq  subset reads number
|  -o,--out TEXT              | output directory
|  -a,--amark INT in [1 - 11] | analysis marker range
|  -i,--imark INT in [1 - 11] | initial analysis marker
|  -e,--emark INT in [1 - 11] | end analysis marker
|  -q,--queue TEXT            | queue to run tasks
|  -c,--ctd                   | continue from last failure
|  -l,--loc                   | run in localhost
|  -g,--gen                   | generate sjms, not run tasks
|  -u,--update Needs: --ctd   | update command to execute
|  -n,--noclean               | not cleanup intermediate files
Installation

1. clone repo  
`git clone https://github.com/vanNul/pipe`

2. compile  
`cd pipe`  
`./autogen.sh`  
`./configure --prefix=/path/to/install/dir/`  
`make`  
`make install`  

3. execute  
`/path/to/install/dir/bin/fusepipe`  

PS  
this is a sjm based fusion pipeline  

|Marker|Analysis         |Software |Version
|------|-----------------|---------|----------- 
|1     |split read       |spliter  |0.0.0       
|2     |cutadapter and qc|fqtool   |0.0.0       
|3     |filter rrna      |filter   |0.0.0       
|4     |downsample fastq |seqtk    |1.3-r106    
|5     |genome slignment |bwa      |0.7.17-r1188
|6     |markdup          |duplexer |0.0.0       
|7     |bam QC           |bamqc    |0.0.0       
|8     |fusion calling   |fusionMap|10.0.1.29   
|9     |express quant    |kallisto |0.45.1      
|10    |report           |genrpt   |0.0.0       
|11    |cleanup          |rm       |8.4         
