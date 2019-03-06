program: fusepipe  
version: 0.0.0  
updated: 22:36:22 Mar  6 2019  
Usage: fusepipe [OPTIONS]  

|  Options                    | Explanations
|-----------------------------|---------------------------------
|  -h,--help                  | Print this help message and exit
|  -s,--slist FILE REQUIRED   | sample list file
|  -r,--regf FILE REQUIRED    | bed region file
|  -v,--vread TEXT            | fastq  subset reads number
|  -o,--out TEXT              | output directory
|  -a,--amark INT in [1 - 10] | analysis marker range
|  -i,--imark INT in [1 - 10] | initial analysis marker
|  -e,--emark INT in [1 - 10] | end analysis marker
|  -q,--queue TEXT            | queue to run tasks
|  -c,--ctd                   | continue from last failure
|  -l,--loc                   | run in localhost
|  -g,--gen                   | generate sjms, not run tasks

Installation

1. clone repo  
`git clone https://github.com/vanNul/pipe`

2. compile  
`cd pipe`
`./autogen.sh`
`./configure --prefix=/path/to/install/dir/`
`make`
`make install`

PS  
this is a sjm based fusion pipeline  
the relation between analysis marker and content is listed below  

|Marker   |Analysis                     
|---------|----------------------------
|1        | cutadaptor by fastp         
|2        | split read by splitr        
|3        | filter rrna by filtdb       
|4        | downsample fastq by seqtk   
|5        | genome slignment by bwa     
|6        | markdup by mkdp             
|7        | bam QC by bamqc             
|8        | fusion calling by fusionMap 
|9        | express quant by kallisto   
|10       | report by genrep            
