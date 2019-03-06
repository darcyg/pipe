program: fusepipe  
version: 0.0.0  
updated: 22:36:22 Mar  6 2019  
Usage: fusepipe [OPTIONS]  

|Options                         | Explanations
|--------------------------------|---------------------------------
| -h,--help                      |Print this help message and exit
| -s,--sample_list FILE REQUIRED |sample list file
| -r,--region_file FILE REQUIRED |region bed file
| -v,--reads_count TEXT          |down sample fastq reads number
| -o,--out_dir TEXT              |output directory
| -a,--ana_mark INT in [1 - 10]  |...analysis stage marker
| -i,--ini_mark INT in [1 - 10]  |initial analysis stage marker
| -e,--end_mark INT in [1 - 10]  |end analysis stage marker
| -c,--continue                  |continue run from last failure
| -l,--loc                       |Running in localhost
| -g,--gensjm                    |Only generate sjms, not run pipeline

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
this is a sjm based fusion pipeline~