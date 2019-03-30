#include <string>

/** struct to store input/output diecttory Options */
struct IODirectoryOptions{
    std::string out_dir = "./";           ///< output directory of analysis
    std::string bin_dir = "./";           ///< binary directory of pipeline
    std::string db_dir = "./";            ///< database directory of pipeline
    std::string sjm_dir = "00.sjm";       ///< directory to store sjm control files
    std::string cut_dir = "01.splitread"; ///< directory to store library splitted
    std::string spl_dir = "02.cutadapter";///< directory to store library with adapter cut
    std::string fil_dir = "03.filtdb";    ///< directory to store library filtered by databases such as rRNA
    std::string dfq_dir = "04.seqtk";     ///< directory to store library downsampled to predefined volume
    std::string aln_dir = "05.alignment"; ///< directory to store alignment results of library
    std::string mkd_dir = "06.markdup";   ///< directory to store markdup results of library
    std::string bqc_dir = "07.bamqc";     ///< directory to store alignment qc results of library
    std::string fus_dir = "08.fusion";    ///< directory to store fusion analysis results of library
    std::string exp_dir = "09.express";   ///< directory to store expression quantification of library
    std::string rep_dir = "report";       ///< directory to store final analysis report of pipeline
    std::string log_dir = "log";          ///< directory to store analysis success/fail status of pipeline
};

/** struct to store pipeline execution control options */
struct PipeControlOptions{
    const int minstage = 1;      ///< minimum stages this pipeline supports
    const int maxstage = 10;     ///< maximum stages this pipeline supports
    std::string sample_list;     ///< sample list path
    std::string ref;             ///< reference fasta file path
    std::string reg;             ///< region bed file path
    std::string gset;            ///< gene set tsv file path
    std::string queue = "";      ///< sge queue name
    std::vector<int> ana_marker; ///< analysis marker number
    int ini_marker = 1;          ///< starting analysis stage marker
    int end_marker = 10;         ///< ending analysis stage marker
    std::string dfq_vol = "VOL"; ///< fastq reads number to downsample to
    bool rerun = false;          ///< resume running from last failing if true
    bool local = false;          ///< running pipeline in local server if true
    bool gensjm = false;         ///< only generate sjm control files without running pipeline if true
    bool update = false;         ///< update sjm control file before running pipeline if true
};

/** struct to store various  options of pipeline */
struct Options{
    IODirectoryOptions ioOpt;         ///< IODirectoryOptions object
    PipeControlOptions clOpt;         ///< PipeControlOptions object
    std::string version = "0.0.0";    ///< pipeline version
    void updateOptions();
};
