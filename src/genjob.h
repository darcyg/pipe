#ifndef GENJOB_H
#define GENJOB_H

#include <string>
#include "job.h"
#include "options.h"

/** Class to store fusionmap command templates */
class FusionMapTemplate{
    public:
        std::string fq1;           ///< read1 file path
        std::string fq2;           ///< read2 file path
        std::string outdir;        ///< output path
        std::string pre;           ///< prefix of output file names
        
        FusionMapTemplate(const std::string& r1, const std::string r2, const std::string& otd, const std::string& pr){
            fq1 = r1;
            fq2 = r2;
            outdir = otd;
            pre = pr;
        }
        
        friend std::ostream& operator<<(std::ostream& os, const FusionMapTemplate& f){
            os << "<Files>\n";
            os << f.fq1 << "\n" << f.fq2 << "\n\n";
            os << "<Options>\n";
            os << "MonoPath=/share/public/software/mono/mono/bin/mono\n";
            os << "PairedEnd=True\nRnaMode=True\nThreadNumber=8\nGzip=False\nOutputFusionReads=True\n\n";
            os << "<Output>\n";
            os << "TempPath=" << f.outdir << "/" << f.pre << "_tmp\n";
            os << "OutputPath=" << f.outdir << "/" << f.pre << "\n";
            os << "OutputName=" << f.pre << "\n";
            return os;
        }
};

/** Class to generate a job */
class GenJob{
    public:
        Options* mOpt;      ///< pointer to Options object
        std::string bam;    ///< bam file
        std::string lib1;   ///< fastq read1 file path
        std::string lib2;   ///< fastq read2 file path
    public:
        /** Construct a GenJob object
         * @param opt pointer to Options object
         */
        GenJob(Options* opt);
        
        /** Destroy a GenJob object */
        ~GenJob();

        /** set library of GenJob 
         * @param l1 read1 file path
         * @param l2 read2 file path
         */
        void setLib(const std::string& l1, const std::string& l2);
        
        /** set bam of GenJob
         * @param bam bam file path
         */
        void setBam(const std::string& b);

        /** generate fastp Job
         * @param j pointer to Job
         */
        void genFastpJob(Job* j);
        
        /** generate splitr Job
         * @param conf barcode configure file path
         * @param j pointer to Job
         */
        void genSplitrJob(const std::string& conf, Job* j);
        
        /** generate seqtk Job
         * @param j pointer to Job
         */
        void genSeqtkJob(Job* j);
        
        /** generate filtdb Job
         * @param j pointer to Job
         */
        void genFiltdbJob(Job* j);
        
        /** generate alignment Job
         * @param j pointer to Job
         */ 
        void genAlignJob(Job* j);
        
        /** generate markdup Job
         * @param j pointer to Job
         */
        void genMkdupJob(Job* j);
        
        /** generate bamqc Job
         * @param j pointer to Job
         */
        void genBamqcJob(Job* j);
        
        /** generate fusionmap Job
         * @param j pointer to Job
         */
        void genFusionJob(Job* j);
        
        /** generate kallisto Job
         * @param j pointer to Job
         */
        void genExpressJob(Job* j);
        
        /** generate genrpt Job
         * @param j pointer to Job
         */
        void genReportJob(Job* j);
        
        /** get minimum reads number of fastq in split libraries from a common library
         * @return minimum reads number of fastq in split libraries
         */
        size_t getMinFqVol();
};

#endif
