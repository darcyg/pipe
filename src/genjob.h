#ifndef GENJOB_H
#define GENJOB_H

#include <string>
#include "options.h"
#include "job.h"

class GenJob{
    public:
        Options* mOpt;
        std::string lib1;
        std::string lib2;
    public:
        GenJob(Options* opt);
        ~GenJob();

        void setLib(const std::string& l1, const std::string& l2);
        void genFastpJob(Job* j);
        void genSplitrJob(const std::string& conf, Job* j);
        void genSeqtkJob(Job* j);
        void genFiltdbJob(Job* j);
        void genAlignJob(Job* j);
        void genMkdupJob(const std::string& bam, Job* j);
        void genBamqcJob(const std::string& bam, Job* j);
        void genFusionJob(Job* j);
        void genExpressJob(Job* j);
        void genReportJob(Job* j);
};

#endif
