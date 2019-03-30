#include <string>
#include <vector>
#include "task.h"

/** struct to store information of a Task */
struct runFile{
    std::string sjmCMD;       ///< sjm execution command
    std::string failMarkFile; ///< file to touch if this task failed
    std::string goodMarkFile; ///< file to touch if this task succeeded
    std::string logFile;      ///< logging file of sjm execution of this task
    int retValue = 0;         ///< return value of this execution 
};

/** Pipeline class hold a series of Tasks */
class Pipeline{
    public:
        std::vector<std::vector<std::vector<runFile*>>> pipelist; ///< stages in pipeline[i] will execute before pipeline[i+1]\n
                                                                  ///< tasks in pipeline[i][j] will execute before pipeline[i][j+1]\n
                                                                  ///< all tasks in pipelist[i][j] will execute parallely 
        std::string failMarkFile;                                 ///< file to touch if this pipeline failed
        std::string goodMarkFile;                                 ///< file to touch if this pipeline succeeded
        int retValue = 0;                                         ///< return value of this pipeline
        bool forceUpdateSJM = false;                              ///< force update sjm file of Task if true

    public:
        /** default constructor of Pipeline */
        Pipeline() = default;

        /** construct a Pipeline which have n stages
         * @param n stages of Pipeline
         * @param fmkf fail marker file of Pipeline
         * @param smkf success marker file of Pipeline
         */
        Pipeline(int n, const std::string& fmkf, const std::string& smkf);

        /** destroy Pipeline */
        ~Pipeline();

        /** add a runFile to pipeline
         * @param r pointer to a runFile
         * @param s stage number of pipeline to add this runFile to
         * @param t substage number of pipeline to add this runFile to
         */
        void addRunFile(runFile* r, int s, int t);

        /** prepare Pipeline to execute */
        void preparePipeline();
        
        /** runing a Task
         * @param r pointer to runFile object
         */
        int runTask(runFile* r);

        /** running a stage of a pipeline
         * @param s stage of pipeline to run
         */
        int runStage(int s);
        
        /** execute Pipeline */
        int runPipeline();
};
