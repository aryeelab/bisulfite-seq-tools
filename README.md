# methylation-wdl

### Build docker image
cd Docker_Build
docker build -t methy .

### Run WDL in Cromwell
#### Note: On Macs this can be installed with `brew install cromwell`
cromwell run methylation_pipeline.wdl methylation_pipeline_testsample.json 
cromwell run methylation_pipeline.wdl methylation_pipeline_testsample2.json 


cromwell run methylation_pipeline_aggregation.wdl methylation_pipeline_aggregation.json 


/Users/maryee/Dropbox/projects/methylation-wdl/cromwell-executions/methpipe/2c284299-5803-4854-ba10-83d69de61352/call-create_rda/execution/testsample_01.rda
/Users/maryee/Dropbox/projects/methylation-wdl/cromwell-executions/methpipe/4d93cd3e-cdd4-472c-8f45-322369ba4f8e/call-create_rda/execution/testsample_02.rda