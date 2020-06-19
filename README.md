# nf-rmvp

Nextflow version of [rMVP](https://github.com/xiaolei-lab/rMVP/)

- [nf-rmvp](#nf-rmvp)
  - [INPUT](#input)
  - [EXAPMLE](#exapmle)

## INPUT

| Param | Default | Notes |
| --- | --- | --- |
|genotype  | file(params.genotype)  | VCF Genotype file (no missing)  
|phenotype | file(params.phenotype) | Phenotype file, sep is '\t'
|methods   | "FarmCPU" | Support GLM, MLM and FarmCPU models, multiple models are separated by commas(eg. "GLM,MLM")
|ncpus     | "auto" | Number of CPU cores used, the default is the maximum number of cores - 1
|threshold | "0.05" | Bonfferoni correction threshold, 0.05 means the threshold line is -log10(0.05 / number of markers)
|impute    | "TRUE" | Whether to automatically imputation, if not missing will automatically skip(WARN: only applicable when there is very little data missing)
|outdir    | "result" | Output file directory
|plot_type | "jpg" | The format of the exported picture can be "jpg", "tiff" or "pdf"

## EXAPMLE

```shell
> make test
nextflow run GWAS.nf --genotype example/mvp.vcf --phenotype example/mvp.phe --methods GLM,MLM
N E X T F L O W  ~  version 20.04.1
Launching `GWAS.nf` [extravagant_rubens] - revision: 2c541b76a6
=============================================
G W A S - N F v0.1
=============================================
genotype  : /Users/hyacz/Code/GWAS-pipeline/example/mvp.vcf
phenotype : /Users/hyacz/Code/GWAS-pipeline/example/mvp.phe
methods   : DataflowQueue(queue=[DataflowVariable(value=GLM), DataflowVariable(value=MLM), DataflowVariable(value=groovyx.gpars.dataflow.operator.PoisonPill@20a7953c)])
ncpus     : parallel::detectCores() - 1
threshold : 0.05
outdir    : result
executor >  local (19)
[15/e04ed9] process > DataFormatConvert           [100%] 1 of 1 ✔
[dd/3224e8] process > RunAnalysis (4)             [100%] 4 of 4 ✔
[e6/ad9aaa] process > FilterSignalSNP (4)         [100%] 4 of 4 ✔
[28/014783] process > PlotSnpDensity              [100%] 1 of 1 ✔
[39/b3707d] process > plot_phenotype_distribution [100%] 1 of 1 ✔
[1c/6c942d] process > PlotQQ (4)                  [100%] 4 of 4 ✔
[54/11fd83] process > PlotManhattan (4)           [100%] 4 of 4 ✔


> tree result
result
├── MVP.SNP_Density.jpg
├── Phe_Distribution.T1.jpg
├── Phe_Distribution.T2.jpg
├── QQplot.T1.GLM.jpg
├── QQplot.T1.MLM.jpg
├── QQplot.T2.GLM.jpg
├── QQplot.T2.MLM.jpg
├── Rectangular-Manhattan.T1.GLM.jpg
├── Rectangular-Manhattan.T1.MLM.jpg
├── Rectangular-Manhattan.T2.GLM.jpg
├── Rectangular-Manhattan.T2.MLM.jpg
├── T1.GLM.csv
├── T1.GLM_signals.csv
├── T1.MLM.csv
├── T1.MLM_signals.csv
├── T2.GLM.csv
├── T2.GLM_signals.csv
├── T2.MLM.csv
└── T2.MLM_signals.csv
```
