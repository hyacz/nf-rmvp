.PHONY: build clean test

build:
	docker build . -t rmvp

test:
	nextflow run GWAS.nf --genotype /data/GWAS-pipeline/example/mvp.vcf --phenotype /data/GWAS-pipeline/example/mvp.phe --methods GLM,MLM

clean:
	rm -rf .nextflow/ .nextflow.log* work/ result/ report/
