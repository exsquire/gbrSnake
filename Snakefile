singularity: "gbrs.simg"
import os
import pandas as pd
from snakemake.utils import validate
##### load config and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="config/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"]).set_index("run", drop = False)
validate(samples, "config/schemas/samples.schema.yaml")
tags = list(set(samples.tag))
#For unique transprob files after pooling
tagSamples = samples.drop_duplicates(subset='tag')
##### Rules #############################

rule all:
	input:
		expand([
		"bam/{run}/{run}.bam",
		"emase/{run}/{run}.h5",
		"output/emase_compress/{run}.h5c"],
		run = samples.run.tolist()),
		expand([
		"output/pooled_emase/{tag}.pooled",
		"output/multiway/{tag}/gbrs.quantified.multiway.genes.alignment_counts",
                "output/multiway/{tag}/gbrs.quantified.multiway.genes.expected_read_counts",
                "output/multiway/{tag}/gbrs.quantified.multiway.genes.tpm",
                "output/multiway/{tag}/gbrs.quantified.multiway.isoforms.alignment_counts",
                "output/multiway/{tag}/gbrs.quantified.multiway.isoforms.expected_read_counts",
		"output/reconstruct/{tag}/gbrs.reconstructed.genoprobs.npz",
		"output/reconstruct/{tag}/gbrs.reconstructed.genotypes.npz",
		"output/reconstruct/{tag}/gbrs.reconstructed.genotypes.tsv",
                "output/diploid/{tag}/gbrs.quantified.diploid.genes.alignment_counts",
		"output/diploid/{tag}/gbrs.quantified.diploid.genes.expected_read_counts",
		"output/diploid/{tag}/gbrs.quantified.diploid.genes.tpm",
		"output/diploid/{tag}/gbrs.quantified.diploid.isoforms.alignment_counts",
		"output/diploid/{tag}/gbrs.quantified.diploid.isoforms.expected_read_counts",
                "output/diploid/{tag}/gbrs.quantified.diploid.isoforms.tpm"],
		tag = tags)

rule align_fastq:
	input:
		"input/{run}.fastq.gz"
	output:
		temp("bam/{run}/{run}.bam")
	priority: 7
	shell:
		"set +o pipefail; "
		"export PATH=/opt/conda/envs/gbrs/bin:$PATH; "
		"bowtie -q -a --best --strata --sam -v 3 /usr/share/gbrs/R84-REL1505/transcripts "
		"{input} | samtools view -bS - > {output}"

rule bam_to_emase:
	input:
		"bam/{run}/{run}.bam"
	output:
		temp("emase/{run}/{run}.h5")
	priority: 6
	shell:
		"set +o pipefail; "
		"export PATH=/opt/conda/envs/gbrs/bin:$PATH; "
		"gbrs bam2emase -i {input} -m /usr/share/gbrs/R84-REL1505/ref.transcripts.info "
		"-s A,B,C,D,E,F,G,H -o {output}"

rule compress_emase:
	input:
		"emase/{run}/{run}.h5"
	output:
		"output/emase_compress/{run}.h5c"
	priority: 5
	shell:
		"set +o pipefail; "
		". /opt/conda/etc/profile.d/conda.sh; "
		"conda activate gbrs; "
		"gbrs compress -i {input} -o {output}"

rule pool_emase:
	input:
		lambda wildcards: list("output/emase_compress/"+samples.run.loc[(samples.tag == wildcards.tag)]+".h5c")
	output:
		"output/pooled_emase/{tag}.pooled"
	priority: 4
	params:
		pool=lambda wildcards, input: ",".join(input)
	shell:
		"set +o pipefail; "
		". /opt/conda/etc/profile.d/conda.sh; "
		"conda activate gbrs; "
		"gbrs compress -i {params.pool} -o {output}"

rule quantify_multiway:
	input:
		emase = "output/pooled_emase/{tag}.pooled"
	output:
		"output/multiway/{tag}/gbrs.quantified.multiway.genes.alignment_counts",
		"output/multiway/{tag}/gbrs.quantified.multiway.genes.expected_read_counts",
		"output/multiway/{tag}/gbrs.quantified.multiway.genes.tpm",
		"output/multiway/{tag}/gbrs.quantified.multiway.isoforms.alignment_counts",
		"output/multiway/{tag}/gbrs.quantified.multiway.isoforms.expected_read_counts",
		"output/multiway/{tag}/gbrs.quantified.multiway.isoforms.tpm"
	priority: 3
	params:
		g2tRef = "/usr/share/gbrs/R84-REL1505/ref.gene2transcripts.tsv",
		hybTarg = "/usr/share/gbrs/R84-REL1505/gbrs.hybridized.targets.info",
		method = "2"
	shell:
		"set +o pipefail; "
		". /opt/conda/etc/profile.d/conda.sh; "
		"conda activate gbrs; "
		"gbrs quantify -i {input.emase} -g {params.g2tRef} "
		"-L {params.hybTarg} -M {params.method} "
		"--report-alignment-counts; "
		"mv gbrs.quantified.multiway* output/multiway/{wildcards.tag}/"

rule reconstruct_genome:
	input:
		counts = "output/multiway/{tag}/gbrs.quantified.multiway.genes.tpm"
	output:
		"output/reconstruct/{tag}/gbrs.reconstructed.genoprobs.npz",
		"output/reconstruct/{tag}/gbrs.reconstructed.genotypes.npz",
		"output/reconstruct/{tag}/gbrs.reconstructed.genotypes.tsv"
	priority: 2
	params:
		transprob = lambda wildcards: tagSamples.loc[(tagSamples.tag == wildcards.tag),"transprob"].item(),
		avecs = "/usr/share/gbrs/R84-REL1505/avecs.npz",
		genPos = "/usr/share/gbrs/R84-REL1505/ref.gene_pos.ordered.npz"
	shell:
		"set +o pipefail; "
		". /opt/conda/etc/profile.d/conda.sh; "
		"conda activate gbrs; "
		"echo {wildcards.tag} {input.counts} {params.transprob} > "
		"gbrs.reconstructed.transprob.log; "
		"gbrs reconstruct -e {input.counts} -t /usr/share/gbrs/R84-REL1505/{params.transprob} "
		"-x {params.avecs} "
		"-g {params.genPos}; "
		"mv gbrs.reconstructed* output/reconstruct/{wildcards.tag}/"

rule quantify_diploid:
	input:
		emase = "output/pooled_emase/{tag}.pooled",
		geno = "output/reconstruct/{tag}/gbrs.reconstructed.genotypes.tsv"
	output:
		"output/diploid/{tag}/gbrs.quantified.diploid.genes.alignment_counts",
		"output/diploid/{tag}/gbrs.quantified.diploid.genes.expected_read_counts",
		"output/diploid/{tag}/gbrs.quantified.diploid.genes.tpm",
		"output/diploid/{tag}/gbrs.quantified.diploid.isoforms.alignment_counts",
		"output/diploid/{tag}/gbrs.quantified.diploid.isoforms.expected_read_counts",
		"output/diploid/{tag}/gbrs.quantified.diploid.isoforms.tpm"
	priority: 1
	params:
		g2tRef = "/usr/share/gbrs/R84-REL1505/ref.gene2transcripts.tsv",
		hybTarg = "/usr/share/gbrs/R84-REL1505/gbrs.hybridized.targets.info",
		method = "2"
	shell:
		"set +o pipefail; "
		". /opt/conda/etc/profile.d/conda.sh; "
		"conda activate gbrs; "
		"gbrs quantify -i {input.emase} -G {input.geno} -g {params.g2tRef} "
		"-L {params.hybTarg} -M {params.method} --report-alignment-counts; "
		"mv gbrs.quantified.diploid* output/diploid/{wildcards.tag}/"
