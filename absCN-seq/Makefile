export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash

all: 592_rem_rel.copynumber.called.segmented.absCNseq.cn 592_rem_rel.copynumber.called.segmented.absCNseq.snv

592_rem_rel.copynumber.called: /data/knoll/VarScan2/592_rem_rel/592_rem_rel.copynumber
	java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar copyCaller \
		/data/knoll/VarScan2/592_rem_rel/592_rem_rel.copynumber \
		--output-file 592_rem_rel.copynumber.called \
		--recenter-down 0.54 \
		--min-coverage 30
		
592_rem_rel.copynumber.called.segmented.absCNseq.cn: 592_rem_rel.copynumber.called /mnt/projects/hdall/scripts/absCN-seq/segment.R
	Rscript /mnt/projects/hdall/scripts/absCN-seq/segment.R
	mv 592_rem_rel.copynumber.called.segmented.part 592_rem_rel.copynumber.called.segmented
	mv $@.part $@
	
592_rem_rel.copynumber.called.segmented.absCNseq.snv: /data/knoll/VarScan2/592_rem_rel/592_rem_rel_somatic_snp.snp
	cat $< | perl -lane '$$F[0] =~ s/chr([\dM])/$$1/; print "$$F[0]\t$$F[1]\t".($$F[10]=~/([\d\.]+)/ ? $$1/100 : $$F[10]) if ($$F[12] eq "somatic_status" or ($$F[12] eq "Somatic" and $$F[14] < 1e-5))' > $@.part
	mv $@.part $@
	
# r <- run.absCNSeq(seg.fn="/mnt/projects/hdall/results/absCN-seq/592_rem_rel.copynumber.called.segmented.absCNseq.cn", snv.fn="/mnt/projects/hdall/results/absCN-seq/592_rem_rel.copynumber.called.segmented.absCNseq.snv", res.dir="/mnt/projects/hdall/results/absCN-seq", smp.name="592_rel", seq.type="WES")
