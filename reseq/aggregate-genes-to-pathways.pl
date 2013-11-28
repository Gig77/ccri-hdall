use warnings FATAL => qw( all );
use strict;
use Getopt::Long;

my ($level);
GetOptions
(
	"level=s" => \$level # coarse or fine 
);

die if ($level ne 'coarse' and $level ne 'fine');

while (<>)
{
	# skip header lines
	if (/^#/)
	{
		print $_;
		next;
	}

	if ($level eq "fine")
	{
		s/^(PAX5|PRDM5|PRDM9|PRDM13|PRDM14|PRDM6|PRDM8|ETV1|PAX7|OLIG2|BCL3|HSPA2|ZNF516|TLX3|E2F8|ETV6|IRF2BPL|IKZF1|IKZF2)\t/TF\t/;		
		s/^(SETDB1|MLL2|TET2|MLL3|WHSC1|KDM5A)\t/methylation\t/;
		s/^(ATM|TRRAP|CREBBP)\t/acetylation\t/;
		s/^(FBXL7|NDC80|NUMA1|CDC27|RB1)\t/cellcycle\t/;
		s/^(MYOM2|MYH2|OBSL1|NEFH|PDE4DIP|MYO19)\t/cytoskeleton\t/;
		s/^(ERBB4|ATXN2)\t/signaling_EGFR\t/;
		s/^(PCDHA4|GPR98|TBXA2R|UBXN11)\t/signaling_GPCR\t/;
		s/^(MAP3K3|FN1|GNAQ|IGF2BP3)\t/signaling_MAPK\t/;
		s/^(IKBKB|ZBP1)\t/signaling_NFkB\t/;
		s/^(NOTCH1|GATA1|MESP2)\t/signaling_NOTCH\t/;
		s/^(PIK3CB|RPS6KB2|RICTOR|SYK|PTCH1|MAGI1|USP9X|MTMR10|MMP9|PI4KA)\t/signaling_OTHER\t/;
		s/^(PTPN11|KRAS|NRAS|RAPGEF2|FLT3|NF1|PTEN)\t/signaling_RAS\t/;
		s/^(OBSCN|KNDC1|CDC42EP1)\t/signaling_RHO\t/;
		s/^(KREMEN1|CTNNB1|MCC)\t/signaling_WNT\t/;
		s/^(CRLF2|IL7R|JAK1|JAK2|JAK3)\t/signaling_JAK_STAT\t/;
		s/^(PRPS1L1|CTBS|CYP11B1)\t/metabolism\t/;
		s/^(FANCD2|ERCC4|TP53)\t/dnarepair\t/;
		s/^(DNTT|GCNT2|GDPD2)\t/enzymes\t/;
		s/^(FRG1|SF3B1)\t/spliceosome\t/;
	}
	else
	{
		s/^(ERBB4|ATXN2|PCDHA4|GPR98|TBXA2R|UBXN11|MAP3K3|FN1|GNAQ|IGF2BP3|IKBKB|ZBP1|NOTCH1|GATA1|MESP2|PIK3CB|RPS6KB2|RICTOR|SYK|PTCH1|MAGI1|USP9X|MTMR10|MMP9|PI4KA|PTEN|PTPN11|KRAS|NRAS|RAPGEF2|FLT3|NF1|OBSCN|KNDC1|CDC42EP1|KREMEN1|CTNNB1|MCC|CRLF2|IL7R|JAK1|JAK2|JAK3)\t/signaling\t/;
		s/^(PAX5|PRDM5|PRDM9|PRDM13|PRDM14|PRDM6|PRDM8|ETV1|PAX7|OLIG2|BCL3|HSPA2|ZNF516|TLX3|E2F8|ETV6|IRF2BPL|IKZF1|IKZF2)\t/TF\t/;		
		s/^(ATM|TRRAP|HIST1H3F|SETDB1|MLL2|TET2|MLL3|WHSC1|KDM5A|CREBBP)\t/chromatin\t/;
		s/^(MYOM2|MYH2|OBSL1|NEFH|PDE4DIP|MYO19)\t/cytoskeleton\t/;
		s/^(FBXL7|NDC80|NUMA1|CDC27|RB1)\t/cellcycle\t/;
		s/^(PRPS1L1|CTBS|CYP11B1)\t/metabolism\t/;
		s/^(FANCD2|ERCC4|TP53)\t/dnarepair\t/;
		s/^(DNTT|GCNT2|GDPD2)\t/enzymes\t/;
		s/^(FRG1|SF3B1)\t/spliceosome\t/;
	}
	
	print $_;		
}
