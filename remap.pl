use strict;
use warnings FATAL => qw( all );

my %remap = (
	'SNRK-AS1' => 'SNRK',
	'MST1P9' => 'MST1L',
	'FAM48B1' =>'SUPT20HL1',
	'HEATR7B2' => 'MROH2B',
	'ODZ4' => 'TENM4',
	'ODZ1' => 'TENM1',
	'USP17' => 'USP17L15',
	'CCDC165' => 'SOGA2',
	'C3orf15' => 'MAATS1',
	'C10orf140' => 'SKIDA1',
	'C19orf51' => 'DNAAF3',
	'PRIC285' => 'HELZ2',
	'FAM125B' => 'MVB12B',
	'BOD1L' => 'BOD1L1',
	'KDM4DL' => 'KDM4E',
	'FAM123B' => 'AMER1',
	'ZNF815' => 'ZNF815P',
	'LOC100192378' => 'ZFHX4-AS1',
	'HNRNPA1' => 'HNRNPA1P10',
	'CCDC13' => 'CCDC13-AS1',
	'MAGEA12' => 'CSAG4',
	'TXNDC5' => 'BLOC1S5-TXNDC5',
	'PGM5' => 'PGM5-AS1',
	'FLJ45983' => 'GATA3-AS1',
	'LOC643955' => 'ZNF733P',
	'LOC650623' => 'BEND3P3',
	'FAM59A' => 'GAREM',
	'C9orf7' => 'CACFD1',
	'KIAA0090' => 'EMC1',
	'KLHDC5' => 'KLHL42',
	'C13orf30' => 'FAM216B',
	'FAM164C' => 'ZC2HC1C',
	'OBFC2A' => 'NABP1',
	'C20orf54' => 'SLC52A3',
	'OTX2' => 'OTX2-AS1',
	'FLJ43860' => 'MROH5',
	'C17orf57' => 'EFCAB13',
	'C3orf19' => 'CCDC174',
	'C9orf68' => 'SPATA6L',
	'C3orf32' => 'SSUH2',
	'AHCTF1' => 'AHCTF1P1',
	'EIF2C1' => 'AGO1',
	'KIAA1340' => 'KLHL42',
	'C5orf44' => 'TRAPPC13',
	'C1orf96' => 'CCSAP',
	'FAM190A' => 'CCSER1',
	'NELF' => 'NSMF',
	'HEATR8' => 'MROH7',
	'MCART3P' => 'SLC25A51P1',
	'C8orf80' => 'NUGGC',
	'FAM75A6' => 'SPATA31A6',
	'LINC00598' => 'TTL',
	'CORO7-PAM16' => 'CORO7',
	
	'LOC84856' => 'LINC00839',
	'AK127292' => 'ANKRD20A19P',
	'LOC149837' => 'LINC00654',
	'LOC100132526' => 'FGD5P1',
	'LOC338739' => 'CSTF3-AS1',
	'TRAPPC2' => 'TRAPPC2P1',
	'AK302238' => 'GOLGA6L7P',
	'LOC154822' => 'LINC00689',
	'LOC126536' => 'LINC00661',
	'LOC646851' => 'FAM227A',
	'LOC649305' => 'LOC729732',
	'LOC146336' => 'SSTR5-AS1'

	# yet unmapped	
	# 'USP17L1P'
	# 'RNA45S5'
	# 'PTGER4P2-CDK2AP2P2'
);

while(<>)
{
	while ( my ($k, $v) = each(%remap) )
	{
		s/\t$k\t/\t$v\t/g;
		s/^$k\t/$v\t/g;
	}
	
	print $_;
}