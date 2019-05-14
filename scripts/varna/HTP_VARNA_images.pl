#!usr/bin/perl -w

# Read Vienna style RNAfold bracket structures and feed them to VARNA for 2D structure generation
# MUST HAVE VARNA FILE IN WORKING DIRECTORY: can be downloaded here: http://varna.lri.fr/index.php?lang=en&page=downloads&css=varna
# Usage: perl HTP_VARNA_images.pl infile.txt > VARNA_cmdline.txt

$infile1 = $ARGV[0];

open (INFILE1, "$infile1") || die "Can't open the infile!\n";
my @File = <INFILE1>;
close(INFILE1) || die "can't close file";
my $N = @File;

for ($i = 0; $i < $N; $i += 3) {
    my $Name = $File[$i];
	chomp $Name;
	$Name =~ s/\R//;
	$Name =~ s/>//;
    my $Sequence = $File[$i + 1];
	chomp $Sequence;
    my $StructureEnergy = $File[$i + 2];
	my @StrEn = split (/\s/, $StructureEnergy); 
	my $Structure = $StrEn[0];
	chomp $Structure;
	
	#Generate VARNA input
    my $VARNA = 'java -cp /mnt/research/germs/shane/transActRNA/scripts/varna/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd  -algorithm radiate -baseOutline "#FFFFFF" -baseInner "#FFFFFF" -sequenceDBN "' . $Sequence . '" -structureDBN "' . $Structure . '" -o "' . $Name . '.jpg"';
	
	#Print the VARNA command line to faciliate manipulations
	print "$VARNA\n";
	
    # Run VARNA script to generate figures
    system $VARNA;
}
