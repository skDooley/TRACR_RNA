#!usr/bin/perl -w

# Read FASTA file and fold sequences.
# Usage: perl HTP_RNAfold.pl infile.fasta > outfile.txt

$infile1 = $ARGV[0];
open (INFILE1, "$infile1") || die "Can't open the infile!\n";
my @BigFasta = <INFILE1>;
close(INFILE1) || die "can't close file";
my $N = @BigFasta;

#Put FASTA seqs into a Hash

my %FASTA = ();
my $CurrentSeq = "";
my @SeqNames = ();
for ($i = 0; $i < $N; $i++) {
    $Line = $BigFasta[$i];
    chomp $Line;
    #Place title and seq into their respective places in the array FASTA
    if (substr($Line, 0, 1) eq ">") {
        my $SeqName = $Line;
        chomp $SeqName;
        $SeqName =~ s/>//g;
        $SeqName =~ s/:/-/g;
        $CurrentSeq = $SeqName;
		push (@SeqNames, $SeqName);
        }

    if ( $Line =~ m/(A|G|C|U|T|Y|R|K|M|B|D|H|V|N|S|W)/g && substr($Line, 0, 1) ne ">") {
    $FASTA{$CurrentSeq} .= $Line;
    }
}

foreach my $SequenceName (@SeqNames) {
    chomp $SequenceName; 
    print ">$SequenceName\n";	
    my $Sequence = $FASTA{$SequenceName};
    chomp $Sequence;
    $Sequence =~ s/\s+//g;
    $Sequence =~ s/-//g;
    $Sequence = uc $Sequence;
    my @Out = `echo "$Sequence" | RNAfold `;
	print @Out;
    }
    
`rm *.ps`;
