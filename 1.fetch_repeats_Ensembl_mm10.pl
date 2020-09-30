#!/usr/bin/env perl

use warnings;

# Simply obtain Ensembl Repeat locations

use strict;
#use lib "/homes/nenbar/local/src/ensembl/modules";
#use lib "/homes/nenbar/local/src/bioperl-live";
use Bio::EnsEMBL::Registry;

$| = 1;

my $ensembl    = 89;
my $species    = "mouse";
my $slice_size = 1000000;   # max slice size to process at a time
my $chr_rename    = 0;  # if true, add 'chr' to numeric chromosome names
my $strand_rename = 1;  # if true, change 1/-1 into +/- (other -> *)

my $out_path = "API_extracted";
mkdir $out_path unless -e $out_path;


# Initialize EnsEMBL Registry
my $registry = 'Bio::EnsEMBL::Registry';
$registry ->load_registry_from_url("mysql://anonymous\@ensembldb.ensembl.org/$ensembl");

my $sa_core = $registry->get_adaptor( $species, 'Core', 'Slice' );
my $ra_core = $registry->get_adaptor( $species, 'Core', 'RepeatFeature' );

my @chrs = @{$sa_core -> fetch_all('toplevel')};

####### Repeats ########
print "Getting Repeats from Ensembl...\n";

my $out_file = "$out_path/$species-$ensembl.repeats.tab.gz";
if (-s $out_file) {
   warn "$out_file already exists, will be over-written!\n";
}

open (OUT, "| gzip -9 > $out_file") || die "Can't gzip $out_file!\n";
print OUT "chr\tstart\tend\tstrand\ttype\tanalysis\tscore\tclass\tclass_desc\n";

my $rep_count = 0;
my %types;
foreach my $chr_slice (@chrs) {
   unless ($chr_slice) {
      warn "Can't fetch next slice!\n";
      next;
   }
   my $chr_len = $chr_slice -> seq_region_length();
	my $num_chr = $chr_slice -> seq_region_name();
#   next unless $num_chr eq "Y";
	print "Chromosome $num_chr (length = $chr_len)\n";

   # Prepare for Ensembl extraction   
   my $chr = $num_chr;
   if ($chr_rename) {
      $chr = "chr$chr" if $chr =~ /^\d+/;  # In case simple number
   }
   
   my %used;
   # Get all repeats
   for (my $slice_start = 1; $slice_start < $chr_len; $slice_start += $slice_size) {
      my %others;
      
      my $slice_end = $slice_start + $slice_size;
      $slice_end = $chr_len if $slice_end > $chr_len;

      my $sub_slice = $chr_slice -> sub_Slice($slice_start,$slice_end);
      next unless $sub_slice; # It will sometimes return undef
      
      my @repeats = @{ $ra_core-> fetch_all_by_Slice($sub_slice) };
      while (my $repeat = shift @repeats) {
         my $analysis = $repeat -> analysis()->logic_name();
         $types{$analysis}++;
#         next unless $analysis eq "tRNAscan";
         my $strand = $repeat -> strand || "NA";
         my $score = $repeat -> score() || "NA";
         my $start = $repeat -> start() || "NA";
         my $end   = $repeat -> end() || "NA";
         next if ($start eq "NA" or $end eq "NA");
         $start += $slice_start-1;
         $end   += $slice_start-1;
         my $disp_id = $repeat -> display_id() || "NA";
         
         # Consensus has more information?
         my $rep_cons = $repeat -> repeat_consensus() || "NA";
#         my $r_c_name = $rep_cons -> name() || "NA";
#         my $r_c_desc = $rep_cons -> desc() || "NA";
#         my $r_c_consensus = $rep_cons -> repeat_consensus() || "NA";
         my $r_c_class = $rep_cons -> repeat_class() || "NA";
         my $r_c_type = $rep_cons -> repeat_type() || "NA";
         

         # Strand is always more useful as +/-?
         if ($strand_rename) {
            if ($strand eq 1) {
               $strand = "+";
            } elsif ($strand eq -1) {
               $strand = "-";
            } else {
               $strand = "*";
            }
         }

         my $line = "$chr\t$start\t$end\t$strand\t$disp_id\t$analysis\t$score\t$r_c_class\t$r_c_type";
         next if $used{$line};
         $rep_count++;
         
#         print "$line";<STDIN>;
         print OUT "$line\n";
         $used{$line} = 1;
      }
      
   }
   undef %used;
   
}
close OUT;

foreach my $analysis (sort keys %types) {
   print "$analysis -> $types{$analysis}\n";
}
print "Total: $rep_count\n";
