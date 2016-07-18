#!/usr/bin/perl

use strict;
use Math::Round;


if ($#ARGV < 1) {
   print "parameter mismatch\nTo run type this command:\nperl $0 fastahack reference input_pos_file output_file human_gff_file\n\n";

   print " first argument = full path to fastahack\n"; 
   print " second argument = full path to reference genome\n"; 
   print " third argument = input file with arbitrary number of columns, but 1st col=chromosome name and 2nd col=position\n"; 
   print " fourth argument = output file with three columns: chromosome name, position of the center nucleotide, and the thre-nucleotide context for that position\n";
   print " fifth argument = full path to human gff file\n\n\n"; 
   exit 1;
}


my $Fastahack=$ARGV[0];
my $Reference=$ARGV[1];
open(InputPositions,             '<', $ARGV[2]) || die("Could not open file!");
open(OutputTrinucleotideContext, '>', $ARGV[3]) || die("Could not open file!");
open(HumanGFF,                   '<', $ARGV[4]) || die("Could not open file!");



################ read in one coordinate at a time and execute fastahack on it

# reading the header
my $head = <InputPositions>;
$head =~ s/\n|\r//;
print OutputTrinucleotideContext "$head\tContext\n";
my $gffHead = <HumanGFF>;
chomp $gffHead;

# creating trinucleotide context data hash, insertion and deletion counts
my %trinucleotide_context_data;
my %context_tally_across_mutated_to;
my %gff_hash;
my $gffMatch;
my %location;
# my %genotype_hash;
my %insertion_hash;
my %deletion_hash;   
my $insertion_total;
my $deletion_total;
my $zygotes_total;
my %annotation_hash;
my $annotation_total;
my %exonic_consequence_hash;
my $intronic;
my $exonic;
my $intergenic;

# reading the positional information
my $line_count = 1;
while (<InputPositions>) {
   $_ =~ s/\n|\r//;
   #print "$_\n";
   my @line = split('\t', $_);

   # getting the chromosome and coordinate fields from input file
   # fastahack will need to the chromosome and coordinate to read the information from the reference
   my $chromosome = $line[0];
   my $coordinate = $line[1];

   # get coordinates of first and last character in the context
   my $start_region = $coordinate - 1;
   my $end_region = $coordinate + 1;

   # if the coordinate is the very first letter on the chromosome, then do not read before that position
   # the context becomes 2 letter code, as opposed to a trinucleotide
   if ( $start_region == 0 ) {
      $start_region = 1;
      $end_region = 2;
   }

   #print "$Fastahack -r $chromosome:$start_region..$end_region $Reference\n";
   my $context = `$Fastahack -r $chromosome:$start_region..$end_region $Reference`;
   
   # capitalize context letters
   $context = uc($context);

   #### IF USING CONTROLLED DATA, split germline column into germline allele and mutated_to allele
   # my @germline = split ('/', $line[6]);

   # if germline allele does not equal reference allele, print "start_region germline allele end_region"
   # specifically, replace the middle letter of the context with the germline allele
   #print "$germline[0], $germline[1]\n";
   # if ($germline[0] ne $germline[1]) {
      # print "germline/reference mismatch, line number $line_count\n";
      # if ($coordinate != 1) {
         # substr($context,1,1)= $germline[1];   
      # }
      # else {
        #  substr($context,0,1)= $germline[1];
      # }
   # }

   print OutputTrinucleotideContext "$_\t$context";
   

   ###############################
   # new section: forming the data structure
   ###############################

   # to create N_N contexts for data structure, context_code is defined as the trinucleotide context with a blank middle allele
   my $context_code=$context;
   $context_code =~ s/\n|\r//;
   substr($context_code,1,1) = "_";
   
   # create variables for mutated_from and  mutated_to nucleotides
   my $mutated_from = $line[3];
   my $mutated_to = $line[4];

   # creating genotype variable from column 10 of VCF
   my $genotype = $line[9];

   # incrementing each genotype
   # $genotype_hash{$genotype} = $genotype_hash{$genotype} + 1;

   # splitting heterozygosity by comma, defining heterozygosity total
   my @zygotes = split (',', $mutated_to);
   my $zygotes_length = scalar(@zygotes);

   # identify heterozygosity, choose one at random to use. Count heterozygosity instances
   if ($zygotes_length > 1) {
      my $zygotesRand = $zygotes_length*rand();
      my $zygotesRound = round($zygotesRand) - 1;
      $zygotes_total = $zygotes_total + 1;

      $mutated_to = $zygotes[$zygotesRound];
      # print "@zygotes\t$mutated_to\n";
   }  
   # print "@zygotes\t$mutated_to\n";
   
   # my $round_rand_test = round(rand());
   # print $round_rand_test;

   # define length of insertions and deletions
   my $insertion_length;
   my $deletion_length;
   if ($mutated_from eq "-") {
      $insertion_length = length( $mutated_to );
   }
   else {
      $insertion_length = length( $mutated_to ) - 1;
   }
   if ($mutated_to eq "-") {
      $deletion_length = length( $mutated_from );
   }
   else {
      $deletion_length = length( $mutated_from ) - 1;
   }

   # context_codes are totalled
   $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to} = $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to} + 1; 
   $context_tally_across_mutated_to{$context_code}{$mutated_from} = $context_tally_across_mutated_to{$context_code}{$mutated_from} + 1; 

   # insertion and deletion lengths are totalled
   if ($insertion_length > $deletion_length) {
      $insertion_hash{$insertion_length} = $insertion_hash{$insertion_length} + 1;
   }
   if ($deletion_length > $insertion_length) {
      $deletion_hash{$deletion_length} = $deletion_hash{$deletion_length} + 1;
   }

   # total insertions and deletions
   if ($insertion_length != $deletion_length) {
      if ($insertion_length > $deletion_length) {
         $insertion_total = $insertion_total + 1;
      }
      elsif ($deletion_length > $insertion_length) {
         $deletion_total = $deletion_total + 1;
      } 
   }

   # Find variant annotation and exonic consequence in ANNOVAR outfile
   my $annotation = $line[7];
   if ( $annotation =~ /Func.refGene=(.{1,30});Gene\.refGene/ ) {
      # print "$1\n";
      $annotation_hash{$1}++;
      $annotation_total++;
   }   
   if ( $annotation =~ /ExonicFunc.refGene=(.{1,30});AAChange\.refGene/ ) {
      # print "$1\n";
      $exonic_consequence_hash{$1}++;
   } 
   if ( $annotation =~ /Func.refGene=.{0,15}intronic\;/ ) {
      $intronic++;
   }
   if ( $annotation !~ /Func.refGene=ncRNA_exonic/ ) {
      if ( $annotation =~ /Func.refGene=.{0,15}exonic\;/ ) {
         $exonic++;
      }
   }
   if ( $annotation =~ /Func.refGene=.{0,15}intergenic\;/ ) {
      $intergenic++;
   }
   elsif ( $annotation =~ /Func.refGene=.{0,15}ncRNA_splicing\;/ ) {
      $intergenic++;
   }
   elsif ( $annotation =~ /Func.refGene=.{0,15}upstream\;/ ) {
      $intergenic++;
   }
   elsif ( $annotation =~ /Func.refGene=.{0,15}downstream\;/ ) {
      $intergenic++;
   }
   
   $location{$coordinate}++;
   # Reading input gff file, incrementing gff variant region hash
   # while (<HumanGFF>) {
      # $_ =~ s/\n|\r//;
      # my @line = split('\t', $_);
      # my $region_name = "$line[3]-$line[4]";
      # if ($coordinate >= $line[3] && $coordinate <= $line[4]) {
         # $gff_hash{$region_name}++;
         # $gffMatch++;
         # print "$coordinate $region_name\n";
      # }
   # }
   #print "$region_name, $gff_hash{$region_name}\n";
   

   # to keep track of progress
   # 1000000 for LARGE dbsnp vcfs, 10000 for smaller vcf/tsv tumor mutation files
   unless ($line_count%10000) {
      print "processed $line_count lines\n";
   }
   $line_count++; 
} 
# end working through the input file

# print total number of mutations
my $mutation_total = $line_count;
print "Number of Mutations -- $mutation_total\n";


################### Reading the input gff and creating custom BED file ####################

my $gffBED = "vars.bed";
open(my $bed_handle, '>', $gffBED) || die("Could not open file!");

# Print BED file Header
print $bed_handle "START\tEND\tVariant_Frequency\n";

# Reading input gff file, incrementing gff variant region hash
while (<HumanGFF>) {
   $_ =~ s/\n|\r//;
   my @line = split('\t', $_);
   my $region_name = "$line[3]-$line[4]";
   my $region_length = $line[4] - $line[3];
   my $region_freq = 0;
   foreach my $coordinate (sort(keys %location)) {
      if ($coordinate >= $line[3] && $coordinate <= $line[4]) {
         $gff_hash{$region_name}++;
         $gffMatch++;
         # print "$coordinate $region_name\n";
      }
   }
   if ($gff_hash{$region_name} == 0) {
      print $bed_handle "$line[3]\t$line[4]\t$region_freq\n";
   }
   if ($gff_hash{$region_name} > 0) {
      $region_freq = $gff_hash{$region_name} / $region_length;
      print $bed_handle "$line[3]\t$line[4]\t$region_freq\n";
      print "Region $region_name variant frequency -- $region_freq\n";
      print "Total variants in region $region_name -- $gff_hash{$region_name}\n";
   }
}
   #print "$region_name, $gff_hash{$region_name}\n";

print "GFF Match -- $gffMatch\n";


######################### open files for writing ##########################


# my $genotype_name = "zygosity.prob";
# open(my $genotype_handle, '>', $genotype_name) || die("Could not open file!");
   
my $insertion_file_name = "SSM_insLength.prob";
open(my $insertion_prob_handle, '>', $insertion_file_name) || die("Could not open file!");

my $deletion_file_name = "SSM_delLength.prob";
open(my $deletion_prob_handle, '>', $deletion_file_name) || die("Could not open file!");

my $overall_file_name = "SSM_overall.prob";
open(my $overall_prob_handle, '>', $overall_file_name) || die("Could not open file!");

my $heterozygosity_file_name = "heterozygosity.prob";
open(my $heterozygosity_prob_handle, '>', $heterozygosity_file_name) || die("Could not open file!");

my $annotation_file_name = "annofreq.prob";
open(my $annotation_handle, '>', $annotation_file_name) || die("Could not open file!");

my $exonic_con_file_name = "exonic_consequences.prob";
open(my $exonic_con_handle, '>', $exonic_con_file_name) || die("Could not open file!");

my $intronic_file_name = "intronic_vars.prob";
open(my $intronic_handle, '>', $intronic_file_name) || die("Could not open file!");

my $exonic_file_name = "exonic_vars.prob";
open(my $exonic_handle, '>', $exonic_file_name) || die ("Could not open file!");

my $intergenic_file_name = "intergenic_vars.prob";
open(my $intergenic_handle, '>', $intergenic_file_name) || die ("Could not open file!");


######################### Calculate frequency models ####################### 


# calculate zygosity ratio frequency, print to file
# foreach my $genotype (sort(keys %genotype_hash)) {
   # my $zygosity_frequency;
   # $zygosity_frequency = $genotype_hash{$genotype}/$mutation_total;
   # print $genotype_handle "$genotype\t$zygosity_frequency\n";
   # print "Genotype, $genotype -- $genotype_hash{$genotype}\n";
# }

# print annotation and exonic consequence frequencies
foreach $1 (sort(keys %annotation_hash)) {
   my $annotation_frequency;
   $annotation_frequency = $annotation_hash{$1}/$mutation_total;
   print "$1 -- $annotation_hash{$1}, $annotation_frequency\n";
   print $annotation_handle "$1\t$annotation_frequency\n";
}
foreach $1 (sort(keys %exonic_consequence_hash)) {
   my $exonic_con_freq;
   if ( $1 ne "." ) {
      $exonic_con_freq = $exonic_consequence_hash{$1}/$mutation_total;
      print "Exonic Consequence: $1 -- $exonic_consequence_hash{$1}, $exonic_con_freq\n";
      print $exonic_con_handle "$1\t$exonic_con_freq\n";
   }
}

# Calculating exonic, intronic, and intergenic frequencies, printing to files
my $intronic_freq;
my $exonic_freq;
my $intergenic_freq;
$intronic_freq = $intronic/$mutation_total;
$exonic_freq = $exonic/$mutation_total;
$intergenic_freq = $intergenic/$mutation_total;
print $intronic_handle "$intronic_freq\n";
print $exonic_handle "$exonic_freq\n";
print $intergenic_handle "$intergenic_freq\n";

print "Intronic -- $intronic\nExonic -- $exonic\nIntergenic -- $intergenic\n";
#print "Total Annotations -- $annotation_total\n";

# print overall likelihood file headers
print $overall_prob_handle "mutation_type\tprobability\n";

# print insertions and deletion probabilities out of all mutations
my $insertion_prob_all = $insertion_total / $mutation_total;
my $deletion_prob_all = $deletion_total / $mutation_total;
print $overall_prob_handle "insertion\t$insertion_prob_all\ndeletion\t$deletion_prob_all\n";
# print $overall_prob_handle "Deletion Probability -- $deletion_prob_all\n";

# print InDel totals
print "Insertions $insertion_total\n";
print "Deletions $deletion_total\n";

# print insertion and deletion headers
print $insertion_prob_handle "insertion_length\tprobability\n";
print $deletion_prob_handle "deletion_length\tprobability\n";

# calculate InDel length totals and probability out of total number of insertions/deletions. Print probabilities to file.
foreach my $insertion_length (sort(keys %insertion_hash)) {
   my $insertion_probability;
   $insertion_probability = $insertion_hash{$insertion_length}/$insertion_total;
   print $insertion_prob_handle "$insertion_length\t$insertion_probability\n";
   # print "Insertion, $insertion_length, total , $insertion_hash{$insertion_length}\n";
}
foreach my $deletion_length (sort(keys %deletion_hash)) {
   my $deletion_probability;
   $deletion_probability = $deletion_hash{$deletion_length}/$deletion_total;
   print $deletion_prob_handle "$deletion_length\t$deletion_probability\n";
   # print "Deletion, $deletion_length, total, $deletion_hash{$deletion_length}\n";
}
 
# print heterozygosity frequency to file
my $zygote_frequency = $zygotes_total / $mutation_total;
print $heterozygosity_prob_handle "$zygote_frequency\n";

print "heterozygous alleles -- $zygotes_total\n";


# define nucleotide array
my @nucleotides = ("A", "C", "G", "T");
   
foreach my $nt1 (@nucleotides) {
   foreach my $nt3 (@nucleotides) {

      # define the output file name and open it for writing
      my $trinucleotide_SNP_probability_file_name = "Context".$nt1."-".$nt3.".trinuc";
      open(my $trinuc_prob_handle, '>', $trinucleotide_SNP_probability_file_name) || die("Could not open file!");


      # print trinucleotide contexts and corresponding totals for every mutated_to nucleotide
      my $context_code=$nt1."_".$nt3;

         #foreach my $mutated_from_nucl_key (keys %{ $trinucleotide_context_data{$context_code} }) {
         foreach my $mutated_from (@nucleotides) {
            # define the "mutated_to" keys in trinuc context hash
            # my $mutated_to_nucl_key;

            # the sum is only across mutated_to, and will be redefined for each mutated_from
            my $context_sum_across_mutated_to = 0;
            my $context_sum_across_indel = 0;

            # print "\nRaw counts for mutated_from $mutated_from \n";


            # foreach $mutated_to_nucl_key (keys %{ $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key} }) {
            foreach my $mutated_to (@nucleotides) {
               my $mutated_from_length = length( $mutated_from );
               my $mutated_to_length = length( $mutated_to );
               if ( $mutated_from_length == 1 ) {
                  if ( $mutated_from ne "-" ) {
                     if ( $mutated_to_length == 1 ) {
                        if ( $mutated_to ne "-" ) {
                           # print "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key -- $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key}\n";
                           $context_sum_across_mutated_to = $context_sum_across_mutated_to + $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to};
                        }# end if statement
                        else {
                           $context_sum_across_indel = $context_sum_across_indel + $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to};
                        }# end else statement
                     }# end if statement
                     else {
                        $context_sum_across_indel = $context_sum_across_indel + $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to};
                     }# end else statement
                  }# end if statement
                  else {
                     $context_sum_across_indel = $context_sum_across_indel + $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to};
                  }# end else statement
               }# end if statement
               else {
                  $context_sum_across_indel = $context_sum_across_indel + $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to};
               }# end else statement
               # print "$context_code, $mutated_from, $mutated_to-- $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to}\n";
            }# end of loop over mutated_to

            # print "\nProbabilities for mutated_from $mutated_from:\n";


            foreach my $mutated_to (@nucleotides) {
            #foreach $mutated_to_nucl_key (keys %{ $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key} }) {
               my $mutated_from_length = length( $mutated_from);
               my $mutated_to_length = length( $mutated_to);
               if ( $mutated_from_length == 1 ) {
                  if ( $mutated_from ne "-" ) {
                     if ( $mutated_to_length == 1 ) {
                        if ( $mutated_to ne "-" ) {
                           my $SNP_probability;
                           if ( $context_sum_across_mutated_to == 0 ) {
                              $SNP_probability = 0;
                           }
                           else {
                              $SNP_probability = $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to}/$context_sum_across_mutated_to;
                           }
                           if ( $mutated_to eq "T" ) {
                              print $trinuc_prob_handle "$SNP_probability";
                           }
                           else {
                              # print "$context_code, $mutated_from, $mutated_to, context_sum_across_mutated_to=$context_sum_across_mutated_to -- $SNP_probability\n";
                              print $trinuc_prob_handle "$SNP_probability\t";
                           }
                        }# end of if statement
                        else {
                           my $indel_probability = $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to}/$context_sum_across_indel;
                           # print $indel_prob_handle "$context_code, $mutated_from, $mutated_to, context_sum_across_indel=$context_sum_across_indel -- $indel_probability\n";
                        }# end else statement
                     }# end of if statement
                     else {
                        # my $indel_probability = $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to}/$context_sum_across_indel;
                        # print $indel_prob_handle "$context_code, $mutated_from, $mutated_to, context_sum_across_indel=$context_sum_across_indel -- $indel_probability\n";
                     }# end else statement
                  }# end of if statement
                  else {
                     # my $indel_probability = $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to}/$context_sum_across_indel;
                     # print $indel_prob_handle "$context_code, $mutated_from, $mutated_to, context_sum_across_indel=$context_sum_across_indel -- $indel_probability\n";
                  }# end else statement
               }# end of if statement
               else {
                  my $indel_probability;
                  if ( $context_sum_across_indel = 0 ) {
                     $indel_probability = 0;
                  }
                  else {
                     # $indel_probability = $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to}/$context_sum_across_indel;
                     # print $indel_prob_handle "$context_code, $mutated_from, $mutated_to, context_sum_across_indel=$context_sum_across_indel -- $indel_probability\n";
                  }
               }# end else statement
            }# end of loop over mutated_to
            print $trinuc_prob_handle "\n";

         }# end of loop over mutated_from

     # print "\n\n";
     

  }# end loop over nt3
}# end loop over nt1




