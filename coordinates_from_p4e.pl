#!/usr/bin/perl
# coordinates_from_p4e.pl     Emanuel Heitlinger  

use warnings;      # avoid d'oh! bugs
use strict;        # avoid d'oh! bugs
$|=1;
use Data::Dumper;  # easy printing, sometimes only used during coding


## use Bio::Tools::pSW;
use Bio::Perl;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast; 
use Getopt::Long;  # support options+arguments
use Pod::Usage;    # avoid redundant &Usage()

my ($opt_help, $opt_man);
my ($p4e_fasta, $raw_fasta, $visual);

GetOptions(
	   'p4e_fasta=s'      => \$p4e_fasta,
           'raw_fasta=s' => \$raw_fasta,
           'visual=s'  => \$visual,
	   'help!'     => \$opt_help,                                # for all the standard 
	   'man!'      => \$opt_man,                                 # command-line  args
	  ) or pod2usage(-verbose => 1) && exit;                    # or exit printing usage
pod2usage(-verbose => 1) && exit if defined $opt_help;        # obligatory arguments wrongly
pod2usage(-verbose => 2) && exit if defined $opt_man;

# check existence of obligatory commmand-line arguments (optional)
pod2usage(-verbose => 1) && exit unless $p4e_fasta and $raw_fasta;

## set the blast for later use filters have to be off!!!
my @params = (program  => 'blastn',
              F => 'F');
my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);


my %trans = %{fastafile2hash($p4e_fasta)};
my %raw = %{fastafile2hash($raw_fasta)};


open VIS, ">$visual" or die $! if $visual;

for my $contig ( keys %raw) { 
  
  ## only if sequence exists in translated to allow giving raw
  ## sequence which contains sequences not used in p4e
  if (exists ($trans{$contig})) {
    
    my $raw_seq =  $raw{$contig}{seq};
    my $query = Bio::Seq->new(-id  => "translation",
                              -seq =>  $trans{$contig}{seq});
    my $subject = Bio::Seq->new(-id  => "raw",
                                -seq =>   $raw_seq); 
    # execute the blast command with this line
    my $blast_report = $factory->bl2seq ($query,  $subject);

    # just one result in a bl2seq report first hsp for this
    my $hsp_obj = $blast_report->next_result->next_hit->next_hsp; 
  
    my $translated_seq = $trans{$contig}{seq};

    my $blast_query = $hsp_obj->query_string;
    my $hom_string = $hsp_obj->homology_string;
    my $blast_hit= $hsp_obj->hit_string;

    my $hit_start = $hsp_obj->start('hit');
    my $query_start = $hsp_obj->start('query');
    my $hit_end = $hsp_obj->end('hit');
    my $query_end = $hsp_obj->end('query');
    
    my $indicator = ">";

    if ($hsp_obj->strand('hit') == -1){
      $translated_seq = revcom($translated_seq)->seq;
      $blast_query = revcom($blast_query)->seq;
      $hom_string = reverse($hom_string);
      $blast_hit = revcom($blast_hit)->seq;
      $indicator = "<";
    }

    my $replacement_string = $blast_query;
    $replacement_string =~ s/-//g;

    ## X'es introduced by prot4est at the beginning or the end of a
    ## orf indicate: XX = one base is explicit in the raw sequence; X
    ## = two bases are explicit in the raw sequence! ie. the X is a
    ## placeholder.  -> substract the codon containing the placeholder
    ## from the orf
    
    my $start_minus = 0;
    if ($translated_seq =~ /(^\w?X+)/){ 
       $start_minus = length($1);
       $start_minus =~ tr/12/21/;
    }
  
    my $stop_minus  = 0;
    if ($translated_seq =~ /(X+\w?$)/){ 
      $stop_minus = length($1);
      $stop_minus =~ tr/12/21/;
    } 

    ## Get the real number of gaps introduced relative to the raw
    ## sequence
    my $blast_query_copy = $blast_query;
    my $raw_correct = 0;
    $raw_correct = $blast_query_copy =~ s/r|y|s|w|k|m|b|d|h|v|n|x//ig;
    my $raw_real =0;
    my $raw_seq_copy = $raw_seq;
    $raw_real = $raw_seq_copy =~ s/r|y|s|w|k|m|b|d|h|v|n|x//ig;
    $raw_correct = $raw_correct - $raw_real;
    
    ## the start of the orf
    my $start_relative = $hit_start + $start_minus - 1;
    ## the complicated end:

    my $stop_relative = $start_relative + length($blast_query) - $stop_minus - $raw_correct ;

    my $before = "";
    $before = substr($raw_seq, 0, $start_relative) if ($start_relative>0);


    $replacement_string = substr($replacement_string, 
                                    $start_minus, length($replacement_string));
    $replacement_string = substr($replacement_string, 
                                    0, length($replacement_string)-$stop_minus);
    
    my $after = "";
    $after = substr( $raw_seq, 
                     $stop_relative,
                     length($raw_seq) - $stop_relative 
                   ) 
      if $stop_relative  < length($raw_seq); # and $hsp_obj->strand('hit') == 1;


    ## In longest-orf predictions prot4 est often leaves out one or
    ## two bases at the end. Thes have to be resubstitutet from the
    ## raw sequence
    my $p4estupidness_patch = "";
    
    my $p4ewrongness=length($replacement_string)%3;
    $p4ewrongness =~ tr/12/21/;
    
    if ($p4ewrongness > 0 and $hsp_obj->strand('hit') == 1) {

      if (length($after) >= $p4ewrongness){ #p4e is wrong
        $p4estupidness_patch = substr($after, 0, $p4ewrongness);
        $after = substr($after, $p4ewrongness, length($after)-$p4ewrongness) 
      }

      ## but sometimes it did just complete the codon, when 3rd
      ## positon was completely degraded -> add an n
      elsif (length($after) < $p4ewrongness){
        $p4estupidness_patch = $after."n";
        $after = "";
      }
    }

    if ($p4ewrongness > 0 and $hsp_obj->strand('hit') == -1) { #p4e is wrong
      if (length($before) >= $p4ewrongness){
        $p4estupidness_patch = substr($before, length($before)-$p4ewrongness, $p4ewrongness);
        $before = substr($before, 0, length($before)-$p4ewrongness)
    }

      ## 3rd positon was completely degraded -> add an n
      elsif (length($before) < $p4ewrongness){
        $replacement_string = "n".$replacement_string;
        $before = "";
      }
    }
    

    my $imputed_seq = uc($before).
      lc($replacement_string).
        lc($p4estupidness_patch).
          uc($after);
    
    my $orf_line="ORF     " .
      " " x $start_relative. 
        "|". "$indicator" x ($stop_relative-$start_relative). "|". "\n"; 
    
    print VIS "--$contig", "-" x length($raw_seq), "\n",
      "p4eNuc   ", " " x $start_relative, $translated_seq, "\n",
        "blst  qe ",  " " x $start_relative, $blast_query,    "\n",
          "         ",  " " x $start_relative, $hom_string, "\n",
            "blst hit ",  " " x $start_relative, $blast_hit,      "\n\n",
              "raw all  ", $raw_seq, "\n",
                "imputed  ", $imputed_seq, "\n",
                  "$orf_line\n",
                    "perc ide " , $hsp_obj-> percent_identity, "\n",
                      "strnd h  ", $hsp_obj->strand('hit'), "\n",
                        "hit start:", $hit_start, "\n", 
                          "start-minus: $start_minus\n",
                            "stop-minus: $stop_minus\n",
                              "start-relative: $start_relative\n",
                                "stop-relative: $stop_relative\n",
                                  "3 mod rep:", length($replacement_string)%3, "\n",
                                    "3 mod corrected orf:", 
                                      length($replacement_string.$p4estupidness_patch)%3, "\n",
                                        "after:", $after, "\n",
                                          "before:", $before, "\n",
                                            "p4ewrongness:", $p4ewrongness, "\n", 
                                              "\n\n" if $visual;
    
    print ">$contig $trans{$contig}{desc}\n",
      $imputed_seq, "\n";
  }

  else {  
    print ">", $contig, " no-prediction\n",
      $raw{$contig}{seq}, "\n";
  }
  
}

sub fastafile2hash{
  my $fastafile = shift @_;
  my %sequences;
  open FA, "<$fastafile" or die $!;
  while (<FA>){
    next unless /^>(\S+)(.*)/;
    chomp($sequences{$1}{seq} = <FA>);
    $sequences{$1}{desc} = $2;
  }
  return \%sequences;
}


sub right{
  my $skalar=shift;
  my $zeichenanzahl=shift;
  return substr($skalar,length($skalar)-$zeichenanzahl,$zeichenanzahl);
}

__END__

=head1 NAME

coordinates_from_p4e.pl

=head1 SYNOPSIS

coordinates_from_p4e.pl -p prot4estNuc.fsa -r raw_fasta.fsa > imputed_sequences.fsa

=head1 DESCRIPTION

Print a fasta file with p4e predicted ORFs in lower-case.

=head1 ARGUMENTS


 --p4e_fasta       prot4est's Nuc.fsa file (obligatory)
 --raw_fasta       raw fasta file you used to predict your sequences from can also be a bigger file containing only a subset of sequences for which proteins were predicted (obligatory)
 --visual          output file to print visual representation of what this  script is doing (optional)
 --help      print Options and Arguments 
 --man       print complete man page

=head1 AUTHOR

Emanuel Heitlinger, emanuelheitlinger@gmail.com

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Emanuel Heitlinger

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
