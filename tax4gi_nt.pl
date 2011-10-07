#!/usr/bin/perl
# coordinates_from_p4e.pl     Emanuel Heitlinger  

use warnings;      # avoid d'oh! bugs
use strict;        # avoid d'oh! bugs
$|=1;
use Data::Dumper;  # easy printing, sometimes only used during coding


use Bio::Perl;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;
use Bio::LITE::Taxonomy::NCBI::Gi2taxid; # qw/new_dict/;

my $taxDB = Bio::LITE::Taxonomy::NCBI-> new (
                                              db=>"NCBI",
                                              names=>
                                              "/home/ele/db/blastdb/taxonomy/names.dmp",
                                              nodes=>
                                              "/home/ele/db/blastdb/taxonomy/nodes.dmp",
                                              dict=>"/home/ele/db/blastdb/taxonomy/gi_taxid_nucl.bin"
                                             );

while (<>) {
  chomp($_);
  my $taxid = $taxDB->get_taxid($_);

  my $family = $taxDB->get_term_at_level($taxid,"family");
  my $phylum = $taxDB->get_term_at_level($taxid,"phylum");
  my $kingdom = $taxDB->get_term_at_level($taxid,"kingdom");

  print "$_,$taxid,$family,$phylum,$kingdom\n";
}


__END__

=head1 NAME

tax4gi.pl

=head1 SYNOPSIS

tax4gi.pl gi-file

=head1 DESCRIPTION

print the taxon number for a gi-number

=head1 ARGUMENTS


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
