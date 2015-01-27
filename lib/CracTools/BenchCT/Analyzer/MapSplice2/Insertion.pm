use strict;
use warnings;
package CracTools::BenchCT::Analyzer::MapSplice2::Insertion;
# ABSTRACT: Analyze deletion bed files (as tophat produce)
#
use parent 'CracTools::BenchCT::Analyzer::BED';

use CracTools::Utils;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'insertion') {
    return 1;
  }
  return 0;
}

sub _init {
  my $self = shift;
  my %args = @_;
  my $insertion_file = $args{file};
  my $insertions_it = CracTools::Utils::getFileIterator(file => $insertion_file,
    parsing_method => \&parseMapSplice2InsertionLine,
    header_regex => '^#'
  );
  while (my $insertions = $insertions_it->()) {
    # We can have multiple insertions per line that start at the same position
    foreach my $ins (@{$insertions->{insertions}}) {
      my $true_insertion = $self->checker->isTrueInsertion($insertions->{chr},$insertions->{start},length($ins->{sequence}));
      if($true_insertion) {
        $self->getStats('insertion')->addTruePositive($true_insertion);
      } else {
        $self->getStats('insertion')->addFalsePositive();
      }
    }
  }
}


=head2 parseMapSplice2ChimeraLine

  chr20~chr17 49411710 59445688 JUNC_1215 354 ++ 255,0,0 2 94,94,201,174, 0,18446744073699517733, 3.886652 6 GTAG 0 4 1.596045 50 50 0 354 0 208 156 52 0 208 146 282 49411616 59445782 59445681,180M12006N66M8046N47M|59445681,180M20118N47M|59445681,180M23477N74M|59445681,180M24387N8M| 49411511,205M| 0 0 0.304892 0.182286 293 188 205 205 doner_exact_matched acceptor_exact_matched CCTGACCCCCGAGCCTGGGGCCGAG AGGGTCACGCTCCTGTCAAAGGTAC 1 from_fusion fusion +,+ BCAS4, BCAS3, 

=cut

sub parseMapSplice2InsertionLine {
  my $line = shift;
  my @fields = split("\t",$line);
  my @insertions = split(",",$fields[29]);
  my @insertions_hash;
  foreach my $ins (@insertions) {
    my($sequence,$cover) = split("-",$ins);
    push(@insertions_hash,{
        cover => $cover,
        sequence => $sequence,
      }
    );
  }
  return {
    chr => $fields[0],
    start => $fields[1],
    insertions => \@insertions_hash,
    
  };
}

1;
