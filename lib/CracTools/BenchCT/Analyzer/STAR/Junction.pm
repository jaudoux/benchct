use strict;
use warnings;
package CracTools::BenchCT::Analyzer::STAR::Junction;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'splice') {
    return 1;
  }
  return 0;
}

sub _init {
  my $self = shift;
  my %args = @_;
  my $junction_file = $args{file};
  my $splice_it = CracTools::Utils::getFileIterator(file => $junction_file,
    parsing_method => \&parseSTARJunctionLine,
    header_regex => '^#',
  );
  while (my $splice = $splice_it->()) {
    my $true_splice = $self->checker->isTrueSplice(
      $splice->{chr},
      $splice->{start},
      $splice->{end} - $splice->{start},
      $splice->{strand},
    );
    
    if($true_splice) {
      $self->getStats('splice')->addTruePositive(id => $true_splice);
    } else {
      $self->getStats('splice')->addFalsePositive();
    }
  }
}

=head2 parseSTARJunctionLine

=cut

sub parseSTARJunctionLine {
  my $line = shift;
  my($chr,$start,$end,$strand,$intron_motif,$annotated,$uniq_cover,$multi_cover,$max_overhang) = split("\t",$line);
  return {
    chr           => $chr,
    start         => $start-1, # Because STAR junctions file is 1-based
    end           => $end-1, # Because STAR junctions file is 1-based
    strand        => $strand == 1? 1 : -1,
    intron_motif  => $intron_motif,
    annotated     => $annotated,
    uniq_cover    => $uniq_cover,
    multi_cover   => $multi_cover,
    max_overhang  => $max_overhang,
  };
}

1;
