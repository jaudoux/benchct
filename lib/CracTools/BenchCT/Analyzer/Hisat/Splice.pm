use strict;
use warnings;
package CracTools::BenchCT::Analyzer::Hisat::Splice;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;
#use Data::Dumper;

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
    parsing_method => \&parseHisatJunctionLine,
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
      $self->getStats('splice')->addTruePositive($true_splice);
    } else {
      $self->getStats('splice')->addFalsePositive();
    }
  }
}

sub parseHisatJunctionLine {
  my $line = shift;
  my($chr,$start,$end,$strand) = split("\t",$line);
  return {
    chr           => $chr,
    start         => $start-1, # Because STAR junctions file is 1-based
    end           => $end-1, # Because STAR junctions file is 1-based
    strand        => CracTools::Utils::convertStrand($strand),
  };
}

1;
