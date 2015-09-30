use strict;
use warnings;
package CracTools::BenchCT::Events::Exon;
# ABSTRACT: A collection of events of type 'exon' on the genome

use parent 'CracTools::BenchCT::Events::Interval';

use CracTools::Interval::Query;

sub printEvent {
  my ($self,$fh,$id) = @_;
  my $exon = $self->getEvent($id);
  print $fh join("\t",$id,split('-',$exon)),"\n";
}

sub addExon {
  my $self = shift;
  return $self->addInterval(@_);
}

sub isTrueExon {
  my $self = shift;
  my ($chr,$start,$end,$strand) = @_; 
  return $self->isTrueInterval($chr,$start,$end,$strand);
}

1;
