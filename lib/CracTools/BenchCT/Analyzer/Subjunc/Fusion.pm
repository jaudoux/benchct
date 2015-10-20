use strict;
use warnings;
package CracTools::BenchCT::Analyzer::Subjunc::Fusion;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
 
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;
#use Data::Dumper;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'chimera') {
    return 1;
  }
  return 0;
}

sub _init {
  my $self = shift;
  my %args = @_;
  my $chimera_file = $args{file};
  my $chim_it = CracTools::Utils::getFileIterator(file => $chimera_file,
    parsing_method => \&parseSubjuncFusionLine,
    header_regex => '^#'
  );
  while (my $chimera = $chim_it->()) {
    my $true_chimera = $self->checker->isTrueChimera($chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2});
    
    if($true_chimera) {
      $self->getStats('chimera')->addTruePositive(id => $true_chimera);
    } else {
      #print STDERR Dumper($chimera);
      $self->getStats('chimera')->addFalsePositive(id => join(";",$chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2}));
    }
  }
}

=head2 parseSubjuncFusionLine

=cut

sub parseSubjuncFusionLine {
  my $line = shift;
  my($chr1,$pos1,$chr2,$pos2,$same_strand,$support) = split("\t",$line);
  my $strand1 = 1;
  my $strand2 = $same_strand eq 'Yes'? 1 : -1;
  return {
    chr1 => $chr1,
    pos1 => $pos1,
    strand1 => $strand1,
    chr2 => $chr2,
    pos2 => $pos2,
    strand2 => $strand2,
    support => $support,
  };
}

1;
