use strict;
use warnings;
package CracTools::BenchCT::Analyzer::CRAC::Chimera;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
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
  my $min_score = $args{options}->{min_score};
  my $chim_it = CracTools::Utils::getFileIterator(file => $chimera_file,
    parsing_method => \&parseCRACChimeraLine,
    header_regex => '^#'
  );
  while (my $chimera = $chim_it->()) {

    if(!defined $min_score || !defined $chimera->{score} || $chimera->{score} >= $min_score) {
      my $true_chimera = $self->checker->isTrueChimera($chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2});
      
      if($true_chimera) {
        $self->getStats('chimera')->addTruePositive(id => $true_chimera);
      } else {
        #print STDERR Dumper($chimera);
        $self->getStats('chimera')->addFalsePositive(out_string => join("\t",($chimera->{chr1},$chimera->{pos1},$chimera->{strand1},$chimera->{chr2},$chimera->{pos2},$chimera->{strand2},$chimera->{score}),join(",",@{$chimera->{ids}})));
      }
    }
  }
}

=head2 parseCRACChimeraLine

=cut

sub parseCRACChimeraLine {
  my $line = shift;
  my($chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$score,$ids,$nb) = split("\t",$line);
  my @ids = split(",",$ids);
  return {
    chr1 => $chr1,
    pos1 => $pos1,
    strand1 => CracTools::Utils::convertStrand($strand1),
    chr2 => $chr2,
    pos2 => $pos2,
    strand2 => CracTools::Utils::convertStrand($strand2),
    score => $score eq 'N/A'? undef : $score,
    ids => \@ids,
    nb => $nb,
  };
}

1;
