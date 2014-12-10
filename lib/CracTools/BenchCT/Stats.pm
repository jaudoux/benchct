use strict;
use warnings;
package CracTools::BenchCT::Stats;
# ABSTRACT: Store statistics about FP/TP/FN on independant elements

=head2 new

=cut

sub new {
  my $class = shift;
  my %args = @_;

  my $self = bless {
    nb_elements => $args{nb_elements}, # This is equal to : number of true positives + number of false negatives
    true_positives => 0,
    false_positives => 0,
  }, $class;

  return $self;
}

sub addTruePositive {
  my $self = shift;
  $self->{true_positives}++;
}

sub addFalsePositive {
  my $self = shift;
  $self->{false_positives}++;
}

sub nbElements {
  my $self = shift;
  return $self->{nb_elements};
}

sub nbTruePositives {
  my $self = shift;
  return $self->{true_positives};
}

sub nbFalsePositives {
  my $self = shift;
  return $self->{false_positives};
}

sub nbFalseNegatives {
  my $self = shift;
  return $self->nbElements - $self->nbTruePositives;
}

sub getSensitivity {
  my $self = shift;
  return 0 if($self->nbTruePositives == 0 && $self->nbFalseNegatives == 0);
  return $self->nbTruePositives / ($self->nbTruePositives + $self->nbFalseNegatives);
}

sub getAccuracy {
  my $self = shift;
  return 0 if($self->nbTruePositives == 0 && $self->nbFalsePositives == 0);
  return $self->nbTruePositives / ($self->nbTruePositives + $self->nbFalsePositives);
}

sub getGain {
  my $self = shift;
  return 0 if($self->nbTruePositives == 0 && $self->nbFalseNegatives == 0);
  return ($self->nbTruePositives - $self->nbFalsePositives) / ($self->nbTruePositives + $self->nbFalseNegatives);
}

sub print {
  my $self = shift;
  my $fh = shift;
  print $fh "NB_ELEMENTS: ".$self->nbElements."\n";
  print $fh "NB_TRUE_POS: ".$self->nbTruePositives."\n";
  print $fh "NB_FALSE_POS: ".$self->nbFalsePositives."\n";
  print $fh "NB_FALSE_NEG: ".$self->nbFalseNegatives."\n";

  print $fh "Sensitivity: ".$self->getSensitivity."\n";
  print $fh "Accuracy   : ".$self->getAccuracy."\n";
}

1;
