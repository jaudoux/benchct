use strict;
use warnings;
package CracTools::BenchCT::Stats;
# ABSTRACT: Store statistics about FP/TP/FN on independant elements

=head2 new

=cut

use CracTools::BitVector;
use CracTools::Utils;

sub new {
  my $class = shift;
  my %args = @_;

  my $self = bless {
    nb_elements => $args{nb_elements}, # This is equal to : number of true positives + number of false negatives
    bitvector => CracTools::BitVector->new($args{nb_elements}),
    true_positives => 0,
    false_positives => 0,
    false_positives_fh => defined $args{false_positives_file}? CracTools::Utils::getWritingFileHandle($args{false_positives_file}) : undef,
  },$class;

  return $self;
}

sub addTruePositive {
  my $self = shift;
  my $id = shift;
  # If we have already seen this event we do not count it
  if(defined $id) {
    if($self->{bitvector}->get($id-1) == 0) {
      $self->{bitvector}->set($id-1);
      $self->{true_positives}++;
    }
  } else {
    $self->{true_positives}++;
  }
}

sub addFalsePositive {
  my $self = shift;
  my $false_positive = shift;
  my $fh = $self->getFalsePositivesFileHandle;
  if(defined $false_positive && defined $fh) {
    chomp $false_positive;
    print $fh $false_positive,"\n";
  }
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

sub getFalsePositivesFileHandle {
  my $self = shift;
  return $self->{false_positives_fh};
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
