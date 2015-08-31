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
    print_header => $args{print_header},
    print_element => $args{print_element},
    true_positives_nb => 0,
    false_positives_nb => 0,
    false_positives_hash => {}, # Hash to store false positives element
    false_positives_fh => _getFilehandleIfDef($args{false_positives_file}),
    true_positives_fh => _getFilehandleIfDef($args{true_positives_file}),
    false_negatives_fh => _getFilehandleIfDef($args{false_negatives_file}),

  },$class;

  $self->_init();

  return $self;
}

=head2 addTruePositive

  [id] String (optional)
       This is the id of the true positive in order to count it only once.
       This is the value that has been returned by the CracTools::BenchCT::Checker
       object.
  [out_string] String (optional)
       This string will be printed on the output file that contains false positive
       information

=cut

sub addTruePositive {
  my $self = shift;
  my %args = @_;
  my $id = $args{id};
  my $out_string = $args{out_string};
  # If we have already seen this event we do not count it
  if(defined $id) {
    if($self->bitvector->get($id-1) == 0) {
      $self->bitvector->set($id-1);
      $self->{true_positives_nb}++;
      # If we have an output stream and an output string we print it
      _printOutputString($self->getTruePositivesFileHandle,$out_string);
    }
  } else {
    $self->{true_positives_nb}++;
    # If we have an output stream and an output string we print it
    _printOutputString($self->getTruePositivesFileHandle,$out_string);
  }
}

=head2 addFalsePositive 

  [id] String (optional)
       This is the id of the false positive in order to count it only once.
       This is the value that has been returned by the CracTools::BenchCT::Checker
       object. If there is no redundant FP in the file processed by the analyzer
       do not indicate an id, otherwise it will take memory for nothing.
  [out_string] String (optional)
       This string will be printed on the output file that contains false positive
       information

=cut

sub addFalsePositive {
  my $self = shift;
  my %args = @_;
  my $id = $args{id};
  my $out_string = $args{out_string};
  my $fh = $self->getFalsePositivesFileHandle;
  if(defined $id) {
    # If we have already seen that FP, we do not count it anymore
    if(!defined $self->getFalsePositivesHash->{$id}) {
      $self->getFalsePositivesHash->{$id} = 1;
      $self->{false_positives_nb}++;
      # If we have an output stream and an output string we print it
      _printOutputString($self->getFalsePositivesFileHandle,$out_string);
    } else {
      $self->getFalsePositivesHash->{$id}++;
    }
  } else {
    $self->{false_positives_nb}++;
    # If we have an output stream and an output string we print it
    _printOutputString($self->getFalsePositivesFileHandle,$out_string);
  }
}

sub bitvector {
  my $self = shift;
  return $self->{bitvector};
}

sub nbElements {
  my $self = shift;
  return $self->{nb_elements};
}

sub nbTruePositives {
  my $self = shift;
  return $self->{true_positives_nb};
}

sub nbFalsePositives {
  my $self = shift;
  return $self->{false_positives_nb};
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

sub getTruePositivesFileHandle {
  my $self = shift;
  return $self->{true_positives_fh};
}

sub getFalseNegativesFileHandle {
  my $self = shift;
  return $self->{false_negatives_fh};
}

sub getFalsePositivesHash {
  my $self = shift;
  return $self->{false_positives_hash};
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

sub closeOutputs {
  my $self = shift;

  # First we print FN before closing this stream
  if(defined $self->getFalseNegativesFileHandle) {
    my $fh = $self->getFalseNegativesFileHandle;
    for(my $i = 0; $i < $self->bitvector->length; $i++) {
      if($self->bitvector->get($i) == 0) {
        if(defined $self->{print_element}) {
          $self->{print_element}->($fh,$i);
        } else {
          print $fh $i,"\n";
        }
      }
    }
  }
  
  close $self->getFalsePositivesFileHandle() if defined $self->getFalsePositivesFileHandle();
  delete $self->{false_positives_fh};
  close $self->getFalseNegativesFileHandle() if defined $self->getFalseNegativesFileHandle();
  delete $self->{false_negatives_fh};
  close $self->getTruePositivesFileHandle() if defined $self->getTruePositivesFileHandle();
  delete $self->{true_positives_fh};

  delete $self->{print_header};
  delete $self->{print_element};
}

=head1 PRIVATE METHODS

=cut

sub _init {
  my $self = shift;
  # print header ins output files, if they exists
  $self->{print_header}->($self->getFalseNegativesFileHandle) if defined $self->getFalseNegativesFileHandle && defined $self->{print_header};
  $self->{print_header}->($self->getTruePositivesFileHandle) if defined $self->getTruePositivesFileHandle && defined $self->{print_header};
}

sub _getFilehandleIfDef {
  my $file = shift;
  return defined $file? CracTools::Utils::getWritingFileHandle($file) : undef,
}

sub _printOutputString {
  my $fh = shift;
  my $out_string = shift;
  if (defined $fh && defined $out_string) {
    chomp $out_string;
    print $fh $out_string,"\n";
  }
}

1;
