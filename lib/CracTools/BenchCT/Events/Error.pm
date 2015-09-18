use strict;
use warnings;
package CracTools::BenchCT::Events::Error;
# ABSTRACT: A collection of read error positions

use parent 'CracTools::BenchCT::Events';

use Carp;
use CracTools::BitVector;

=head2 new

=cut

sub new {
  my $class = shift;

  # Call parent constructor
  my $self  = $class->SUPER::new(@_);

  # Get args
  my %args = @_;

  my $nb_reads      = $args{nb_reads};
  my $max_length    = $args{max_length};
  my $sampling_rate = $args{sampling_rate};

  croak "Missing 'nb_reads' argument"   unless defined $nb_reads; 
  croak "Missing 'max_length' argument" unless defined $max_length; 

  $self->{nb_reads}   = $nb_reads;
  $self->{max_length} = $max_length;

  # Hash that contains for each read $id its alignement
  $self->{errors_bv}  = CracTools::BitVector->new($nb_reads*$max_length);
   
  # Since rank operation has no strong implementation, we have to cheat
  $self->{sampling_rate}        = defined $sampling_rate? $sampling_rate : $max_length*10;
  $self->{sampled_rank}         = [];
  $self->{sampled_rank_status}  = 0;

  for(my $i = 0; $i < int($self->{errors_bv}->length/$self->{sampling_rate}); $i++) {
    $self->{sampled_rank}->[$i] = 0;
  }
  
  return $self;
}

sub addError {
  my $self = shift;
  my ($read_id,$error_pos) = @_;
  
  my $bv_pos = $self->_getErrorBVPos($read_id,$error_pos);

  if(!defined $bv_pos) {
    croak "Error position ($error_pos) for read $read_id is greater that the maximum read length (".$self->{max_length}.")";
  }

  # Check if the error is not already set
  if(!$self->{errors_bv}->get($bv_pos)) {

    # We update the sampled rank array
    $self->{sampled_rank}->[int($bv_pos/$self->{sampling_rate})]++;
    # TODO UPDATE SAMPLED RANK
    #$self->{sampled_rank_status} = 0;

    # Increment the number of evevent
    my $id = $self->addEvent();
    #my $id = $self->addEvent($bv_pos);

    $self->{errors_bv}->set($bv_pos);
    return $id;
  }
}

sub isTrueError {
  my $self = shift;
  my ($read_id,$error_pos) = @_;

  my $bv_pos = $self->_getErrorBVPos($read_id,$error_pos);

  # If the given error position is greater that the maximum read_lenght,
  # we print an error and return false
  if(!defined $bv_pos) {
    carp "Error position ($error_pos) for read $read_id is greater that the maximum read length (".$self->{max_length}.")";
    return 0;
  }

  # If there is an error at this position, we return its ID plus 1
  if($self->{errors_bv}->get($bv_pos)) {
    return $self->_getErrorId($bv_pos) + 1;
  }

  # If we have found no valid block, then we return false
  return 0;
}

# Given a read id (starting at 0) and an error position, return its position
# in the bitvector
sub _getErrorBVPos {
  my $self = shift;
  my ($read_id,$error_pos) = @_;

  #if($read_id >= $self->nbEvents()) {
  #  croak "Read ID is not valid";
  #}

  if($error_pos >= $self->{max_length}) {
    return undef;
  } else {
    return $read_id*$self->{max_length} + $error_pos;
  }
}

# Given and error position in the bitvector, return a unique id
# from 0 to nb_errors - 1
sub _getErrorId {
  my $self      = shift;
  my $bv_pos    = shift;
  my $error_id  = 0;

  if($self->{sampled_rank_status} == 0) {
    for(my $i = 1; $i < @{$self->{sampled_rank}}; $i++) {
      $self->{sampled_rank}->[$i] += $self->{sampled_rank}->[$i-1];
    }
    $self->{sampled_rank_status} = 1;
  }

  # Look over the sampled_rank array and sums up
  my $sampling_group = int($bv_pos/$self->{sampling_rate});
  if($sampling_group > 0) {
    $error_id += $self->{sampled_rank}->[$sampling_group-1];
  }
  
  # Now we campute the last positions
  my $sampling_bound = $sampling_group * $self->{sampling_rate} - 1;
  while($bv_pos > $sampling_bound) {
    $error_id++;
    $bv_pos = $self->{errors_bv}->prev($bv_pos-1);
  }

  return $error_id -1;
}

1;
