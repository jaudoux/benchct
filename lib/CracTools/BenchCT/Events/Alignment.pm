use strict;
use warnings;
package CracTools::BenchCT::Events::Alignment;
# ABSTRACT: A collection of reads alignments on the genome

use parent 'CracTools::BenchCT::Events';

=head2 new

=cut

sub new {
  my $class = shift;

  # Call parent constructor
  my $self  = $class->SUPER::new(@_);

  # Get args
  my %args = @_;

  # Hash that contains for each read $id its alignement
  #$self->{reads_id} = {};
  
  return $self;
}

sub addAlignment {
  my $self = shift;
  my $bed_line = shift; 

  # We create a string that hold the alignment to be as "light" as possible
  my $alignment_string = join(",",
    map { join(":", 
      $_->{chr}, 
      $_->{strand}, 
      $_->{block_start}, 
      $_->{size}, 
      $_->{ref_start},
    ) } @{$bed_line->{blocks}}
  );

  # We place the alignement string in a new event
  my $id = $self->addEvent($alignment_string);

  # We create a new entry to the conversion hash
  #$self->{reads_id}->{$bed_line->{name}} = $id;
  return $id;
}

sub isGoodAlignment {
  my $self = shift;
  my ($read_id,$ref_name,$pos_start,$ref_start,$strand) = @_;

  # Retrieve the alignment string that correspond to this read
  #my $event_id = $self->{reads_id}->{$read_name};
  

  # Retrieve the string that contains alignment informations
  my $alignment_string = $self->getEvent($read_id);
   
  # If no event id exists for this read name, we return false
  return 0 if !defined $alignment_string;
  
  foreach my $block (split ",", $alignment_string) {
    my ($b_chr,$b_strand,$b_start,$b_size,$b_ref_start) = split ":", $block;

    next if ($b_start + $b_size) < $pos_start;
    
    # Check if the read has been mapped to the right reference
    next if $b_chr ne $ref_name;
     
    # Check if protocol is stranded and the read is mapped to the wrong strand
    next if $b_strand != $strand;

    # Difference between the first pos of the block to the
    # pos where the alignement is started
    my $delta = $pos_start - $b_start;

    # Check if the position if right within a THRESHOLD_MAPPING window
    if(abs($b_ref_start + $delta - $ref_start) <= $self->threshold) {
      return $read_id + 1;
    }
  }

  # If we have found no valid block, then we return false
  return 0;
}

1;
