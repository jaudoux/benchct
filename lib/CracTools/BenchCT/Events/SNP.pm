use strict;
use warnings;
package CracTools::BenchCT::Events::SNP;
# ABSTRACT: A collection of events of type 'SNP' (snp,indel) on the genome

use parent 'CracTools::BenchCT::Events';

use CracTools::GenomeMask;

=head2 new

=cut

sub new {
  my $class = shift;

  # Call parent constructor
  my $self  = $class->SUPER::new(@_);

  # Get args
  my %args = @_;

  #my $crac_index_conf = $args{crac_index_conf};

  #print STDERR "[Events::SNP] Creating GenomeMask over the genome\n" if $self->verbose;

  #$self->{genome_mask} = CracTools::GenomeMask->new(crac_index_conf => $crac_index_conf,
  #  verbose => $self->{verbose},
  #);

  $self->{interval_query} = CracTools::Interval::Query->new();
  
  return $self;
}

#sub genomeMask {
#  my $self = shift;
#  return $self->{genome_mask};
#}

sub intervalQuery {
  my $self = shift;
  return $self->{interval_query};
}

=head2 addMutation

Add a new event to the collection

=cut

sub addMutation {
  my $self  = shift;
  my ($chr,$pos,$new_nuc) = @_;
  # First we check if this event is already registered
  my $id = $self->isTrueMutation($chr,$pos,$new_nuc,0);
  if($id) {
    return $id - 1;
  } else {
    $id = $self->addEvent($new_nuc);
    $self->intervalQuery->addInterval($chr,
      $pos,
      $pos,
      undef, # A mutation does not have a strand....
      $id,
    );
    return $id;
  }
}

sub isTrueMutation {
  my $self = shift;
  my ($chr,$pos,$nuc) = @_;
  return $self->_foundSNP($chr,$pos,$nuc,$self->threshold);
}

sub _foundSNP {
  my $self = shift;
  my ($chr,$pos,$nuc,$threshold) = @_;
  my @snps = @{$self->intervalQuery->fetchByRegion($chr,$pos - $threshold,$pos + $threshold)};
  # We return true if we have found a matching snp that have the same length
  foreach my $snp_id (@snps) {
    my $snp_nuc = $self->getEvent($snp_id);
    if($snp_nuc eq $nuc) {
      return $snp_id + 1;
    }
  }
  return 0;
}

1;
