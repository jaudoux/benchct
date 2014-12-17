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

  my $crac_index_conf = $args{crac_index_conf};

  print STDERR "[Events::SNP] Creating GenomeMask over the genome\n" if $self->verbose;

  $self->{genome_mask} = CracTools::GenomeMask->new(crac_index_conf => $crac_index_conf,
    verbose => $self->{verbose},
  );
  
  return $self;
}

sub genomeMask {
  my $self = shift;
  return $self->{genome_mask};
}

=head2 addMutation

Add a new event to the collection

=cut

sub addMutation {
  my $self  = shift;
  my $info_line = shift;
  $self->addEvent();
  $self->genomeMask->setPos($info_line->{chr},
    $info_line->{old_pos}, # TODO +1 ?
  );
}

sub isTrueMutation {
  my $self = shift;
  my ($chr,$pos) = @_;
  # Set search bounds with the threshold
  my $start = $pos - $self->threshold;
  my $end = $pos + $self->threshold;
  # Adjust bounds if we have gone too far
  $start = 0 if $start < 0;
  $end = $self->genomeMask->getChrLength($chr) - 1 if $end > $self->genomeMask->getChrLength($chr) - 1;
  return $self->genomeMask->getNbBitsSetInRegion($chr,$pos - $self->threshold,$pos + $self->threshold);
}

1;
