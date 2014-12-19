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
  my $info_line = shift;
  my ($old_nuc,$new_nuc) = $info_line->{mutation} =~ /(\S)\s->\s(\S)/;
  my $id = $self->addEvent($new_nuc);
  $self->intervalQuery->addInterval($info_line->{chr},
    $info_line->{old_pos},
    $info_line->{old_pos},
    undef, # A mutation does not have a strand....
    $id,
  );
}

sub isTrueMutation {
  my $self = shift;
  my ($chr,$pos,$nuc) = @_;
  my @snps = @{$self->intervalQuery->fetchByRegion($chr,$pos - $self->threshold,$pos + $self->threshold)};
  # We return true if we have found a matching snp that have the same length
  foreach my $snp (@snps) {
    my $snp_nuc = $self->getEvent($snp);
    if($snp_nuc eq $nuc) {
      return $snp + 1;
    }
  }
  return 0;
}

1;
