use strict;
use warnings;
package CracTools::BenchCT::Stats::Events;
# ABSTRACT: Store statistics about FP/TP/FN on events and keep tracks of added elements
#

use parent 'CracTools::BenchCT::Stats';

use CracTools::GenomeMask;

=head2 new

=cut

sub new {
  my $class = shift;

  # Call parent constructor
  my $self  = $class->SUPER::new(@_);

  $self->{genome_mask} = CracTools::GenomeMask->new();

  return $self;
}

sub addTruePositive {
  my $self = shift;
}

1;
