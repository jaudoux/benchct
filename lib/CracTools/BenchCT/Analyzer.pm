use strict;
use warnings;
package CracTools::BenchCT::Analyzer;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use CracTools::BenchCT::Stats;

=head2 new

=cut

sub new {
  my $class = shift;
  my %args = @_;

  my $self = bless {
    checker => $args{checker},
  }, $class;


  # If we are looking at reads alignement
  if($args{check_reads_mapping} && $self->canCheckMapping) {
    $self->{mapping_stats} = CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbReads());
  }

  # If we are looking at reads alignement
  if($args{check_errors} && $self->canCheckErrors) {
    $self->{errors_stats} = CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbErrors());
  }

  # If we are looking at reads alignement
  if($args{check_snps} && $self->canCheckSnps) {
    $self->{snps_stats} = CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbEvents('snp'));
  }

  # If we are looking at reads alignement
  if($args{check_splices} && $self->canCheckSplices) {
    $self->{splices_stats} = CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbEvents('splice'));
  }

  # If we are looking at reads alignement
  if($args{check_deletions} && $self->canCheckDeletions) {
    $self->{deletions_stats} = CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbEvents('deletion'));
  }

  # If we are looking at reads alignement
  if($args{check_insertions} && $self->canCheckInsertions) {
    $self->{insertions_stats} = CracTools::BenchCT::Stats->new(nb_elements => $self->checker->nbEvents('insertion'));
  }

  $self->_init(@_);

  return $self;
}

sub _init {
  my $self = shift;
}

sub checker {
  my $self = shift;
  return $self->{checker};
}

sub canCheckMapping {
  my $self = shift;
  return 0;
}

sub canCheckErrors {
  my $self = shift;
  return 0;
}

sub canCheckSnps {
  my $self = shift;
  return 0;
}

sub canCheckInsertions {
  my $self = shift;
  return 0;
}

sub canCheckDeletions {
  my $self = shift;
  return 0;
}

sub canCheckSplices {
  my $self = shift;
  return 0;
}

sub mappingStats {
  my $self = shift;
  return $self->{mapping_stats};
}

sub errorsStats {
  my $self = shift;
  return $self->{errors_stats};
}

sub snpsStats {
  my $self = shift;
  return $self->{snps_stats};
}

sub deletionsStats {
  my $self = shift;
  return $self->{deletions_stats};
}

sub insertionsStats {
  my $self = shift;
  return $self->{insertions_stats};
}

sub splicesStats {
  my $self = shift;
  return $self->{splices_stats};
}



1;
