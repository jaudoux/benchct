use strict;
use warnings;
package CracTools::BenchCT::Analyzer::VCF;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
#
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;


sub canCheckSnps {
  my $self = shift;
  return 1;
}

sub canCheckInsertions {
  my $self = shift;
  return 1;
}

sub canCheckDeletions {
  my $self = shift;
  return 1;
}

sub _init {
  my $self = shift;
  my %args = @_;

  my %stats_objects = ( snp => $self->snpsStats,
    deletion => $self->deletionsStats,
    insertion => $self->insertionsStats,
  );

  $self->{stats_objects} = \%stats_objects;

  my $vcf_file = $args{file};
  my $vcf_it = CracTools::Utils::vcfFileIterator($vcf_file);
  while (my $vcf_line = $vcf_it->()) {
    # This is a hook line for subclasses
    $self->_processLine($vcf_line);
  }
}

sub getStats {
  my $self = shift;
  my $type = shift;
  return $self->{stats_objects}->{$type};
}

# We do nothing here... childs will.
sub _processLine {
  my $self = shift;
  my $vcf_line = shift;
  my $ref_length = length $vcf_line->{ref};
  foreach my $alt (@{$vcf_line->{alt}}) {
    my $alt_length = length $alt;
    my ($type,$pos,$length);
    
    # First we try to determinate the type of the mutation
    # It is a substitution
    if($ref_length == $alt_length) {
      $type = 'snp';
      # We shift the pos if the reference has more than one nucleotide
      $pos = $vcf_line->{pos} + ($ref_length - 1);
      $length = 1;
    # This is a deletion
    } elsif($ref_length > $alt_length) {
      $type = 'deletion';
      $length = $ref_length - $alt_length;
      $pos = $vcf_line->{pos} + ($alt_length - 1);
    # This is an insertion
    } else {
      $type = 'insertion';
      $length = $alt_length - $ref_length;
      $pos = $vcf_line->{pos} + ($ref_length - 1);
    }

    # Now we check the validity of the event
    if($self->checker->isTrueMutation($type,$vcf_line->{chr},$pos,$length)) {
      $self->getStats($type)->addTruePositive();
    } else {
      $self->getStats($type)->addFalsePositive();
    }
  }
}

1;
