use strict;
use warnings;
package CracTools::BenchCT::Analyzer::VCF;
# ABSTRACT: Analyze results from a software experiment over the simulated data handled by BenchCT
 
use parent 'CracTools::BenchCT::Analyzer';

use CracTools::Utils;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type =~ /^(snp|insertion|deletion)$/) {
    return 1;
  }
  return 0;
}

sub _init {
  my $self = shift;
  my %args = @_;

  my $vcf_file = $args{file};
  my $vcf_it = CracTools::Utils::vcfFileIterator($vcf_file);
  while (my $vcf_line = $vcf_it->()) {
    # This is a hook line for subclasses
    $self->_processLine($vcf_line,$args{options});
  }
}

# We do nothing here... childs will.
sub _processLine {
  my $self = shift;
  my $vcf_line = shift;
  my $options = shift;

  # Options
  my $qual = $vcf_line->{qual};
  return 0 if defined $options->{min_qual} && $qual ne '.' && $qual < $options->{min_qual};

  my $ref_length = length $vcf_line->{ref};
  foreach my $alt (@{$vcf_line->{alt}}) {
    my $alt_length = length $alt;
    my ($type,$pos,$length);
    my $true_mutation;
    
    # First we try to determinate the type of the mutation
    # It is a substitution
    if($ref_length == $alt_length) {
      $type = 'snp';
      # We shift the pos if the reference has more than one nucleotide
      $pos = $vcf_line->{pos} + ($ref_length - 1);
      next if $alt =~ /N/i;
      $true_mutation = $self->checker->isTrueSNP($vcf_line->{chr},$pos,$alt);
    # This is a deletion
    } elsif($ref_length > $alt_length) {
      $type = 'deletion';
      $length = $ref_length - $alt_length;
      $pos = $vcf_line->{pos} + ($alt_length - 1);
      $true_mutation = $self->checker->isTrueDeletion($vcf_line->{chr},$pos,$length);
    # This is an insertion
    } else {
      $type = 'insertion';
      $length = $alt_length - $ref_length;
      $pos = $vcf_line->{pos} + ($ref_length - 1);
      $true_mutation = $self->checker->isTrueInsertion($vcf_line->{chr},$pos,$length);
    }

    # Pos -1 because checker is 0 based and VCF is 1-based
    $pos--;

    # Now we check the validity of the event
    if(defined $self->getStats($type)) {
      if($true_mutation) {
        $self->getStats($type)->addTruePositive(id => $true_mutation, out_string =>
          join("\t",
            $vcf_line->{chr},
            $vcf_line->{pos},
            $vcf_line->{ref},
            $alt,
          ),
        );
      } else {
        $self->getStats($type)->addFalsePositive(out_string => join("\t",
            $vcf_line->{chr},
            $vcf_line->{pos},
            $vcf_line->{id},
            $vcf_line->{ref},
            $alt,
            $vcf_line->{qual},
            $vcf_line->{filter},
            defined $vcf_line->{info}? join(";",map { defined $vcf_line->{info}->{$_}? "$_=".$vcf_line->{info}->{$_} : "$_="} keys %{$vcf_line->{info}}) : "",
          ),
        );
      }
    }
  }
}

1;
