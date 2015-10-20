use strict;
use warnings;
package CracTools::BenchCT::Analyzer::BED::Junction;
# ABSTRACT: Analyze junction bed fils (as tophat produce)
 
use parent 'CracTools::BenchCT::Analyzer::BED';

use CracTools::Utils;
#use Data::Dumper;

sub canCheck {
  my $self = shift;
  my $event_type = shift;
  if($self->SUPER::canCheck($event_type) ||  $event_type eq 'splice') {
    return 1;
  }
  return 0;
}

sub _processLine {
  my $self = shift;
  my $bed_line = shift;
  # Check if the bed line has blocks
  # otherwise start and end pos are the splice positions
  if(@{$bed_line->{blocks}}) {
    # Loop over blocks
    for(my $i=1; $i < @{$bed_line->{blocks}}; $i++) {

      # Extract splice coordinates
      my $chr = $bed_line->{chr};
      my $strand = CracTools::Utils::convertStrand($bed_line->{strand});
      my $start = $bed_line->{blocks}[$i-1]->{ref_end};
      my $end = $bed_line->{blocks}[$i]->{ref_start};
      my $length = $end - $start + 1;
      my $true_splice = $self->checker->isTrueSplice($chr,$start,$length,$strand);
      
      if($true_splice) {
        $self->getStats('splice')->addTruePositive(id => $true_splice);
      } else {
        #print STDERR Dumper($bed_line);
        $self->getStats('splice')->addFalsePositive(out_string => join("\t",
            $chr,
            $start,
            $end,
            $strand,
            $bed_line->{name},
        ));
      }
    }
  } else {
    my $length      = $bed_line->{end} - $bed_line->{start} + 1;
    my $true_splice = $self->checker->isTrueSplice(
      $bed_line->{chr},
      $bed_line->{start},
      $length,
      CracTools::Utils::convertStrand($bed_line->{strand}),
    );
    
    if($true_splice) {
      $self->getStats('splice')->addTruePositive(id => $true_splice);
    } else {
      #print STDERR Dumper($bed_line);
      $self->getStats('splice')->addFalsePositive(out_string => join("\t",
          $bed_line->{chr},
          $bed_line->{start},
          $bed_line->{end},
          $bed_line->{strand},
          $bed_line->{name},
      ));
    }
  }
}

1;
