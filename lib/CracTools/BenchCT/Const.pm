use strict;
use warnings;
package CracTools::BenchCT::Const;
# ABSTRACT: BenchCT constants


# Tolerance on mapping position
our $THRESHOLD_MAPPING  = 5;
our $THRESHOLD_SNP      = 5;
our $THRESHOLD_INS      = 5;
our $THRESHOLD_DEL      = 5;
our $THRESHOLD_CHIMERA  = 20;
our $THRESHOLD_ERR      = 5;
our $THRESHOLD_SPLICE   = 5;

our $NB_DECIMALS        = 4;

1;
