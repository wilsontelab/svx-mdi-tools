use strict;
use warnings;

use vars qw(%options);
  
push @{$options{consensus}}, (
    {           
        long        => 'downsample-n',
        short       => 'D',
        type        => 'int',            
        required    => 0,
        default     => 11,
        message     => "use at most this many read pairs during consensus making",
    },
    {           
        long        => 'min-score-factor',
        short       => 's',
        type        => 'dbl',            
        required    => 0,
        default     => 0.7,
        message     => "MIN_SCORE_FACTOR * length = minimum acceptable SW score",
    },
    {           
        long        => 'consensus-factor',
        short       => 'C',
        type        => 'dbl',            
        required    => 0,
        default     => 0.667,
        message     => "fraction of bases that must have a common value to be a consensus",
    },
);

