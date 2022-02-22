use strict;
use warnings;

# utility functions available to perl-based pipeline scripts

use vars qw($action);
sub getActionName { $action || "unknown" }

#----------------------------------------------------------
# error handling
#----------------------------------------------------------
sub throwError {
    my ($msg) = @_;
    die getActionName()." error: $msg\n";
}

#----------------------------------------------------------
# facilitate conversion of environment variables to named variables
#----------------------------------------------------------
sub fillEnvVar {
    my ($var, $varName, $allowNull, $default) = @_;
    my $val = $ENV{$varName};
    $val and $val eq 'null' and $val = undef;
    defined $val or $allowNull or die getActionName()." error: missing environment variable: $varName\n";
    !defined $val and defined $default and $val = $default;
    $$var = $val;
}

#----------------------------------------------------------
# summary count reporting to STDERR and log file
#----------------------------------------------------------
my $countFile;
sub resetCountFile { # called by parent to initialize log
    my $action_ = getActionName();
    $countFile = "$ENV{LOG_FILE_PREFIX}.$action_.log.txt";
    unlink $countFile;
}
sub printCount { # used by parent to print to report and log file
    $ENV{SUPPRESS_PRINT_COUNT} and return;
    my ($val, $name, $msg) = @_;
    my $action_ = getActionName();
    open my $outH, ">>", $countFile or die "$action_ error: could not open $countFile $!\n";
    my $report = join("\t", "$action_", $val || 0, $name, $msg)."\n";
    print STDERR $report;
    print $outH $report;
    close $outH;
}

#----------------------------------------------------------
# enable multi-threading
#----------------------------------------------------------
# https://docstore.mik.ua/orelly/perl/cookbook/ch16_11.htm
our (@readH, @writeH, %pids);
sub launchChildThreads { # initialize children in parent
    my ($childSub) = @_; # where parent distributes elements from one input data stream
    foreach my $childN(1..$ENV{N_CPU}) {
        pipe($readH[$childN], $writeH[$childN]); # establish the read/write pipe for sending data to child
        my $pid = fork();
        die "Failed to fork: $!\n" unless defined $pid;
        if (!$pid) {
            close $writeH[$childN]; # close all instances of all handles
            &$childSub($childN);
            close $readH[$childN];  # not strictly necessary, but for clarity
            exit;
        }
        close $readH[$childN];
        $pids{$pid} = $childN;
    }    
}
sub finishChildThreads { # called when parent has finished its work
    my ($suppressFileHClose) = @_;
    unless($suppressFileHClose){ # true if caller has already closed the handles
        foreach my $childN(1..$ENV{N_CPU}) { # close the children write handles to flush all data
            close $writeH[$childN];
        }        
    }
    while(scalar(keys %pids)){ # wait for all children to finish their work
        my $pid = wait();
        delete $pids{$pid}; 
    } 
}
sub runChildThreads_array {      # initialize children in parent
    my ($childSub, @array) = @_; # where every child operates on a single input value (e.g. chrom)
    my $nThreads = $ENV{N_CPU};
    $nThreads > @array and $nThreads = @array;
    my $lastChildI;
    foreach my $childN(1..$nThreads) {
        my $childI = $childN - 1;
        forkChild($childSub, $childI);
        $lastChildI = $childI;
    }
    while(scalar(keys %pids)){ # wait for all children to finish their work
        my $pid = wait();
        delete $pids{$pid};
        if($lastChildI < $#array){ # and launch a new process for next item, until done
            $lastChildI++;
            forkChild($childSub, $lastChildI);
        }
    }
    sub forkChild {
        my ($childSub, $childI) = @_;
        my $pid = fork();
        die "Failed to fork: $!\n" unless defined $pid;
        if (!$pid) {
            &$childSub($childI);
            exit;
        }
        $pids{$pid} = $childI;        
    }
}

1;
