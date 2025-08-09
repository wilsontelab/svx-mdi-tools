#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use File::Basename;

#========================================================================
# create an alias, i.e., named shortcut, to this suite's 'run' utility
# can have multiple aliases to the same script, but only one alias of a given name
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
my $usage = "perl alias.pl ALIAS_NAME\nadds an alias to './run', called ALIAS NAME, to ~/.bashrc";
#========================================================================

#========================================================================
# main execution block
#------------------------------------------------------------------------
my $MDI_CENTRIC   = "mdi-centric";
my $SUITE_CENTRIC = "suite-centric";
my $SUITE_MODE    = $MDI_CENTRIC;
my $MDI_DIR       = "../../..";
if( !-f "$MDI_DIR/mdi" or !-d "$MDI_DIR/frameworks" or !-d "$MDI_DIR/suites" ) {
    $SUITE_MODE = $SUITE_CENTRIC; # i.e., this is the top-level of a single-suite installation
}
if( $SUITE_MODE eq $MDI_CENTRIC ) {
    print STDERR "\nNothing to do.\n\n";
    print STDERR "This copy of the tool suite is part of an MDI installation in directory:\n    $MDI_DIR\n\n";
    print STDERR "This alias.pl script is only useful in the top directory of a single-suite installation.\n\n";
    exit 1;
}

# parse the options and apply defaults
my ($alias) = @ARGV;
$alias or throwError("missing alias name");
my $bashrc = "$ENV{HOME}/.bashrc";
my $suiteDir = $ENV{PWD};
my $aliasCommand = "alias $alias=\"$suiteDir/run\"";
my $outLine = "$aliasCommand # written by MDI alias\n";

# get user permission to modify their profile
getPermissionGeneral(
    "The following line:\n".
    "    $outLine".  
    "will be written to file:\n".
    "    $bashrc\n"
);

# check the profile file path
if(!-f $bashrc){
    open HANDLE, ">>$bashrc" or die "touch failed: $bashrc: $!\n";
    close HANDLE;
} 
-f $bashrc or throwError("file not found:\n    $bashrc");

# collect the contents of the current file as an array of lines
my ($replaced, @profile);    
open my $inH, "<", $bashrc or die "could not read file: $bashrc: $!\n";
while (my $line = <$inH>){
    $line =~ m/\n/ or $line = "$line\n"; # guard against incomplete final line
    $line eq $outLine and exit; # nothing to do, exit quietly
    if($line =~ m/^alias\s+$alias=/){ 
        getPermissionGeneral(
            "Alias '$alias' already exists and will be overwritten from:\n". 
            "    $line".
            "to:\n".
            "    $outLine"
        );
        push @profile, $outLine; 
        $replaced = 1;
    } else {
        push @profile, $line;
    }
}
close $inH;
$replaced or push @profile, $outLine;  

# print the new file
my $buFile = "$bashrc.mdiAliasBackup";
-e $buFile or copy($bashrc, $buFile);
open my $outH, ">", $bashrc or die "could not write file: $bashrc: $!\n";
print $outH join("", @profile);
close $outH;

# give user feedback on how to activate the alias in their current shell
print "\ndone\n";
print "\nEither log out and back in or execute the following\n";
print   "command to activate the alias in this shell also:\n\n$aliasCommand\n\n";

#========================================================================

#========================================================================
# support functions
#------------------------------------------------------------------------
sub throwError {
    my ($message) = @_;
    print "\n$message\nusage: $usage\n\n";
}
sub getPermissionGeneral {
    my ($msg) = @_;
    my $leader = "-" x 80;
    print "\n$leader\n$msg\n";  
    print "Continue? [y|N]: ";
    my $permission = <STDIN>;
    chomp $permission;
    $permission = "\U$permission";
    ($permission eq 'YES' or $permission eq 'Y') and return 1;
    print "aborting with no action taken\n";
    exit 1;
}
#========================================================================
