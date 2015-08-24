#!/usr/bin/perl

use strict;

my $type = shift(@ARGV);
$::hb = "/acrm/usr/local/share/hb";

if(($type eq "N") || ($type eq "n"))
{
    $::chb = "$::hb/bin/checkhbond_Ndonor -c 0.5 -m $::hb/data/hbmatricesS35_SCMC.dat -n $::hb/data/hbmatricesS35_N.dat";
}
elsif(($type eq "O") || ($type eq "o"))
{
    $::chb = "$::hb/bin/checkhbond_Oacceptor -c 0.5 -m $::hb/data/hbmatricesS35_SCMC.dat -n $::hb/data/hbmatricesS35_O.dat";
}
elsif(($type eq "S") || ($type eq "s"))
{
    $::chb = "$::hb/bin/checkhbond -c 0.5 -m $::hb/data/hbmatricesS35.dat";
}
else
{
    die "Must specify type: calce.pl [N|O|S] datafile";
}

my($file);

$file = "";
while(<>)
{
    chomp;
    if(/^INFO: Processing file (.*)/)
    {
        $file = $1;
    }
    elsif(/^WARNING:/ || /^INFO: / || /^Warning/ || /Unable/)
    {
        ;
    }
    else                        # It's an HBond!
    {
        Process($file, $_);
    }
}

sub Process
{
    my($file, $hbond) = @_;
    my($chain, $res1, $res2, $energy, $type);
    my($name) = $file;
    $name =~ s/\/.*\///;        # Remove up to last /

    # We can't process those with a file length of 4 since they are
    # single chains and we don't know which chain
    if(length($name) > 4)
    {
        # If it's a split domain file, we need to find the chain
        if(!($name =~ /\.ent/))
        {
            $chain = substr($name, 4, 1);
            $chain = "" if($chain eq "0");
        }
        ($res1, $res2, $type) = GetResidues($hbond, $chain);
        $energy = RunHBond($file, $res1, $res2, $type);
    }
}

sub RunHBond
{
    my($file, $res1, $res2, $type) = @_;
    my($results);

#    print STDERR "$::chb $file $res1 $res2 $type\n";
    $results = `$::chb $file $res1 $res2 $type`;
    print $results;
}

sub GetResidues
{
    my($hbond, $chain) = @_;
    my($res1, $res2, $type, @fields);

    @fields = split(/\s+/, $hbond);
    $res1 = $chain . $fields[1];
    $type = $fields[4];
    $res2 = $chain . $fields[5];

    return($res1, $res2, $type);
}
