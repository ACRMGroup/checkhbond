#!/usr/bin/perl -s

use strict;

my $type = shift(@ARGV);
$::hb = "/acrm/usr/local/share/hb";
$::ehb = "/home/bsm/martin/bin/ehb3";
$::pdbhadd = "/home/bsm/martin/bin/pdbhadd";
$::hstrip = "/home/bsm/martin/bin/hstrip";

if(($type eq "N") || ($type eq "n"))
{
    $::chb = "$::hb/bin/checkhbond_Ndonor -c 0.25 -m $::hb/data/hbmatricesS35_SCMC.dat -n $::hb/data/hbmatricesS35_N.dat";
    die "Option N not supported";
}
elsif(($type eq "O") || ($type eq "o"))
{
    $::chb = "$::hb/bin/checkhbond_Oacceptor -c 0.25 -m $::hb/data/hbmatricesS35_SCMC.dat -n $::hb/data/hbmatricesS35_O.dat";
    die "Option O not supported";
}
elsif(($type eq "S") || ($type eq "s"))
{
    $::chb = "$::hb/bin/checkhbond -c 0.25 -m $::hb/data/hbmatricesS35.dat";
}
else
{
    die "Must specify type: calce.pl [N|O|S] datafile";
}

$::ehb .= " -r" if(defined($::r));
$::ehb .= " -o" if(defined($::o));

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
    my($chain, $res1, $res2, $energy, $type, $enthalpy);
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
        if($energy ne "ERROR")
        {
            ($res1, $res2) = GetResAtom($hbond, $chain);
            $enthalpy = RunECalc($file, $res1, $res2);
            
            print "$enthalpy $energy\n";
        }
    }
}

sub RunHBond
{
    my($file, $res1, $res2, $type) = @_;
    my($results);

#    print STDERR "$::chb $file $res1 $res2 $type\n";
    $results = `$::chb $file $res1 $res2 $type 2>/dev/null`;
    if($results =~ /Pseudoenergy/)
    {
        $results =~ /bond:\s+(.*)/;
        $results = $1;
    }
    else
    {
        $results = "ERROR";
    }
    return $results;
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


sub GetResAtom
{
    my($hb, $chain) = @_;
    my($res1,$res2);
    my(@fields) = split;

    $res1 = $chain . $fields[1] . "." . $fields[2];
    $res2 = $chain . $fields[5] . "." . $fields[6];

    return($res1, $res2);
}

sub RunECalc
{
    my($file, $res1, $res2) = @_;
    my($enthalpy);
    my($tfile) = "/tmp/ecalc_$$.pdh";

    `$::hstrip $file | $::pdbhadd > $tfile 2>/dev/null`;

    $enthalpy = `$::ehb $tfile $res1 $res2 2>/dev/null`;
    chomp $enthalpy;
    return($enthalpy);
}
