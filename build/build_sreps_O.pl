#!/acrm/usr/local/bin/perl
#*************************************************************************
#
#   Program:    
#   File:       
#   
#   Version:    
#   Date:       
#   Function:   
#   
#   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2005
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   Phone:      +44 (0)171 679 7034
#   EMail:      andrew@bioinf.org.uk
#               martin@biochem.ucl.ac.uk
#   Web:        http://www.bioinf.org.uk/
#               
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************
use strict;

$::cathfile = "/acrm/data/cath/latest/CathDomainList.S35";
$::datadir  = "../data";
$::hbfile   = "$::datadir/hbmatricesS35_O.dat";
$::hbstem   = "$::datadir/hbmatricesS35_O_%s.dat";
$::hmatprog = "../bin/hydrogen_matrices_Oacceptor";

#*************************************************************************
my($version, $cversion, $matfile);

$version = GetCATHVersion($::cathfile);
$cversion = GetHBVersion($::hbfile);

if($version eq $cversion)
{
    print "Matrix file is up to date\n";
}
else
{
    $matfile = BuildMatrixFile($::cathfile, $version, $::hbstem);
    InstallFile($matfile, $::hbfile);
}

#*************************************************************************
sub InstallFile
{
    my($matfile, $hbfile) = @_;
#    print "ln -sf $matfile $hbfile\n";
    `ln -sf $matfile $hbfile`;
}

#*************************************************************************
sub BuildMatrixFile
{
    my($cathfile, $version, $hbstem) = @_;
    my($hbfile, $log);

    $hbfile = sprintf($hbstem, $version);
    system("$::hmatprog $cathfile $hbfile");

    return($hbfile);
}

#*************************************************************************
sub GetCATHVersion
{
    my($cathfile) = @_;
    my($line, @fields, $fullname, $version);

    if(-e $cathfile)
    {
        $line = `ls -l $cathfile`;
        @fields = split(/\s+/,$line);
        $fullname = $fields[10];
        $fullname =~ /\/(.*)\//;
        $version = $1;
    }
    else
    {
        die "Cath file does not exist: $cathfile";
    }

    return($version);
}

#*************************************************************************
sub GetHBVersion
{
    my($hbfile) = @_;
    my($line, @fields, $fullname, $version);

    if(-e $hbfile)
    {
        $line = `ls -l $hbfile`;
        @fields = split(/\s+/,$line);
        $fullname = $fields[10];
        $fullname =~ /_(.*).dat/;
        $version = $1;
    }
    else
    {
        $version = 0;
    }

    return($version);
}
