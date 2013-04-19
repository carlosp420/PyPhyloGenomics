#!/usr/bin/env python

import urllib2
import urllib

import sys

#Usage:  MUSCLE.designPrimers("alignment.fas")


"""
parameters needed, default between brackets:
	fasta_file     :
	tm             : [55]
	min_amplength  : [100]
	max_amplength  : [500]
	mode           : [primers]
	gencode        : [universal]
	clustype       : [dna]
"""

def designPrimers(aln, temp, min_length, max_length, gencode, clustype):
	url = 'http://floresta.eead.csic.es/primers4clades/primers4clades.cgi'
	params = {
		'sequencefile':		aln,
		'tm':				temp,
		'min_amplength':	min_length,
		'max_amplength':	max_length,
		'mode':				'primers',
		'gencode':			gencode,
		'clustype':			clustype
			};
	data = urllib.urlencode(params);
	req = urllib2.urlopen(url, data);
	response = req.read();
	print response
"""
############
$ch = curl_init();

$file = $argv[1];

$data = array('sequencefile' => '@' . $argv[1], 
			  'tm' => '55',
			  'min_amplength' => '100', # may need to change later. Minimum accepted amplicon length
			  'max_amplength' => '500', # may need to change later. Minimum accepted amplicon length
			  'mode' => 'primers',
			  'gencode' => 'universal',
			  'clustype' => 'dna'
			  
			  );

curl_setopt($ch, CURLOPT_URL, 'http://floresta.eead.csic.es/primers4clades/primers4clades.cgi');
curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
curl_setopt($ch, CURLOPT_POST, true);
curl_setopt($ch, CURLOPT_POSTFIELDS, $data);

$result = curl_exec($ch);
$result = explode("\n", $result);

$primers = array();
foreach( $result as $line ) {
	$pattern = '/.+degen_corr/';
	preg_match($pattern, $line, $matches);
	if( count($matches) > 0 ) {
		$primers[] = $matches[0];
	}
}


$primer_1 = explode(" ", $primers[0]);
$primer_2 = explode(" ", $primers[1]);

$primer_F = $primer_1[0];
$primer_R = $primer_2[0];

echo $primer_F;
echo "\n";
echo $primer_R;
curl_close($ch);
"""
