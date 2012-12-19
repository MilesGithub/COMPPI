<?php

namespace Comppi\BuildBundle\Service\DatabaseProvider\Parser\Map;

class UniprotFullname extends AbstractMapParser
{
    protected $altNameBlacklist = array(
        'Fragment',
        'Fragments'
    );

    protected $headerCount = 1;

    protected static $parsableFileNames = array(
        'uniprot_hs_fullname.tab',
        'uniprot_dm_fullname.tab',
        'uniprot_sc_fullname.tab',
        'uniprot_ce_fullname.tab'
    );

    protected $recordReady = array();

    protected function readRecord() {
        if (!empty($this->recordReady)) {
            $this->currentRecord = array_shift($this->recordReady);
            return;
        }

        $line = $this->readLine();

        if ($line === false) {
            // EOF
            return;
        }

        $recordArray = explode("\t", $line);
        $this->checkRecordFieldCount($recordArray, 4);

        $firstParenPos = strpos($recordArray[3], '(');

        if ($firstParenPos !== false) {
            // alt name found, strip "full name" (the first one)
            $fullName = substr($recordArray[3], 0, $firstParenPos - 1);

            // extract alt names
            $altNameString = substr($recordArray[3], $firstParenPos+1, -1);
            $altNames = explode(') (', $altNameString);

            foreach ($altNames as $altName) {
                if (in_array($altName, $this->altNameBlacklist)) {
                    continue;
                }

                $this->recordReady[] = array(
                    'namingConventionA' => 'UniProtKB-AC',
                    'namingConventionB'	=> 'UniProtAlt',
                    'proteinNameA'	=> $recordArray[0],
                    'proteinNameB'	=> $altName
                );
            }
        } else {
            $fullName = $recordArray[3];
        }

        $this->recordReady[] = array(
            'namingConventionA' => 'UniProtKB-AC',
            'namingConventionB'	=> 'UniProtFull',
            'proteinNameA'	=> $recordArray[0],
            'proteinNameB'	=> $fullName
        );


        if ($recordArray[2] == 'reviewed') {
            $this->currentRecord = array(
                'namingConventionA' => 'UniProtKB-AC',
                'namingConventionB'	=> 'UniProtKB/Swiss-Prot',
                'proteinNameA'	=> $recordArray[0],
                'proteinNameB'	=> $recordArray[0]
            );
        } else {
            $this->currentRecord = array(
                'namingConventionA' => 'UniProtKB-AC',
                'namingConventionB'	=> 'UniProtKB/TrEmbl',
                'proteinNameA'	=> $recordArray[0],
                'proteinNameB'	=> $recordArray[0]
            );
        }
    }
}