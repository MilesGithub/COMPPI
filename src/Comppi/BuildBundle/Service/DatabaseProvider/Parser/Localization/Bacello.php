<?php

namespace Comppi\BuildBundle\Service\DatabaseProvider\Parser\Localization;

class Bacello extends AbstractLocalizationParser
{
    protected static $parsableFileNames = array(
        'pred_sce'
    );
    
    protected $localizationToGoCode = array (
        'Cytoplasm' => 'GO:0005737',
        'Secretory' => 'secretory_pathway',
        'Mitochondrion' => 'GO:0005739',
        'Nucleus' => 'GO:0005634'
    );
    
    protected $hasHeader = false;
    
    protected function readRecord() {
        $line = $this->readLine();
        if ($line === false) {
            // EOF
            return;
        }
        
        $recordArray = preg_split('/ +/', $line);
        $this->checkRecordFieldCount($recordArray, 2);
        
        $this->currentRecord = array(
            'proteinId' => $recordArray[0],
            'namingConvention' => 'EnsemblPeptideId',
            'localization' => $this->getGoCodeByLocalizationName($recordArray[1]),
            'pubmedId' => 16873501,
            'experimentalSystemType' => 'SVM decision tree (predicted)'
        );
    }
}