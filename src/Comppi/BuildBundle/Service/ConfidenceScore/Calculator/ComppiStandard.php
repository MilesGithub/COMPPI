<?

namespace Comppi\BuildBundle\Service\ConfidenceScore\Calculator;

use Doctrine\DBAL\Types\IntegerType as IntegerParameter;

class ComppiStandard implements CalculatorInterface
{
    private $localizationTranslator;
    private $conn;
    
    public function __construct($id) {
        $this->id = $id;
		ProteinScoreCalculator::$classWeights[0]=0.95;
		ProteinScoreCalculator::$classWeights[1]=0.5;
		ProteinScoreCalculator::$classWeights[2]=0.5;
    }

    public function setLocalizationTranslator($translator) {
        $this->localizationTranslator = $translator;
    }

    public function calculate(\Doctrine\DBAL\Connection $connection) 
    {    
        $this->conn=$connection; //needed for CalculateLinkConfidence;
        $locQuery=$connection->prepare('SELECT proteinId,localizationId,specieId,confidenceType,protLocId,systemTypeId FROM ProtLocToSystemType LEFT JOIN ProteinToLocalization on ProtLocToSystemType.protLocId=ProteinToLocalization.id LEFT JOIN Protein on Protein.id=ProteinToLocalization.proteinId LEFT JOIN SystemType on SystemType.id=ProtLocToSystemType.systemTypeId WHERE proteinId=?');
        echo (" CompPPI Standard calculator init...\n");
        $this->initCalculation($connection);
        //for each protein
        ProteinScoreCalculator::$fullProteinList=array();
        echo (" Adding proteins from database...\n");
        $protRes=$connection->query("SELECT DISTINCT id FROM Protein");
        $protIDs=$protRes->fetchAll(\Doctrine\ORM\AbstractQuery::HYDRATE_ARRAY);
        $num=0;
        foreach($protIDs as $protein)
        {
         $num++;
         $locQuery->bindValue(1,$protein['id'], IntegerParameter::INTEGER);
         $locQuery->execute();
         $localizations=$locQuery->fetchAll(\Doctrine\ORM\AbstractQuery::HYDRATE_ARRAY);
         $proteinCalc=new ProteinScoreCalculator();
         if(!$proteinCalc->addEntries($localizations))
            {echo("Warning: No localization information for protein #".$protein['id']."!\n"); continue;}
         ProteinScoreCalculator::$fullProteinList[$protein['id']]=$proteinCalc;
         if($num%10000==0) echo($num."/".count($protIDs)." proteins added.\n");
        }
        
        echo("Calculating localization evidence scores...\n");
        foreach(ProteinScoreCalculator::$fullProteinList as $protein)
         $protein->calculateMyLocationScores();
        
        $insert = $connection->prepare(
            'INSERT INTO ConfidenceScore(interactionId, calculatorId, score)' .
            ' VALUES(?, ?, ?)'
        );
        
        $insert->bindValue(2, $this->id, IntegerParameter::INTEGER);

        //no cross-species links
        $interactionSelect = $connection->prepare(
            'SELECT Interaction.id as id FROM Interaction ORDER BY id ASC LIMIT ?, ?'
        );

        $interactionOffset = 0;
        $blockSize = 500;

        $interactionSelect->bindValue(1, $interactionOffset, IntegerParameter::INTEGER);
        $interactionSelect->bindValue(2, $blockSize, IntegerParameter::INTEGER);

        $interactionSelect->execute();
        
        echo("Calculating and inserting confidence scores...\n");
        
        while ($interactionSelect->rowCount() > 0) {
            $interactions = $interactionSelect->fetchAll(\Doctrine\ORM\AbstractQuery::HYDRATE_ARRAY);

            $connection->beginTransaction();

            foreach ($interactions as $interaction) {

				echo("DEBUG: Memory usage (step 1):" . memory_get_usage()."\n");
				$score = 0;
                $score = $this->CalculateLinkConfidence($interaction['id']);
				echo("DEBUG: Memory usage (step 2):" . memory_get_usage()."\n");
                $insert->bindValue(1, $interaction['id']);
                $insert->bindValue(3, $score);
                $insert->execute();
				echo("DEBUG: Memory usage (step 3):" . memory_get_usage()."\n\n");
            }

            $connection->commit();

            // advance cursor
            $interactionOffset += $blockSize;

            $interactionSelect->closeCursor();
            $interactionSelect->bindValue(1, $interactionOffset, IntegerParameter::INTEGER);
            $interactionSelect->execute();
        }
        
        echo("Confidence calculation complete.\n");
    }

    public function getName() {
        return "ComPPI Standard";
    }

    private function initCalculation(\Doctrine\DBAL\Connection $connection) {
        // @TODO this hack is required here because of a PDO bug
        // https://bugs.php.net/bug.php?id=44639
        $connection->getWrappedConnection()->setAttribute(\PDO::ATTR_EMULATE_PREPARES, false);
        
        ProteinScoreCalculator::$localizationTranslator=$this->localizationTranslator;
        ProteinScoreCalculator::$compartments=array_keys($this->localizationTranslator->getLargelocs());
        }
    
    private function CalculateLinkConfidence($linkID)
    {
     $interactionQuery=$this->conn->prepare('SELECT actorAId,actorBId,specieId FROM Interaction LEFT JOIN Protein on Interaction.actorAId=Protein.id WHERE Interaction.id=?');
     $interactionQuery->bindValue(1,$linkID,IntegerParameter::INTEGER);
     $interactionQuery->execute();
     $interactionResult=$interactionQuery->fetchAll(\Doctrine\ORM\AbstractQuery::HYDRATE_ARRAY);
     if((!isset(ProteinScoreCalculator::$fullProteinList[$interactionResult[0]['actorAId']])) or (!isset(ProteinScoreCalculator::$fullProteinList[$interactionResult[0]['actorBId']])))
     {echo("Protein for interaction ".$linkID." not found in database (or no localization info), skipping.\n");return 0;}
     $startLocScores=ProteinScoreCalculator::$fullProteinList[$interactionResult[0]['actorAId']]->getMyLocationScores();
     $endLocScores=ProteinScoreCalculator::$fullProteinList[$interactionResult[0]['actorBId']]->getMyLocationScores();
     $linkLocalizationScores=array();
     $value=1;     
     foreach (ProteinScoreCalculator::$compartments as $compartment)
	 {
        $value*=(1-$startLocScores[$compartment]*$endLocScores[$compartment]);
		#echo("DEBUG: compartment {$compartment}: startLocScore:{$startLocScores[$compartment]} -- endLocScore:{$endLocScores[$compartment]} --> value (current round): {$value}\n");
	}

     $LocalizationConfidence=1-$value;
	 
     //currently no interaction confidence limitation
     $InteractionConfidence=1;
     #echo("DEBUG: Confidence score for link {$linkID} is: {$LocalizationConfidence}\n");
     return $LocalizationConfidence*$InteractionConfidence;
    }

}

class ProteinScoreCalculator
{
  public $id;
  
  public static $compartments=array();
  public static $classWeights;
  public static $fullProteinList;
  
  private $connection;
  
  private $entries; 
  public $locationScores;
    
  public static $localizationTranslator;
  
  public function addEntries($SqlAssocArray)
  {
   $rowcount=0;  
   foreach ($SqlAssocArray as $row)
   {
    try{
    $compartment=ProteinScoreCalculator::$localizationTranslator->getLargelocById($row['localizationId']);	
    }
    catch(\InvalidArgumentException $e)
	{
    echo("Warning: could not determine compartment (protein {$row['proteinId']}, loc {$row['localizationId']}, evidence {$row['protLocId']})!\n");
    continue;
    }	
	$rowcount++;
    if (!isset($this->entries[$compartment][$row['confidenceType']])) $this->entries[$compartment][$row['confidenceType']]=array();
	if (!isset($this->entries[$compartment][$row['confidenceType']][$row['systemTypeId']])) $this->entries[$compartment][$row['confidenceType']][$row['systemTypeId']]=1; //only 1 per experimental system type
	//same experimental system type is likely from the same experiment in different DBs
   }
   if($rowcount==0) //no localization information, empty set received
    return FALSE;
   #echo("DEBUG: Inserted {$rowcount} entries.\n");
   return TRUE;
  }

  public function calculateMyLocationScores()
  {
   $this->locationScores=array();
   foreach(ProteinScoreCalculator::$compartments as $compartment)
    {     
     if(!isset($this->entries[$compartment]))
     {$this->locationScores[$compartment]=0;continue;}
     
     $score=1;
     foreach($this->entries[$compartment] as $entryClass => $entryNum)
        $score*=pow( 1 - ProteinScoreCalculator::$classWeights[$entryClass] , count($entryNum) );
     $this->locationScores[$compartment]=1 - $score;
    }
  }

  public function getMyLocationScores()
  {
   return $this->locationScores;
  }
}
