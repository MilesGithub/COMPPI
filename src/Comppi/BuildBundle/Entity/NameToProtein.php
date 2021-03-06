<?php

namespace Comppi\BuildBundle\Entity;

use Doctrine\ORM\Mapping as ORM;

/**
 * @ORM\Entity
 * @ORM\Table(indexes={@ORM\index(name="search_idx", columns={"name", "namingConvention", "specieId"}), @ORM\index(name="synonym_idx", columns={"proteinId"})})
 */
class NameToProtein
{
    /**
     * @ORM\Id
     * @ORM\Column(type="integer")
     * @ORM\GeneratedValue(strategy="AUTO")
     */
    protected $id;

    /**
     * @ORM\Column(type="integer")
     */
    protected $specieId;

    /**
     * @ORM\Column(type="string", length="255")
     */
    protected $namingConvention;

    /**
     * @ORM\Column(type="string", length="255")
     */
    protected $name;

    /**
     * @ORM\Column(type="integer")
     * @ORM\ManyToOne(targetEntity="Protein")
     */
    protected $proteinId;

}