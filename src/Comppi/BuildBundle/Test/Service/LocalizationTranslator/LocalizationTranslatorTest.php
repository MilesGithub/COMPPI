<?php

namespace Comppi\BuildBundle\Test\Service\LocalizationTranslator;

use Comppi\BuildBundle\Service\LocalizationTranslator\LocalizationTranslator;

class LocalizationTranslatorTest extends \PHPUnit_Framework_TestCase
{
    /**
     * @var Comppi\BuildBundle\Service\LocalizationTranslator\LocalizationTranslator
     */
    protected $translator;

    protected function setUp() {
        $loctreeFile = __DIR__ . DIRECTORY_SEPARATOR . 'loctree.textile';
        $largelocFile = __DIR__ . DIRECTORY_SEPARATOR . 'largelocs.yml';
        $this->translator = new LocalizationTranslator($loctreeFile, $largelocFile);
    }

    public function testGetIdByLocalization() {
        $this->assertEquals(
            1,    // expected
            $this->translator->getIdByLocalization('GO:0005575')    // actual
        );

        $this->assertEquals(
            2,
            $this->translator->getIdByLocalization('GO:0032991')
        );

        $this->assertEquals(
            5,
            $this->translator->getIdByLocalization('GO:0030686')
        );

        $this->assertEquals(
            9,
            $this->translator->getIdByLocalization('GO:0043234')
        );

        $this->assertEquals(
            10,
            $this->translator->getIdByLocalization('GO:0000151')
        );

        $this->assertEquals(
            23,
            $this->translator->getIdByLocalization('secretory_pathway')
        );
    }

    public function testGetLocalizationById() {
        $this->assertEquals(
            'GO:0030686',
            $this->translator->getLocalizationById(5)
        );

        $this->assertEquals(
            'GO:0005575',
            $this->translator->getLocalizationById(1)
        );

        $this->assertEquals(
            'GO:0000151',
            $this->translator->getLocalizationById(10)
        );

        $this->assertEquals(
            'secretory_pathway',
            $this->translator->getLocalizationById(23)
        );
    }

    public function testGetSecondaryIdByLocalization() {
        $this->assertEquals(
            6,
            $this->translator->getSecondaryIdByLocalization('GO:0030686')
        );

        $this->assertEquals(
            7,
            $this->translator->getSecondaryIdByLocalization('GO:0030684')
        );

        $this->assertEquals(
            8,
            $this->translator->getSecondaryIdByLocalization('GO:0030529')
        );

        $this->assertEquals(
            22,
            $this->translator->getSecondaryIdByLocalization('GO:0005575')
        );

        $this->assertEquals(
            24,
            $this->translator->getSecondaryIdByLocalization('secretory_pathway')
        );
    }

    public function testGetHumanReadableLocalizationById() {
        $this->assertEquals(
            'protein complex',
            $this->translator->getHumanReadableLocalizationById(9)
        );

        $this->assertEquals(
            'ubiquitin ligase complex',
            $this->translator->getHumanReadableLocalizationById(10)
        );

        $this->assertEquals(
            'secretory pathway',
            $this->translator->getHumanReadableLocalizationById(23)
        );
    }

    public function testGetLargelocs() {
        $this->assertEquals(
            array(
                'cellular-component' => array(1),
                'macromolecular-complex' =>
                    array(2, 3, 4, 5, 9, 10, 11, 12, 14, 18)
            ),
            $this->translator->getLargelocs()
        );
    }

    public function testGetLargelocById() {
        $this->assertEquals(
            'cellular-component',
            $this->translator->getLargelocById(1)
        );

        $this->assertEquals(
            'macromolecular-complex',
            $this->translator->getLargelocById(2)
        );

        $this->assertEquals(
            'macromolecular-complex',
            $this->translator->getLargelocById(3)
        );

        $this->assertEquals(
            'macromolecular-complex',
            $this->translator->getLargelocById(4)
        );

        $this->assertEquals(
            'macromolecular-complex',
            $this->translator->getLargelocById(5)
        );

        $this->assertEquals(
            'macromolecular-complex',
            $this->translator->getLargelocById(10)
        );
    }

    /**
     * @expectedException \InvalidArgumentException
     */
    public function testInvalidgetLargelocById() {
        $this->translator->getLargelocById(23);
    }
}