SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL';

DROP SCHEMA IF EXISTS `proteindatabase` ;
CREATE SCHEMA IF NOT EXISTS `proteindatabase` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci ;
USE `proteindatabase` ;

-- -----------------------------------------------------
-- Table `proteindatabase`.`goinfo`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`goinfo` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`goinfo` (
  `GOTermID` INT NOT NULL ,
  `GOName` LONGTEXT NULL ,
  `GOType` VARCHAR(45) NULL ,
  `GOPaths` LONGTEXT NULL ,
  `LevelOne` LONGTEXT NULL ,
  `LevelTwo` LONGTEXT NULL ,
  PRIMARY KEY (`GOTermID`) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`proteininfo`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`proteininfo` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`proteininfo` (
  `UPAccession` VARCHAR(10) NOT NULL ,
  `ProteinName` VARCHAR(45) NULL ,
  `A` FLOAT NULL DEFAULT 0.0 ,
  `C` FLOAT NULL DEFAULT 0.0 ,
  `D` FLOAT NULL DEFAULT 0.0 ,
  `E` FLOAT NULL DEFAULT 0.0 ,
  `F` FLOAT NULL DEFAULT 0.0 ,
  `G` FLOAT NULL DEFAULT 0.0 ,
  `H` FLOAT NULL DEFAULT 0.0 ,
  `I` FLOAT NULL DEFAULT 0.0 ,
  `K` FLOAT NULL DEFAULT 0.0 ,
  `L` FLOAT NULL DEFAULT 0.0 ,
  `M` FLOAT NULL DEFAULT 0.0 ,
  `P` FLOAT NULL DEFAULT 0.0 ,
  `N` FLOAT NULL DEFAULT 0.0 ,
  `Q` FLOAT NULL DEFAULT 0.0 ,
  `R` FLOAT NULL DEFAULT 0.0 ,
  `S` FLOAT NULL DEFAULT 0.0 ,
  `T` FLOAT NULL DEFAULT 0.0 ,
  `V` FLOAT NULL DEFAULT 0.0 ,
  `W` FLOAT NULL DEFAULT 0.0 ,
  `Y` FLOAT NULL DEFAULT 0.0 ,
  `NegativelyCharged` FLOAT NULL DEFAULT 0.0 ,
  `PositivelyCharged` FLOAT NULL DEFAULT 0.0 ,
  `Basic` FLOAT NULL DEFAULT 0.0 ,
  `Charged` FLOAT NULL DEFAULT 0.0 ,
  `Polar` FLOAT NULL DEFAULT 0.0 ,
  `NonPolar` FLOAT NULL DEFAULT 0.0 ,
  `Aromatic` FLOAT NULL DEFAULT 0.0 ,
  `Aliphatic` FLOAT NULL DEFAULT 0.0 ,
  `Small` FLOAT NULL DEFAULT 0.0 ,
  `Tiny` FLOAT NULL DEFAULT 0.0 ,
  `PESTMotif` INT NULL DEFAULT 0 ,
  `LowComplexity` INT NULL DEFAULT 0 ,
  `Hydrophobicity` FLOAT NULL DEFAULT 0.0 ,
  `Isoelectric` FLOAT NULL DEFAULT 0.0 ,
  `ModeOfAction` VARCHAR(45) NULL DEFAULT 'NA' ,
  `ECNumber` LONGTEXT NULL ,
  `OGlycosylation` LONGTEXT NULL ,
  `NGlycosylation` LONGTEXT NULL ,
  `Phosphoserine` LONGTEXT NULL ,
  `Phosphothreonine` LONGTEXT NULL ,
  `Phosphotyrosine` LONGTEXT NULL ,
  `SubcellularLocation` LONGTEXT NULL ,
  `TopologicalDomain` LONGTEXT NULL ,
  `PredictedSubcellularLocation` LONGTEXT NULL ,
  `SignalPeptide` LONGTEXT NULL ,
  `TransmembraneHelices` LONGTEXT NULL ,
  `Turns` LONGTEXT NULL ,
  `AlphaHelices` LONGTEXT NULL ,
  `BetaStrands` LONGTEXT NULL ,
  `PredictedAlphaHelices` LONGTEXT NULL ,
  `PredictedBetaSheets` LONGTEXT NULL ,
  `Isoforms` LONGTEXT NULL ,
  `Target` VARCHAR(1) NULL DEFAULT 'N' ,
  `Sequence` LONGTEXT NULL ,
  PRIMARY KEY (`UPAccession`) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`entrezgene`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`entrezgene` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`entrezgene` (
  `GeneID` VARCHAR(45) NOT NULL ,
  `Cancer_Ovarian` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Melanoma` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Breast` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Lung` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Pancreatic` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Liver` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Colon` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Prostate` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Testicular` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Oesophageal` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Stomach` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Heart` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Oral_Cavity` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Leukaemia` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Lymphoma` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Intestinal` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_Brain` VARCHAR(1) NULL DEFAULT 'N' ,
  `Cancer_All` VARCHAR(1) NULL DEFAULT 'N' ,
  PRIMARY KEY (`GeneID`) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`ensemblgene`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`ensemblgene` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`ensemblgene` (
  `EnsemblGeneID` VARCHAR(45) NOT NULL ,
  `NumberTranscripts` INT NULL DEFAULT 0 ,
  `ProteinCodingTranscripts` INT NULL DEFAULT 0 ,
  `RetainedIntronTranscripts` INT NULL DEFAULT 0 ,
  `ProcessedTranscripts` INT NULL DEFAULT 0 ,
  `NonsenseMediatedDecayTranscripts` INT NULL DEFAULT 0 ,
  PRIMARY KEY (`EnsemblGeneID`) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`germvariants`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`germvariants` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`germvariants` (
  `EnsemblTranscriptID` VARCHAR(45) NOT NULL ,
  `VariantID` VARCHAR(45) NOT NULL ,
  `EnsemblGeneID` VARCHAR(45) NOT NULL ,
  `InitialAminoAcid` VARCHAR(45) NULL DEFAULT 'NA' ,
  `FinalAminoAcid` VARCHAR(45) NULL DEFAULT 'NA' ,
  `CDSStart` INT NULL DEFAULT -1 ,
  `CDSEnd` INT NULL DEFAULT -1 ,
  `3Untranslated` INT NULL DEFAULT 0 ,
  `5Untranslated` INT NULL DEFAULT 0 ,
  `CodingUnknown` INT NULL DEFAULT 0 ,
  `ComplexInDel` INT NULL DEFAULT 0 ,
  `Downstream` INT NULL DEFAULT 0 ,
  `EssentialSpliceSite` INT NULL DEFAULT 0 ,
  `FrameshiftCoding` INT NULL DEFAULT 0 ,
  `Intergenic` INT NULL DEFAULT 0 ,
  `Intronic` INT NULL DEFAULT 0 ,
  `NMDTranscript` INT NULL DEFAULT 0 ,
  `NonSynonymousCoding` INT NULL DEFAULT 0 ,
  `PartialCodon` INT NULL DEFAULT 0 ,
  `RegulatoryRegion` INT NULL DEFAULT 0 ,
  `SpliceSite` INT NULL DEFAULT 0 ,
  `StopGained` INT NULL DEFAULT 0 ,
  `StopLost` INT NULL DEFAULT 0 ,
  `SynonymousCoding` INT NULL DEFAULT 0 ,
  `TranscriptionFactorBindingMotif` INT NULL DEFAULT 0 ,
  `Upstream` INT NULL DEFAULT 0 ,
  `WithinMatureMIRNA` INT NULL DEFAULT 0 ,
  `WithinNonCodingGene` INT NULL DEFAULT 0 ,
  INDEX `EnsemblGeneIDGermVar` (`EnsemblGeneID` ASC) ,
  PRIMARY KEY (`EnsemblTranscriptID`, `VariantID`) ,
  CONSTRAINT `EnsemblGeneIDGermVar`
    FOREIGN KEY (`EnsemblGeneID` )
    REFERENCES `proteindatabase`.`ensemblgene` (`EnsemblGeneID` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`unigene`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`unigene` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`unigene` (
  `UniGeneID` INT NOT NULL ,
  `DS_Embryoid_Body` INT NULL DEFAULT 0 ,
  `DS_Blastocyst` INT NULL DEFAULT 0 ,
  `DS_Fetus` INT NULL DEFAULT 0 ,
  `DS_Neonate` INT NULL DEFAULT 0 ,
  `DS_Infant` INT NULL DEFAULT 0 ,
  `DS_Juvenile` INT NULL DEFAULT 0 ,
  `DS_Adult` INT NULL DEFAULT 0 ,
  `HS_Adrenal_Tumor` INT NULL DEFAULT 0 ,
  `HS_Bladder_Carcinoma` INT NULL DEFAULT 0 ,
  `HS_Breast_Mammary_Gland_Tumor` INT NULL DEFAULT 0 ,
  `HS_Cervical_Tumor` INT NULL DEFAULT 0 ,
  `HS_Chondrosarcoma` INT NULL DEFAULT 0 ,
  `HS_Colorectal_Tumor` INT NULL DEFAULT 0 ,
  `HS_Esophageal_Tumor` INT NULL DEFAULT 0 ,
  `HS_Gastrointestinal_Tumor` INT NULL DEFAULT 0 ,
  `HS_Germ_Cell_Tumor` INT NULL DEFAULT 0 ,
  `HS_Glioma` INT NULL DEFAULT 0 ,
  `HS_Head_And_Neck_Tumor` INT NULL DEFAULT 0 ,
  `HS_Kidney_Tumor` INT NULL DEFAULT 0 ,
  `HS_Leukemia_Tumor` INT NULL DEFAULT 0 ,
  `HS_Liver_Tumor` INT NULL DEFAULT 0 ,
  `HS_Lung_Tumor` INT NULL DEFAULT 0 ,
  `HS_Lymphoma` INT NULL DEFAULT 0 ,
  `HS_Non_neoplasia` INT NULL DEFAULT 0 ,
  `HS_Normal` INT NULL DEFAULT 0 ,
  `HS_Ovarian_Tumor` INT NULL DEFAULT 0 ,
  `HS_Pancreatic_Tumor` INT NULL DEFAULT 0 ,
  `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT NULL DEFAULT 0 ,
  `HS_Prostate_Cancer` INT NULL DEFAULT 0 ,
  `HS_Retinoblastoma` INT NULL DEFAULT 0 ,
  `HS_Skin_Tumor` INT NULL DEFAULT 0 ,
  `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT NULL DEFAULT 0 ,
  `HS_Uterine_Tumor` INT NULL DEFAULT 0 ,
  `BS_Adipose_Tissue` INT NULL DEFAULT 0 ,
  `BS_Adrenal_Gland` INT NULL DEFAULT 0 ,
  `BS_Ascites` INT NULL DEFAULT 0 ,
  `BS_Bladder` INT NULL DEFAULT 0 ,
  `BS_Blood` INT NULL DEFAULT 0 ,
  `BS_Bone` INT NULL DEFAULT 0 ,
  `BS_Bone_Marrow` INT NULL DEFAULT 0 ,
  `BS_Brain` INT NULL DEFAULT 0 ,
  `BS_Cervix` INT NULL DEFAULT 0 ,
  `BS_Connective_Tissue` INT NULL DEFAULT 0 ,
  `BS_Ear` INT NULL DEFAULT 0 ,
  `BS_Embryonic_Tissue` INT NULL DEFAULT 0 ,
  `BS_Esophagus` INT NULL DEFAULT 0 ,
  `BS_Eye` INT NULL DEFAULT 0 ,
  `BS_Heart` INT NULL DEFAULT 0 ,
  `BS_Intestine` INT NULL DEFAULT 0 ,
  `BS_Kidney` INT NULL DEFAULT 0 ,
  `BS_Larynx` INT NULL DEFAULT 0 ,
  `BS_Liver` INT NULL DEFAULT 0 ,
  `BS_Lung` INT NULL DEFAULT 0 ,
  `BS_Lymph` INT NULL DEFAULT 0 ,
  `BS_Lymph_Node` INT NULL DEFAULT 0 ,
  `BS_Mammary_Gland` INT NULL DEFAULT 0 ,
  `BS_Mouth` INT NULL DEFAULT 0 ,
  `BS_Muscle` INT NULL DEFAULT 0 ,
  `BS_Nerve` INT NULL DEFAULT 0 ,
  `BS_Ovary` INT NULL DEFAULT 0 ,
  `BS_Pancreas` INT NULL DEFAULT 0 ,
  `BS_Parathyroid` INT NULL DEFAULT 0 ,
  `BS_Pharynx` INT NULL DEFAULT 0 ,
  `BS_Pituitary_Gland` INT NULL DEFAULT 0 ,
  `BS_Placenta` INT NULL DEFAULT 0 ,
  `BS_Prostate` INT NULL DEFAULT 0 ,
  `BS_Salivary_Gland` INT NULL DEFAULT 0 ,
  `BS_Skin` INT NULL DEFAULT 0 ,
  `BS_Spleen` INT NULL DEFAULT 0 ,
  `BS_Stomach` INT NULL DEFAULT 0 ,
  `BS_Testis` INT NULL DEFAULT 0 ,
  `BS_Thymus` INT NULL DEFAULT 0 ,
  `BS_Thyroid` INT NULL DEFAULT 0 ,
  `BS_Tonsil` INT NULL DEFAULT 0 ,
  `BS_Trachea` INT NULL DEFAULT 0 ,
  `BS_Umbilical_Cord` INT NULL DEFAULT 0 ,
  `BS_Uterus` INT NULL DEFAULT 0 ,
  `BS_Vascular` INT NULL DEFAULT 0 ,
  PRIMARY KEY (`UniGeneID`) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`uniprot2go`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`uniprot2go` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`uniprot2go` (
  `UPAccession` VARCHAR(45) NOT NULL ,
  `GOTermID` INT NOT NULL ,
  PRIMARY KEY (`UPAccession`, `GOTermID`) ,
  INDEX `UPAccessionGO` (`UPAccession` ASC) ,
  INDEX `GOTermIDUP` (`GOTermID` ASC) ,
  CONSTRAINT `UPAccessionGO`
    FOREIGN KEY (`UPAccession` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `GOTermIDUP`
    FOREIGN KEY (`GOTermID` )
    REFERENCES `proteindatabase`.`goinfo` (`GOTermID` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`uniprot2entrez`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`uniprot2entrez` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`uniprot2entrez` (
  `UPAccession` VARCHAR(45) NOT NULL ,
  `GeneID` VARCHAR(45) NOT NULL ,
  PRIMARY KEY (`UPAccession`, `GeneID`) ,
  INDEX `UPAccessionEntrez` (`UPAccession` ASC) ,
  INDEX `GeneIDUP` (`GeneID` ASC) ,
  CONSTRAINT `UPAccessionEntrez`
    FOREIGN KEY (`UPAccession` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `GeneIDUP`
    FOREIGN KEY (`GeneID` )
    REFERENCES `proteindatabase`.`entrezgene` (`GeneID` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`uniprot2unigene`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`uniprot2unigene` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`uniprot2unigene` (
  `UPAccession` VARCHAR(45) NOT NULL ,
  `UniGeneID` INT NOT NULL ,
  PRIMARY KEY (`UPAccession`, `UniGeneID`) ,
  INDEX `UPAccessionUG` (`UPAccession` ASC) ,
  INDEX `UniGeneIDUP` (`UniGeneID` ASC) ,
  CONSTRAINT `UPAccessionUG`
    FOREIGN KEY (`UPAccession` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `UniGeneIDUP`
    FOREIGN KEY (`UniGeneID` )
    REFERENCES `proteindatabase`.`unigene` (`UniGeneID` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`uniprot2ensembl`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`uniprot2ensembl` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`uniprot2ensembl` (
  `UPAccession` VARCHAR(45) NOT NULL ,
  `EnsemblGeneID` VARCHAR(45) NOT NULL ,
  `EnsemblTranscriptID` VARCHAR(45) NOT NULL ,
  `EnsemblProteinID` VARCHAR(45) NOT NULL ,
  PRIMARY KEY (`UPAccession`, `EnsemblGeneID`, `EnsemblTranscriptID`, `EnsemblProteinID`) ,
  INDEX `UPAccessionEnsembl` (`UPAccession` ASC) ,
  INDEX `EnsemblGeneIDUP` (`EnsemblGeneID` ASC) ,
  CONSTRAINT `UPAccessionEnsembl`
    FOREIGN KEY (`UPAccession` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `EnsemblGeneIDUP`
    FOREIGN KEY (`EnsemblGeneID` )
    REFERENCES `proteindatabase`.`ensemblgene` (`EnsemblGeneID` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`unigenetotals`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`unigenetotals` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`unigenetotals` (
  `StageStateSite` VARCHAR(100) NOT NULL ,
  `Total` INT NULL DEFAULT 0 ,
  PRIMARY KEY (`StageStateSite`) )
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`blastresults`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`blastresults` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`blastresults` (
  `ProteinA` VARCHAR(10) NOT NULL ,
  `ProteinB` VARCHAR(10) NOT NULL ,
  `Similarity` FLOAT NULL ,
  `Length` INT NULL ,
  `EValue` FLOAT NULL ,
  PRIMARY KEY (`ProteinA`, `ProteinB`) ,
  INDEX `QueryUP` (`ProteinA` ASC) ,
  INDEX `HitUP` (`ProteinB` ASC) ,
  CONSTRAINT `QueryUP`
    FOREIGN KEY (`ProteinA` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `HitUP`
    FOREIGN KEY (`ProteinB` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`nonredundant`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`nonredundant` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`nonredundant` (
  `UPAccession` VARCHAR(10) NOT NULL ,
  `AllTargetPositive` VARCHAR(1) NULL DEFAULT 'N' ,
  `AllTargetNegative` VARCHAR(1) NULL DEFAULT 'N' ,
  `GPCRTargetPositive` VARCHAR(1) NULL DEFAULT 'N' ,
  `GPCRTargetNegative` VARCHAR(1) NULL DEFAULT 'N' ,
  `IonChannelTargetPositive` VARCHAR(1) NULL DEFAULT 'N' ,
  `IonChannelTargetNegative` VARCHAR(1) NULL DEFAULT 'N' ,
  `KinaseTargetPositive` VARCHAR(1) NULL DEFAULT 'N' ,
  `KinaseTargetNegative` VARCHAR(1) NULL DEFAULT 'N' ,
  `ProteaseTargetPositive` VARCHAR(1) NULL DEFAULT 'N' ,
  `ProteaseTargetNegative` VARCHAR(1) NULL DEFAULT 'N' ,
  `CancerTargetPositive` VARCHAR(1) NULL DEFAULT 'N' ,
  `CancerTargetNegative` VARCHAR(1) NULL DEFAULT 'N' ,
  `CancerTypePositive` VARCHAR(1) NULL DEFAULT 'N' ,
  `CancerTypeNegative` VARCHAR(1) NULL DEFAULT 'N' ,
  `CancerProteinPositive` VARCHAR(1) NULL DEFAULT 'N' ,
  `CancerProteinNegative` VARCHAR(1) NULL DEFAULT 'N' ,
  PRIMARY KEY (`UPAccession`) ,
  INDEX `UPAccessionNR` (`UPAccession` ASC) ,
  CONSTRAINT `UPAccessionNR`
    FOREIGN KEY (`UPAccession` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`cosmicvariants`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`cosmicvariants` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`cosmicvariants` (
  `EnsemblTranscriptID` VARCHAR(45) NOT NULL ,
  `VariantID` VARCHAR(45) NOT NULL ,
  `EnsemblGeneID` VARCHAR(45) NOT NULL ,
  `InitialAminoAcid` VARCHAR(45) NULL DEFAULT 'NA' ,
  `FinalAminoAcid` VARCHAR(45) NULL DEFAULT 'NA' ,
  `CDSStart` INT NULL DEFAULT -1 ,
  `CDSEnd` INT NULL DEFAULT -1 ,
  `3Untranslated` INT NULL DEFAULT 0 ,
  `5Untranslated` INT NULL DEFAULT 0 ,
  `CodingUnknown` INT NULL DEFAULT 0 ,
  `ComplexInDel` INT NULL DEFAULT 0 ,
  `Downstream` INT NULL DEFAULT 0 ,
  `EssentialSpliceSite` INT NULL DEFAULT 0 ,
  `FrameshiftCoding` INT NULL DEFAULT 0 ,
  `Intergenic` INT NULL DEFAULT 0 ,
  `Intronic` INT NULL DEFAULT 0 ,
  `NMDTranscript` INT NULL DEFAULT 0 ,
  `NonSynonymousCoding` INT NULL DEFAULT 0 ,
  `PartialCodon` INT NULL DEFAULT 0 ,
  `RegulatoryRegion` INT NULL DEFAULT 0 ,
  `SpliceSite` INT NULL DEFAULT 0 ,
  `StopGained` INT NULL DEFAULT 0 ,
  `StopLost` INT NULL DEFAULT 0 ,
  `SynonymousCoding` INT NULL DEFAULT 0 ,
  `TranscriptionFactorBindingMotif` INT NULL DEFAULT 0 ,
  `Upstream` INT NULL DEFAULT 0 ,
  `WithinMatureMIRNA` INT NULL DEFAULT 0 ,
  `WithinNonCodingGene` INT NULL DEFAULT 0 ,
  INDEX `EnsemblGeneIDCosmic` (`EnsemblGeneID` ASC) ,
  PRIMARY KEY (`EnsemblTranscriptID`, `VariantID`) ,
  CONSTRAINT `EnsemblGeneIDCosmic`
    FOREIGN KEY (`EnsemblGeneID` )
    REFERENCES `proteindatabase`.`ensemblgene` (`EnsemblGeneID` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`homologs`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`homologs` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`homologs` (
  `HumanGene` VARCHAR(45) NOT NULL ,
  `HomologGene` VARCHAR(45) NOT NULL ,
  `HomologSpecies` VARCHAR(255) NULL ,
  `HomologyType` VARCHAR(45) NULL ,
  `Ancestor` VARCHAR(255) NULL ,
  `dN` FLOAT NULL DEFAULT -1 ,
  `dS` FLOAT NULL DEFAULT -1 ,
  `PeptideAligned` FLOAT NULL DEFAULT 0.0 ,
  `Identity` FLOAT NULL DEFAULT 0.0 ,
  `Positivity` FLOAT NULL DEFAULT 0.0 ,
  PRIMARY KEY (`HumanGene`, `HomologGene`) ,
  INDEX `HomologEnsemblHumanGene` (`HumanGene` ASC) ,
  CONSTRAINT `HomologEnsemblHumanGene`
    FOREIGN KEY (`HumanGene` )
    REFERENCES `proteindatabase`.`ensemblgene` (`EnsemblGeneID` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`ppi`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`ppi` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`ppi` (
  `PPIProteinOne` VARCHAR(10) NOT NULL ,
  `PPIProteinTwo` VARCHAR(10) NOT NULL ,
  `IsoformID` VARCHAR(45) NOT NULL DEFAULT 'No Isoform' ,
  `OrganismsDiffer` VARCHAR(10) NULL DEFAULT false ,
  `NumExperiments` INT NULL DEFAULT 0 ,
  PRIMARY KEY (`PPIProteinOne`, `PPIProteinTwo`, `IsoformID`) ,
  INDEX `ProteinOneUP` (`PPIProteinOne` ASC) ,
  CONSTRAINT `ProteinOneUP`
    FOREIGN KEY (`PPIProteinOne` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB
COMMENT = 'PPIProteinTwo: not foreign key as may not be human.';


-- -----------------------------------------------------
-- Table `proteindatabase`.`drugs`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`drugs` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`drugs` (
  `UPAccession` VARCHAR(10) NOT NULL ,
  `DrugID` VARCHAR(45) NOT NULL ,
  `DrugName` VARCHAR(255) NOT NULL ,
  `KiValue` FLOAT NULL DEFAULT -1 ,
  `KdValue` FLOAT NULL DEFAULT -1 ,
  PRIMARY KEY (`UPAccession`, `DrugID`, `DrugName`) ,
  INDEX `UPAccDrug` (`UPAccession` ASC) ,
  CONSTRAINT `UPAccDrug`
    FOREIGN KEY (`UPAccession` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB
COMMENT = 'The Kd and Ki values are in nM.';


-- -----------------------------------------------------
-- Table `proteindatabase`.`pathways`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`pathways` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`pathways` (
  `UPAccession` VARCHAR(10) NOT NULL ,
  `NumberOfPathways` INT NULL ,
  PRIMARY KEY (`UPAccession`) ,
  INDEX `UPAccPathway` (`UPAccession` ASC) ,
  CONSTRAINT `UPAccPathway`
    FOREIGN KEY (`UPAccession` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `proteindatabase`.`stability`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `proteindatabase`.`stability` ;

CREATE  TABLE IF NOT EXISTS `proteindatabase`.`stability` (
  `UPAccession` VARCHAR(10) NOT NULL ,
  `HalfLife` FLOAT NULL DEFAULT 0.0 ,
  `InstabilityIndex` FLOAT NULL DEFAULT 0.0 ,
  PRIMARY KEY (`UPAccession`) ,
  INDEX `UPAccStability` (`UPAccession` ASC) ,
  CONSTRAINT `UPAccStability`
    FOREIGN KEY (`UPAccession` )
    REFERENCES `proteindatabase`.`proteininfo` (`UPAccession` )
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`all_all_targ_r_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`all_all_targ_r_p` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`all_all_targ_r_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`all_all_targ_r_n` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_gpcr_targ_r_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_gpcr_targ_r_p` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_gpcr_targ_r_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_gpcr_targ_r_n` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_kinase_targ_r_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_kinase_targ_r_p` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_protease_targ_r_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_protease_targ_r_p` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_ionchannel_targ_r_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_ionchannel_targ_r_p` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_kinase_targ_r_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_kinase_targ_r_n` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_protease_targ_r_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_protease_targ_r_n` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_ionchannel_targ_r_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_ionchannel_targ_r_n` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_targ_r_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_targ_r_n` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_targ_r_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_targ_r_p` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_transcripts`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_transcripts` (`id` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_ppi`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_ppi` (`UPAccession` INT, `BinaryPPI` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_paralogs`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_paralogs` (`UPAccession` INT, `Paralogs` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_germvariants`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_germvariants` (`UPAccession` INT, `VariantID` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_expression`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_expression` (`UPAccession` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `Total` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`all_all_targ_nr_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`all_all_targ_nr_p` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_genetrans`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_genetrans` (`UPAccession` INT, `EnsemblGeneID` INT, `EnsemblTranscriptID` INT, `ProteinCodingTranscripts` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_variantsnumber`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_variantsnumber` (`UPAccession` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Total` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`all_all_targ_nr_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`all_all_targ_nr_n` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_targ_nr_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_targ_nr_p` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_targ_nr_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_targ_nr_n` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_gpcr_targ_nr_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_gpcr_targ_nr_p` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_gpcr_targ_nr_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_gpcr_targ_nr_n` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_ionchannel_targ_nr_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_ionchannel_targ_nr_p` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_ionchannel_targ_nr_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_ionchannel_targ_nr_n` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_kinase_targ_nr_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_kinase_targ_nr_p` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_kinase_targ_nr_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_kinase_targ_nr_n` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_protease_targ_nr_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_protease_targ_nr_p` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`type_protease_targ_nr_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`type_protease_targ_nr_n` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_ovarian`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_ovarian` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_melanoma`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_melanoma` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_breast`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_breast` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_lung`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_lung` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_pancreatic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_pancreatic` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_liver`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_liver` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_colon`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_colon` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_prostate`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_prostate` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_testicular`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_testicular` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_oesophageal`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_oesophageal` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_stomach`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_stomach` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_heart`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_heart` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_oral_cavity`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_oral_cavity` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_leukaemia`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_leukaemia` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_lymphoma`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_lymphoma` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_intestinal`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_intestinal` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`upacc_cancer_brain`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`upacc_cancer_brain` (`UPAccession` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_prot_r_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_prot_r_p` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_prot_r_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_prot_r_n` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_prot_nr_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_prot_nr_p` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_prot_nr_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_prot_nr_n` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_type_r_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_type_r_p` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_type_r_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_type_r_n` (`UPAccession` INT, `Sequence` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_type_nr_n`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_type_nr_n` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- Placeholder table for view `proteindatabase`.`ill_cancer_type_nr_p`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `proteindatabase`.`ill_cancer_type_nr_p` (`UPAccession` INT, `A` INT, `C` INT, `D` INT, `E` INT, `F` INT, `G` INT, `H` INT, `I` INT, `K` INT, `L` INT, `M` INT, `P` INT, `N` INT, `Q` INT, `R` INT, `S` INT, `T` INT, `V` INT, `W` INT, `Y` INT, `NegativelyCharged` INT, `PositivelyCharged` INT, `Basic` INT, `Charged` INT, `Polar` INT, `NonPolar` INT, `Aromatic` INT, `Aliphatic` INT, `Small` INT, `Tiny` INT, `PESTMotif` INT, `LowComplexity` INT, `Hydrophobicity` INT, `Isoelectric` INT, `ECNumber` INT, `OGlycosylation` INT, `NGlycosylation` INT, `Phosphoserine` INT, `Phosphothreonine` INT, `Phosphotyrosine` INT, `SubcellularLocation` INT, `TopologicalDomain` INT, `PredictedSubcellularLocation` INT, `SignalPeptide` INT, `TransmembraneHelices` INT, `AlphaHelices` INT, `BetaStrands` INT, `PredictedAlphaHelices` INT, `PredictedBetaSheets` INT, `Sequence` INT, `3Untranslated` INT, `5Untranslated` INT, `NonSynonymousCoding` INT, `SynonymousCoding` INT, `Paralogs` INT, `BinaryPPI` INT, `AlternativeTranscripts` INT, `DS_Embryoid_Body` INT, `DS_Blastocyst` INT, `DS_Fetus` INT, `DS_Neonate` INT, `DS_Infant` INT, `DS_Juvenile` INT, `DS_Adult` INT, `HS_Adrenal_Tumor` INT, `HS_Bladder_Carcinoma` INT, `HS_Breast_Mammary_Gland_Tumor` INT, `HS_Cervical_Tumor` INT, `HS_Chondrosarcoma` INT, `HS_Colorectal_Tumor` INT, `HS_Esophageal_Tumor` INT, `HS_Gastrointestinal_Tumor` INT, `HS_Germ_Cell_Tumor` INT, `HS_Glioma` INT, `HS_Head_And_Neck_Tumor` INT, `HS_Kidney_Tumor` INT, `HS_Leukemia_Tumor` INT, `HS_Liver_Tumor` INT, `HS_Lung_Tumor` INT, `HS_Lymphoma` INT, `HS_Non_neoplasia` INT, `HS_Normal` INT, `HS_Ovarian_Tumor` INT, `HS_Pancreatic_Tumor` INT, `HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS` INT, `HS_Prostate_Cancer` INT, `HS_Retinoblastoma` INT, `HS_Skin_Tumor` INT, `HS_Soft_Tissue_Muscle_Tissue_Tumor` INT, `HS_Uterine_Tumor` INT, `BS_Adipose_Tissue` INT, `BS_Adrenal_Gland` INT, `BS_Ascites` INT, `BS_Bladder` INT, `BS_Blood` INT, `BS_Bone` INT, `BS_Bone_Marrow` INT, `BS_Brain` INT, `BS_Cervix` INT, `BS_Connective_Tissue` INT, `BS_Ear` INT, `BS_Embryonic_Tissue` INT, `BS_Esophagus` INT, `BS_Eye` INT, `BS_Heart` INT, `BS_Intestine` INT, `BS_Kidney` INT, `BS_Larynx` INT, `BS_Liver` INT, `BS_Lung` INT, `BS_Lymph` INT, `BS_Lymph_Node` INT, `BS_Mammary_Gland` INT, `BS_Mouth` INT, `BS_Muscle` INT, `BS_Nerve` INT, `BS_Ovary` INT, `BS_Pancreas` INT, `BS_Parathyroid` INT, `BS_Pharynx` INT, `BS_Pituitary_Gland` INT, `BS_Placenta` INT, `BS_Prostate` INT, `BS_Salivary_Gland` INT, `BS_Skin` INT, `BS_Spleen` INT, `BS_Stomach` INT, `BS_Testis` INT, `BS_Thymus` INT, `BS_Thyroid` INT, `BS_Tonsil` INT, `BS_Trachea` INT, `BS_Umbilical_Cord` INT, `BS_Uterus` INT, `BS_Vascular` INT, `HalfLife` INT, `InstabilityIndex` INT);

-- -----------------------------------------------------
-- View `proteindatabase`.`all_all_targ_r_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`all_all_targ_r_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`all_all_targ_r_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`all_all_targ_r_p` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`all_all_targ_r_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`all_all_targ_r_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`all_all_targ_r_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`all_all_targ_r_n` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "N";

-- -----------------------------------------------------
-- View `proteindatabase`.`type_gpcr_targ_r_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_gpcr_targ_r_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_gpcr_targ_r_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_gpcr_targ_r_p` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "Y" AND `ModeOfAction` = "G-protein coupled receptor";

-- -----------------------------------------------------
-- View `proteindatabase`.`type_gpcr_targ_r_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_gpcr_targ_r_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_gpcr_targ_r_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_gpcr_targ_r_n` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "N" AND `ModeOfAction` = "G-protein coupled receptor";

-- -----------------------------------------------------
-- View `proteindatabase`.`type_kinase_targ_r_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_kinase_targ_r_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_kinase_targ_r_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_kinase_targ_r_p` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "Y" AND `ModeOfAction` = "Kinase";

-- -----------------------------------------------------
-- View `proteindatabase`.`type_protease_targ_r_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_protease_targ_r_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_protease_targ_r_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_protease_targ_r_p` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "Y" AND `ModeOfAction` = "Protease";

-- -----------------------------------------------------
-- View `proteindatabase`.`type_ionchannel_targ_r_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_ionchannel_targ_r_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_ionchannel_targ_r_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_ionchannel_targ_r_p` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "Y" AND `ModeOfAction` = "Ion Channel";

-- -----------------------------------------------------
-- View `proteindatabase`.`type_kinase_targ_r_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_kinase_targ_r_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_kinase_targ_r_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_kinase_targ_r_n` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "N" AND `ModeOfAction` = "Kinase";

-- -----------------------------------------------------
-- View `proteindatabase`.`type_protease_targ_r_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_protease_targ_r_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_protease_targ_r_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_protease_targ_r_n` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "N" AND `ModeOfAction` = "Protease";

-- -----------------------------------------------------
-- View `proteindatabase`.`type_ionchannel_targ_r_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_ionchannel_targ_r_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_ionchannel_targ_r_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_ionchannel_targ_r_n` AS
SELECT `UPAccession`, `Sequence` FROM `proteindatabase`.`proteininfo`
    WHERE `Target` = "N" AND `ModeOfAction` = "Ion Channel";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_All` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_targ_r_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_targ_r_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_targ_r_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_targ_r_n` AS
    SELECT
        `proteininfo`.`UPAccession`, `Sequence`
    FROM
        `proteindatabase`.`proteininfo`
        INNER JOIN
        `proteindatabase`.`upacc_cancer`
        ON `proteininfo`.`UPAccession` = `upacc_cancer`.`UPAccession`
    WHERE
        `proteininfo`.`Target` = "N";

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_targ_r_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_targ_r_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_targ_r_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_targ_r_p` AS
    SELECT
        `proteininfo`.`UPAccession`, `Sequence`
    FROM
        `proteindatabase`.`proteininfo`
        INNER JOIN
        `proteindatabase`.`upacc_cancer`
        ON `proteininfo`.`UPAccession` = `upacc_cancer`.`UPAccession`
    WHERE
        `proteininfo`.`Target` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_transcripts`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_transcripts` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_transcripts`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_transcripts` AS
    SELECT
        maxi.*
    FROM
            upacc_genetrans AS maxi
        LEFT OUTER JOIN
            upacc_genetrans AS checker
        ON
            checker.UPAccession = maxi.UPAccession AND
            checker.ProteinCodingTranscripts > maxi.ProteinCodingTranscripts
    WHERE
        checker.ProteinCodingTranscripts IS NULL;

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_ppi`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_ppi` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_ppi`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_ppi` AS
    SELECT
        prot.UPAccession,
        IFNULL((SELECT COUNT(*) FROM ppi WHERE prot.UPAccession = ppi.PPIProteinOne AND ppi.OrganismsDiffer = 'false'), 0) As BinaryPPI
    FROM
        proteininfo AS prot
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_paralogs`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_paralogs` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_paralogs`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_paralogs` AS
    SELECT
        trans.UPAccession,
        IFNULL((SELECT COUNT(DISTINCT hom.HomologGene) FROM homologs AS hom WHERE trans.EnsemblGeneID = hom.HumanGene AND hom.HomologyType = "within_species_paralog"), 0) AS Paralogs
    FROM
        upacc_transcripts AS trans
    GROUP BY
        trans.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_germvariants`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_germvariants` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_germvariants`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_germvariants` AS
    SELECT
        DISTINCT
        trans.UPAccession,
        germ.VariantID,
        germ.3Untranslated,
        germ.5Untranslated,
        germ.NonSynonymousCoding,
        germ.SynonymousCoding
    FROM
        upacc_transcripts AS trans,
        germvariants as germ
    WHERE
        trans.EnsemblTranscriptID = germ.EnsemblTranscriptID;

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_expression`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_expression` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_expression`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_expression` AS
    SELECT
        prot.UPAccession,
        IFNULL((SELECT SUM(unigene.DS_Embryoid_Body) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS DS_Embryoid_Body,
        IFNULL((SELECT SUM(unigene.DS_Blastocyst) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS DS_Blastocyst,
        IFNULL((SELECT SUM(unigene.DS_Fetus) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS DS_Fetus,
        IFNULL((SELECT SUM(unigene.DS_Neonate) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS DS_Neonate,
        IFNULL((SELECT SUM(unigene.DS_Infant) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS DS_Infant,
        IFNULL((SELECT SUM(unigene.DS_Juvenile) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS DS_Juvenile,
        IFNULL((SELECT SUM(unigene.DS_Adult) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS DS_Adult,
        IFNULL((SELECT SUM(unigene.HS_Adrenal_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Adrenal_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Bladder_Carcinoma) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Bladder_Carcinoma,
        IFNULL((SELECT SUM(unigene.HS_Breast_Mammary_Gland_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Cervical_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Cervical_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Chondrosarcoma) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Chondrosarcoma,
        IFNULL((SELECT SUM(unigene.HS_Colorectal_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Colorectal_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Esophageal_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Esophageal_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Gastrointestinal_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Germ_Cell_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Germ_Cell_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Glioma) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Glioma,
        IFNULL((SELECT SUM(unigene.HS_Head_And_Neck_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Kidney_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Kidney_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Leukemia_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Leukemia_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Liver_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Liver_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Lung_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Lung_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Lymphoma) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Lymphoma,
        IFNULL((SELECT SUM(unigene.HS_Non_neoplasia) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Non_neoplasia,
        IFNULL((SELECT SUM(unigene.HS_Normal) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Normal,
        IFNULL((SELECT SUM(unigene.HS_Ovarian_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Ovarian_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Pancreatic_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Pancreatic_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL((SELECT SUM(unigene.HS_Prostate_Cancer) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Prostate_Cancer,
        IFNULL((SELECT SUM(unigene.HS_Retinoblastoma) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Retinoblastoma,
        IFNULL((SELECT SUM(unigene.HS_Skin_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Skin_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Soft_Tissue_Muscle_Tissue_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL((SELECT SUM(unigene.HS_Uterine_Tumor) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS HS_Uterine_Tumor,
        IFNULL((SELECT SUM(unigene.BS_Adipose_Tissue) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Adipose_Tissue,
        IFNULL((SELECT SUM(unigene.BS_Adrenal_Gland) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Adrenal_Gland,
        IFNULL((SELECT SUM(unigene.BS_Ascites) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Ascites,
        IFNULL((SELECT SUM(unigene.BS_Bladder) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Bladder,
        IFNULL((SELECT SUM(unigene.BS_Blood) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Blood,
        IFNULL((SELECT SUM(unigene.BS_Bone) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Bone,
        IFNULL((SELECT SUM(unigene.BS_Bone_Marrow) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Bone_Marrow,
        IFNULL((SELECT SUM(unigene.BS_Brain) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Brain,
        IFNULL((SELECT SUM(unigene.BS_Cervix) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Cervix,
        IFNULL((SELECT SUM(unigene.BS_Connective_Tissue) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Connective_Tissue,
        IFNULL((SELECT SUM(unigene.BS_Ear) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Ear,
        IFNULL((SELECT SUM(unigene.BS_Embryonic_Tissue) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Embryonic_Tissue,
        IFNULL((SELECT SUM(unigene.BS_Esophagus) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Esophagus,
        IFNULL((SELECT SUM(unigene.BS_Eye) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Eye,
        IFNULL((SELECT SUM(unigene.BS_Heart) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Heart,
        IFNULL((SELECT SUM(unigene.BS_Intestine) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Intestine,
        IFNULL((SELECT SUM(unigene.BS_Kidney) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Kidney,
        IFNULL((SELECT SUM(unigene.BS_Larynx) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Larynx,
        IFNULL((SELECT SUM(unigene.BS_Liver) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Liver,
        IFNULL((SELECT SUM(unigene.BS_Lung) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Lung,
        IFNULL((SELECT SUM(unigene.BS_Lymph) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Lymph,
        IFNULL((SELECT SUM(unigene.BS_Lymph_Node) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Lymph_Node,
        IFNULL((SELECT SUM(unigene.BS_Mammary_Gland) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Mammary_Gland,
        IFNULL((SELECT SUM(unigene.BS_Mouth) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Mouth,
        IFNULL((SELECT SUM(unigene.BS_Muscle) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Muscle,
        IFNULL((SELECT SUM(unigene.BS_Nerve) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Nerve,
        IFNULL((SELECT SUM(unigene.BS_Ovary) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Ovary,
        IFNULL((SELECT SUM(unigene.BS_Pancreas) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Pancreas,
        IFNULL((SELECT SUM(unigene.BS_Parathyroid) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Parathyroid,
        IFNULL((SELECT SUM(unigene.BS_Pharynx) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Pharynx,
        IFNULL((SELECT SUM(unigene.BS_Pituitary_Gland) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Pituitary_Gland,
        IFNULL((SELECT SUM(unigene.BS_Placenta) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Placenta,
        IFNULL((SELECT SUM(unigene.BS_Prostate) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Prostate,
        IFNULL((SELECT SUM(unigene.BS_Salivary_Gland) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Salivary_Gland,
        IFNULL((SELECT SUM(unigene.BS_Skin) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Skin,
        IFNULL((SELECT SUM(unigene.BS_Spleen) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Spleen,
        IFNULL((SELECT SUM(unigene.BS_Stomach) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Stomach,
        IFNULL((SELECT SUM(unigene.BS_Testis) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Testis,
        IFNULL((SELECT SUM(unigene.BS_Thymus) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Thymus,
        IFNULL((SELECT SUM(unigene.BS_Thyroid) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Thyroid,
        IFNULL((SELECT SUM(unigene.BS_Tonsil) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Tonsil,
        IFNULL((SELECT SUM(unigene.BS_Trachea) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Trachea,
        IFNULL((SELECT SUM(unigene.BS_Umbilical_Cord) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Umbilical_Cord,
        IFNULL((SELECT SUM(unigene.BS_Uterus) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Uterus,
        IFNULL((SELECT SUM(unigene.BS_Vascular) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS BS_Vascular,
        IFNULL((SELECT SUM(unigene.DS_Embryoid_Body) + SUM(unigene.DS_Blastocyst) + SUM(unigene.DS_Fetus) + SUM(unigene.DS_Neonate) + SUM(unigene.DS_Infant) + SUM(unigene.DS_Juvenile) + SUM(unigene.DS_Adult) + SUM(unigene.HS_Adrenal_Tumor) + SUM(unigene.HS_Bladder_Carcinoma) + SUM(unigene.HS_Breast_Mammary_Gland_Tumor) + SUM(unigene.HS_Cervical_Tumor) + SUM(unigene.HS_Chondrosarcoma) + SUM(unigene.HS_Colorectal_Tumor) + SUM(unigene.HS_Esophageal_Tumor) + SUM(unigene.HS_Gastrointestinal_Tumor) + SUM(unigene.HS_Germ_Cell_Tumor) + SUM(unigene.HS_Glioma) + SUM(unigene.HS_Head_And_Neck_Tumor) + SUM(unigene.HS_Kidney_Tumor) + SUM(unigene.HS_Leukemia_Tumor) + SUM(unigene.HS_Liver_Tumor) + SUM(unigene.HS_Lung_Tumor) + SUM(unigene.HS_Lymphoma) + SUM(unigene.HS_Non_neoplasia) + SUM(unigene.HS_Normal) + SUM(unigene.HS_Ovarian_Tumor) + SUM(unigene.HS_Pancreatic_Tumor) + SUM(unigene.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS) + SUM(unigene.HS_Prostate_Cancer) + SUM(unigene.HS_Retinoblastoma) + SUM(unigene.HS_Skin_Tumor) + SUM(unigene.HS_Soft_Tissue_Muscle_Tissue_Tumor) + SUM(unigene.HS_Uterine_Tumor) + SUM(unigene.BS_Adipose_Tissue) + SUM(unigene.BS_Adrenal_Gland) + SUM(unigene.BS_Ascites) + SUM(unigene.BS_Bladder) + SUM(unigene.BS_Blood) + SUM(unigene.BS_Bone) + SUM(unigene.BS_Bone_Marrow) + SUM(unigene.BS_Brain) + SUM(unigene.BS_Cervix) + SUM(unigene.BS_Connective_Tissue) + SUM(unigene.BS_Ear) + SUM(unigene.BS_Embryonic_Tissue) + SUM(unigene.BS_Esophagus) + SUM(unigene.BS_Eye) + SUM(unigene.BS_Heart) + SUM(unigene.BS_Intestine) + SUM(unigene.BS_Kidney) + SUM(unigene.BS_Larynx) + SUM(unigene.BS_Liver) + SUM(unigene.BS_Lung) + SUM(unigene.BS_Lymph) + SUM(unigene.BS_Lymph_Node) + SUM(unigene.BS_Mammary_Gland) + SUM(unigene.BS_Mouth) + SUM(unigene.BS_Muscle) + SUM(unigene.BS_Nerve) + SUM(unigene.BS_Ovary) + SUM(unigene.BS_Pancreas) + SUM(unigene.BS_Parathyroid) + SUM(unigene.BS_Pharynx) + SUM(unigene.BS_Pituitary_Gland) + SUM(unigene.BS_Placenta) + SUM(unigene.BS_Prostate) + SUM(unigene.BS_Salivary_Gland) + SUM(unigene.BS_Skin) + SUM(unigene.BS_Spleen) + SUM(unigene.BS_Stomach) + SUM(unigene.BS_Testis) + SUM(unigene.BS_Thymus) + SUM(unigene.BS_Thyroid) + SUM(unigene.BS_Tonsil) + SUM(unigene.BS_Trachea) + SUM(unigene.BS_Umbilical_Cord) + SUM(unigene.BS_Uterus) + SUM(unigene.BS_Vascular) FROM uniprot2unigene AS u2u, unigene WHERE prot.UPAccession = u2u.UPAccession AND u2u.UniGeneID = unigene.UniGeneID), 0) AS Total
    FROM
        proteininfo AS prot
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`all_all_targ_nr_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`all_all_targ_nr_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`all_all_targ_nr_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`all_all_targ_nr_p` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.AllTargetPositive = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_genetrans`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_genetrans` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_genetrans`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_genetrans` AS
    SELECT
        UPAccession,
        ensemblgene.EnsemblGeneID,
        uniprot2ensembl.EnsemblTranscriptID,
        ensemblgene.ProteinCodingTranscripts
    FROM
        uniprot2ensembl
        JOIN ensemblgene
        ON uniprot2ensembl.EnsemblGeneID = ensemblgene.EnsemblGeneID;

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_variantsnumber`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_variantsnumber` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_variantsnumber`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_variantsnumber` AS
    SELECT
        trans.UPAccession,
        IFNULL((SELECT SUM(3Untranslated) FROM upacc_germvariants AS germ WHERE trans.UPAccession = germ.UPAccession), 0) AS 3Untranslated,
        IFNULL((SELECT SUM(5Untranslated) FROM upacc_germvariants AS germ WHERE trans.UPAccession = germ.UPAccession), 0) AS 5Untranslated,
        IFNULL((SELECT SUM(NonSynonymousCoding) FROM upacc_germvariants AS germ WHERE trans.UPAccession = germ.UPAccession), 0) AS NonSynonymousCoding,
        IFNULL((SELECT SUM(SynonymousCoding) FROM upacc_germvariants AS germ WHERE trans.UPAccession = germ.UPAccession), 0) AS SynonymousCoding,
        IFNULL((SELECT SUM(3Untranslated) + SUM(5Untranslated) + SUM(NonSynonymousCoding) + SUM(SynonymousCoding) FROM upacc_germvariants AS germ WHERE trans.UPAccession = germ.UPAccession), 0) AS Total
    FROM
        upacc_transcripts AS trans
    GROUP BY
        trans.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`all_all_targ_nr_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`all_all_targ_nr_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`all_all_targ_nr_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`all_all_targ_nr_n` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.AllTargetNegative = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_targ_nr_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_targ_nr_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_targ_nr_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_targ_nr_p` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.CancerTargetPositive = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_targ_nr_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_targ_nr_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_targ_nr_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_targ_nr_n` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.CancerTargetNegative = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`type_gpcr_targ_nr_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_gpcr_targ_nr_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_gpcr_targ_nr_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_gpcr_targ_nr_p` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.GPCRTargetPositive = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`type_gpcr_targ_nr_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_gpcr_targ_nr_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_gpcr_targ_nr_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_gpcr_targ_nr_n` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.GPCRTargetNegative = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`type_ionchannel_targ_nr_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_ionchannel_targ_nr_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_ionchannel_targ_nr_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_ionchannel_targ_nr_p` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.IonChannelTargetPositive = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`type_ionchannel_targ_nr_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_ionchannel_targ_nr_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_ionchannel_targ_nr_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_ionchannel_targ_nr_n` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.IonChannelTargetNegative = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`type_kinase_targ_nr_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_kinase_targ_nr_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_kinase_targ_nr_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_kinase_targ_nr_p` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.KinaseTargetPositive = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`type_kinase_targ_nr_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_kinase_targ_nr_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_kinase_targ_nr_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_kinase_targ_nr_n` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.KinaseTargetNegative = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`type_protease_targ_nr_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_protease_targ_nr_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_protease_targ_nr_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_protease_targ_nr_p` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.ProteaseTargetPositive = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`type_protease_targ_nr_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`type_protease_targ_nr_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`type_protease_targ_nr_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`type_protease_targ_nr_n` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.ProteaseTargetNegative = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_ovarian`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_ovarian` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_ovarian`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_ovarian` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Ovarian` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_melanoma`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_melanoma` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_melanoma`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_melanoma` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Melanoma` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_breast`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_breast` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_breast`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_breast` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Breast` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_lung`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_lung` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_lung`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_lung` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Lung` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_pancreatic`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_pancreatic` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_pancreatic`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_pancreatic` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Pancreatic` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_liver`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_liver` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_liver`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_liver` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Liver` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_colon`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_colon` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_colon`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_colon` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Colon` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_prostate`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_prostate` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_prostate`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_prostate` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Prostate` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_testicular`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_testicular` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_testicular`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_testicular` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Testicular` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_oesophageal`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_oesophageal` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_oesophageal`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_oesophageal` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Oesophageal` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_stomach`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_stomach` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_stomach`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_stomach` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Stomach` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_heart`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_heart` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_heart`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_heart` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Heart` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_oral_cavity`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_oral_cavity` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_oral_cavity`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_oral_cavity` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Oral_Cavity` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_leukaemia`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_leukaemia` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_leukaemia`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_leukaemia` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Leukaemia` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_lymphoma`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_lymphoma` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_lymphoma`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_lymphoma` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Lymphoma` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_intestinal`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_intestinal` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_intestinal`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_intestinal` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Intestinal` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`upacc_cancer_brain`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`upacc_cancer_brain` ;
DROP TABLE IF EXISTS `proteindatabase`.`upacc_cancer_brain`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`upacc_cancer_brain` AS
    SELECT
        DISTINCT `UPAccession`
    FROM
        `proteindatabase`.`uniprot2entrez` AS UP2Entrez
        INNER JOIN
        `proteindatabase`.`entrezgene` AS entrez
        ON UP2Entrez.`GeneID` = entrez.`GeneID`
    WHERE
        entrez.`Cancer_Brain` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_prot_r_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_prot_r_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_prot_r_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_prot_r_p` AS
    SELECT
        `proteininfo`.`UPAccession`, `Sequence`
    FROM
        `proteindatabase`.`proteininfo`
        INNER JOIN
        `proteindatabase`.`upacc_cancer`
        ON `proteininfo`.`UPAccession` = `upacc_cancer`.`UPAccession`;

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_prot_r_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_prot_r_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_prot_r_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_prot_r_n` AS
    SELECT
        prot.`UPAccession`, prot.`Sequence`
    FROM
        `proteindatabase`.`proteininfo` AS prot LEFT JOIN `proteindatabase`.`upacc_cancer` AS cancer
    ON
        prot.`UPAccession` = cancer.`UPAccession`
    WHERE
        cancer.`UPAccession` IS NULL;

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_prot_nr_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_prot_nr_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_prot_nr_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_prot_nr_p` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.CancerProteinPositive = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_prot_nr_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_prot_nr_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_prot_nr_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_prot_nr_n` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.CancerProteinNegative = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_type_r_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_type_r_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_type_r_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_type_r_p` AS
    SELECT
        `proteininfo`.`UPAccession`, `Sequence`
    FROM
        `proteindatabase`.`proteininfo`
        INNER JOIN
        `proteindatabase`.`upacc_cancer`
        ON `proteininfo`.`UPAccession` = `upacc_cancer`.`UPAccession`
    WHERE
        `proteininfo`.`Target` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_type_r_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_type_r_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_type_r_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_type_r_n` AS
    SELECT
        prot.`UPAccession`, prot.`Sequence`
    FROM
        `proteindatabase`.`proteininfo` AS prot LEFT JOIN `proteindatabase`.`upacc_cancer` AS cancer
    ON
        prot.`UPAccession` = cancer.`UPAccession`
    WHERE
        cancer.`UPAccession` IS NULL AND
        prot.`Target` = "Y";

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_type_nr_n`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_type_nr_n` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_type_nr_n`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_type_nr_n` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.CancerTypeNegative = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;

-- -----------------------------------------------------
-- View `proteindatabase`.`ill_cancer_type_nr_p`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `proteindatabase`.`ill_cancer_type_nr_p` ;
DROP TABLE IF EXISTS `proteindatabase`.`ill_cancer_type_nr_p`;
USE `proteindatabase`;
CREATE  OR REPLACE VIEW `proteindatabase`.`ill_cancer_type_nr_p` AS
    SELECT
        DISTINCT
        prot.UPAccession,
        prot.A,
        prot.C,
        prot.D,
        prot.E,
        prot.F,
        prot.G,
        prot.H,
        prot.I,
        prot.K,
        prot.L,
        prot.M,
        prot.P,
        prot.N,
        prot.Q,
        prot.R,
        prot.S,
        prot.T,
        prot.V,
        prot.W,
        prot.Y,
        prot.NegativelyCharged,
        prot.PositivelyCharged,
        prot.Basic,
        prot.Charged,
        prot.Polar,
        prot.NonPolar,
        prot.Aromatic,
        prot.Aliphatic,
        prot.Small,
        prot.Tiny,
        prot.PESTMotif,
        prot.LowComplexity,
        prot.Hydrophobicity,
        prot.Isoelectric,
        prot.ECNumber,
        prot.OGlycosylation,
        prot.NGlycosylation,
        prot.Phosphoserine,
        prot.Phosphothreonine,
        prot.Phosphotyrosine,
        prot.SubcellularLocation,
        prot.TopologicalDomain,
        prot.PredictedSubcellularLocation,
        prot.SignalPeptide,
        prot.TransmembraneHelices,
        prot.AlphaHelices,
        prot.BetaStrands,
        prot.PredictedAlphaHelices,
        prot.PredictedBetaSheets,
        prot.Sequence,
        IFNULL(gv.3Untranslated / NULLIF(gv.Total, 0), 0) AS 3Untranslated,
        IFNULL(gv.5Untranslated / NULLIF(gv.Total, 0), 0) AS 5Untranslated,
        IFNULL(gv.NonSynonymousCoding / NULLIF(gv.Total, 0), 0) AS NonSynonymousCoding,
        IFNULL(gv.SynonymousCoding / NULLIF(gv.Total, 0), 0) AS SynonymousCoding,
        IFNULL(para.Paralogs, 0) AS Paralogs,
        IFNULL(ppi.BinaryPPI, 0) AS BinaryPPI,
        IFNULL(trans.ProteinCodingTranscripts, 1) AS AlternativeTranscripts,
        IFNULL(exp.DS_Embryoid_Body / NULLIF(exp.Total, 0), 0) AS DS_Embryoid_Body,
        IFNULL(exp.DS_Blastocyst / NULLIF(exp.Total, 0), 0) AS DS_Blastocyst,
        IFNULL(exp.DS_Fetus / NULLIF(exp.Total, 0), 0) AS DS_Fetus,
        IFNULL(exp.DS_Neonate / NULLIF(exp.Total, 0), 0) AS DS_Neonate,
        IFNULL(exp.DS_Infant / NULLIF(exp.Total, 0), 0) AS DS_Infant,
        IFNULL(exp.DS_Juvenile / NULLIF(exp.Total, 0), 0) AS DS_Juvenile,
        IFNULL(exp.DS_Adult / NULLIF(exp.Total, 0), 0) AS DS_Adult,
        IFNULL(exp.HS_Adrenal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Adrenal_Tumor,
        IFNULL(exp.HS_Bladder_Carcinoma / NULLIF(exp.Total, 0), 0) AS HS_Bladder_Carcinoma,
        IFNULL(exp.HS_Breast_Mammary_Gland_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Breast_Mammary_Gland_Tumor,
        IFNULL(exp.HS_Cervical_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Cervical_Tumor,
        IFNULL(exp.HS_Chondrosarcoma / NULLIF(exp.Total, 0), 0) AS HS_Chondrosarcoma,
        IFNULL(exp.HS_Colorectal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Colorectal_Tumor,
        IFNULL(exp.HS_Esophageal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Esophageal_Tumor,
        IFNULL(exp.HS_Gastrointestinal_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Gastrointestinal_Tumor,
        IFNULL(exp.HS_Germ_Cell_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Germ_Cell_Tumor,
        IFNULL(exp.HS_Glioma / NULLIF(exp.Total, 0), 0) AS HS_Glioma,
        IFNULL(exp.HS_Head_And_Neck_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Head_And_Neck_Tumor,
        IFNULL(exp.HS_Kidney_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Kidney_Tumor,
        IFNULL(exp.HS_Leukemia_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Leukemia_Tumor,
        IFNULL(exp.HS_Liver_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Liver_Tumor,
        IFNULL(exp.HS_Lung_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Lung_Tumor,
        IFNULL(exp.HS_Lymphoma / NULLIF(exp.Total, 0), 0) AS HS_Lymphoma,
        IFNULL(exp.HS_Non_neoplasia / NULLIF(exp.Total, 0), 0) AS HS_Non_neoplasia,
        IFNULL(exp.HS_Normal / NULLIF(exp.Total, 0), 0) AS HS_Normal,
        IFNULL(exp.HS_Ovarian_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Ovarian_Tumor,
        IFNULL(exp.HS_Pancreatic_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Pancreatic_Tumor,
        IFNULL(exp.HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS / NULLIF(exp.Total, 0), 0) AS HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS,
        IFNULL(exp.HS_Prostate_Cancer / NULLIF(exp.Total, 0), 0) AS HS_Prostate_Cancer,
        IFNULL(exp.HS_Retinoblastoma / NULLIF(exp.Total, 0), 0) AS HS_Retinoblastoma,
        IFNULL(exp.HS_Skin_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Skin_Tumor,
        IFNULL(exp.HS_Soft_Tissue_Muscle_Tissue_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Soft_Tissue_Muscle_Tissue_Tumor,
        IFNULL(exp.HS_Uterine_Tumor / NULLIF(exp.Total, 0), 0) AS HS_Uterine_Tumor,
        IFNULL(exp.BS_Adipose_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Adipose_Tissue,
        IFNULL(exp.BS_Adrenal_Gland / NULLIF(exp.Total, 0), 0) AS BS_Adrenal_Gland,
        IFNULL(exp.BS_Ascites / NULLIF(exp.Total, 0), 0) AS BS_Ascites,
        IFNULL(exp.BS_Bladder / NULLIF(exp.Total, 0), 0) AS BS_Bladder,
        IFNULL(exp.BS_Blood / NULLIF(exp.Total, 0), 0) AS BS_Blood,
        IFNULL(exp.BS_Bone / NULLIF(exp.Total, 0), 0) AS BS_Bone,
        IFNULL(exp.BS_Bone_Marrow / NULLIF(exp.Total, 0), 0) AS BS_Bone_Marrow,
        IFNULL(exp.BS_Brain / NULLIF(exp.Total, 0), 0) AS BS_Brain,
        IFNULL(exp.BS_Cervix / NULLIF(exp.Total, 0), 0) AS BS_Cervix,
        IFNULL(exp.BS_Connective_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Connective_Tissue,
        IFNULL(exp.BS_Ear / NULLIF(exp.Total, 0), 0) AS BS_Ear,
        IFNULL(exp.BS_Embryonic_Tissue / NULLIF(exp.Total, 0), 0) AS BS_Embryonic_Tissue,
        IFNULL(exp.BS_Esophagus / NULLIF(exp.Total, 0), 0) AS BS_Esophagus,
        IFNULL(exp.BS_Eye / NULLIF(exp.Total, 0), 0) AS BS_Eye,
        IFNULL(exp.BS_Heart / NULLIF(exp.Total, 0), 0) AS BS_Heart,
        IFNULL(exp.BS_Intestine / NULLIF(exp.Total, 0), 0) AS BS_Intestine,
        IFNULL(exp.BS_Kidney / NULLIF(exp.Total, 0), 0) AS BS_Kidney,
        IFNULL(exp.BS_Larynx / NULLIF(exp.Total, 0), 0) AS BS_Larynx,
        IFNULL(exp.BS_Liver / NULLIF(exp.Total, 0), 0) AS BS_Liver,
        IFNULL(exp.BS_Lung / NULLIF(exp.Total, 0), 0) AS BS_Lung,
        IFNULL(exp.BS_Lymph / NULLIF(exp.Total, 0), 0) AS BS_Lymph,
        IFNULL(exp.BS_Lymph_Node / NULLIF(exp.Total, 0), 0) AS BS_Lymph_Node,
        IFNULL(exp.BS_Mammary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Mammary_Gland,
        IFNULL(exp.BS_Mouth / NULLIF(exp.Total, 0), 0) AS BS_Mouth,
        IFNULL(exp.BS_Muscle / NULLIF(exp.Total, 0), 0) AS BS_Muscle,
        IFNULL(exp.BS_Nerve / NULLIF(exp.Total, 0), 0) AS BS_Nerve,
        IFNULL(exp.BS_Ovary / NULLIF(exp.Total, 0), 0) AS BS_Ovary,
        IFNULL(exp.BS_Pancreas / NULLIF(exp.Total, 0), 0) AS BS_Pancreas,
        IFNULL(exp.BS_Parathyroid / NULLIF(exp.Total, 0), 0) AS BS_Parathyroid,
        IFNULL(exp.BS_Pharynx / NULLIF(exp.Total, 0), 0) AS BS_Pharynx,
        IFNULL(exp.BS_Pituitary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Pituitary_Gland,
        IFNULL(exp.BS_Placenta / NULLIF(exp.Total, 0), 0) AS BS_Placenta,
        IFNULL(exp.BS_Prostate / NULLIF(exp.Total, 0), 0) AS BS_Prostate,
        IFNULL(exp.BS_Salivary_Gland / NULLIF(exp.Total, 0), 0) AS BS_Salivary_Gland,
        IFNULL(exp.BS_Skin / NULLIF(exp.Total, 0), 0) AS BS_Skin,
        IFNULL(exp.BS_Spleen / NULLIF(exp.Total, 0), 0) AS BS_Spleen,
        IFNULL(exp.BS_Stomach / NULLIF(exp.Total, 0), 0) AS BS_Stomach,
        IFNULL(exp.BS_Testis / NULLIF(exp.Total, 0), 0) AS BS_Testis,
        IFNULL(exp.BS_Thymus / NULLIF(exp.Total, 0), 0) AS BS_Thymus,
        IFNULL(exp.BS_Thyroid / NULLIF(exp.Total, 0), 0) AS BS_Thyroid,
        IFNULL(exp.BS_Tonsil / NULLIF(exp.Total, 0), 0) AS BS_Tonsil,
        IFNULL(exp.BS_Trachea / NULLIF(exp.Total, 0), 0) AS BS_Trachea,
        IFNULL(exp.BS_Umbilical_Cord / NULLIF(exp.Total, 0), 0) AS BS_Umbilical_Cord,
        IFNULL(exp.BS_Uterus / NULLIF(exp.Total, 0), 0) AS BS_Uterus,
        IFNULL(exp.BS_Vascular / NULLIF(exp.Total, 0), 0) AS BS_Vascular,
        stab.HalfLife,
        stab.InstabilityIndex
    FROM
        proteininfo AS prot,
        nonredundant AS nr,
        stability AS stab,
        upacc_expression AS exp,
        upacc_variantsnumber AS gv,
        upacc_paralogs AS para,
        upacc_ppi AS ppi,
        upacc_transcripts AS trans
    WHERE
        prot.UPAccession = nr.UPAccession AND
        nr.CancerTypePositive = "Y" AND
        prot.UPAccession = stab.UPAccession AND
        prot.UPAccession = exp.UPAccession AND
        prot.UPAccession = gv.UPAccession AND
        prot.UPAccession = para.UPAccession AND
        prot.UPAccession = ppi.UPAccession AND
        prot.UPAccession = trans.UPAccession
    GROUP BY
        prot.UPAccession;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;
