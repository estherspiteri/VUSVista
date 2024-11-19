import { IPhenotype } from "./phenotype.model";

export interface IVusSample {
  id: string;
  hgvs: string;
  noOfVariants: number;
  consequence: string;
}

export interface INotVusSample {
  id: string;
  noOfVariants: number;
}

export interface IVus {
  id: number;
  chromosome: string;
  chromosomePosition: number;
  classification: string;
  clinvarClassification?: string;
  clinvarClassificationLastEval?: string;
  clinvarClassificationReviewStatus?: string;
  clinvarCanonicalSpdi?: string;
  clinvarId?: number;
  clinvarVariationId?: string;
  clinvarErrorMsg?: string;
  gene: string;
  altAllele: string;
  refAllele: string;
  rsid?: string;
  rsidDbsnpVerified: boolean;
  rsidDbsnpErrorMsgs: string;
  type: string;
  numHeterozygous?: number;
  numHomozygous?: number;
  samples: IVusSample[];
  notVusSamples: INotVusSample[];
  phenotypes: IPhenotype[];
  acmgRuleIds: number[];
  numOfPublications: number;
}
