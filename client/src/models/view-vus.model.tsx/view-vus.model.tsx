export interface IVus {
  chromosome: string;
  chromosomePosition: string;
  classification: string;
  clinvarClassification?: string;
  clinvarClassificationLastEval?: string;
  clinvarClassificationReviewStatus?: string;
  clinvarCanonicalSpdi?: string;
  clinvarErrorMsg?: string;
  gene: string;
  genotype: string;
  observedAllele: string;
  refAllele: string;
  rsid?: string;
  rsidDbsnpVerified: boolean;
  rsidDbsnpErrorMsgs: string;
  type: string; //TODO: change to enum?
  vusID: number;
}
