export interface IVus {
  id: number;
  chromosome: string;
  chromosomePosition: string;
  classification: string;
  clinvarClassification?: string;
  clinvarClassificationLastEval?: string;
  clinvarClassificationReviewStatus?: string;
  clinvarCanonicalSpdi?: string;
  clinvarUid?: string;
  clinvarErrorMsg?: string;
  gene: string;
  altAllele: string;
  refAllele: string;
  rsid?: string;
  rsidDbsnpVerified: boolean;
  rsidDbsnpErrorMsgs: string;
  type: string; //TODO: change to enum?
  numHeterozygous?: number;
  numHomozygous?: number;
}
