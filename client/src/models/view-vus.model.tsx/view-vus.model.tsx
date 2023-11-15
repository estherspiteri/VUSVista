export interface IVus {
  chromosome: string;
  chromosomePosition: string;
  classification: string;
  clinvarClassification?: string;
  clinvarErrorMsg?: string;
  gene: string;
  genotype: string;
  observedAllele: string;
  refAllele: string;
  rsid?: string;
  rsidDbsnpVerified: boolean;
  type: string; //TODO: change to enum?
  vusID: number;
}
