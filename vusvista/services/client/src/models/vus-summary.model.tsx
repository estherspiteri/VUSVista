export interface IVUSSummary {
  id: number;
  chromosome: string;
  chromosomePosition: number;
  gene: string;
  altAllele: string;
  refAllele: string;
  rsid?: string;
  rsidDbsnpVerified?: boolean;
  classification?: string;
  isFoundInClinvar?: boolean;
  rsidReviewRequired?: boolean;
}
