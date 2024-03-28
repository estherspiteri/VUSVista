export interface IVUSSummary {
  id: number;
  chromosome: string;
  chromosomePosition: string;
  gene: string;
  altAllele: string;
  refAllele: string;
  rsid?: string;
  rsidDbsnpVerified?: boolean;
}
