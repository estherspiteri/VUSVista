import { IPhenotype } from "./phenotype.model";

export interface IAcmgRuleUpload {
  id: number;
  name: string;
}

export interface IVusUpload {
  chromosome: string;
  chromosomePosition: string;
  classification: string;
  gene: string;
  geneId: number;
  altAllele: string;
  refAllele: string;
  genotype: string;
  type: string; //TODO: change to enum?
  samples: string[];
  phenotypes: IPhenotype[];
  acmgRules: IAcmgRuleUpload[];
  hgvs: string;
}
