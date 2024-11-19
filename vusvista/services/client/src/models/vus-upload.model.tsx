import { IPhenotype } from "./phenotype.model";

export interface IAcmgRuleUpload {
  id: number;
  name: string;
}

export interface IVusUpload {
  chromosome: string;
  chromosomePosition: string;
  gene: string;
  geneId: number;
  altAllele?: string | null;
  refAllele?: string | null;
  genotype: string;
  type: string;
  samples: string[];
  phenotypes: IPhenotype[];
  acmgRules: IAcmgRuleUpload[];
  hgvs: string;
  rsid: string;
  literatureLinks: string;
}
